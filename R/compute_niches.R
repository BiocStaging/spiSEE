#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom S4Vectors metadata metadata<-
NULL

# Spatial graph, cellular niches and spatial domains.
#
# The k-NN graph is built once on the whole dataset (stored as a colPair,
# so subsetting the object subsets the graph) and used for two clusterings:
#
#  1. Cellular niches ("neighbourhoods"): each cell is described by the
#     cell-type composition of its k nearest spatial neighbours; k-means on
#     those composition vectors gives recurrent micro-environments.
#     Approach of Schurch et al. 2020 (Cell, CODEX "cellular neighborhoods")
#     and imcRtools::aggregateNeighbors + k-means (OSTA neighbourhood chapter).
#  2. Spatial domains: expression PCs are concatenated with their
#     neighbourhood average (a simplified BANKSY-style augmentation, Singhal
#     et al. 2024 Nat Genet) and clustered with a graph method. lambda = 0
#     gives ordinary expression clusters, so the two can be compared directly.
#
# Everything heavy is a plain function that works offline on the full object
# (see scripts/precompute_niches.R); SpatialNichePanel only displays the
# stored result, and recomputes on the current selection on request.


# ---- matrix-level helpers -------------------------------------------------

# Directed k-NN adjacency (n x n, dgCMatrix, x = 1): row i has a 1 in the
# columns of its k nearest neighbours (self excluded).
.knnAdjacency <- function(coords, k = 10L,
                            BNPARAM = BiocNeighbors::KmknnParam()) {
    coords <- as.matrix(coords)
    n <- nrow(coords)
    if (n < 2L) stop("need at least 2 cells to build a neighbour graph")
    k <- min(as.integer(k), n - 1L)
    nn <- BiocNeighbors::findKNN(coords, k = k, BNPARAM = BNPARAM,
                                get.distance = FALSE)
    Matrix::sparseMatrix(i = rep(seq_len(n), times = k),
        j = as.vector(nn$index),
        x = 1, dims = c(n, n),
        dimnames = list(rownames(coords), rownames(coords)))
}

# Row-normalise so each row sums to 1 (rows with no neighbours stay 0).
.rowNormalise <- function(adj) {
    rs <- Matrix::rowSums(adj)
    rs[rs == 0] <- 1
    Matrix::Diagonal(x = 1 / rs) %*% adj
}

# Cell x label matrix of neighbour label proportions.
.neighbourComposition <- function(adj, labels, include_self = FALSE) {
    labels <- droplevels(factor(labels))
    if (anyNA(labels)) stop("labels contain NA; drop unlabelled cells first")
    if (include_self) adj <- adj + Matrix::Diagonal(nrow(adj))
    onehot <- Matrix::sparseMatrix(i = seq_along(labels),
                                    j = as.integer(labels),
        x = 1, dims = c(length(labels), nlevels(labels)))
    comp <- as.matrix(.rowNormalise(adj) %*% onehot)
    dimnames(comp) <- list(rownames(adj), levels(labels))
    comp
}

# Shannon entropy (bits) of each row of a proportion matrix: 0 = all
# neighbours share one label, log2(L) = perfectly mixed.
.shannonEntropy <- function(p) {
    lp <- ifelse(p > 0, log2(p), 0)
    -rowSums(p * lp)
}

# log2 (observed + 0.5) / (expected + 0.5) of a contingency table, expected
# under independence of rows and columns.
.log2Enrichment <- function(groups, labels) {
    tab <- table(groups, labels)
    expected <- outer(rowSums(tab), colSums(tab)) / sum(tab)
    out <- log2((tab + 0.5) / (expected + 0.5))
    unclass(out)
}

# Name each cluster after its two most over-represented labels, so the
# legend reads "N3: Astro/Oligo" instead of "3".
.labelClusters <- function(cl, labels, prefix = "N") {
    enr <- .log2Enrichment(cl, labels)
    tops <- apply(enr, 1, function(r) {
        top <- names(sort(r, decreasing = TRUE))
        paste(top[seq_len(min(2, length(r)))], collapse = "/")
    })
    ids <- paste0(prefix, rownames(enr))
    stats::setNames(paste0(ids, ": ", tops), rownames(enr))
}

.safeKmeans <- function(x, centers, seed = 1L) {
    n_distinct <- nrow(unique(round(x, 6)))
    centers <- min(as.integer(centers), n_distinct)
    if (centers < 2L)
        stop("not enough distinct neighbourhoods to form 2 clusters")
    # withr::with_seed: reproducible without changing the user's random state
    big <- nrow(x) > 5e4
    withr::with_seed(seed,
        stats::kmeans(x, centers = centers,
            nstart = if (big) 3 else 10, iter.max = if (big) 300 else 100,
            algorithm = if (big) "MacQueen" else "Hartigan-Wong")$cluster)
}

# Concatenate an embedding with its neighbourhood mean, weighted so lambda
# is the share of variance given to the neighbourhood part (BANKSY-style).
.augmentEmbedding <- function(emb, adj, lambda = 0.2) {
    stopifnot(lambda >= 0, lambda <= 1)
    if (lambda == 0) return(emb)
    nbr <- as.matrix(.rowNormalise(adj) %*% emb)
    cbind(sqrt(1 - lambda) * emb, sqrt(lambda) * nbr)
}

# Graph clustering of the rows of an embedding: k-NN graph in feature space,
# symmetrised, then Louvain modularity with a resolution parameter.
.graphCluster <- function(emb, k = 20L, resolution = 1, seed = 1L,
                            BNPARAM = BiocNeighbors::KmknnParam()) {
    adj <- .knnAdjacency(emb, k = k, BNPARAM = BNPARAM)
    sym <- (adj + Matrix::t(adj)) > 0
    g <- igraph::graph_from_adjacency_matrix(sym * 1, mode = "undirected",
        weighted = NULL, diag = FALSE)
    withr::with_seed(seed,
                    igraph::membership(igraph::cluster_louvain(g,
                        resolution = resolution)))
}

# ---- argument checks ------------------------------------------------------

.checkCount <- function(x, name, min = 1) {
    if (length(x) != 1 || is.na(x) || !is.numeric(x) || x < min ||
        x != round(x))
        stop("'", name, "' must be a whole number >= ", min, call. = FALSE)
    as.integer(x)
}

.checkFraction <- function(x, name) {
    if (length(x) != 1 || is.na(x) || !is.numeric(x) || x < 0 || x > 1)
        stop("'", name, "' must be a number between 0 and 1", call. = FALSE)
    x
}

# ---- SpatialExperiment-level API ------------------------------------------

#' Build the spatial k-NN graph once and store it on the object
#'
#' Computes a k-nearest-neighbour graph on `spatialCoords(spe)` and stores it
#' as a column pairing (`colPair(spe, name)`). Because column pairings are
#' subset together with the cells, any region selected later reuses this
#' graph instead of rebuilding it.
#'
#' @param spe A `SpatialExperiment`.
#' @param k Number of neighbours per cell. The default, 20, is that of
#'   [computeNiches()], so the graph is reused there.
#' @param name Name of the stored column pairing.
#' @param BNPARAM A `BiocNeighborParam`; use `BiocNeighbors::AnnoyParam()`
#'   for very large datasets.
#' @return `spe` with `colPair(spe, name)` set.
#' @seealso [computeNiches()], [computeSpatialDomains()].
#' @examples
#' # two zones of cell types on a 20 x 20 grid
#' xy <- as.matrix(expand.grid(x = 1:20, y = 1:20))
#' type <- ifelse(xy[, "x"] <= 10,
#'     c("Tumor", "T.cells")[(seq_len(400) %% 5 == 0) + 1],
#'     c("Fibroblast", "T.cells")[(seq_len(400) %% 4 == 0) + 1])
#' counts <- matrix(seq_len(400 * 10) %% 7, 10, 400,
#'     dimnames = list(paste0("g", 1:10), paste0("c", 1:400)))
#' spe <- SpatialExperiment::SpatialExperiment(
#'     assays = list(counts = counts, logcounts = log1p(counts)),
#'     colData = S4Vectors::DataFrame(Level1 = type), spatialCoords = xy)
#' spe <- buildSpatialGraph(spe, k = 8)
#' SingleCellExperiment::colPairNames(spe)
#' @export
buildSpatialGraph <- function(spe, k = 20L, name = "spatial_knn",
                                BNPARAM = BiocNeighbors::KmknnParam()) {
    k <- .checkCount(k, "k")
    adj <- .knnAdjacency(SpatialExperiment::spatialCoords(spe), k = k,
                        BNPARAM = BNPARAM)
    SingleCellExperiment::colPair(spe, name) <- adj
    S4Vectors::metadata(spe)$spiSEE_graph <- list(name = name, k = k)
    spe
}

# Use the stored graph when it exists and matches k; build otherwise.
.getSpatialAdjacency <- function(spe, k, name = "spatial_knn") {
    meta <- S4Vectors::metadata(spe)$spiSEE_graph
    if (name %in% SingleCellExperiment::colPairNames(spe) &&
        !is.null(meta) && identical(as.integer(meta$k), as.integer(k))) {
        adj <- SingleCellExperiment::colPair(spe, name, asSparse = TRUE)
        return(methods::as(adj, "CsparseMatrix"))
    }
    .knnAdjacency(SpatialExperiment::spatialCoords(spe), k = k)
}

#' Cellular niches from the spatial neighbour graph
#'
#' For every cell, computes the cell-type composition of its `k` nearest
#' spatial neighbours, then groups cells with similar surroundings by
#' k-means. Adds to `colData`: `niche` (factor, named after the two most
#' enriched cell types) and `niche_entropy` (Shannon entropy of the
#' neighbourhood, bits). Stores the composition as `reducedDim(spe,
#' "NicheComposition")` and the niche x cell-type log2 enrichment in
#' `metadata(spe)$spiSEE_niches`.
#'
#' @param spe A `SpatialExperiment`.
#' @param label_col `colData` column with cell-type labels.
#' @param k Neighbours per cell (default 20, as [buildSpatialGraph()] and
#'   the niche panel).
#' @param n_niches Number of niches (k-means centres).
#' @param include_self Count the cell itself in its neighbourhood.
#' @param seed Random seed for k-means.
#' @return `spe` with niche annotations added. Cells with an `NA` label get
#'   `NA` niche.
#' @seealso [buildSpatialGraph()], [SpatialNichePanel()].
#' @examples
#' # two zones of cell types on a 20 x 20 grid
#' xy <- as.matrix(expand.grid(x = 1:20, y = 1:20))
#' type <- ifelse(xy[, "x"] <= 10,
#'     c("Tumor", "T.cells")[(seq_len(400) %% 5 == 0) + 1],
#'     c("Fibroblast", "T.cells")[(seq_len(400) %% 4 == 0) + 1])
#' counts <- matrix(seq_len(400 * 10) %% 7, 10, 400,
#'     dimnames = list(paste0("g", 1:10), paste0("c", 1:400)))
#' spe <- SpatialExperiment::SpatialExperiment(
#'     assays = list(counts = counts, logcounts = log1p(counts)),
#'     colData = S4Vectors::DataFrame(Level1 = type), spatialCoords = xy)
#' spe <- computeNiches(spe, label_col = "Level1", k = 8, n_niches = 2)
#' table(spe$niche, spe$Level1)
#' @export
computeNiches <- function(spe, label_col, k = 20L, n_niches = 6L,
                            include_self = FALSE, seed = 1L) {
    cd <- SummarizedExperiment::colData(spe)
    if (!is.character(label_col) || length(label_col) != 1 ||
        !label_col %in% colnames(cd))
        stop("'label_col' must be the name of one colData column")
    k <- .checkCount(k, "k")
    n_niches <- .checkCount(n_niches, "n_niches", min = 2)
    keep <- !is.na(cd[[label_col]])
    sub <- spe[, keep]
    adj <- .getSpatialAdjacency(sub, k)
    labels <- cd[[label_col]][keep]
    comp <- .neighbourComposition(adj, labels, include_self = include_self)
    cl <- .safeKmeans(comp, n_niches, seed = seed)
    nice <- .labelClusters(cl, labels)

    niche <- rep(NA_character_, ncol(spe))
    niche[keep] <- nice[as.character(cl)]
    entropy <- rep(NA_real_, ncol(spe))
    entropy[keep] <- .shannonEntropy(comp)
    full_comp <- matrix(NA_real_, ncol(spe), ncol(comp),
                        dimnames = list(colnames(spe), colnames(comp)))
    full_comp[keep, ] <- comp

    SummarizedExperiment::colData(spe)$niche <- factor(niche,
        levels = unname(nice[order(as.integer(names(nice)))]))
    SummarizedExperiment::colData(spe)$niche_entropy <- entropy
    SingleCellExperiment::reducedDim(spe, "NicheComposition") <- full_comp
    S4Vectors::metadata(spe)$spiSEE_niches <- list(
        label_col = label_col, k = k, n_niches = length(nice),
        enrichment = .log2Enrichment(niche[keep], labels))
    spe
}

#' Spatial domains: clustering on neighbourhood-augmented expression
#'
#' Clusters cells on their expression PCs concatenated with the average PCs
#' of their spatial neighbours (weight `lambda`; a simplified version of the
#' BANKSY idea). With `lambda = 0` this is plain expression clustering, which
#' is also returned (`expr_cluster`) so the two can be compared side by side.
#'
#' @param spe A `SpatialExperiment` with a `logcounts` assay (or a
#'   `reducedDim` named by `use_dimred`).
#' @param k_spatial Spatial neighbours used for the neighbourhood average
#'   (default 20, as [buildSpatialGraph()], whose graph is then reused).
#' @param lambda Share (0-1) of the neighbourhood part.
#' @param use_dimred Name of an existing `reducedDim` to use; if absent,
#'   PCA is computed on `logcounts` with `scater::calculatePCA()`.
#' @param n_pcs Number of PCs when computing PCA.
#' @param k_graph Neighbours for the feature-space clustering graph.
#' @param resolution Louvain resolution (higher = more domains).
#' @param seed Random seed.
#' @return `spe` with `colData` columns `spatial_domain` and `expr_cluster`.
#'   The parameters are stored in `metadata(spe)$spiSEE_domains`;
#'   *Compute* in [SpatialNichePanel()] reuses `resolution` and `k_graph`.
#' @seealso [buildSpatialGraph()], [SpatialNichePanel()].
#' @examples
#' # two zones of cell types on a 20 x 20 grid
#' xy <- as.matrix(expand.grid(x = 1:20, y = 1:20))
#' type <- ifelse(xy[, "x"] <= 10,
#'     c("Tumor", "T.cells")[(seq_len(400) %% 5 == 0) + 1],
#'     c("Fibroblast", "T.cells")[(seq_len(400) %% 4 == 0) + 1])
#' counts <- matrix(seq_len(400 * 10) %% 7, 10, 400,
#'     dimnames = list(paste0("g", 1:10), paste0("c", 1:400)))
#' spe <- SpatialExperiment::SpatialExperiment(
#'     assays = list(counts = counts, logcounts = log1p(counts)),
#'     colData = S4Vectors::DataFrame(Level1 = type), spatialCoords = xy)
#' spe <- computeSpatialDomains(spe, k_spatial = 8, n_pcs = 3, k_graph = 10)
#' table(spe$spatial_domain)
#' @export
computeSpatialDomains <- function(spe, k_spatial = 20L, lambda = 0.2,
                                    use_dimred = "PCA", n_pcs = 20L,
                                    k_graph = 20L, resolution = 0.8,
                                    seed = 1L) {
    k_spatial <- .checkCount(k_spatial, "k_spatial")
    k_graph <- .checkCount(k_graph, "k_graph")
    n_pcs <- .checkCount(n_pcs, "n_pcs")
    lambda <- .checkFraction(lambda, "lambda")
    if (use_dimred %in% SingleCellExperiment::reducedDimNames(spe)) {
        emb <- SingleCellExperiment::reducedDim(spe, use_dimred)
        emb <- emb[, seq_len(min(n_pcs, ncol(emb))), drop = FALSE]
    } else {
        emb <- withr::with_seed(seed,
            scater::calculatePCA(spe, ncomponents = min(n_pcs, nrow(spe) - 1L)))
    }
    adj <- .getSpatialAdjacency(spe, k_spatial)
    dom <- .graphCluster(.augmentEmbedding(emb, adj, lambda), k = k_graph,
        resolution = resolution, seed = seed)
    ecl <- .graphCluster(emb, k = k_graph, resolution = resolution, seed = seed)
    SummarizedExperiment::colData(spe)$spatial_domain <- factor(paste0("D",
        dom), levels = paste0("D", sort(unique(dom))))
    SummarizedExperiment::colData(spe)$expr_cluster   <- factor(paste0("C",
        ecl), levels = paste0("C", sort(unique(ecl))))
    S4Vectors::metadata(spe)$spiSEE_domains <- list(
        k_spatial = k_spatial, lambda = lambda, k_graph = k_graph,
        resolution = resolution)
    spe
}
