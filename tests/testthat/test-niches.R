# Matrix-level helpers: no SpatialExperiment needed.

.twoZoneData <- function(n = 400, seed = 1) {
    # Left half mostly type A, right half mostly type B, with a mixed strip.
    set.seed(seed)
    xy <- cbind(x = stats::runif(n, 0, 100), y = stats::runif(n, 0, 100))
    p_a <- ifelse(xy[, "x"] < 40, 0.95, ifelse(xy[, "x"] > 60, 0.05, 0.5))
    lab <- factor(ifelse(stats::runif(n) < p_a, "A", "B"))
    rownames(xy) <- paste0("cell", seq_len(n))
    list(xy = xy, lab = lab)
}

test_that(".knnAdjacency gives k neighbours per cell and excludes self", {
    d <- .twoZoneData()
    adj <- .knnAdjacency(d$xy, k = 7)
    expect_equal(dim(adj), c(400L, 400L))
    expect_true(all(Matrix::rowSums(adj) == 7))
    expect_true(all(Matrix::diag(adj) == 0))
})

test_that(".neighbourComposition rows are proportions", {
    d <- .twoZoneData()
    comp <- .neighbourComposition(.knnAdjacency(d$xy, 10), d$lab)
    expect_equal(colnames(comp), c("A", "B"))
    expect_equal(unname(rowSums(comp)), rep(1, 400))
    # cells deep in the left zone see mostly A
    left <- d$xy[, "x"] < 30
    expect_gt(mean(comp[left, "A"]), 0.8)
})

test_that(".shannonEntropy is 0 for pure and 1 bit for 50/50", {
    p <- rbind(c(1, 0), c(0.5, 0.5))
    expect_equal(.shannonEntropy(p), c(0, 1))
})

test_that("niche k-means separates the two zones", {
    d <- .twoZoneData()
    comp <- .neighbourComposition(.knnAdjacency(d$xy, 15), d$lab)
    cl <- .safeKmeans(comp, 2)
    left <- d$xy[, "x"] < 30
    right <- d$xy[, "x"] > 70
    maj <- function(v) names(which.max(table(v)))
    expect_false(maj(cl[left]) == maj(cl[right]))
    nice <- .labelClusters(cl, d$lab)
    expect_match(nice, "^N[12]: ")
})

test_that(".log2Enrichment: ~0 if independent, > 0 on an aligned diagonal", {
    g <- rep(c("x", "y"), each = 50)
    e <- .log2Enrichment(g, g)
    expect_gt(e["x", "x"], 0)
    expect_lt(e["x", "y"], 0)
})

test_that("lambda 0 is plain expression; domains find the zones", {
    d <- .twoZoneData(n = 600)
    set.seed(2)
    # Expression signal follows the zone weakly; neighbourhood averaging
    # should sharpen it.
    emb <- cbind(PC1 = (d$xy[, "x"] > 50) * 1 + stats::rnorm(600, sd = 0.8),
        PC2 = stats::rnorm(600))
    adj <- .knnAdjacency(d$xy, 15)
    expect_identical(.augmentEmbedding(emb, adj, 0), emb)
    aug <- .augmentEmbedding(emb, adj, 0.5)
    expect_equal(ncol(aug), 4L)
    dom <- .graphCluster(aug, k = 20, resolution = 0.3)
    expect_gte(length(unique(dom)), 2L)
    # purity of the left/right split is better with spatial augmentation
    purity <- function(cl) {
        tab <- table(cl, d$xy[, "x"] > 50)
        sum(apply(tab, 1, max)) / sum(tab)
    }
    expect_gt(purity(.safeKmeans(aug, 2)), purity(.safeKmeans(emb, 2)))
})

# SpatialExperiment-level API (needs SpatialExperiment; runs on your machine).
test_that("graph and niches work on an SPE and survive subsetting", {
    skip_if_not_installed("SpatialExperiment")
    spe <- .makeTestSPE(n_cells = 300)
    spe <- buildSpatialGraph(spe, k = 10)
    expect_true("spatial_knn" %in% SingleCellExperiment::colPairNames(spe))
    spe <- computeNiches(spe, label_col = "SingleR_label", k = 10, n_niches = 3)
    expect_true(all(c("niche", "niche_entropy") %in% colnames(colData(spe))))
    expect_equal(nrow(SingleCellExperiment::reducedDim(spe,
        "NicheComposition")), 300L)
    sub <- .subsetBySelection(spe, "1,2")
    expect_true(all(SingleCellExperiment::colPair(sub,
        "spatial_knn")@to <= ncol(sub)))
})

test_that("SpatialNichePanel constructor defaults", {
    p <- SpatialNichePanel()
    expect_identical(p@PlotType, "niche_map")
    # NA until resolved against the data (stored niches, else 20 / 6)
    expect_identical(p@K, NA_integer_)
    expect_identical(p@NNiches, NA_integer_)
    expect_identical(p@Selection, "")
})

test_that("niche results are reused when only the plot type changes", {
    skip_if_not_installed("SpatialExperiment")
    spe <- .makeTestSPE(n_cells = 300)
    calls <- 0
    real <- .nicheData
    local_mocked_bindings(.nicheData = function(x, se) {
        calls <<- calls + 1
        real(x, se)
    })
    cache <- new.env()
    p <- SpatialNichePanel(LabelCol = "SingleR_label", K = 10L, NNiches = 3L,
                           ComputeTrigger = "t1")
    out1 <- .nicheOutput(p, spe, cache, note = "Computed on 300 cells")
    p@PlotType <- "niche_heatmap"
    out2 <- .nicheOutput(p, spe, cache)
    p@PlotType <- "niche_bar"
    out3 <- .nicheOutput(p, spe, cache)
    expect_equal(calls, 1)
    expect_s3_class(out1, "ggplot")
    expect_identical(out1$labels$caption, "Computed on 300 cells")
    expect_s3_class(out2, "pheatmap")
})

test_that("niches are reproducible; the random state is left alone", {
    set.seed(42)
    before <- .Random.seed
    x <- matrix(c(rnorm(50), rnorm(50, 5)), ncol = 2)
    set.seed(42)
    invisible(runif(1))
    state <- .Random.seed
    a <- .safeKmeans(x, 2, seed = 7)
    expect_identical(.Random.seed, state)
    expect_identical(a, .safeKmeans(x, 2, seed = 7))
})

test_that("the niche panel opens on the stored niches' parameters", {
    spe <- .makeTestSPE(n_cells = 300)
    spe$cell_type <- sample(c("A", "B", "C"), ncol(spe), replace = TRUE)
    spe <- computeNiches(spe, label_col = "cell_type", k = 7, n_niches = 3)
    x <- .refineParameters(SpatialNichePanel(), spe)
    expect_identical(x@K, 7L)
    expect_identical(x@NNiches, 3L)
    expect_identical(x@LabelCol, "cell_type")
    y <- .refineParameters(SpatialNichePanel(K = 10L), spe)
    expect_identical(y@K, 10L)
    z <- .refineParameters(SpatialNichePanel(), .makeTestSPE())
    expect_identical(c(z@K, z@NNiches), c(20L, 6L))
})

test_that("without a stored result, nothing is computed before Compute", {
    spe <- .makeTestSPE(n_cells = 300)
    calls <- 0
    local_mocked_bindings(
        computeNiches = function(...) { calls <<- calls + 1; stop("called") },
        computeSpatialDomains = function(...) {
            calls <<- calls + 1
            stop("called")
        })
    p <- SpatialNichePanel(LabelCol = "SingleR_label", K = 10L, NNiches = 3L)
    out <- .nicheOutput(p, spe)
    p@PlotType <- "domain_map"
    out2 <- .nicheOutput(p, spe)
    expect_equal(calls, 0)
    expect_identical(out$labels$title, .PRESS_COMPUTE)
    expect_identical(out2$labels$title, .PRESS_COMPUTE)
})

test_that("niche and graph functions check their arguments", {
    spe <- .makeTestSPE(n_cells = 60)
    expect_error(computeNiches(spe, "SingleR_label", k = NA), "'k'")
    expect_error(computeNiches(spe, c("a", "b")), "label_col")
    expect_error(computeNiches(spe, "SingleR_label", n_niches = 1),
                    "n_niches")
    expect_error(buildSpatialGraph(spe, k = -1), "'k'")
    expect_error(computeSpatialDomains(spe, lambda = 2), "lambda")
})
