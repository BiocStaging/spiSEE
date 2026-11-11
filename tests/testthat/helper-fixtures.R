# Small synthetic CosMx-shaped SpatialExperiment for fast, dependency-light
# smoke tests -- real CosMx data is a ~450MB RDS living outside the repo,
# far too heavy to ship or load in a test suite. This fixture is just big
# enough (a handful of FOVs, a few dozen cells) to exercise the actual code
# paths that don't need real segmentation polygons (.subsetBySelection, panel
# constructors, the live compute functions' graceful-fallback branches).
# Tests that need a real neighbour graph (e.g. the "graph_plot" spatial
# stat) skip unless the real dataset is present locally -- see
# test-panels.R.
.makeTestSPE <- function(n_cells = 40, n_genes = 10, n_fovs = 4, seed = 1) {
    set.seed(seed)
    counts <- matrix(rpois(n_cells * n_genes, lambda = 3), nrow = n_genes,
        dimnames = list(paste0("Gene", seq_len(n_genes)),
            paste0("cell", seq_len(n_cells))))
    coords <- cbind(x = stats::runif(n_cells, 0, 100), y = stats::runif(n_cells, 0, 100))
    cd <- S4Vectors::DataFrame(
        fov           = factor(sample(seq_len(n_fovs), n_cells, replace = TRUE)),
        total         = rpois(n_cells, 50),
        detected      = rpois(n_cells, 8),
        Area_um       = stats::runif(n_cells, 10, 40),
        SingleR_label = factor(sample(c("A", "B", "C"), n_cells, replace = TRUE))
    )
    spe <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts, logcounts = log1p(counts)),
        colData = cd, spatialCoords = coords)
    S4Vectors::metadata(spe)$gene_ranks <- rownames(spe)
    spe
}

# Path to the real analysis-ready dataset, if this happens to be running on
# a machine that has it (see R/launch.R / run_app.R). Tests gated on this
# give real regression coverage locally without breaking CI or a fresh
# checkout that hasn't fetched the (external, ~450MB) data. Two sets of
# relative candidates: devtools::test()/R CMD check both run with the
# working directory at tests/testthat/, but a plain source()'d test file
# (or devtools::load_all() + manual testthat::test_that()) runs from the
# package root -- cover both rather than assume one.
.realCosmxPath <- function() {
    candidates <- c(
        "../../../cos_isee_ready.RDS", "../../../data-raw/cos_isee_ready.RDS",
        "../cos_isee_ready.RDS",       "data-raw/cos_isee_ready.RDS",
        "../../../../cos_isee_ready.RDS"
    )
    found <- candidates[file.exists(candidates)]
    if (length(found) == 0) NA_character_ else found[1]
}

# Small synthetic Xenium-shaped SpatialExperiment, same reasoning as
# .makeTestSPE() above but for .subsetBySelection()/Xenium panel constructors.
.makeTestXeniumSPE <- function(n_cells = 40, n_genes = 10, seed = 1) {
    set.seed(seed)
    counts <- matrix(rpois(n_cells * n_genes, lambda = 3), nrow = n_genes,
        dimnames = list(paste0("Gene", seq_len(n_genes)),
            paste0("cell", seq_len(n_cells))))
    coords <- cbind(x = stats::runif(n_cells, 0, 100), y = stats::runif(n_cells, 0, 100))
    cd <- S4Vectors::DataFrame(
        transcript_counts = rpois(n_cells, 20),
        total_counts       = rpois(n_cells, 25),
        cell_area           = stats::runif(n_cells, 10, 40),
        nucleus_area         = stats::runif(n_cells, 5, 20),
        nucleus_count        = sample(0:1, n_cells, replace = TRUE)
    )
    SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts, logcounts = log1p(counts)),
        colData = cd, spatialCoords = coords)
}

# Path to the small, real, committed Xenium test bundle -- ships with the
# repo (data-raw/xenium_test_data, ~17MB), so tests gated on this always
# run rather than skip -- see test-xenium-panels.R.
.xeniumTestDataPath <- function() {
    candidates <- c("../../data-raw/xenium_test_data", "data-raw/xenium_test_data",
        "../../../iseeSpatial/data-raw/xenium_test_data", "../iseeSpatial/data-raw/xenium_test_data")
    found <- candidates[dir.exists(candidates)]
    if (length(found) == 0) NA_character_ else found[1]
}

# Small synthetic Visium HD-shaped object (columns as produced by the OSTA
# workflow: um2, nuc_um2, sum, detected, Level1/Level0).
.makeTestVisiumHDSPE <- function(n_cells = 60, n_genes = 12, seed = 1) {
    set.seed(seed)
    counts <- matrix(rpois(n_cells * n_genes, lambda = 2), nrow = n_genes,
        dimnames = list(paste0("GENE", seq_len(n_genes)), paste0("cellid_", seq_len(n_cells))))
    coords <- cbind(x = stats::runif(n_cells, 0, 1000), y = stats::runif(n_cells, 0, 1000))
    cd <- S4Vectors::DataFrame(
        um2 = 4 * rpois(n_cells, 30), nuc = stats::runif(n_cells),
        sum = colSums(counts), detected = colSums(counts > 0),
        Level1 = factor(sample(c("Tumor", "T.cells", "Fibroblast"), n_cells, replace = TRUE)))
    cd$nuc_um2 <- cd$um2 * cd$nuc
    SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts, logcounts = log1p(counts)), colData = cd, spatialCoords = coords)
}
