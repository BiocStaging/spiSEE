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
  left <- d$xy[, "x"] < 30; right <- d$xy[, "x"] > 70
  maj <- function(v) names(which.max(table(v)))
  expect_false(maj(cl[left]) == maj(cl[right]))
  nice <- .labelClusters(cl, d$lab)
  expect_match(nice, "^N[12]: ")
})

test_that(".log2Enrichment is ~0 under independence and positive on the diagonal when aligned", {
  g <- rep(c("x", "y"), each = 50)
  e <- .log2Enrichment(g, g)
  expect_gt(e["x", "x"], 0); expect_lt(e["x", "y"], 0)
})

test_that(".augmentEmbedding with lambda 0 is the identity; spatial clustering finds the zones", {
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
test_that("buildSpatialGraph/computeNiches work on a SpatialExperiment and survive subsetting", {
  skip_if_not_installed("SpatialExperiment")
  spe <- .makeTestSPE(n_cells = 300)
  spe <- buildSpatialGraph(spe, k = 10)
  expect_true("spatial_knn" %in% SingleCellExperiment::colPairNames(spe))
  spe <- computeNiches(spe, label_col = "SingleR_label", k = 10, n_niches = 3)
  expect_true(all(c("niche", "niche_entropy") %in% colnames(colData(spe))))
  expect_equal(nrow(SingleCellExperiment::reducedDim(spe, "NicheComposition")), 300L)
  sub <- .subsetBySelection(spe, "1,2")
  expect_true(all(SingleCellExperiment::colPair(sub, "spatial_knn")@to <= ncol(sub)))
})

test_that("SpatialNichePanel constructor defaults", {
  p <- SpatialNichePanel()
  expect_identical(p@PlotType, "niche_map")
  expect_identical(p@K, 20L)
  expect_identical(p@NNiches, 6L)
  expect_identical(p@Selection, "")
})
