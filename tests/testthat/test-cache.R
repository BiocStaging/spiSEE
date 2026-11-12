test_that("dataset fingerprint separates datasets and ignores metadata", {
    a <- .makeTestSPE(n_cells = 40, seed = 1)
    b <- .makeTestSPE(n_cells = 40, seed = 2)
    expect_false(identical(.datasetFingerprint(a), .datasetFingerprint(b)))
    a2 <- a
    S4Vectors::metadata(a2)$anything <- 1
    expect_identical(.datasetFingerprint(a), .datasetFingerprint(a2))
    a3 <- a
    colData(a3)$total[1] <- colData(a3)$total[1] + 1
    expect_false(identical(.datasetFingerprint(a), .datasetFingerprint(a3)))
})

test_that("disk cache: folder, size limit, survives a new cache object, clear",
          {
    dir <- withr::local_tempdir()
    local_mocked_bindings(.isDevLoad = function() FALSE)
    c1 <- .configurePlotCache(dir, size = 5e6)
    expect_s3_class(c1, "cachem")
    expect_identical(c1$info()$max_size, 5e6)
    c1$set("k1", 1:10)
    c2 <- .configurePlotCache(dir)            # e.g. after an app restart
    expect_identical(c2$get("k1"), 1:10)
    clearSpatialCache(dir)
    expect_true(cachem::is.key_missing(c2$get("k1")))
    expect_identical(.configurePlotCache(FALSE), "app")
})

test_that("under load_all the cache stays in memory", {
    local_mocked_bindings(.isDevLoad = function() TRUE)
    expect_identical(.configurePlotCache(TRUE), "app")
})

test_that("default image cache folder is separate from the map backgrounds", {
    expect_identical(basename(spatialCacheDir()), "plots")
    expect_identical(dirname(spatialCacheDir()),
                     tools::R_user_dir("iseeSpatial", "cache"))
})

test_that("Compute keys depend on what is computed, not on when", {
    a <- SpatialStatsPanel(PlotType = "localmoran_plot", Genes = "Gene1",
                           K = 6L, Selection = "1,2")
    b <- a
    expect_identical(.computeKey(a), .computeKey(b))
    b@K <- 10L
    expect_false(identical(.computeKey(a), .computeKey(b)))
    c <- a
    c@Selection <- "1,3"
    expect_false(identical(.computeKey(a), .computeKey(c)))
    # display-only slots do not change the key
    d <- a
    d@PanelWidth <- 6L
    expect_identical(.computeKey(a), .computeKey(d))
    expect_match(.computeKey(SVGPanel()), "^k[0-9a-f]+$")
})

test_that("a trigger set in a layout becomes the key of its parameters", {
    spe <- .makeTestSPE()
    x <- .refineParameters(SpatialStatsPanel(Genes = "Gene1",
                                             ComputeTrigger = "launch"), spe)
    y <- .refineParameters(SpatialStatsPanel(Genes = "Gene2",
                                             ComputeTrigger = "launch"), spe)
    expect_false(identical(x@ComputeTrigger, y@ComputeTrigger))
    expect_identical(x@ComputeTrigger, .computeKey(x))
    z <- .refineParameters(SpatialStatsPanel(), spe)
    expect_identical(z@ComputeTrigger, "")
})

test_that("computed results are reused from the cache", {
    old <- .CACHE$cache
    on.exit(.CACHE$cache <- old)
    .configurePlotCache(withr::local_tempdir())
    n <- 0
    f <- function() { n <<- n + 1; data.frame(gene = "A", moran_i = 0.5) }
    key <- paste0("test", as.integer(Sys.time()))
    r1 <- .cachedResult(key, f())
    r2 <- .cachedResult(key, f())
    expect_identical(r1, r2)
    expect_identical(n, 1)
})
