test_that("dataset fingerprint separates datasets and ignores metadata", {
  a <- .makeTestSPE(n_cells = 40, seed = 1)
  b <- .makeTestSPE(n_cells = 40, seed = 2)
  expect_false(identical(.datasetFingerprint(a), .datasetFingerprint(b)))
  a2 <- a; S4Vectors::metadata(a2)$anything <- 1
  expect_identical(.datasetFingerprint(a), .datasetFingerprint(a2))
  a3 <- a; colData(a3)$total[1] <- colData(a3)$total[1] + 1
  expect_false(identical(.datasetFingerprint(a), .datasetFingerprint(a3)))
})

test_that("disk cache: folder, size limit, survives a new cache object, clear", {
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
