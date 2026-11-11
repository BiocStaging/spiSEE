test_that("platform is detected from colData and can be overridden", {
    expect_identical(getPlatform(.makeTestSPE()), "cosmx")
    expect_identical(getPlatform(.makeTestXeniumSPE()), "xenium")
    expect_identical(getPlatform(.makeTestVisiumHDSPE()), "visiumhd")
    spe <- setPlatform(.makeTestSPE(), "generic")
    expect_identical(getPlatform(spe), "generic")
    expect_error(setPlatform(spe, "merfish"))
})

test_that("selection unit follows the platform", {
    expect_identical(.selectionUnit(.makeTestSPE()), "fov")
    expect_identical(.selectionUnit(.makeTestXeniumSPE()), "cell")
    expect_identical(.selectionUnit(.makeTestVisiumHDSPE()), "cell")
})

test_that(".subsetBySelection uses FOVs on CosMx and cell ids elsewhere", {
    spe <- .makeTestSPE()
    sub <- .subsetBySelection(spe, "1,2")
    expect_true(all(as.character(colData(sub)$fov) %in% c("1", "2")))
    expect_identical(.subsetBySelection(spe, ""), spe)
    x <- .makeTestXeniumSPE()
    ids <- colnames(x)[1:5]
    expect_setequal(colnames(.subsetBySelection(x,
                                                paste(ids, collapse = ","))),
                                                ids)
})

test_that("QC metrics and label columns put platform columns first", {
    expect_identical(.qcMetrics(.makeTestSPE())[1], "total")
    expect_identical(.qcMetrics(.makeTestXeniumSPE())[1], "transcript_counts")
    expect_identical(.qcMetrics(.makeTestVisiumHDSPE())[1], "sum")
    expect_identical(.defaultLabelCol(.makeTestSPE()), "SingleR_label")
    expect_identical(.defaultLabelCol(.makeTestVisiumHDSPE()), "Level1")
})

test_that(".allRegionIDs returns sorted FOVs or cell ids", {
    expect_identical(.allRegionIDs(.makeTestSPE(n_fovs = 4)), as.character(1:4))
    x <- .makeTestXeniumSPE()
    expect_identical(.allRegionIDs(x), colnames(x))
})
