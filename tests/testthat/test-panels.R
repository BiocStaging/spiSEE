test_that("panel constructors set their defaults", {
  expect_identical(SpatialMapPanel()@PendingIDs, "")
  p <- SpatialPlotPanel(); expect_identical(p@ColorBy, ""); expect_identical(p@Selection, "")
  expect_identical(SpatialQCPanel()@UseFences, TRUE)
  s <- SpatialStatsPanel(); expect_identical(s@PlotType, "graph_plot"); expect_identical(s@K, 6L)
  v <- SVGPanel(K = 12); expect_identical(v@K, 12L); expect_identical(v@TopN, 200L)
  n <- SpatialNichePanel(K = 30); expect_identical(n@K, 30L)
})

test_that("every panel reports a name and a valid hex colour", {
  panels <- list(SpatialMapPanel(), SpatialPlotPanel(), FOVZoomPanel(), SpatialQCPanel(),
                 FOVStatsPanel(), SpatialStatsPanel(), SVGPanel(), SpatialNichePanel())
  for (p in panels) {
    expect_type(.fullName(p), "character")
    expect_match(.panelColor(p), "^#[0-9A-Fa-f]{6}$")
  }
})

test_that("all selection-following panels share the Selection slot", {
  for (p in list(SpatialPlotPanel(), FOVZoomPanel(), SpatialQCPanel(), FOVStatsPanel(),
                 SpatialStatsPanel(), SVGPanel(), SpatialNichePanel()))
    expect_true(is(p, "SpatialPanelBase"))
  expect_false(is(SpatialMapPanel(), "SpatialPanelBase"))
})

test_that("deprecated constructors return the new classes with the selection mapped", {
  expect_s4_class(CosMxSpatialPanel(SelectedFOVs = "1,2"), "SpatialPlotPanel")
  expect_identical(CosMxSpatialPanel(SelectedFOVs = "1,2")@Selection, "1,2")
  expect_identical(CosMxGenePanel(Gene = "Gene1")@ColorBy, "Gene1")
  expect_s4_class(XeniumQCPanel(), "SpatialQCPanel")
  expect_identical(XeniumSVGPanel(LinkedPanel = "SpatialPlotPanel1")@LinkedPanel, "SpatialPlotPanel1")
  expect_s4_class(CosMxFOVMapPanel(), "SpatialMapPanel")
})

test_that("default layouts adapt to the platform", {
  cos_qc <- .defaultPanels(.makeTestSPE(), "qc")
  expect_true(any(vapply(cos_qc, is, logical(1), "FOVZoomPanel")))
  xen_qc <- .defaultPanels(.makeTestXeniumSPE(), "qc")
  expect_false(any(vapply(xen_qc, is, logical(1), "FOVZoomPanel")))
  an <- .defaultPanels(.makeTestVisiumHDSPE(), "analysis")
  expect_true(any(vapply(an, is, logical(1), "SpatialNichePanel")))
  expect_identical(an[[2]]@ColorBy, "Level1")
})

test_that("plot panels render on synthetic data without the app", {
  skip_if_not_installed("SpaceTrooper")
  spe <- .makeTestVisiumHDSPE()
  expect_s3_class(.generateOutput(SpatialPlotPanel(ColorBy = "sum"), spe), "ggplot")
  expect_s3_class(.generateOutput(SpatialPlotPanel(ColorBy = "GENE1"), spe), "ggplot")
  p <- .generateOutput(SpatialPlotPanel(UsePolygons = TRUE), spe)
  expect_s3_class(p, "ggplot")    # graceful message, no polygons loaded
  expect_s3_class(.generateOutput(FOVZoomPanel(), spe), "ggplot")  # not a tile platform
})

test_that("stats degrade gracefully on too few cells and honour k", {
  spe <- .makeTestSPE(n_cells = 3)
  expect_s3_class(.computeSpatialStat(spe, "graph_plot", k = 6), "ggplot")
  tbl <- .computeSVGRanking(spe)
  expect_equal(nrow(tbl), 0); expect_named(tbl, c("gene", "moran_i"))
  expect_s3_class(.computeSpatialStat(.makeTestSPE(), "joincount_plot", group_col = NA), "ggplot")
})

# ── Integration on real data, when present locally ──
test_that("stats and SVG ranking run on the real CosMx data", {
  path <- .realCosmxPath()
  skip_if(is.na(path), "real CosMx dataset not present on this machine")
  spe <- setPlatform(readRDS(path), "cosmx")
  sub <- .subsetBySelection(spe, "72,73,74")
  expect_s3_class(.computeSpatialStat(sub, "graph_plot", k = 10), "ggplot")
  tbl <- .computeSVGRanking(sub, top_n = 20, k = 10)
  expect_gt(nrow(tbl), 0); expect_true(all(diff(tbl$moran_i) <= 0))
})

test_that("Xenium test bundle: loader, k-NN graph and known top SVG", {
  path <- .xeniumTestDataPath()
  skip_if(is.na(path), "bundled Xenium test data not found")
  skip_if_not_installed("SpatialExperimentIO")
  spe <- load_xenium(path)
  expect_identical(getPlatform(spe), "xenium")
  expect_equal(ncol(spe), 36)
  expect_s3_class(.computeSpatialStat(spe, "graph_plot"), "ggplot")
  tbl <- .computeSVGRanking(spe, top_n = 50)
  expect_identical(tbl$gene[1], "Defa5")   # verified on this bundle with k = 6
})

test_that("genes with a hyphen in their name can be plotted", {
  skip_if_not_installed("SpaceTrooper")
  spe <- .makeTestVisiumHDSPE()
  rownames(spe)[1] <- "MT-ATP6"
  p <- .generateOutput(SpatialPlotPanel(ColorBy = "MT-ATP6"), spe)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("SVG tables are kept per session, not shared", {
  a <- new.env(); b <- new.env()
  assign("p", data.frame(gene = "A"), envir = .svgTables(a))
  assign("p", data.frame(gene = "B"), envir = .svgTables(b))
  expect_identical(get("p", envir = .svgTables(a))$gene, "A")
  expect_identical(get("p", envir = .svgTables(b))$gene, "B")
})

test_that("the launcher offers the spatial panels in 'Organize panels'", {
  expect_s3_class(iseeSpatial(.makeTestVisiumHDSPE(), mode = "qc"), "shiny.appobj")
  expect_false(any(vapply(.extraPanels(.makeTestVisiumHDSPE()), is, logical(1), "FOVZoomPanel")))
  expect_true(any(vapply(.extraPanels(.makeTestSPE()), is, logical(1), "FOVZoomPanel")))
})

test_that("join counts use binary symmetric weights (F14)", {
  sfe <- .toSFEWithKNNNeighbors(.makeTestSPE(n_cells = 120), k = 6L)
  lw <- .joincountWeights(sfe)
  expect_identical(lw$style, "B")
  m <- spdep::listw2mat(lw)
  expect_true(isSymmetric(unname(m)))
  expect_true(all(m %in% c(0, 1)))
  df <- .joincountDF(sfe, "SingleR_label")
  expect_true(all(df$Variance[!is.na(df$Variance)] >= 0))
})
