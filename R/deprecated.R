# ============================================================
# Old (pre-0.99) constructors and launchers, kept so existing scripts and
# the deployed apps keep working. Each returns the new unified class.
# Remove after the first Bioconductor release.
# ============================================================

#' Deprecated constructors and launchers
#'
#' Platform-specific panels (`CosMx*Panel`, `Xenium*Panel`) were merged into
#' one class per panel. These wrappers return the new classes:
#'
#' | Old | New |
#' |---|---|
#' | `CosMxFOVMapPanel`, `XeniumMapPanel` | [SpatialMapPanel()] |
#' | `CosMxSpatialPanel`, `CosMxGenePanel`, `XeniumSpatialPanel` | [SpatialPlotPanel()] |
#' | `CosMxFOVPanel` | [FOVZoomPanel()] |
#' | `CosMxQCPanel`, `XeniumQCPanel` | [SpatialQCPanel()] |
#' | `CosMxFOVStatsPanel` | [FOVStatsPanel()] |
#' | `CosMxSpatialStatsPanel`, `XeniumSpatialStatsPanel` | [SpatialStatsPanel()] |
#' | `CosMxSVGPanel`, `XeniumSVGPanel` | [SVGPanel()] |
#' | `iseeSpatialAnalysis`, `iseeXeniumAnalysis`, `iseeXeniumQC` | [iseeSpatial()] |
#'
#' `SelectedFOVs` / `SelectedCells` arguments map to `Selection`.
#'
#' @param ... Arguments of the old functions.
#' @param spe A `SpatialExperiment`.
#' @return The corresponding new object.
#' @name iseeSpatial-deprecated
#' @keywords internal
NULL

.sel <- function(SelectedFOVs = "", SelectedCells = "") {
  if (nzchar(SelectedFOVs)) SelectedFOVs else SelectedCells
}

#' @rdname iseeSpatial-deprecated
#' @export
CosMxFOVMapPanel <- function(...) SpatialMapPanel()
#' @rdname iseeSpatial-deprecated
#' @export
XeniumMapPanel <- function(...) SpatialMapPanel()

#' @rdname iseeSpatial-deprecated
#' @export
CosMxSpatialPanel <- function(ColorBy = "SingleR_label", UsePolygons = FALSE, SelectedFOVs = "", ...)
  SpatialPlotPanel(ColorBy = ColorBy, UsePolygons = UsePolygons, Selection = SelectedFOVs)
#' @rdname iseeSpatial-deprecated
#' @export
CosMxGenePanel <- function(Gene = "", UsePolygons = FALSE, SelectedFOVs = "", ...)
  SpatialPlotPanel(ColorBy = Gene, UsePolygons = UsePolygons, Selection = SelectedFOVs)
#' @rdname iseeSpatial-deprecated
#' @export
XeniumSpatialPanel <- function(ColorBy = "", UsePolygons = FALSE, SelectedCells = "", ...)
  SpatialPlotPanel(ColorBy = ColorBy, UsePolygons = UsePolygons, Selection = SelectedCells)

#' @rdname iseeSpatial-deprecated
#' @export
CosMxFOVPanel <- function(ColorBy = "QC_score", SelectedFOVs = "", ...)
  FOVZoomPanel(ColorBy = ColorBy, Selection = SelectedFOVs)

#' @rdname iseeSpatial-deprecated
#' @export
CosMxQCPanel <- function(Metric = "total", UseFences = TRUE, SelectedFOVs = "", ...)
  SpatialQCPanel(Metric = Metric, UseFences = UseFences, Selection = SelectedFOVs)
#' @rdname iseeSpatial-deprecated
#' @export
XeniumQCPanel <- function(Metric = "total_counts", UseFences = TRUE, SelectedCells = "", ...)
  SpatialQCPanel(Metric = Metric, UseFences = UseFences, Selection = SelectedCells)

#' @rdname iseeSpatial-deprecated
#' @export
CosMxFOVStatsPanel <- function(SelectedFOVs = "", ...) FOVStatsPanel(Selection = SelectedFOVs)

#' @rdname iseeSpatial-deprecated
#' @export
CosMxSpatialStatsPanel <- function(SelectedFOVs = "", PlotType = "graph_plot", Genes = character(0), ...)
  SpatialStatsPanel(PlotType = PlotType, Genes = Genes, Selection = SelectedFOVs)
#' @rdname iseeSpatial-deprecated
#' @export
XeniumSpatialStatsPanel <- function(SelectedCells = "", PlotType = "graph_plot", Genes = character(0), ...)
  SpatialStatsPanel(PlotType = PlotType, Genes = Genes, Selection = SelectedCells)

#' @rdname iseeSpatial-deprecated
#' @export
CosMxSVGPanel <- function(SelectedFOVs = "", TopN = 200L, LinkedGenePanel = "", ...)
  SVGPanel(TopN = TopN, LinkedPanel = LinkedGenePanel, Selection = SelectedFOVs)
#' @rdname iseeSpatial-deprecated
#' @export
XeniumSVGPanel <- function(SelectedCells = "", TopN = 200L, LinkedPanel = "", ...)
  SVGPanel(TopN = TopN, LinkedPanel = LinkedPanel, Selection = SelectedCells)

#' @rdname iseeSpatial-deprecated
#' @export
iseeSpatialAnalysis <- function(spe) {
  .Deprecated("iseeSpatial(spe, mode = 'analysis')")
  iseeSpatial(spe, mode = "analysis")
}
#' @rdname iseeSpatial-deprecated
#' @export
iseeXeniumAnalysis <- function(spe) {
  .Deprecated("iseeSpatial(spe, mode = 'analysis')")
  iseeSpatial(spe, mode = "analysis", platform = "xenium")
}
#' @rdname iseeSpatial-deprecated
#' @export
iseeXeniumQC <- function(spe) {
  .Deprecated("iseeSpatial(spe, mode = 'qc')")
  iseeSpatial(spe, mode = "qc", platform = "xenium")
}
