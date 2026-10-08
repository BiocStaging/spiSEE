# ============================================================
# All panel classes in one place (collated first, alphabetically).
#
# Every panel that follows the overview map's selection extends
# SpatialPanelBase and therefore carries a `Selection` slot: a
# comma-separated list of region IDs whose meaning depends on the platform
# (FOV ids on CosMx, cell ids on Xenium / Visium HD; see platform.R). The map
# panel finds receivers by that slot, so panels added at runtime through
# iSEE's "Organize panels" also follow the selection.
# ============================================================

#' Base class for iseeSpatial panels
#'
#' Virtual class extending [iSEE::Panel-class]. It adds a `Selection` slot
#' that the [SpatialMapPanel()] fills when the user presses *Validate*.
#'
#' @section Slots:
#' \describe{
#'   \item{`Selection`}{Comma-separated region identifiers (FOV ids or cell
#'     ids, depending on the platform). Empty string = whole dataset.}
#' }
#' @name SpatialPanelBase-class
#' @aliases SpatialPanelBase-class
#' @export
setClass("SpatialPanelBase", contains = c("Panel", "VIRTUAL"),
  slots = c(Selection = "character"))

setMethod("initialize", "SpatialPanelBase", function(.Object, ...) {
  .Object <- callNextMethod(.Object, ...)
  if (length(.Object@Selection) == 0) .Object@Selection <- ""
  .Object
})

#' @rdname SpatialMapPanel
#' @export
setClass("SpatialMapPanel", contains = "Panel",
  slots = c(PendingIDs = "character", PendingRemoveIDs = "character"))

#' @rdname SpatialPlotPanel
#' @export
setClass("SpatialPlotPanel", contains = "SpatialPanelBase",
  slots = c(ColorBy = "character", UsePolygons = "logical", PointSize = "numeric"))

#' @rdname FOVZoomPanel
#' @export
setClass("FOVZoomPanel", contains = "SpatialPanelBase",
  slots = c(ColorBy = "character"))

#' @rdname SpatialQCPanel
#' @export
setClass("SpatialQCPanel", contains = "SpatialPanelBase",
  slots = c(Metric = "character", UseFences = "logical"))

#' @rdname FOVStatsPanel
#' @export
setClass("FOVStatsPanel", contains = "SpatialPanelBase")

#' @rdname SpatialStatsPanel
#' @export
setClass("SpatialStatsPanel", contains = "SpatialPanelBase",
  slots = c(PlotType = "character", Genes = "character", LabelCol = "character",
            K = "integer", ComputeTrigger = "character"))

#' @rdname SVGPanel
#' @export
setClass("SVGPanel", contains = "SpatialPanelBase",
  slots = c(TopN = "integer", K = "integer", LinkedPanel = "character",
            ComputeTrigger = "character"))

#' @rdname SpatialNichePanel
#' @export
setClass("SpatialNichePanel", contains = "SpatialPanelBase",
  slots = c(PlotType = "character", LabelCol = "character",
            K = "integer", NNiches = "integer", Lambda = "numeric",
            ComputeTrigger = "character"))
