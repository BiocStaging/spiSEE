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

#' @param x A panel instance (the class this page documents).
#' @param se The `SpatialExperiment` being explored.
#' @param select_info A list of selection information, passed by iSEE.
#' @param output,pObjects,rObjects The app's Shiny output list and the
#'   persistent/reactive state environments iSEE passes to every panel
#'   method.
#' @param input,session The current Shiny `input` and `session` objects.
#' @param all_memory,all_contents iSEE's full panel-state list and the
#'   current selection contents of every panel, passed to
#'   `.generateOutput()` for export/snapshot but unused by this package's
#'   own panels (selection is handled client-side, see `SpatialMapPanel`).
#' @rdname SpatialMapPanel
#' @export
setClass("SpatialMapPanel", contains = "Panel",
  # CommittedIDs mirrors what was last sent to the panels (i.e. the browser's
  # `committed` set) -- without it, Validate only ever knew about this round's
  # PendingIDs/PendingRemoveIDs and silently dropped everything committed in
  # earlier rounds (see panel_map.R's Validate observer).
  slots = c(PendingIDs = "character", PendingRemoveIDs = "character",
            CommittedIDs = "character"))

#' @inheritParams SpatialMapPanel
#' @rdname SpatialPlotPanel
#' @export
setClass("SpatialPlotPanel", contains = "SpatialPanelBase",
  slots = c(ColorBy = "character", UsePolygons = "logical", PointSize = "numeric"))

#' @inheritParams SpatialMapPanel
#' @rdname FOVZoomPanel
#' @export
setClass("FOVZoomPanel", contains = "SpatialPanelBase",
  slots = c(ColorBy = "character"))

#' @inheritParams SpatialMapPanel
#' @rdname SpatialQCPanel
#' @export
setClass("SpatialQCPanel", contains = "SpatialPanelBase",
  slots = c(Metric = "character", UseFences = "logical"))

#' @inheritParams SpatialMapPanel
#' @rdname FOVStatsPanel
#' @export
setClass("FOVStatsPanel", contains = "SpatialPanelBase")

#' @inheritParams SpatialMapPanel
#' @rdname SpatialStatsPanel
#' @export
setClass("SpatialStatsPanel", contains = "SpatialPanelBase",
  slots = c(PlotType = "character", Genes = "character", LabelCol = "character",
            K = "integer", ComputeTrigger = "character"))

#' @inheritParams SpatialMapPanel
#' @rdname SVGPanel
#' @export
setClass("SVGPanel", contains = "SpatialPanelBase",
  slots = c(TopN = "integer", K = "integer", LinkedPanel = "character",
            ComputeTrigger = "character"))

#' @inheritParams SpatialMapPanel
#' @rdname SpatialNichePanel
#' @export
setClass("SpatialNichePanel", contains = "SpatialPanelBase",
  slots = c(PlotType = "character", LabelCol = "character",
            K = "integer", NNiches = "integer", Lambda = "numeric",
            ComputeTrigger = "character", UsePolygons = "logical"))
