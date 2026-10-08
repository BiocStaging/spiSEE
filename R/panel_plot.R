# ============================================================
# SpatialPlotPanel - cells on the tissue, coloured by an annotation, a QC
# metric or a gene. Replaces CosMxSpatialPanel, CosMxGenePanel and
# XeniumSpatialPanel (one class, any platform).
# ============================================================

setMethod("initialize", "SpatialPlotPanel", function(.Object, ...) {
  .Object <- callNextMethod(.Object, ...)
  if (length(.Object@ColorBy)     == 0) .Object@ColorBy     <- ""
  if (length(.Object@UsePolygons) == 0) .Object@UsePolygons <- FALSE
  if (length(.Object@PointSize)   == 0) .Object@PointSize   <- 0.5
  .Object
})

#' @export
setMethod(".fullName",   "SpatialPlotPanel", function(x) "Spatial Plot")
#' @rdname SpatialPlotPanel
#' @export
setMethod(".panelColor", "SpatialPlotPanel", function(x) "#EF476F")
#' @rdname SpatialPlotPanel
#' @export
setMethod(".defineOutput", "SpatialPlotPanel", function(x)
  shiny::plotOutput(.getEncodedName(x), height = "600px"))

# Menu content: annotations, then QC metrics, then genes.
.colorByChoices <- function(se) {
  c(.labelColumns(se), grep("^PC", colnames(colData(se)), value = TRUE), .qcMetrics(se), .geneChoices(se))
}

# "" means: first label column, else first QC metric, else first gene.
.resolveColorBy <- function(x, se) {
  cb <- x@ColorBy
  if (nzchar(cb) && (cb %in% colnames(colData(se)) || cb %in% rownames(se))) return(cb)
  cand <- c(.defaultLabelCol(se), .qcMetrics(se), .geneChoices(se))
  cand[!is.na(cand)][1]
}

#' @rdname SpatialPlotPanel
#' @export
setMethod(".defineInterface", "SpatialPlotPanel", function(x, se, select_info) {
  enc <- .getEncodedName(x)
  list(iSEE::collapseBox(
    id = paste0(enc, "_DataBoxOpen"), title = "Data parameters", open = TRUE,
    # filled server-side in .createObservers (gene lists can be very long)
    shiny::selectizeInput(paste0(enc, "_ColorBy"), "Colour by (annotation, QC metric or gene):",
                          choices = .resolveColorBy(x, se), selected = .resolveColorBy(x, se)),
    shiny::checkboxInput(paste0(enc, "_UsePolygons"), "Draw cell polygons (slow on large regions)",
                         value = x@UsePolygons),
    shiny::sliderInput(paste0(enc, "_PointSize"), "Point size:", min = 0.1, max = 3,
                       value = x@PointSize, step = 0.1),
    .selectionHelp(se)
  ))
})

#' @rdname SpatialPlotPanel
#' @export
setMethod(".generateOutput", "SpatialPlotPanel", function(x, se, ...) {
  colorby <- .resolveColorBy(x, se)
  se <- .subsetBySelection(se, x@Selection)
  if (ncol(se) == 0) return(.notAvailablePlot("No cells in this selection"))
  if (identical(colorby, "fov") && "fov" %in% colnames(colData(se))) return(.plotFOVOverview(se))
  se <- .forceCategorical(se, colorby)
  if (!colorby %in% colnames(colData(se)) && colorby %in% rownames(se))
    colData(se)[[colorby]] <- as.numeric(SingleCellExperiment::logcounts(se)[colorby, ])

  args <- list(se, colourBy = colorby)
  # Stored label colours (CosMx: SingleR_label_colors) are reused when present.
  pal_col <- paste0(colorby, "_colors")
  if (pal_col %in% colnames(colData(se))) {
    map <- unique(data.frame(l = as.character(colData(se)[[colorby]]),
                             c = as.character(colData(se)[[pal_col]])))
    args$palette <- stats::setNames(map$c, map$l)
  }
  if (.selectionUnit(se) == "fov") args["sampleId"] <- list(NULL)   # as the CosMx panels always did

  if (x@UsePolygons) {
    # plotPolygons() stops without a "polygons" column, and that error used
    # to hang the panel silently instead of surfacing: check first.
    if (!"polygons" %in% colnames(colData(se)))
      return(.notAvailablePlot("No cell polygons loaded for this dataset: untick 'Draw cell polygons'"))
    do.call(SpaceTrooper::plotPolygons, args)
  } else {
    do.call(SpaceTrooper::plotCentroids, c(args, size = x@PointSize))
  }
})

#' @rdname SpatialPlotPanel
#' @export
setMethod(".renderOutput", "SpatialPlotPanel", function(x, se, output, pObjects, rObjects) {
  .renderCachedSpatialPlot(x, se, output, pObjects, rObjects,
                           c("ColorBy", "UsePolygons", "PointSize", "Selection"))
})

#' @rdname SpatialPlotPanel
#' @export
setMethod(".createObservers", "SpatialPlotPanel", function(x, se, input, session, pObjects, rObjects) {
  callNextMethod(); enc <- .getEncodedName(x)
  .serverSelectize(session, enc, "ColorBy", .colorByChoices(se), .resolveColorBy(x, se))
  .slotObserver(enc, "ColorBy", input, pObjects, rObjects, cast = as.character)
  .slotObserver(enc, "UsePolygons", input, pObjects, rObjects, cast = as.logical)
  .slotObserver(enc, "PointSize", input, pObjects, rObjects, cast = as.numeric)
})

#' Spatial plot of cells coloured by annotation, QC metric or gene
#'
#' @param ColorBy A `colData` column or a gene name. Empty string = first
#'   cell-type column, else first QC metric.
#' @param UsePolygons Draw segmentation polygons instead of centroids
#'   (needs a `polygons` column, as created by SpaceTrooper).
#' @param PointSize Centroid size.
#' @param ... Other slots, e.g. `Selection`.
#' @return A `SpatialPlotPanel` object.
#' @examples
#' SpatialPlotPanel(ColorBy = "total")
#' @export
SpatialPlotPanel <- function(ColorBy = "", UsePolygons = FALSE, PointSize = 0.5, ...)
  new("SpatialPlotPanel", ColorBy = ColorBy, UsePolygons = UsePolygons, PointSize = PointSize, ...)

# ============================================================
# FOVZoomPanel - SpaceTrooper zoom on the selected FOVs (tile platforms).
# ============================================================

setMethod("initialize", "FOVZoomPanel", function(.Object, ...) {
  .Object <- callNextMethod(.Object, ...)
  if (length(.Object@ColorBy) == 0) .Object@ColorBy <- ""
  .Object
})

#' @export
setMethod(".fullName",   "FOVZoomPanel", function(x) "FOV Zoom")
#' @rdname FOVZoomPanel
#' @export
setMethod(".panelColor", "FOVZoomPanel", function(x) "#3A86FF")
#' @rdname FOVZoomPanel
#' @export
setMethod(".defineOutput", "FOVZoomPanel", function(x)
  shiny::plotOutput(.getEncodedName(x), height = "550px"))

#' @rdname FOVZoomPanel
#' @export
setMethod(".defineInterface", "FOVZoomPanel", function(x, se, select_info) {
  enc <- .getEncodedName(x)
  m <- .qcMetrics(se)
  list(iSEE::collapseBox(
    id = paste0(enc, "_DataBoxOpen"), title = "Data parameters", open = TRUE,
    shiny::selectInput(paste0(enc, "_ColorBy"), "Colour polygons by:", choices = m,
                       selected = if (x@ColorBy %in% m) x@ColorBy else m[1]),
    .selectionHelp(se)))
})

#' @rdname FOVZoomPanel
#' @export
setMethod(".generateOutput", "FOVZoomPanel", function(x, se, ...) {
  if (.selectionUnit(se) != "fov")
    return(.notAvailablePlot("FOV zoom needs FOV-based data (CosMx)"))
  fovs <- suppressWarnings(as.integer(.splitIDs(x@Selection)))
  fovs <- fovs[!is.na(fovs)]
  if (length(fovs) == 0) return(.notAvailablePlot("Select FOVs on the map, then press Validate."))
  cb <- if (nzchar(x@ColorBy)) x@ColorBy else .qcMetrics(se)[1]
  SpaceTrooper::plotZoomFovsMap(se, fovs = fovs, colourBy = cb)
})

#' @rdname FOVZoomPanel
#' @export
setMethod(".renderOutput", "FOVZoomPanel", function(x, se, output, pObjects, rObjects)
  .renderCachedSpatialPlot(x, se, output, pObjects, rObjects, c("ColorBy", "Selection")))

#' @rdname FOVZoomPanel
#' @export
setMethod(".createObservers", "FOVZoomPanel", function(x, se, input, session, pObjects, rObjects) {
  callNextMethod(); enc <- .getEncodedName(x)
  .slotObserver(enc, "ColorBy", input, pObjects, rObjects, cast = as.character)
})

#' Zoom on selected FOVs (CosMx)
#'
#' @param ColorBy QC metric used to colour cell polygons.
#' @param ... Other slots, e.g. `Selection`.
#' @return A `FOVZoomPanel` object.
#' @examples
#' FOVZoomPanel(ColorBy = "QC_score")
#' @export
FOVZoomPanel <- function(ColorBy = "", ...) new("FOVZoomPanel", ColorBy = ColorBy, ...)
