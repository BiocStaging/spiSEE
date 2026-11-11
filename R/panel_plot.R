# ============================================================
# SpatialPlotPanel - cells on the tissue, coloured by an annotation, a QC
# metric or a gene. Replaces CosMxSpatialPanel, CosMxGenePanel and
# XeniumSpatialPanel (one class, any platform).
# ============================================================

# Cells with 'polygons' above this count are drawn as centroids unless the
# user ticks "Draw cell polygons" by hand.
.POLYGON_AUTO_MAX <- 20000L

# SpaceTrooper::plotPolygons() runs data.frame(colData(se)), which splits the
# sf 'polygons' column into polygons.* columns and then fails to find
# 'polygons'. Draw the geometry directly instead.
.polygonPlot <- function(se, colorby, palette = NULL) {
    value <- colData(se)[[colorby]]
    df <- sf::st_sf(value = value, geometry = sf::st_geometry(colData(se)[["polygons"]]))
    p <- ggplot2::ggplot(df, ggplot2::aes(fill = value)) +
        ggplot2::geom_sf(colour = NA) +
        ggplot2::theme_minimal() + ggplot2::labs(fill = colorby)
    if (is.numeric(value))
        return(p + ggplot2::scale_fill_viridis_c(option = "D", na.value = "grey85"))
    lv <- levels(factor(value))
    pal <- if (!is.null(palette) && all(lv %in% names(palette))) palette[lv]
    else stats::setNames(.discretePalette(length(lv)), lv)
    p + ggplot2::scale_fill_manual(values = pal, na.value = "grey85")
}

setMethod("initialize", "SpatialPlotPanel", function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    if (length(.Object@ColorBy)     == 0) .Object@ColorBy     <- ""
    # NA = not yet resolved against a dataset; .refineParameters() below
    # turns this into TRUE/FALSE once the data (and whether it has a
    # 'polygons' column) is known, so polygons default on wherever that
    # data exists instead of needing to be switched on by hand every time.
    # A user's own explicit tick/untick is never NA, so it's never touched
    # again once set.
    if (length(.Object@UsePolygons) == 0) .Object@UsePolygons <- NA
    if (length(.Object@PointSize)   == 0) .Object@PointSize   <- 0.5
    .Object
})

#' @rdname SpatialPlotPanel
#' @export
setMethod(".refineParameters", "SpatialPlotPanel", function(x, se) {
    x <- callNextMethod()
    # Default on only when polygons exist and the dataset is small enough to
    # draw them: thousands of polygons per panel freeze the browser.
    if (!is.null(x) && isTRUE(is.na(x@UsePolygons)))
        x@UsePolygons <- "polygons" %in% colnames(colData(se)) && ncol(se) <= .POLYGON_AUTO_MAX
    x
})

#' @rdname SpatialPlotPanel
#' @export
setMethod(".fullName",   "SpatialPlotPanel", function(x) "Spatial Plot")
#' @rdname SpatialPlotPanel
#' @export
setMethod(".panelColor", "SpatialPlotPanel", function(x) "#EF476F")
#' @rdname SpatialPlotPanel
#' @export
setMethod(".defineOutput", "SpatialPlotPanel", function(x)
    .withBusy(shiny::plotOutput(.getEncodedName(x), height = "600px")))

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
            value = isTRUE(x@UsePolygons)),
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
    label <- colorby
    pal_col <- paste0(label, "_colors")
    if (!colorby %in% colnames(colData(se)) && colorby %in% rownames(se))
        colData(se)[[colorby]] <- as.numeric(SingleCellExperiment::logcounts(se)[colorby, ])
    # SpaceTrooper rebuilds a data.frame from colData, which renames "MT-ATP6"
    # to "MT.ATP6" and breaks the lookup: plot a copy under a safe name.
    if (!identical(colorby, make.names(colorby))) {
        colData(se)[[".colour_by"]] <- colData(se)[[colorby]]
        colorby <- ".colour_by"
    }

    args <- list(se, colourBy = colorby)
    # Stored label colours (CosMx: SingleR_label_colors) are reused when present.
    if (pal_col %in% colnames(colData(se))) {
        map <- unique(data.frame(l = as.character(colData(se)[[colorby]]),
            c = as.character(colData(se)[[pal_col]])))
        args$palette <- stats::setNames(map$c, map$l)
    }
    if (.selectionUnit(se) == "fov") args["sampleId"] <- list(NULL)   # as the CosMx panels always did

    if (isTRUE(x@UsePolygons)) {
        # Without a "polygons" column there is nothing to draw: say so rather
        # than hang the panel silently.
        if (!"polygons" %in% colnames(colData(se)))
            return(.notAvailablePlot("No cell polygons loaded for this dataset: untick 'Draw cell polygons'"))
        .polygonPlot(se, colorby, args$palette) + ggplot2::labs(fill = label)
    } else {
        p <- do.call(SpaceTrooper::plotCentroids,
            c(args, size = x@PointSize, .scaleBarArgs(SpaceTrooper::plotCentroids, se)))
        if (!identical(label, colorby)) p <- p + ggplot2::labs(colour = label, fill = label)
        p
    }
})

# A server-side selectize reports "" while its choices are still loading.
# Storing it would change the cache key, so the plot was drawn once with ""
# (resolved to the default column, not the panel's gene) and again when the
# real value arrived. "" is never a menu choice: ignore it.
.castColorBy <- function(v) {
    v <- as.character(v)
    if (identical(v, "")) NA_character_ else v
}

#' @rdname SpatialPlotPanel
#' @export
setMethod(".renderOutput", "SpatialPlotPanel", function(x, se, output, pObjects, rObjects) {
    .renderCachedSpatialPlot(x, se, output, pObjects, rObjects,
        c("ColorBy", "UsePolygons", "PointSize", "Selection"))
})

#' @rdname SpatialPlotPanel
#' @export
setMethod(".createObservers", "SpatialPlotPanel", function(x, se, input, session, pObjects, rObjects) {
    callNextMethod()
    enc <- .getEncodedName(x)
    .serverSelectize(session, enc, "ColorBy", .colorByChoices(se), .resolveColorBy(x, se))
    .slotObserver(enc, "ColorBy", input, pObjects, rObjects, cast = .castColorBy)
    .slotObserver(enc, "UsePolygons", input, pObjects, rObjects, cast = as.logical)
    .slotObserver(enc, "PointSize", input, pObjects, rObjects, cast = as.numeric)
})

#' Spatial plot of cells coloured by annotation, QC metric or gene
#'
#' @param ColorBy A `colData` column or a gene name. Empty string = first
#'   cell-type column, else first QC metric.
#' @param UsePolygons Draw segmentation polygons instead of centroids (needs
#'   a `polygons` column, as created by SpaceTrooper). `NA` (the default)
#'   follows the dataset: on at load time if a `polygons` column exists,
#'   off otherwise.
#' @param PointSize Centroid size.
#' @param ... Other slots, e.g. `Selection`.
#' @return A `SpatialPlotPanel` object.
#' @examples
#' SpatialPlotPanel(ColorBy = "total")
#' @export
SpatialPlotPanel <- function(ColorBy = "", UsePolygons = NA, PointSize = 0.5, ...)
    new("SpatialPlotPanel", ColorBy = ColorBy, UsePolygons = UsePolygons, PointSize = PointSize, ...)

# ============================================================
# FOVZoomPanel - SpaceTrooper zoom on the selected FOVs (tile platforms).
# ============================================================

setMethod("initialize", "FOVZoomPanel", function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    if (length(.Object@ColorBy) == 0) .Object@ColorBy <- ""
    .Object
})

#' @rdname FOVZoomPanel
#' @export
setMethod(".fullName",   "FOVZoomPanel", function(x) "FOV Zoom")
#' @rdname FOVZoomPanel
#' @export
setMethod(".panelColor", "FOVZoomPanel", function(x) "#3A86FF")
#' @rdname FOVZoomPanel
#' @export
setMethod(".defineOutput", "FOVZoomPanel", function(x)
    .withBusy(shiny::plotOutput(.getEncodedName(x), height = "550px")))

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
    ids <- .splitIDs(x@Selection)
    fovs <- as.integer(ids[grepl("^[0-9]+$", ids)])
    if (length(fovs) == 0) return(.notAvailablePlot("Select FOVs on the map, then press Validate."))
    if (!"polygons" %in% colnames(colData(se)))
        return(.notAvailablePlot("FOV zoom draws cell polygons; this dataset has none."))
    cb <- if (nzchar(x@ColorBy)) x@ColorBy else .qcMetrics(se)[1]
    do.call(SpaceTrooper::plotZoomFovsMap,
        c(list(se, fovs = fovs, colourBy = cb),
            .scaleBarArgs(SpaceTrooper::plotZoomFovsMap, se, "scaleBars")))
})

#' @rdname FOVZoomPanel
#' @export
setMethod(".renderOutput", "FOVZoomPanel", function(x, se, output, pObjects, rObjects)
    .renderCachedSpatialPlot(x, se, output, pObjects, rObjects, c("ColorBy", "Selection")))

#' @rdname FOVZoomPanel
#' @export
setMethod(".createObservers", "FOVZoomPanel", function(x, se, input, session, pObjects, rObjects) {
    callNextMethod()
    enc <- .getEncodedName(x)
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
