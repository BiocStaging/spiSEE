#' Static plot of cells on the tissue
#'
#' Draws the cells of a `SpatialExperiment` as their polygons (or as points
#' if there are none), coloured by annotations, QC metrics or genes, with
#' the colours and legends of the app's spatial plots. Several values give
#' one panel each, on a shared colour scale.
#'
#' @param spe A `SpatialExperiment`.
#' @param colour_by One or more `colData` columns or genes (row names);
#'   genes are coloured by `logcounts` (or `counts`). All numeric or all
#'   categorical.
#' @param polygons Draw the cell polygons of a `polygons` column (as created
#'   by SpaceTrooper) if present; `FALSE` draws points.
#' @param palette Named colours for categorical values (default: distinct
#'   colours, see the app).
#' @param ncol Number of columns of panels.
#' @param borders Draw cell borders (polygons only).
#' @param point_size Point size when drawing points.
#' @return A `ggplot`.
#' @examples
#' dir <- system.file("extdata", "CosMx_DBKero_Tiny", package = "SpaceTrooper")
#' if (nzchar(dir)) {
#'     spe <- SpaceTrooper::readCosmxSPE(dir, sampleName = "DBKero")
#'     spe <- SpaceTrooper::spatialPerCellQC(spe)
#'     spe <- SpaceTrooper::readAndAddPolygonsToSPE(spe)
#'     plotSpatialCells(spe, "Area_um")
#' }
#' @export
plotSpatialCells <- function(spe, colour_by, polygons = TRUE, palette = NULL,
                            ncol = NULL, borders = FALSE, point_size = 0.3) {
    stopifnot(methods::is(spe, "SpatialExperiment"), length(colour_by) >= 1)
    vals <- lapply(colour_by, function(v) .cellValues(spe, v))
    num <- vapply(vals, is.numeric, logical(1))
    if (length(unique(num)) > 1)
        stop("'colour_by' mixes numeric and categorical values")
    numeric <- num[1]
    if (!numeric) {
        lv <- unique(unlist(lapply(vals, function(v) levels(factor(v)))))
        vals <- lapply(vals, function(v) factor(as.character(v), levels = lv))
    }
    panel <- factor(rep(colour_by, each = ncol(spe)), levels = colour_by)
    value <- if (numeric) unlist(vals, use.names = FALSE)
    else factor(unlist(lapply(vals, as.character), use.names = FALSE),
                levels = levels(vals[[1]]))
    use_poly <- polygons && "polygons" %in% colnames(colData(spe))
    p <- .cellsBase(spe, panel, value, use_poly, borders, point_size)
    aes_name <- if (use_poly) "fill" else "colour"
    lab <- list(if (length(colour_by) == 1) colour_by else NULL)
    names(lab) <- aes_name
    p <- p + ggplot2::theme_void() + do.call(ggplot2::labs, lab)
    if (length(colour_by) > 1)
        p <- p + ggplot2::facet_wrap(~panel, ncol = ncol)
    if (numeric) {
        sc <- if (aes_name == "fill") ggplot2::scale_fill_viridis_c
        else ggplot2::scale_colour_viridis_c
        return(p + sc(na.value = "grey85"))
    }
    lv <- levels(value)
    pal <- if (!is.null(palette) && all(lv %in% names(palette))) palette[lv]
    else stats::setNames(.discretePalette(length(lv)), lv)
    sc <- if (aes_name == "fill") ggplot2::scale_fill_manual
    else ggplot2::scale_colour_manual
    .compactLegend(p + sc(values = pal, na.value = "grey85"), value)
}

# Values of one colData column, or the expression of one gene.
.cellValues <- function(spe, v) {
    if (v %in% colnames(colData(spe))) {
        x <- colData(spe)[[v]]
        if (is.logical(x) || is.character(x)) x <- factor(x)
        return(x)
    }
    if (v %in% rownames(spe)) {
        a <- intersect(c("logcounts", "counts"),
                        SummarizedExperiment::assayNames(spe))
        if (!length(a)) stop("no logcounts or counts assay for gene ", v)
        return(as.numeric(SummarizedExperiment::assay(spe, a[1])[v, ]))
    }
    stop("'", v, "' is neither a colData column nor a gene")
}

# The cells as polygons (fill) or points (colour), one copy per panel.
.cellsBase <- function(spe, panel, value, use_poly, borders, point_size) {
    n <- nlevels(panel)
    if (use_poly) {
        geo <- sf::st_geometry(colData(spe)[["polygons"]])
        df <- sf::st_sf(panel = panel, value = value, geometry = rep(geo, n))
        return(ggplot2::ggplot(df) +
            ggplot2::geom_sf(ggplot2::aes(fill = value),
                            colour = if (borders) "grey20" else NA,
                            linewidth = 0.05))
    }
    xy <- SpatialExperiment::spatialCoords(spe)
    df <- data.frame(x = rep(xy[, 1], n), y = rep(xy[, 2], n),
                    panel = panel, value = value)
    ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = value)) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::coord_fixed()
}
