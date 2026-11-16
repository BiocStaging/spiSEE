# SpatialQCPanel - histogram of a QC metric over the selection, with
# SpaceTrooper's medcouple outlier fences when `<metric>_outlier_mc` exists.
# FOVStatsPanel - per-FOV boxplots (tile platforms; metadata(se)$fov_stats).

setMethod("initialize", "SpatialQCPanel", function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    if (length(.Object@Metric)    == 0) .Object@Metric    <- ""
    if (length(.Object@UseFences) == 0) .Object@UseFences <- TRUE
    .Object
})

#' @rdname SpatialQCPanel
#' @export
setMethod(".fullName",   "SpatialQCPanel", function(x) "QC Histogram")
#' @rdname SpatialQCPanel
#' @export
setMethod(".panelColor", "SpatialQCPanel", function(x) "#E63946")
#' @rdname SpatialQCPanel
#' @export
setMethod(".defineOutput", "SpatialQCPanel", function(x)
    .withBusy(shiny::plotOutput(.getEncodedName(x), height = "400px")))

.resolveMetric <- function(x, se) {
    m <- .qcMetrics(se)
    if (nzchar(x@Metric) && x@Metric %in% m) x@Metric else m[1]
}

#' @rdname SpatialQCPanel
#' @export
setMethod(".defineInterface", "SpatialQCPanel", function(x, se, select_info) {
    enc <- .getEncodedName(x)
    list(iSEE::collapseBox(
        id = paste0(enc, "_DataBoxOpen"), title = "Data parameters",
        open = TRUE,
        shiny::selectInput(paste0(enc, "_Metric"), "QC metric:",
                            choices = .qcMetrics(se),
            selected = .resolveMetric(x, se)),
        shiny::checkboxInput(paste0(enc, "_UseFences"),
                            "Show outlier fences (medcouple), if computed",
            value = x@UseFences),
        .selectionHelp(se)))
})

#' @rdname SpatialQCPanel
#' @export
setMethod(".generateOutput", "SpatialQCPanel", function(x, se, ...) {
    metric <- .resolveMetric(x, se)
    if (is.na(metric))
        return(.notAvailablePlot("No numeric QC metric in colData"))
    fence_col <- paste0(metric, "_outlier_mc")
    fences <- if (x@UseFences && fence_col %in% colnames(colData(se)))
        fence_col
    se <- .subsetBySelection(se, x@Selection)
    if (ncol(se) == 0) return(.notAvailablePlot("No cells in this selection"))
    p <- SpaceTrooper::plotMetricHist(se, metric = metric, useFences = fences) +
        ggplot2::labs(title = .prettyMetric(metric), x = .prettyMetric(metric))
    if (is.null(fences)) return(p)
    # SpaceTrooper names the legend after the internal column
    # ("Fences transcript_counts_outlier_mc") and shows no counts.
    flagged <- .isFlagged(colData(se)[[fences]])
    p + ggplot2::labs(colour = "Outlier fences (medcouple)",
        subtitle = .flaggedText(flagged))
})

# "transcript_counts" -> "Transcript counts"; "nCount_RNA" stays as is.
.prettyMetric <- function(m) {
    if (!grepl("_", m)) return(m)
    s <- gsub("_", " ", m)
    paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
}

# SpaceTrooper stores medcouple outliers as "HIGH"/"LOW"/"NO" (class
# outlier.filter); scuttle::isOutlier() gives a logical.
.isFlagged <- function(x) {
    v <- unclass(x)
    attributes(v) <- NULL
    if (is.character(v) || is.factor(x)) return(as.character(v) != "NO")
    as.logical(v)
}

.flaggedText <- function(flagged) {
    n <- sum(flagged, na.rm = TRUE)
    tot <- length(flagged)
    sprintf("%s of %s cells flagged as outliers (%.1f%%)",
        format(n, big.mark = ","), format(tot, big.mark = ","),
        if (tot) 100 * n / tot else 0)
}

#' @rdname SpatialQCPanel
#' @export
setMethod(".renderOutput", "SpatialQCPanel",
            function(x, se, output, pObjects, rObjects)
    .renderCachedSpatialPlot(x, se, output, pObjects, rObjects,
                            c("Metric", "UseFences", "Selection")))

#' @rdname SpatialQCPanel
#' @export
setMethod(".createObservers", "SpatialQCPanel",
            function(x, se, input, session, pObjects, rObjects) {
    callNextMethod()
    enc <- .getEncodedName(x)
    .slotObserver(enc, "Metric", input, pObjects, rObjects, cast = as.character)
    .slotObserver(enc, "UseFences", input, pObjects, rObjects,
                    cast = as.logical)
})

#' QC metric histogram for the selected region
#'
#' iSEE panel with the histogram of one QC metric (a numeric `colData`
#' column) over the cells of the region selected on the
#' [SpatialMapPanel()]. When SpaceTrooper's medcouple outlier detection has
#' been run (`SpaceTrooper::computeSpatialOutlier()`, which adds a column
#' `<Metric>_outlier_mc`), the outlier fences are drawn and the number of
#' flagged cells is given.
#'
#' @param Metric A numeric `colData` column. Empty = first known QC metric
#'   of the platform.
#' @param UseFences Show medcouple outlier fences when a column
#'   `<Metric>_outlier_mc` exists (SpaceTrooper).
#' @param ... Other slots, e.g. `Selection`.
#' @return A `SpatialQCPanel` object.
#' @seealso [FOVStatsPanel()] for metrics per FOV, [SpatialPlotPanel()] to
#'   see a metric on the tissue.
#' @examples
#' SpatialQCPanel(Metric = "total")
#' @export
SpatialQCPanel <- function(Metric = "", UseFences = TRUE, ...)
    new("SpatialQCPanel", Metric = Metric, UseFences = UseFences, ...)

# FOVStatsPanel

#' @rdname FOVStatsPanel
#' @export
setMethod(".fullName",   "FOVStatsPanel", function(x) "Per-FOV QC Summary")
#' @rdname FOVStatsPanel
#' @export
setMethod(".panelColor", "FOVStatsPanel", function(x) "#8338EC")
#' @rdname FOVStatsPanel
#' @export
setMethod(".defineOutput", "FOVStatsPanel", function(x)
    .withBusy(shiny::plotOutput(.getEncodedName(x), height = "600px")))
#' @rdname FOVStatsPanel
#' @export
setMethod(".defineInterface", "FOVStatsPanel", function(x, se, select_info) {
    enc <- .getEncodedName(x)
    list(iSEE::collapseBox(id = paste0(enc, "_DataBoxOpen"), title = "Display",
                            open = TRUE,
        shiny::helpText("Boxplots of QC metrics per FOV",
                        "(metadata(se)$fov_stats, else the QC columns).")))
})

# Long table fov / name / value: metadata(se)$fov_stats if given, else the
# first six QC metrics of colData.
.fovLong <- function(se) {
    st <- S4Vectors::metadata(se)$fov_stats
    if (!is.null(st)) return(as.data.frame(st)[, c("fov", "name", "value")])
    m <- utils::head(.qcMetrics(se), 6)
    if (!length(m)) return(NULL)
    cd <- colData(se)
    fov <- as.character(cd[[.platformInfo(se)$region_col]])
    do.call(rbind, lapply(m, function(v) data.frame(fov = fov, name = v,
        value = as.numeric(cd[[v]]))))
}

#' @rdname FOVStatsPanel
#' @export
setMethod(".generateOutput", "FOVStatsPanel", function(x, se, ...) {
    if (.selectionUnit(se) != "fov")
        return(.notAvailablePlot("Per-FOV summaries need FOV-based data"))
    fov_long <- .fovLong(se)
    if (is.null(fov_long))
        return(.notAvailablePlot("No QC metric to summarise per FOV"))
    sel <- .splitIDs(x@Selection)
    if (length(sel)) fov_long <- fov_long[as.character(fov_long$fov) %in% sel,
                                            , drop = FALSE]
    if (!nrow(fov_long))
        return(.notAvailablePlot("No per-FOV values for this selection"))
    u <- unique(as.character(fov_long$fov))
    all_ids <- .allRegionIDs(se)
    fov_long$fov <- factor(fov_long$fov,
        levels = c(intersect(all_ids, u), setdiff(u, all_ids)))
    n <- nlevels(fov_long$fov)
    brk <- levels(fov_long$fov)
    if (n > 15) brk <- brk[unique(c(1, seq(10, n, by = 10)))]
    ggplot2::ggplot(fov_long, ggplot2::aes(fov, value, fill = fov)) +
        ggplot2::facet_wrap(~name, nrow = 3, scales = "free_y") +
        ggplot2::geom_boxplot(linewidth = 0.2, outlier.size = 0.2,
                                show.legend = FALSE) +
        ggplot2::scale_x_discrete(breaks = brk) + ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none",
                        panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::labs(x = "FOV", y = NULL,
                        title = "QC metrics per field of view")
})

#' @rdname FOVStatsPanel
#' @export
setMethod(".renderOutput", "FOVStatsPanel",
            function(x, se, output, pObjects, rObjects)
    .renderCachedSpatialPlot(x, se, output, pObjects, rObjects, "Selection"))

#' @rdname FOVStatsPanel
#' @export
setMethod(".createObservers", "FOVStatsPanel",
    function(x, se, input, session, pObjects, rObjects) callNextMethod())

#' Per-FOV QC summary (CosMx)
#'
#' iSEE panel with one boxplot per field of view (FOV) for each of several
#' QC metrics, on the FOVs selected on the [SpatialMapPanel()] (all FOVs
#' when nothing is selected). Useful to spot FOVs with weak staining or few
#' transcripts before filtering.
#'
#' By default the first six QC metrics of `colData` are shown, on their own
#' scale. A table in `metadata(spe)$fov_stats` replaces them, for instance
#' to show transformed values: a `data.frame` in long format with one row
#' per cell and metric and the columns `fov` (FOV id, as in `spe$fov`),
#' `name` (metric name, one facet each) and `value` (number).
#'
#' @param ... Slots, e.g. `Selection`.
#' @seealso [SpatialQCPanel()] for the distribution of one metric,
#'   [FOVZoomPanel()] for the cells of the selected FOVs.
#' @return A `FOVStatsPanel` object.
#' @examples
#' FOVStatsPanel()
#'
#' # optional: transformed values, one row per cell and metric
#' spe <- SpatialExperiment::SpatialExperiment(
#'     assays = list(counts = matrix(1, 2, 6)),
#'     colData = S4Vectors::DataFrame(fov = rep(1:2, 3), total = 1:6),
#'     spatialCoords = matrix(runif(12), 6, 2,
#'         dimnames = list(NULL, c("x", "y"))))
#' S4Vectors::metadata(spe)$fov_stats <- data.frame(fov = spe$fov,
#'     name = "log10 total", value = log10(spe$total))
#' @export
FOVStatsPanel <- function(...) new("FOVStatsPanel", ...)
