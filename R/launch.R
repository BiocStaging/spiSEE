#' Launch the iseeSpatial app
#'
#' Opens an iSEE app whose panels adapt to the platform of `spe` (CosMx,
#' Xenium, Visium HD, or any `SpatialExperiment`). A selection made on the
#' overview map applies to every other panel.
#'
#' * `mode = "qc"`: map, spatial plot, two QC histograms (+ FOV zoom and
#'   per-FOV summary for CosMx). Meant for unfiltered data.
#' * `mode = "analysis"`: map, two spatial plots (annotation and gene),
#'   spatial statistics, spatially variable genes, niches and domains.
#'
#' @param spe A `SpatialExperiment` (e.g. from [load_cosmx()],
#'   [load_xenium()] or [load_visiumhd()]).
#' @param mode `"analysis"` or `"qc"`.
#' @param platform Optional platform key to override detection
#'   (see [setPlatform()]).
#' @param selection Initial selection (comma-separated FOV or cell ids).
#'   Empty = whole dataset.
#' @param initial Optional list of panels; replaces the default layout.
#' @param cache Where plot images are cached: `TRUE` (default) on disk in
#'   [spatialCacheDir()], kept across restarts; a folder path (used only for
#'   the images: do not share it with other files); or `FALSE`
#'   for memory only (lost when R stops). Under `devtools::load_all()` the
#'   cache is always in memory.
#' @param cache_size Maximum size of the disk cache in bytes (default 1 GB);
#'   the least recently used images are removed first.
#' @param ... Passed to [iSEE::iSEE()].
#' @return A Shiny app object (prints/launches in an interactive session).
#' @examples
#' spe <- SpatialExperiment::SpatialExperiment(
#'     assays = list(counts = matrix(rpois(200, 3), 10, 20,
#'         dimnames = list(paste0("g", 1:10), paste0("c", 1:20)))),
#'     spatialCoords = matrix(runif(40), 20, 2,
#'         dimnames = list(NULL, c("x", "y"))))
#' app <- iseeSpatial(spe, mode = "qc", cache = FALSE)
#' if (interactive()) shiny::runApp(app)
#' @export
iseeSpatial <- function(spe, mode = c("analysis", "qc"), platform = NULL,
                        selection = "", initial = NULL,
                        cache = getOption("iseeSpatial.cache", TRUE),
                        cache_size = getOption("iseeSpatial.cache_size",
                                                1024^3), ...) {
    stopifnot(methods::is(spe, "SpatialExperiment"))
    mode <- match.arg(mode)
    if (!is.null(platform)) spe <- setPlatform(spe, platform)
    if ("counts" %in% SummarizedExperiment::assayNames(spe))
        spe <- .addLogcounts(spe)
    .configurePlotCache(cache, cache_size)
    S4Vectors::metadata(spe)$iseeSpatial_gene_order <- .geneOrder(spe)
    S4Vectors::metadata(spe)$iseeSpatial_fingerprint <- .datasetFingerprint(spe)
    if (is.null(initial)) initial <- .defaultPanels(spe, mode, selection)
    args <- list(...)
    if (is.null(args$extra)) args$extra <- .extraPanels(spe)
    # quote(spe): iSEE() runs deparse(substitute(se)) for its code export; with
    # the object itself in the call, that deparsed the whole dataset at launch.
    do.call(iSEE::iSEE, c(list(quote(spe), initial = initial), args))
}

# Panel types offered in "Organize panels" in addition to the initial ones.
.extraPanels <- function(spe) {
    extra <- list(SpatialPlotPanel(), SpatialQCPanel(), SpatialStatsPanel(),
        SVGPanel(), SpatialNichePanel())
    if (.selectionUnit(spe) == "fov") extra <- c(extra,
        list(FOVZoomPanel(), FOVStatsPanel()))
    extra
}

.defaultPanels <- function(spe, mode, selection = "") {
    tiles <- .selectionUnit(spe) == "fov"
    qc <- .qcMetrics(spe)
    area <- intersect(c("Area_um", "cell_area", "um2"), qc)
    lab <- .defaultLabelCol(spe)
    if (mode == "qc") {
        panels <- list(SpatialMapPanel(CommittedIDs = selection),
            SpatialPlotPanel(ColorBy = if (length(qc)) qc[1] else "",
                            Selection = selection),
            SpatialQCPanel(Metric = if (length(qc)) qc[1] else "",
                            Selection = selection))
        if (length(area))
            panels <- append(panels, list(
                SpatialQCPanel(Metric = area[1], Selection = selection)))
        # SpaceTrooper's quality score: QScore (devel) or QC_score (older).
        score <- intersect(c("QScore", "QC_score"), qc)
        if (tiles) panels <- append(panels, list(FOVZoomPanel(
            ColorBy = if (length(score)) score[1] else "",
            Selection = selection)))
        if (tiles && !is.null(S4Vectors::metadata(spe)$fov_stats))
            panels <- append(panels, list(FOVStatsPanel(Selection = selection)))
        return(panels)
    }
    gene <- .geneChoices(spe)[1]
    list(SpatialMapPanel(CommittedIDs = selection),
        SpatialPlotPanel(ColorBy = if (!is.na(lab)) lab else "",
                        Selection = selection),
        SpatialPlotPanel(ColorBy = gene, Selection = selection),
        SpatialStatsPanel(Genes = gene, LabelCol = lab, Selection = selection),
        SVGPanel(LinkedPanel = "SpatialPlotPanel2", Selection = selection),
        SpatialNichePanel(LabelCol = lab, Selection = selection))
}

# Most expressed first, so the default gene of a panel shows something.
.geneOrder <- function(spe) {
    if (!"logcounts" %in% SummarizedExperiment::assayNames(spe) ||
        nrow(spe) == 0) return(NULL)
    m <- Matrix::rowMeans(SingleCellExperiment::logcounts(spe))
    rownames(spe)[order(m, decreasing = TRUE)]
}
