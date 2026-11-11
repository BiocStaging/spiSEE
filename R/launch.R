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
#'   [spatialCacheDir()], kept across restarts; a folder path; or `FALSE`
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
#'     spatialCoords = matrix(runif(40), 20, 2, dimnames = list(NULL, c("x", "y"))))
#' app <- iseeSpatial(spe, mode = "qc")
#' if (interactive()) shiny::runApp(app)
#' @export
iseeSpatial <- function(spe, mode = c("analysis", "qc"), platform = NULL,
                        selection = "", initial = NULL,
                        cache = getOption("iseeSpatial.cache", TRUE),
                        cache_size = getOption("iseeSpatial.cache_size", 1024^3), ...) {
  stopifnot(methods::is(spe, "SpatialExperiment"))
  mode <- match.arg(mode)
  if (!is.null(platform)) spe <- setPlatform(spe, platform)
  if (!"logcounts" %in% SummarizedExperiment::assayNames(spe) &&
      "counts" %in% SummarizedExperiment::assayNames(spe))
    SingleCellExperiment::logcounts(spe) <- log1p(SingleCellExperiment::counts(spe))
  .configurePlotCache(cache, cache_size)
  S4Vectors::metadata(spe)$iseeSpatial_fingerprint <- .datasetFingerprint(spe)
  if (is.null(initial)) initial <- .defaultPanels(spe, mode, selection)
  args <- list(...)
  if (is.null(args$extra)) args$extra <- .extraPanels(spe)
  do.call(iSEE::iSEE, c(list(spe, initial = initial), args))
}

# Panel types offered in "Organize panels" in addition to the initial ones.
.extraPanels <- function(spe) {
  extra <- list(SpatialPlotPanel(), SpatialQCPanel(), SpatialStatsPanel(),
                SVGPanel(), SpatialNichePanel())
  if (.selectionUnit(spe) == "fov") extra <- c(extra, list(FOVZoomPanel(), FOVStatsPanel()))
  extra
}

.defaultPanels <- function(spe, mode, selection = "") {
  tiles <- .selectionUnit(spe) == "fov"
  qc <- .qcMetrics(spe)
  area <- intersect(c("Area_um", "cell_area", "um2"), qc)
  lab <- .defaultLabelCol(spe)
  if (mode == "qc") {
    panels <- list(SpatialMapPanel(),
                   SpatialPlotPanel(ColorBy = if (length(qc)) qc[1] else "", Selection = selection),
                   SpatialQCPanel(Metric = if (length(qc)) qc[1] else "", Selection = selection))
    if (length(area)) panels <- append(panels, list(SpatialQCPanel(Metric = area[1], Selection = selection)))
    if (tiles) panels <- append(panels, list(FOVZoomPanel(
      ColorBy = if ("QC_score" %in% qc) "QC_score" else "", Selection = selection)))
    if (tiles && !is.null(S4Vectors::metadata(spe)$fov_stats))
      panels <- append(panels, list(FOVStatsPanel(Selection = selection)))
    return(panels)
  }
  gene <- .geneChoices(spe)[1]
  list(SpatialMapPanel(),
       SpatialPlotPanel(ColorBy = if (!is.na(lab)) lab else "", Selection = selection),
       SpatialPlotPanel(ColorBy = gene, Selection = selection),
       SpatialStatsPanel(Genes = gene, LabelCol = lab, Selection = selection),
       SVGPanel(LinkedPanel = "SpatialPlotPanel2", Selection = selection),
       SpatialNichePanel(LabelCol = lab, Selection = selection))
}
