#' Launch the spiSEE app
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
#' @param spe A `SpatialExperiment` (e.g. from [loadCosmx()],
#'   [loadXenium()] or [loadVisiumHD()]).
#' @param mode `"analysis"` or `"qc"`.
#' @param platform Optional platform key to override detection
#'   (see [setPlatform()]).
#' @param selection Initial selection: FOV or cell ids, as a character
#'   vector or one comma-separated string. Empty = whole dataset. Unknown
#'   ids are dropped with a warning.
#' @param sample With several samples (`spe$sample_id`), the one to show:
#'   the app shows one sample at a time (one map, one coordinate system).
#'   Open one app per sample to compare them.
#' @param title Title of the app (browser tab and header).
#' @param initial Optional list of panels; replaces the default layout.
#' @param cache Where plot images are cached: `TRUE` (default) on disk in
#'   [spatialCacheDir()], kept across restarts; a folder path (used only for
#'   the images: do not share it with other files); or `FALSE`
#'   for memory only (lost when R stops).
#' @param cache_size Maximum size of the disk cache in bytes (default 1 GB);
#'   the least recently used images are removed first.
#' @param ... Passed to [iSEE::iSEE()].
#' @return A Shiny app object (prints/launches in an interactive session).
#' @section Optional metadata:
#' Some entries of `metadata(spe)` are used when present:
#' * `gene_ranks`: gene names, in the order of the gene menus and of the
#'   genes scanned by [SVGPanel()] (default: by mean expression).
#' * `fov_stats`: per-FOV values for [FOVStatsPanel()] (see there).
#' * `technology`: as set by SpaceTrooper's readers; gives the unit of the
#'   scale bars of the spatial plots (none without it).
#' * `fov_positions`, `fov_dim`: FOV layout from
#'   `SpaceTrooper::readCosmxSPE()`, drawn on the map (CosMx); without them
#'   the tiles follow the extent of the cells of each FOV.
#' * Niches and domains stored by [computeNiches()] and
#'   [computeSpatialDomains()]: shown by [SpatialNichePanel()] without
#'   recomputing.
#' @seealso [SpatialMapPanel()] and the other panels; [spatialCacheDir()].
#' @examples
#' spe <- SpatialExperiment::SpatialExperiment(
#'     assays = list(counts = matrix(rpois(200, 3), 10, 20,
#'         dimnames = list(paste0("g", 1:10), paste0("c", 1:20)))),
#'     spatialCoords = matrix(runif(40), 20, 2,
#'         dimnames = list(NULL, c("x", "y"))))
#' app <- spiSEE(spe, mode = "qc", cache = FALSE)
#' if (interactive()) shiny::runApp(app)
#' @export
spiSEE <- function(spe, mode = c("analysis", "qc"), platform = NULL,
                        selection = "", sample = NULL, initial = NULL,
                        cache = getOption("spiSEE.cache", TRUE),
                        cache_size = getOption("spiSEE.cache_size",
                                                1024^3),
                        title = "spiSEE", ...) {
    stopifnot(methods::is(spe, "SpatialExperiment"))
    mode <- match.arg(mode)
    spe <- .oneSample(spe, sample)
    if (anyDuplicated(colnames(spe)) || is.null(colnames(spe))) {
        message("Cell names are missing or not unique: made unique.")
        colnames(spe) <- make.unique(if (is.null(colnames(spe)))
            as.character(seq_len(ncol(spe))) else colnames(spe))
    }
    if (!is.null(platform)) spe <- setPlatform(spe, platform)
    selection <- .checkSelection(selection, spe)
    if ("counts" %in% SummarizedExperiment::assayNames(spe))
        spe <- .addLogcounts(spe)
    .configurePlotCache(cache, cache_size)
    S4Vectors::metadata(spe)$spiSEE_gene_order <- .geneOrder(spe)
    S4Vectors::metadata(spe)$spiSEE_fingerprint <- .datasetFingerprint(spe)
    if (is.null(initial)) initial <- .defaultPanels(spe, mode, selection)
    args <- list(...)
    if (is.null(args$extra)) args$extra <- .extraPanels(spe)
    if (is.null(args$appTitle)) args$appTitle <- title
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
        if (tiles)
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

# One sample per app: several samples share FOV numbers and coordinates.
.oneSample <- function(spe, sample) {
    ids <- unique(as.character(spe$sample_id))
    if (!is.null(sample)) {
        if (!sample %in% ids)
            stop("'sample' must be one of: ", paste(ids, collapse = ", "))
        return(spe[, spe$sample_id == sample])
    }
    if (length(ids) > 1)
        stop("'spe' holds ", length(ids), " samples (",
            paste(utils::head(ids, 5), collapse = ", "),
            if (length(ids) > 5) ", ...", "): choose one with 'sample ='.",
            call. = FALSE)
    spe
}

# Selection given at launch: trimmed, checked against the dataset's ids.
.checkSelection <- function(selection, spe) {
    ids <- trimws(unlist(strsplit(as.character(selection), ",")))
    ids <- unique(ids[nzchar(ids)])
    bad <- setdiff(ids, .allRegionIDs(spe))
    if (length(bad))
        warning("Unknown ids dropped from 'selection': ",
            paste(utils::head(bad, 5), collapse = ", "),
            if (length(bad) > 5) ", ...", call. = FALSE)
    paste(setdiff(ids, bad), collapse = ",")
}
