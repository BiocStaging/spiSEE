# Image cache of the plot panels, on disk by default (size-limited, least
# recently used removed first). Keys: dataset fingerprint, package version,
# panel name and the panel's key slots. Under devtools::load_all() a memory
# cache is used, since the code changes without a version bump.

.CACHE <- new.env(parent = emptyenv())

#' Location of the iseeSpatial image cache
#'
#' Plot images are cached on disk by [iseeSpatial()] (argument `cache`).
#' `spatialCacheDir()` gives the default folder for the images (a `plots`
#' subfolder of the package's user cache, next to the map backgrounds);
#' `clearSpatialCache()` deletes every cached image in a folder.
#'
#' @param dir Cache folder.
#' @return `spatialCacheDir()`: the folder path. `clearSpatialCache()`:
#'   `dir`, invisibly.
#' @examples
#' spatialCacheDir()
#' @export
spatialCacheDir <- function() file.path(tools::R_user_dir("iseeSpatial",
    which = "cache"), "plots")

#' @rdname spatialCacheDir
#' @export
clearSpatialCache <- function(dir = spatialCacheDir()) {
    if (dir.exists(dir)) unlink(list.files(dir, full.names = TRUE))
    invisible(dir)
}

.isDevLoad <- function() exists(".__DEVTOOLS__",
                                envir = asNamespace("iseeSpatial"),
                                inherits = FALSE)

# cache: TRUE (default folder), FALSE (memory only, lost when R stops) or a
# folder path. size: maximum size on disk, in bytes.
.configurePlotCache <- function(cache = TRUE, size = 1024^3) {
    if (isTRUE(cache) && .isDevLoad()) cache <- FALSE
    .CACHE$cache <- if (isFALSE(cache)) {
        # shiny's default: memory, shared by the sessions of this R process
        "app"
    } else {
        # Own folder: cachem's size limit counts and evicts every .rds file in
        # it, so it must not share one with the map backgrounds (map_*.rds).
        dir <- if (isTRUE(cache)) spatialCacheDir() else as.character(cache)
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
        cachem::cache_disk(dir, max_size = size, evict = "lru")
    }
    invisible(.CACHE$cache)
}

.plotCache <- function() if (is.null(.CACHE$cache)) "app" else .CACHE$cache

# Cheap identity of the data a plot is drawn from: cell and gene names,
# coordinates, cell annotations (without the polygon geometries) and the
# per-cell sums of the expression assay.
.datasetFingerprint <- function(se) {
    cd <- as.data.frame(colData(se)[,
                                    setdiff(colnames(colData(se)), "polygons"),
                                    drop = FALSE])
    a <- intersect(c("logcounts", "counts"),
                   SummarizedExperiment::assayNames(se))
    sums <- if (length(a)) Matrix::colSums(SummarizedExperiment::assay(se,
        a[1])) else NULL
    digest::digest(list(dimnames(se), SpatialExperiment::spatialCoords(se), cd,
                        sums,
        .selectionUnit(se)), algo = "xxhash64")
}

.fingerprintOf <- function(se) {
    fp <- S4Vectors::metadata(se)$iseeSpatial_fingerprint
    if (is.null(fp)) .datasetFingerprint(se) else fp
}

# Heavy analyses (statistics, gene ranking, niches) are identified by what
# they compute, not by when Compute was pressed: the same analysis on the
# same region gets the same key, so it is served from the cache, also after
# a restart, instead of being recomputed and stored again.
.COMPUTE_SLOTS <- list(
    SpatialStatsPanel = c("PlotType", "Genes", "LabelCol", "K", "Alpha",
                          "Selection"),
    SVGPanel          = c("TopN", "K", "Selection"),
    SpatialNichePanel = c("LabelCol", "K", "NNiches", "Lambda", "Selection"))

.computeKey <- function(x) {
    cl <- class(x)[1]
    vals <- lapply(.COMPUTE_SLOTS[[cl]], function(s) methods::slot(x, s))
    paste0("k", digest::digest(list(cl, vals), algo = "xxhash64"))
}

# A trigger set by hand (e.g. ComputeTrigger = "launch" in a layout) becomes
# the key of the panel's parameters, so two layouts never share a result.
.keyComputeTrigger <- function(x) {
    if (nzchar(x@ComputeTrigger)) x@ComputeTrigger <- .computeKey(x)
    x
}

# Computed results that are not plots (the gene ranking table): in the disk
# cache when there is one, else in memory for this R session.
.resultCache <- function() {
    if (inherits(.CACHE$cache, "cachem")) return(.CACHE$cache)
    if (is.null(.CACHE$mem))
        .CACHE$mem <- cachem::cache_mem(max_size = 256 * 1024^2)
    .CACHE$mem
}

.cachedResult <- function(key, expr) {
    cc <- .resultCache()
    hit <- cc$get(key)
    if (!cachem::is.key_missing(hit)) return(hit)
    val <- force(expr)
    cc$set(key, val)
    val
}
