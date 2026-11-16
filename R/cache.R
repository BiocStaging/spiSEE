# Image cache of the plot panels, on disk by default (size-limited, least
# recently used removed first). Keys: dataset fingerprint, package version,
# panel name and the panel's key slots. Under devtools::load_all() a memory
# cache is used, since the code changes without a version bump.

.CACHE <- new.env(parent = emptyenv())

#' Location of the spiSEE image cache
#'
#' Plot images and map backgrounds are cached on disk by [spiSEE()]
#' (argument `cache`). `spatialCacheDir()` gives the default folder for the
#' images (a `plots` subfolder of the package's user cache, next to the map
#' backgrounds); `clearSpatialCache()` deletes the cached images in a folder
#' and, for the default folder, the map backgrounds too.
#'
#' @param dir Cache folder. Only cache files (`.rds`) in it are deleted.
#' @return `spatialCacheDir()`: the folder path. `clearSpatialCache()`:
#'   `dir`, invisibly.
#' @seealso [spiSEE()], whose `cache` argument chooses the folder.
#' @examples
#' spatialCacheDir()
#' @export
spatialCacheDir <- function() file.path(.userCacheDir(), "plots")

.userCacheDir <- function() tools::R_user_dir("spiSEE", which = "cache")

#' @rdname spatialCacheDir
#' @export
clearSpatialCache <- function(dir = spatialCacheDir()) {
    rds <- function(d) list.files(d, pattern = "[.]rds$", full.names = TRUE)
    if (dir.exists(dir)) unlink(rds(dir))
    if (identical(normalizePath(dir, mustWork = FALSE),
                    normalizePath(spatialCacheDir(), mustWork = FALSE)))
        unlink(grep("/map_[^/]*[.]rds$", rds(.userCacheDir()), value = TRUE))
    if (!is.null(.CACHE$maps)) rm(list = ls(.CACHE$maps), envir = .CACHE$maps)
    invisible(dir)
}

.isDevLoad <- function() exists(".__DEVTOOLS__",
                                envir = asNamespace("spiSEE"),
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

# Map backgrounds follow the plot cache: on disk next to it when it is on
# disk, else in memory for this R session (nothing written with
# cache = FALSE).
.mapStore <- function() {
    if (inherits(.CACHE$cache, "cachem")) return(.userCacheDir())
    if (is.null(.CACHE$maps)) .CACHE$maps <- new.env(parent = emptyenv())
    .CACHE$maps
}

# Identity of the data the plots and analyses are made from: names,
# coordinates, cell annotations and polygons, reduced dimensions and the
# full content of every assay (not a summary such as per-cell sums, which
# two different matrices can share). Computed once, at launch.
# metadata() is included (gene_ranks, fov_stats, technology... change the
# plots), except the fingerprint itself. Assays kept on disk (DelayedArray,
# e.g. HDF5) are identified by their file and its size and modification
# time.
.datasetFingerprint <- function(se) {
    cd <- colData(se)
    poly <- if ("polygons" %in% colnames(cd))
        sf::st_as_binary(sf::st_geometry(cd[["polygons"]]))
    cd <- as.data.frame(cd[, setdiff(colnames(cd), "polygons"), drop = FALSE])
    assays <- lapply(SummarizedExperiment::assays(se), .assayIdentity)
    md <- S4Vectors::metadata(se)
    md <- md[setdiff(names(md), "spiSEE_fingerprint")]
    digest::digest(list(dimnames(se), SpatialExperiment::spatialCoords(se),
        cd, poly, assays, as.list(SingleCellExperiment::reducedDims(se)),
        md, .selectionUnit(se)), algo = "xxhash64")
}

# A DelayedArray serialises as its seed (file path, dataset) and delayed
# operations, without reading the data: the file's size and modification
# time stand for its content.
.assayIdentity <- function(a) {
    if (!methods::is(a, "DelayedArray")) return(a)
    f <- tryCatch(unlist(BiocGenerics::path(a)), error = function(e) NULL)
    list(a, if (length(f)) file.info(f)[, c("size", "mtime")])
}

.fingerprintOf <- function(se) {
    fp <- S4Vectors::metadata(se)$spiSEE_fingerprint
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
