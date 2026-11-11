# ============================================================
# Image cache of the plot panels (.renderCachedSpatialPlot).
#
# By default images are kept on disk and survive app restarts (on the demo
# server an app stops after 10 minutes idle). The cache has a size limit;
# when it is full the least recently used images are removed first, so the
# images seen at every start and after every Clear (no selection) stay.
#
# Keys = dataset fingerprint + package version + panel name + the panel's
# key slots. The fingerprint keeps two datasets with the same panel settings
# apart; the version drops images drawn by older plot code. Under
# devtools::load_all() the code changes without a version bump, so a memory
# cache (cleared when R stops) is used instead.
# ============================================================

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
spatialCacheDir <- function() file.path(tools::R_user_dir("iseeSpatial", which = "cache"), "plots")

#' @rdname spatialCacheDir
#' @export
clearSpatialCache <- function(dir = spatialCacheDir()) {
  if (dir.exists(dir)) unlink(list.files(dir, full.names = TRUE))
  invisible(dir)
}

.isDevLoad <- function() exists(".__DEVTOOLS__", envir = asNamespace("iseeSpatial"), inherits = FALSE)

# cache: TRUE (default folder), FALSE (memory only, lost when R stops) or a
# folder path. size: maximum size on disk, in bytes.
.configurePlotCache <- function(cache = TRUE, size = 1024^3) {
  if (isTRUE(cache) && .isDevLoad()) cache <- FALSE
  .CACHE$cache <- if (isFALSE(cache)) {
    "app"   # shiny's default: memory, shared by the sessions of this R process
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
  cd <- as.data.frame(colData(se)[, setdiff(colnames(colData(se)), "polygons"), drop = FALSE])
  a <- intersect(c("logcounts", "counts"), SummarizedExperiment::assayNames(se))
  sums <- if (length(a)) Matrix::colSums(SummarizedExperiment::assay(se, a[1])) else NULL
  digest::digest(list(dimnames(se), SpatialExperiment::spatialCoords(se), cd, sums,
                      .selectionUnit(se)), algo = "xxhash64")
}

.fingerprintOf <- function(se) {
  fp <- S4Vectors::metadata(se)$iseeSpatial_fingerprint
  if (is.null(fp)) .datasetFingerprint(se) else fp
}
