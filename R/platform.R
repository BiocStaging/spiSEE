# ============================================================
# Platform adapter.
#
# The only things that really differ between platforms are a handful of
# facts: how a selection is expressed (FOV tiles vs single cells), which
# colData columns hold QC metrics and cell-type labels, and which metric to
# show first. They live here, in one table, instead of being hard-coded in
# parallel CosMx*/Xenium* panel classes.
#
# Adding a platform = adding one entry to .PLATFORMS (and, if needed, a
# loader in load.R). No new panel class.
# ============================================================

.PLATFORMS <- list(
  cosmx = list(
    name = "CosMx", unit = "fov", region_col = "fov",
    qc_metrics = c("total", "detected", "Area_um", "QC_score", "QScore", "log2CountArea",
                   "Mean.DAPI", "Mean.GFAP", "dist_border"),
    label_cols = c("SingleR_label"),
    detect = c("fov", "Area_um")),
  xenium = list(
    name = "Xenium", unit = "cell", region_col = NA_character_,
    qc_metrics = c("transcript_counts", "total_counts", "cell_area",
                   "nucleus_area", "nucleus_count"),
    label_cols = c("Level1", "Level0", "Level2", "Leiden"),
    detect = c("transcript_counts", "cell_area")),
  visiumhd = list(
    name = "Visium HD", unit = "cell", region_col = NA_character_,
    qc_metrics = c("sum", "detected", "subset.proportion.mt", "log_sum",
                   "um2", "nuc_um2", "nuc"),
    label_cols = c("Level1", "Level0"),
    detect = c("um2", "nuc_um2")),
  generic = list(
    name = "Spatial", unit = "cell", region_col = NA_character_,
    qc_metrics = character(0), label_cols = character(0), detect = character(0))
)

#' Set or detect the spatial platform of a dataset
#'
#' iseeSpatial adapts its panels to the platform: CosMx data are selected
#' by field of view (FOV), Xenium and Visium HD data by individual cells.
#' The platform is stored in `metadata(spe)$iseeSpatial_platform`; when it
#' is absent it is guessed from characteristic `colData` columns.
#'
#' @param spe A `SpatialExperiment`.
#' @param platform One of `"cosmx"`, `"xenium"`, `"visiumhd"`, `"generic"`.
#' @return `setPlatform()` returns `spe` with the platform recorded;
#'   `getPlatform()` returns the platform key (character).
#' @examples
#' spe <- SpatialExperiment::SpatialExperiment(
#'     assays = list(counts = matrix(1, 2, 3)),
#'     spatialCoords = matrix(runif(6), 3, 2,
#'         dimnames = list(NULL, c("x", "y"))))
#' getPlatform(spe)
#' spe <- setPlatform(spe, "xenium")
#' getPlatform(spe)
#' @export
setPlatform <- function(spe, platform = c("cosmx", "xenium", "visiumhd", "generic")) {
  platform <- match.arg(platform)
  S4Vectors::metadata(spe)$iseeSpatial_platform <- platform
  spe
}

#' @rdname setPlatform
#' @export
getPlatform <- function(spe) {
  p <- S4Vectors::metadata(spe)$iseeSpatial_platform
  if (!is.null(p) && p %in% names(.PLATFORMS)) return(p)
  cn <- colnames(SummarizedExperiment::colData(spe))
  for (k in setdiff(names(.PLATFORMS), "generic")) {
    if (all(.PLATFORMS[[k]]$detect %in% cn)) return(k)
  }
  "generic"
}

.platformInfo <- function(se) .PLATFORMS[[getPlatform(se)]]

# "fov" (tiles; CosMx-style map) or "cell" (point map).
.selectionUnit <- function(se) {
  info <- .platformInfo(se)
  if (identical(info$unit, "fov") && info$region_col %in% colnames(colData(se))) "fov" else "cell"
}

# Numeric colData columns usable as QC metrics: platform ones first, then any
# other numeric column (so unknown platforms still get something sensible).
.qcMetrics <- function(se) {
  cd <- colData(se)
  num <- colnames(cd)[vapply(cd, function(v) is.numeric(v) && !is.factor(v), logical(1))]
  num <- setdiff(num, c("niche_entropy"))
  known <- intersect(.platformInfo(se)$qc_metrics, num)
  c(known, setdiff(num, c(known, grep("_outlier", num, value = TRUE))))
}

# Categorical colData columns usable as cell-type labels: factor/character
# with 2-60 distinct values. Platform label columns first.
.labelColumns <- function(se) {
  cd <- colData(se)
  ok <- vapply(colnames(cd), function(nm) {
    v <- cd[[nm]]
    (is.factor(v) || is.character(v) || is.logical(v)) &&
      length(unique(stats::na.omit(v))) %in% 2:60
  }, logical(1))
  cand <- setdiff(colnames(cd)[ok], c("sample_id", "polygons"))
  known <- intersect(c(.platformInfo(se)$label_cols, "cell_type", "celltype"), cand)
  c(known, setdiff(cand, known))
}

.defaultLabelCol <- function(se) {
  cand <- setdiff(.labelColumns(se), c("niche", "spatial_domain", "expr_cluster"))
  if (length(cand)) cand[1] else NA_character_
}

# Genes ordered for menus: precomputed ranking if present, else alphabetical.
.geneChoices <- function(se) {
  r <- S4Vectors::metadata(se)$gene_ranks
  if (!is.null(r)) c(intersect(r, rownames(se)), setdiff(sort(rownames(se)), r)) else sort(rownames(se))
}
