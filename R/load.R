# Loaders. Each returns a SpatialExperiment with a `logcounts` assay and
# the platform recorded (setPlatform), ready for iseeSpatial().
# Reader packages are in Suggests: only the loader you use must be installed.

.needs <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE))
        stop("Package '", pkg,
             "' is needed for this loader: BiocManager::install('", pkg, "')",
             call. = FALSE)
}

#' Load a prepared CosMx dataset
#'
#' Reads a `SpatialExperiment` saved as RDS (for instance by the CosMx
#' preparation script, `data-raw/DATASET.R`) and records the platform.
#'
#' @param rds_path Path to the `.RDS` file.
#' @return A `SpatialExperiment`.
#' @examples
#' f <- "cos_isee_ready.RDS"   # a CosMx SpatialExperiment saved with saveRDS()
#' if (file.exists(f)) spe <- load_cosmx(f)
#' @export
load_cosmx <- function(rds_path) {
    spe <- readRDS(rds_path)
    spe <- setPlatform(spe, "cosmx")
    message(sprintf("CosMx: %d cells, %d genes, %d FOVs", ncol(spe), nrow(spe),
        length(unique(colData(spe)$fov))))
    spe
}

#' Load a 10x Xenium output directory
#'
#' Reads the bundle with `SpatialExperimentIO::readXeniumSXE()`, adds
#' `logcounts` (log1p of counts) and, when `cell_boundaries.parquet` is
#' present, cell polygons via SpaceTrooper. Missing polygons are not an
#' error: the panels work on centroids.
#'
#' @param xenium_dir Path to an unzipped Xenium output directory.
#' @return A `SpatialExperiment`.
#' @examples
#' d <- "xenium_output"        # a Xenium output folder (outs/)
#' if (dir.exists(d)) spe <- load_xenium(d)
#' @export
load_xenium <- function(xenium_dir) {
    .needs("SpatialExperimentIO")
    # Bundles shared without transcripts.parquet (the largest file) are
    # common; readXeniumSXE() fails on them unless parquet paths are off.
    tx <- file.exists(file.path(xenium_dir, "transcripts.parquet"))
    spe <- SpatialExperimentIO::readXeniumSXE(xenium_dir,
                                              addParquetPaths = tx)
    SingleCellExperiment::logcounts(spe) <- log1p(SingleCellExperiment::counts(spe))
    bnd <- file.path(xenium_dir, "cell_boundaries.parquet")
    if (file.exists(bnd)) {
        # readAndAddPolygonsToSPE() needs metadata that readXeniumSXE()
        # does not set, so read and add the polygons in two steps.
        pol <- tryCatch(SpaceTrooper::readPolygonsXenium(bnd, type = "parquet"),
            error = function(e) {
                message("  (cell boundaries not read: ", conditionMessage(e),
                        ")")
                NULL
            })
        if (!is.null(pol)) spe <- SpaceTrooper::addPolygonsToSPE(spe, pol,
            polygonsCol = "polygons")
    }
    spe <- setPlatform(spe, "xenium")
    message(sprintf("Xenium: %d cells, %d genes", ncol(spe), nrow(spe)))
    spe
}

#' Load 10x Visium HD cell-segmented outputs
#'
#' Reads the `segmented_outputs` folder of a Visium HD run (Space Ranger
#' cell segmentation) with `VisiumIO::TENxVisiumHD()`, as in the OSTA book
#' chapter "Visium HD (segmented)". Gene symbols become unique row names,
#' a `logcounts` assay is added if missing, and the platform is recorded.
#' For the full workflow (QC, annotation, niches) see
#' `system.file("scripts", "visiumhd_osta_workflow.R",
#' package = "iseeSpatial")`.
#'
#' @param segmented_outputs Path to the `segmented_outputs` directory.
#' @param images Image resolution to load (passed to `TENxVisiumHD`).
#' @return A `SpatialExperiment`.
#' @examples
#' d <- "VisiumHD/segmented_outputs"   # Space Ranger cell segmentation output
#' if (dir.exists(d)) spe <- load_visiumhd(d)
#' @export
load_visiumhd <- function(segmented_outputs, images = "lowres") {
    .needs("VisiumIO")
    spe <- BiocIO::import(VisiumIO::TENxVisiumHD(format = "h5", images = images,
        segmented_outputs = segmented_outputs))
    rd <- SummarizedExperiment::rowData(spe)
    if ("Symbol" %in% colnames(rd)) rownames(spe) <- make.unique(rd$Symbol)
    if (!"logcounts" %in% SummarizedExperiment::assayNames(spe))
        SingleCellExperiment::logcounts(spe) <- log1p(SingleCellExperiment::counts(spe))
    spe <- setPlatform(spe, "visiumhd")
    message(sprintf("Visium HD: %d cells, %d genes", ncol(spe), nrow(spe)))
    spe
}
