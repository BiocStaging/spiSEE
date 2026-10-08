# ============================================================
# iseeSpatial — Xenium Data Preparation Script
# ============================================================
# Xenium counterpart of data-raw/DATASET.R: prepares real Xenium data for
# the Xenium apps and saves the same two-RDS split CosMx's pipeline
# produces upstream of DATASET.R (see load_cosmx.R) -- an unfiltered
# snapshot for the QC app, an analysis-ready one for the Analysis app.
#
# Computation is adapted from Xenium_OSTA_workflow.R (the OSTA book's
# Xenium chapter, reproduced there almost verbatim) -- this script
# restructures the same steps into DATASET.R's steps-then-save shape
# instead of ending in saved PNGs, and adds the QC-app snapshot
# Xenium_OSTA_workflow.R has no reason to produce itself.
#
# NOT YET RUN as of writing. Two real costs before this can be:
#   - Downloads the real "Xenium_HumanColon_Oliveira" dataset *and* a
#     ~280k-cell scRNA-seq reference ("Chromium_HumanColon_Oliveira",
#     same source) for SingleR annotation.
#   - SingleR against that reference is the workflow's own documented
#     slow step (Xenium_OSTA_workflow.R checkpoints it separately, to
#     /tmp/xenium_osta_sub_annotated.rds, for exactly this reason).
# Run once, from the repo root, when ready to pay that cost:
#   Rscript data-raw/XENIUM_DATASET.R
# ============================================================

library(BiocParallel)
library(DropletUtils)
library(OSTA.data)
library(scater)
library(scrapper)
library(SingleR)
library(SpaceTrooper)
library(SpatialExperiment)
library(SpatialExperimentIO)
library(SingleCellExperiment)

set.seed(112358)
# Policlinico has 192 cores and is otherwise near-idle -- 32 is generous
# without hogging a shared server other people are actively using.
bp <- MulticoreParam(32)

cat("\n=====================================================\n")
cat("   iseeSpatial Xenium Data Preparation\n")
cat("=====================================================\n\n")

# ============================================================
# 1. Load raw Xenium data
# ============================================================
cat("Step 1/7 — Loading Xenium Human Colon data...\n")
id <- "Xenium_HumanColon_Oliveira"
repo_wd <- getwd()  # OSTA.data_load() setwd()s into its download tempdir
pa <- OSTA.data_load(id, mol = FALSE)  # and never restores it -- see below
setwd(repo_wd)      # so every relative path after this stays repo-relative
dir.create(td <- tempfile())
unzip(pa, exdir = td)
# See load_xenium.R for why addParquetPaths=FALSE: book uses
# readXeniumSXE(td, addTx=FALSE), which doesn't exist in this installed
# SpatialExperimentIO version -- addParquetPaths=FALSE achieves the same
# thing (skip transcripts.parquet, not present since mol=FALSE above).
spe <- readXeniumSXE(td, addParquetPaths = FALSE)
cat(sprintf("  -> %d cells, %d genes\n", ncol(spe), nrow(spe)))

# Cell boundary polygons, when the download includes them -- same approach
# as load_xenium.R's handling of the small bundled test dataset. Missing/
# unreadable boundaries are non-fatal (UsePolygons defaults to FALSE
# everywhere); this just means the real dataset gets to use them too
# when they're there, instead of only ever being available for the tiny
# test bundle.
boundaries_file <- file.path(td, "cell_boundaries.parquet")
if (file.exists(boundaries_file)) {
  cat("  found cell_boundaries.parquet -- attaching polygons...\n")
  polygons <- tryCatch(
    SpaceTrooper::readPolygonsXenium(boundaries_file, type = "parquet"),
    error = function(e) {
      cat("  (couldn't read cell boundaries: ", conditionMessage(e), ")\n")
      NULL
    })
  if (!is.null(polygons)) {
    spe <- SpaceTrooper::addPolygonsToSPE(spe, polygons, polygonsCol = "polygons")
    cat("  -> polygons attached\n")
  }
} else {
  cat("  no cell_boundaries.parquet in this download -- polygons unavailable, centroids only.\n")
}

# ============================================================
# 2. QC metrics, outlier fences + unfiltered snapshot (feeds the QC app)
# ============================================================
cat("\nStep 2/7 — Computing QC metrics and outlier fences...\n")
qc  <- computeRnaQcMetrics(counts(spe), subsets = list())
thr <- suggestRnaQcThresholds(qc)
spe$sum      <- qc$sum
spe$detected <- qc$detected
spe$keep     <- filterRnaQcMetrics(thr, qc)
logcounts(spe) <- log1p(counts(spe))  # cheap default for the QC app's raw view

# Same medcouple-based fence CosMx's Import_data.R computes
# (SpaceTrooper::computeSpatialOutlier(method="mc")) -- generic, not
# CosMx-specific, so it applies to Xenium's QC metrics unchanged. Feeds
# XeniumQCPanel's "Show outlier fences" checkbox (xenium_qc_panel.R).
for (metric in c("transcript_counts", "total_counts", "cell_area")) {
  spe <- computeSpatialOutlier(spe, computeBy = metric, method = "mc")
}
cat(sprintf("  -> %.1f%% of cells pass QC\n", 100 * mean(spe$keep)))

cat("\nSaving unfiltered snapshot (QC app)...\n")
saveRDS(spe, "../xenium_unfiltered.RDS")
if (dir.exists("data-raw")) saveRDS(spe, "data-raw/xenium_unfiltered.RDS")

# ============================================================
# 3. Filter + crop to the book's demo region
# ============================================================
cat("\nStep 3/7 — Filtering and cropping...\n")
spe <- spe[, spe$keep]
# Same crop box as the book / Xenium_OSTA_workflow.R.
box <- list(xmin = 2e3, xmax = 5e3, ymin = 1e3, ymax = 4e3)
xy  <- spatialCoords(spe)
sub <- spe[, xy[, 1] > box$xmin & xy[, 1] < box$xmax &
               xy[, 2] > box$ymin & xy[, 2] < box$ymax]
cat(sprintf("  -> %d cells in the demo region\n", ncol(sub)))

# ============================================================
# 4. Normalize + PCA
# ============================================================
cat("\nStep 4/7 — Normalizing and running PCA...\n")
sfs <- sub$cell_area / median(sub$cell_area)
logcounts(sub) <- normalizeCounts(counts(sub), size.factors = sfs)
pca <- runPca(as.matrix(logcounts(sub)[rownames(sub), ]), number = 25)
reducedDim(sub, "PCA") <- t(pca$components)

# ============================================================
# 5. Leiden clustering
# ============================================================
cat("\nStep 5/7 — Clustering (Leiden)...\n")
snn <- buildSnnGraph(pca$components, weight.scheme = "jaccard", num.threads = 32)
sub$Leiden <- factor(clusterGraph(snn, method = "leiden", leiden.resolution = 0.5)$membership)
cat(sprintf("  -> %d clusters\n", nlevels(sub$Leiden)))

# ============================================================
# 6. Cell-type annotation (SingleR vs. matched scRNA-seq reference)
# ============================================================
# The slow step -- see Xenium_OSTA_workflow.R's own checkpointing of
# this exact step for why (SingleR against ~280k reference cells).
cat("\nStep 6/7 — Annotating cell types (SingleR)...\n")
ref_id <- "Chromium_HumanColon_Oliveira"
ref_pa <- OSTA.data_load(ref_id)
dir.create(ref_td <- tempfile())
unzip(ref_pa, exdir = ref_td)

ref_sce <- read10xCounts(list.files(ref_td, "h5$", full.names = TRUE))
ref_cd  <- read.csv(list.files(ref_td, "cell_meta", full.names = TRUE))
colData(ref_sce) <- cbind(colData(ref_sce), ref_cd[, -1])
ref_sce <- ref_sce[, ref_sce$QCFilter == "Keep" & grepl("P2", ref_sce$Patient)]
assay(ref_sce) <- as(assay(ref_sce), "dgCMatrix")
logcounts(ref_sce) <- normalizeCounts(assay(ref_sce), size.factors = rep(1, ncol(ref_sce)))
ref_sce <- ref_sce[rowData(ref_sce)$ID %in% rowData(sub)$ID, ]
rownames(ref_sce) <- rowData(ref_sce)$Symbol

res <- SingleR(
  test = sub, ref = ref_sce, labels = ref_sce$Level2,
  de.method = "wilcox", aggr.ref = TRUE, BPPARAM = bp)
sub$Level2 <- factor(res$pruned.labels)

idx <- match(sub$Level2, ref_sce$Level2)
sub$Level1 <- factor(ref_sce$Level1[idx])

lab <- list(
  tum = c("Tumor"),
  epi = c("Intestinal Epithelial"),
  imm = c("B cells", "T cells", "Myeloid"),
  str = c("Endothelial", "Fibroblast", "Smooth Muscle"))
idx2     <- match(sub$Level1, unlist(lab))
lab_flat <- rep.int(names(lab), lengths(lab))
sub$Level0 <- factor(lab_flat[idx2])
cat(sprintf("  -> %d Level1 cell types\n", nlevels(sub$Level1)))

# ============================================================
# 7. Save analysis-ready snapshot (Analysis app)
# ============================================================
cat("\nStep 7/7 — Saving...\n")
output_path <- "../xenium_isee_ready.RDS"
saveRDS(sub, output_path)
if (dir.exists("data-raw")) saveRDS(sub, "data-raw/xenium_isee_ready.RDS")

cat("\n=====================================================\n")
cat("   Preparation Complete\n")
cat("=====================================================\n\n")
cat("Summary:\n")
cat("  Unfiltered cells: ", ncol(spe), "\n")
cat("  Analysis cells:   ", ncol(sub), "\n")
cat("  Genes:            ", nrow(sub), "\n")
cat("  Level1 types:     ", nlevels(sub$Level1), "\n\n")
cat("Launch with:\n")
cat("  library(iseeSpatial)\n")
cat("  iseeXeniumQC(readRDS('data-raw/xenium_unfiltered.RDS'))\n")
cat("  iseeXeniumAnalysis(readRDS('data-raw/xenium_isee_ready.RDS'))\n\n")
