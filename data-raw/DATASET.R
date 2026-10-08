# ============================================================
# iseeCosmos — Data Preparation Script
# ============================================================
# Run this to prepare CosMx data for iseeCosmos
# ============================================================

library(SpatialExperiment)
library(SingleCellExperiment)
library(scater)
library(scRNAseq)
library(SingleR)
library(tidyr)
library(dplyr)
library(sf)
library(SpatialFeatureExperiment)
library(spdep)
library(Voyager)
library(spatialFDA)
library(ggplot2)

set.seed(112358)

cat("\n=====================================================\n")
cat("   iseeCosmos Data Preparation\n")
cat("=====================================================\n\n")

# ============================================================
# 1. Load filtered data
# ============================================================
cat("Step 1/9 — Loading data...\n")
input_path <- "../cos_filtered.RDS"
spe <- readRDS(input_path)
cat(sprintf("  → %d cells, %d genes\n", ncol(spe), nrow(spe)))

# ============================================================
# 2. Fix duplicate column names
# ============================================================
cat("\nStep 2/9 — Fixing duplicate columns...\n")
original_names <- colnames(colData(spe))
colnames(colData(spe)) <- make.unique(colnames(colData(spe)))
dups_fixed <- sum(original_names != colnames(colData(spe)))
if (dups_fixed > 0) {
  cat(sprintf("  → Fixed %d duplicate names\n", dups_fixed))
}

# ============================================================
# 3. Add logcounts (normalized expression)
# ============================================================
cat("\nStep 3/9 — Computing normalized expression...\n")
# Normalize by cell area
logcounts(spe) <- sweep(assay(spe, "counts"), 2, spe$Area_um, `/`)
cat("  → Normalized by Area_um\n")

# ============================================================
# 4. Dimensionality reduction
# ============================================================
cat("\nStep 4/9 — Running PCA and UMAP...\n")
spe <- runPCA(spe, ncomponents = 50)
cat("  → PCA complete (50 components)\n")

spe <- runUMAP(spe, dimred = "PCA")
cat("  → UMAP complete\n")

# Add first 10 PCs to colData
pcs <- reducedDim(spe, "PCA")[, 1:10]
colnames(pcs) <- paste0("PC", 1:10)
colData(spe) <- cbind(colData(spe), pcs)

# ============================================================
# 5. Cell type annotation
# ============================================================
cat("\nStep 5/9 — Annotating cell types with SingleR...\n")
ref <- ZeiselBrainData()
ref <- logNormCounts(ref)

pred <- SingleR(
  test = spe,
  ref = ref,
  labels = ref$level1class,
  de.method = "wilcox"
)

spe$SingleR_label <- factor(pred$pruned.labels)
n_types <- nlevels(spe$SingleR_label)
cat(sprintf("  → Found %d cell types\n", n_types))

# Add colors
pal <- setNames(hcl.colors(n_types, "Spectral"), levels(spe$SingleR_label))
spe$SingleR_label_colors <- pal[as.character(spe$SingleR_label)]

# ============================================================
# 6. Gene ranks
# ============================================================
cat("\nStep 6/9 — Ranking genes by expression...\n")
gene_means <- rowMeans(logcounts(spe))
metadata(spe)$gene_ranks <- names(sort(gene_means, decreasing = TRUE))
cat("  → Top genes:", paste(head(metadata(spe)$gene_ranks, 5), collapse = ", "), "\n")

# ============================================================
# 7. Per-FOV statistics (CRITICAL!)
# ============================================================
cat("\nStep 7/9 — Computing per-FOV statistics...\n")

make_fov_long <- function(spe) {
  df <- data.frame(colData(spe))

  # Metrics to include
  metric_cols <- c("Mean.DAPI", "Mean.GFAP", "detected", "total", "Area_um", "log2CountArea")
  available <- intersect(metric_cols, names(df))

  cat("  → Metrics:", paste(available, collapse = ", "), "\n")

  tidyr::pivot_longer(df, all_of(available)) |>
    dplyr::mutate(
      value = dplyr::case_when(
        grepl("^Mean", name) ~ asinh(value),
        name == "total" ~ log10(value + 1),
        name == "Area_um" ~ log10(value + 1),
        TRUE ~ value
      ),
      fov = factor(as.integer(fov))
    )
}

metadata(spe)$fov_stats <- make_fov_long(spe)
cat(sprintf("  → FOV stats: %d rows\n", nrow(metadata(spe)$fov_stats)))

# ============================================================
# 8. Spatially-variable genes (vignette region)
# ============================================================
# The Spatial Statistics panel itself (graph/Moran's I/join-count/L-function
# plots) is computed live in the app now -- see R/spatial_stats_compute.R,
# driven by whatever FOVs are selected via the overview map, not a fixed
# precomputed region. Only the SVG ranking below still needs a one-time
# precompute here: it's a Moran's I scan over 200 genes, a different order
# of computational cost than a single-graph-type stat, not something to
# redo on every FOV change.
cat("\nStep 8/9 — Computing spatially-variable gene ranking (vignette region)...\n")

# Same FOVs as the Analysis app's default zoom (see R/launch.R) and the
# OSTA-book-style walkthrough in Intermediate_processing.R.
vignette_fovs <- c(72:74, 90:92, 97:99, 114:116)
region_spe <- spe[, as.integer(as.character(spe$fov)) %in% vignette_fovs]
region_spe <- region_spe[, !is.na(region_spe$SingleR_label)]
cat(sprintf("  → Vignette region: %d cells across %d FOVs\n",
            ncol(region_spe), length(vignette_fovs)))

# Convert to SpatialFeatureExperiment for Voyager's neighbor-graph tooling
region_sfe <- as(region_spe, "SpatialFeatureExperiment")
colnames(region_sfe$polygons)[6] <- "geometry"
st_geometry(region_sfe$polygons) <- "geometry"
colGeometry(region_sfe, "cellSeg") <- region_sfe$polygons

colGraph(region_sfe, "poly2nb") <- findSpatialNeighbors(
  region_sfe, type = "cellSeg", method = "poly2nb", style = "W", snap = 50)

# Spatially-variable genes: rank the top 200 most-expressed genes (from the
# gene_ranks computed above) by global Moran's I on this region's poly2nb
# graph, so the Gene Expression panel's gene picker can be driven by spatial
# autocorrelation instead of only average expression.
svg_candidates <- intersect(head(metadata(spe)$gene_ranks, 200), rownames(region_sfe))
moran_df <- calculateMoransI(region_sfe, features = svg_candidates,
                              colGraphName = "poly2nb", colGeometryName = "cellSeg",
                              zero.policy = TRUE)
svg_table <- data.frame(gene = rownames(moran_df), moran_i = moran_df$moran,
                         row.names = NULL)
svg_table <- svg_table[order(-svg_table$moran_i), ]

metadata(spe)$regions <- list(
  vignette = list(
    label     = "Vignette region (FOVs 72-74, 90-92, 97-99, 114-116)",
    svg_table = svg_table
  )
)
cat(sprintf("  → SVG ranking: %d genes ranked by Moran's I\n", nrow(svg_table)))

# ============================================================
# 9. Save
# ============================================================
cat("\nStep 9/9 — Saving...\n")

# Save to main location
output_path <- "../cos_isee_ready.RDS"
saveRDS(spe, output_path)
cat("  → Saved to:", output_path, "\n")

# Also save to package data-raw
if (dir.exists("data-raw")) {
  saveRDS(spe, "data-raw/cos_isee_ready.RDS")
  cat("  → Also saved to: data-raw/cos_isee_ready.RDS\n")
}

cat("\n=====================================================\n")
cat("   ✓ Preparation Complete!\n")
cat("=====================================================\n\n")
cat("Summary:\n")
cat("  Cells:      ", ncol(spe), "\n")
cat("  Genes:      ", nrow(spe), "\n")
cat("  FOVs:       ", length(unique(spe$fov)), "\n")
cat("  Cell types: ", n_types, "\n")
cat("  FOV stats:  ", nrow(metadata(spe)$fov_stats), "rows\n\n")

cat("Launch with:\n")
cat("  library(iseeCosmos)\n")
cat("  spe <- readRDS('", output_path, "')\n", sep = "")
cat("  iseeCosmosAnalysis(spe)\n\n")
