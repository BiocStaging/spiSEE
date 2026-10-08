# ============================================================
# Xenium prototype -- mirrors Intermediate_processing.R's approach,
# but for validating whether iseeSpatial's plotting layer (SpaceTrooper)
# works on Xenium data before touching any package code.
#
# Test data: 10x Genomics public "Xenium Prime Mouse Ileum (tiny)" bundle
# (CC BY 4.0), artificially subset to 3 small patches across 2 FOVs, used
# here purely to confirm the reader + plotting functions work -- not for
# real analysis.
#   https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/resources/xenium-example-data
#   https://cf.10xgenomics.com/samples/xenium/3.0.0/Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny/Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs.zip
# ============================================================

library(SpatialExperimentIO)
library(SpaceTrooper)
library(SingleCellExperiment)
library(SpatialExperiment)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load Xenium output directory into a SpatialExperiment
# ------------------------------------------------------------
xenium_dir <- "data-raw/xenium_test_data"  # unzipped 10x Xenium output bundle
spe <- readXeniumSXE(xenium_dir)

cat(sprintf("Loaded: %d genes x %d cells\n", nrow(spe), ncol(spe)))
cat("colData columns:\n")
print(colnames(colData(spe)))

# KEY FINDING: unlike CosMx, Xenium colData has no `fov` column at all --
# it's continuous imaging, not tiled into discrete FOVs. Confirmed by
# inspecting real output, not just docs.
stopifnot(!"fov" %in% colnames(colData(spe)))

# ------------------------------------------------------------
# 2. Normalize (same approach as Intermediate_processing.R)
# ------------------------------------------------------------
logcounts(spe) <- log1p(counts(spe))

# ------------------------------------------------------------
# 3. Plain, static vignette-style plots -- no Shiny, no interactivity.
#    Saved directly to disk. This dataset is tiny (artificially subset
#    for testing the reader/plotting mechanics, not real analysis), so
#    these are sanity-check plots, not scientifically meaningful ones.
# ------------------------------------------------------------
plot_dir <- "xenium_plots"
dir.create(plot_dir, showWarnings = FALSE)

# QC: total counts per cell
p_qc <- ggplot(as.data.frame(colData(spe)), aes(x = total_counts)) +
  geom_histogram(bins = 15, fill = "#118AB2") +
  theme_minimal() +
  labs(title = "Total transcript counts per cell", x = "Total counts", y = "Cells")
ggsave(file.path(plot_dir, "qc_total_counts.png"), p_qc, width = 6, height = 4, dpi = 150)

# Spatial: cells coloured by total counts (confirms plotCentroids works
# unmodified on Xenium data -- verified via source inspection that it has
# no `fov` dependency, unlike plotCellsFovs/plotZoomFovsMap).
colData(spe)[["total_counts_num"]] <- as.numeric(spe$total_counts)
p_spatial_qc <- SpaceTrooper::plotCentroids(spe, colourBy = "total_counts_num", size = 3)
ggsave(file.path(plot_dir, "spatial_total_counts.png"), p_spatial_qc, width = 6, height = 6, dpi = 150)

# Spatial: cells coloured by the most highly expressed gene
gene <- rownames(spe)[which.max(Matrix::rowSums(counts(spe)))]
colData(spe)[[gene]] <- as.numeric(logcounts(spe)[gene, ])
p_gene <- SpaceTrooper::plotCentroids(spe, colourBy = gene, size = 3)
ggsave(file.path(plot_dir, sprintf("spatial_gene_%s.png", gene)), p_gene, width = 6, height = 6, dpi = 150)

cat(sprintf("\nSaved 3 plots to %s/\n", plot_dir))

# ------------------------------------------------------------
# What's confirmed vs. what needs adaptation for a real Xenium panel:
#
# Works as-is (no `fov` dependency, verified via source + this test):
#   - SpaceTrooper::plotCentroids() / plotPolygons()  -> CosMxSpatialPanel,
#     CosMxGenePanel's underlying plot calls
#   - SpaceTrooper::plotMetricHist()                  -> CosMxQCPanel
#
# Needs adaptation (hard-depend on `fov`, verified via source):
#   - SpaceTrooper::plotCellsFovs() / plotZoomFovsMap() -> everything
#     CosMxFOVMapPanel/CosMxFOVPanel do (the whole click-to-select map
#     built for CosMx). Xenium has no discrete tiles to click on.
#     The natural port: the same plotly box-select mechanism already
#     built for the CosMx map (drag a rectangle to select) works
#     directly on Xenium's continuous coordinates -- no FOV grid
#     needed, just select cells within a dragged region.
# ============================================================
