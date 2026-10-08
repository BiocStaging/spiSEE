#!/usr/bin/env Rscript
# ============================================================
# Visium HD (cell-segmented): full OSTA chapter + app-ready object.
#
# Part 1 reproduces EVERY code chunk of the OSTA book chapter
#   "Workflow: Visium HD (seg)"
#   https://bioconductor.org/books/release/OSTA/pages/seq-workflow-visium-hd-seg.html
# in order (sections 17.2-17.6), with each figure saved to <figdir>/ and each
# printed result written to the console (pipe to a log file to keep it).
# Code is the chapter's; only comments, print() and figure saving were added.
#
# Part 2 (iseeSpatial additions, NOT in the chapter): gene ranking, PCA,
# neighbour graph, niches, spatial domains; saves the filtered object for
# the analysis app and the unfiltered one for the QC app.
#
#   Rscript inst/scripts/visiumhd_osta_workflow.R [out.RDS] [k] [figdir] 2>&1 | tee visiumhd_osta.log
#
# Needs: OSTA.data, VisiumIO, scrapper, SpotSweeper, SingleR, DropletUtils,
# Voyager, SpatialFeatureExperiment, ggspavis, sf, arrow, dplyr, tidyr,
# patchwork, pals, scater, devtools (+ iseeSpatial from this folder).
# Run from the package root. Memory-heavy: prefer the server.
# ============================================================

args   <- commandArgs(trailingOnly = TRUE)
out    <- if (length(args) >= 1) args[1] else "visiumhd_isee_ready.RDS"
k      <- if (length(args) >= 2) as.integer(args[2]) else 20L   # neighbourhood size; adapt to the structures
figdir <- if (length(args) >= 3) args[3] else "osta_visiumhd_figures"
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
fig <- function(p, name, w = 10, h = 6) {
  ggsave(file.path(figdir, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  message("figure: ", name)
}
t_start <- Sys.time()

# ---- 17.1.2 Dependencies ------------------------------------------------
suppressPackageStartupMessages({
  library(sf); library(arrow); library(dplyr); library(tidyr)
  library(scrapper); library(SingleR); library(Voyager); library(ggplot2)
  library(ggspavis); library(VisiumIO); library(OSTA.data); library(patchwork)
  library(SpotSweeper); library(DropletUtils)
  library(SpatialExperiment); library(SpatialFeatureExperiment)
})
# set seed for random number generation in order to make results reproducible
set.seed(777)

# utility for cropping by bounding box
.crop <- \(spe, box) {
  box <- as.list(box)
  xy <- spatialCoords(spe)
  spe[,
      xy[,1] > box$xmin & xy[,1] < box$xmax &
      xy[,2] > box$ymin & xy[,2] < box$ymax ]
}

# ---- 17.2.1 From SPE... -------------------------------------------------
# retrieve data from OSF repo
id <- "VisiumHD_HumanColon_Oliveira"
pa <- OSTA.data_load(id)
dir.create(td <- tempfile())
unzip(pa, exdir = td)

# read into 'SpatialExperiment'
seg <- file.path(td, "segmented_outputs")
spe <- TENxVisiumHD(
  format = "h5",
  images = "lowres",
  segmented_outputs = seg) |>
  import()

# make gene symbols unique
gs <- rowData(spe)$Symbol
rownames(spe) <- make.unique(gs)

# needed for 'ggspavis'
spe$in_tissue <- TRUE

# view bin-to-cell mapping information
print(head(spe$map[[1]]))
# number of segmented 2 µm bins
print(summary(sapply(spe$map, nrow)))
# fraction of 2 µm bins assigned to nucleus
print(summary(sapply(spe$map, \(.) mean(.$in_nucleus))))

xy <- spatialCoords(spe)
xs <- range(xy[, 1])
ys <- range(xy[, 2])
dx <- diff(xs)/4
dy <- diff(ys)/4
box <- list(
  xmin = xs[1]+dx, xmax = xs[2]-dx,
  ymin = ys[1]+dy, ymax = ys[2]-dy)
sub <- .crop(spe, box)

# reverse y-coordinates of bounding box
.box <- box
.box$ymax <- -box$ymin
.box$ymin <- -box$ymax
aes <- list(col = "red", fill = NA, linewidth = 2)
p <- plotCoords(
  spe[, sample(ncol(spe), 5e4)]) +
  do.call(geom_rect, c(.box, aes)) +
  plotCoords(
    sub[, sample(ncol(sub), 5e4)])
fig(p, "17.2.1_crop")

print(round(100*ncol(sub)/ncol(spe), 2))
print(ncol(spe <- sub))

roi <- c(xmin = 54e3, xmax = 58e3, ymin = 12e3, ymax = 16e3)
spe$roi <- colnames(spe) %in% colnames(.crop(spe, roi))
p <- plotCoords(spe, annotate = "roi", point_size = 0) + theme(legend.position = "none")
fig(p, "17.2.1_roi", w = 7, h = 7)

# ---- 17.2.2 ...to SFE ---------------------------------------------------
# convert from SPE to SFE
sfe <- toSpatialFeatureExperiment(spe)
# add cell segmentation boundaries
seg <- metadata(spe)$cellseg
i <- colnames(spe)
j <- match(i, seg$cell_id)
seg <- seg[j, ]; rownames(seg) <- i
colGeometries(sfe) <- list(cellseg = seg)

# ---- 17.3.1 Cell masks --------------------------------------------------
# specify bounding box
box <- c(
  xmin = 56e3, xmax = 58e3,
  ymin = 15e3, ymax = 16e3)
# plot exemplary cell masks
geo <- colGeometry(crop(sfe, box))
png(file.path(figdir, "17.3.1_cell_masks.png"), width = 1500, height = 800, res = 150)
plot(st_geometry(geo), col = rep(colors(), 2))
dev.off(); message("figure: 17.3.1_cell_masks")

# ---- 17.3.2 Bin mapping -------------------------------------------------
# get cell areas = 4x number of 2 µm bins
sfe$um2 <- 4*sapply(sfe$map, nrow)
# get fraction of 2 µm bins that are nuclear
sfe$nuc <- sapply(sfe$map, \(.) mean(.$in_nucleus))
# get nucleus areas = cell area x nuclear fraction
sfe$nuc_um2 <- sfe$um2*sfe$nuc

names(bin) <- bin <- c("square_008um", "square_016um")
ncs <- lapply(bin, \(um) {
  ns <- lapply(spe$map, \(df) unique(df[[um]]))
  as.vector(table(unlist(ns)))
})
print(sapply(ncs, summary))

# ---- 17.4.1 Quality control: non-spatial --------------------------------
mt <- grepl("^MT-", rownames(sfe))
sfe <- quickRnaQc.se(sfe, subsets = list(mt = mt))
sfe$log_sum <- log1p(sfe$sum)
print(head(colData(sfe)[c(
  "um2", "nuc_um2", "nuc",
  "sum", "subset.proportion.mt")]))

p <- plotSpatialFeature(sfe[, sfe$roi],
  colGeometryName = "cellseg",
  features = "log_sum") +
  ggtitle("log-library size") +
  plotSpatialFeature(sfe[, sfe$roi],
    colGeometryName = "cellseg",
    features = "subset.proportion.mt") +
  ggtitle("mitochondrial proportion") +
  plot_layout(nrow = 1) &
  theme(plot.title = element_text(hjust = 0.5)) &
  scale_fill_gradientn(NULL, colors = pals::jet())
fig(p, "17.4.1_qc_roi", w = 12, h = 6)

# ---- 17.4.2 Quality control: spatially-aware ----------------------------
# determine spatial outliers for different metrics
sfe <- localOutliers(sfe,
  workers = 4, metric = "sum",
  log = TRUE, direction = "lower")
sfe <- localOutliers(sfe,
  workers = 4, metric = "detected",
  log = TRUE, direction = "lower")
sfe <- localOutliers(sfe,
  workers = 4, metric = "subset.proportion.mt",
  log = FALSE, direction = "higher")
# get cells x flags matrix
cd <- colData(sfe)
ol <- grep("outliers$", names(cd))
ol <- as.matrix(cd[ol])
# percentage of cells excluded for each
print(round(100*colMeans(ol), 2))
# number of cells kept/removed overall
print(table(sfe$ex <- rowAnys(ol)))

p <- plotCoords(sfe[, order(sfe$ex)],
  annotate = "ex", point_size = 0.2) +
  theme(legend.key.size = ggplot2::unit(0, "pt")) +
  scale_color_manual(values = c("grey90", "red")) +
  guides(col = guide_legend(override.aes = list(size = 2)))
fig(p, "17.4.2_outliers", w = 7, h = 7)

# iseeSpatial: keep the unfiltered object for the QC app (flags visible there).
sfe_unfiltered <- sfe

# apply filtering
sfe <- sfe[, !sfe$ex]

# ---- 17.5 Annotation ----------------------------------------------------
# retrieve single-cell reference from OSF repo
id <- "Chromium_HumanColon_Oliveira"
pa <- OSTA.data_load(id)
dir.create(td <- tempfile())
unzip(pa, exdir = td)

# read into 'SingleCellExperiment'
fs <- list.files(td, full.names = TRUE)
h5 <- grep("h5$", fs, value = TRUE)
sce <- read10xCounts(h5, col.names = TRUE)

# add cell metadata
csv <- grep("csv$", fs, value = TRUE)
cd <- read.csv(csv, row.names = 1)
colData(sce)[names(cd)] <- cd[colnames(sce), ]

# use gene symbols as feature names
rownames(sce) <- make.unique(rowData(sce)$Symbol)

# exclude cells deemed to be of low-quality
sce <- sce[, sce$QCFilter == "Keep"]

# replace problematic characters
print(table(sce$Level1 <- gsub("\\s", "\\.", sce$Level1)))

# (this is done only to keep runtime/memory low)
idx <- split(seq_len(ncol(sce)), sce$Level1)
idx <- lapply(idx, \(.) sample(., min(length(.), 1e3)))
print(ncol(.sce <- sce[, unlist(idx)]))

# log-library size normalization
sfe <- normalizeRnaCounts.se(sfe)
.sce <- normalizeRnaCounts.se(.sce)
# perform label transfer at the Visium HD cell-level,
# using pseudo-bulk Chromium profiles as reference
res <- SingleR(test = sfe, ref = .sce, labels = .sce$Level1, aggr.ref = TRUE)

sfe$Level1 <- factor(res$pruned.labels)

lab <- list(
  tum = c("Tumor"),
  epi = c("Intestinal.Epithelial"),
  imm = c("B.cells", "T.cells", "Myeloid"),
  str = c("Endothelial", "Fibroblast", "Smooth.Muscle"))
idx <- match(sfe$Level1, unlist(lab))
lab <- rep.int(names(lab), sapply(lab, length))
print(table(sfe$Level0 <- factor(lab[idx])))

pal_lv1 <- unname(pals::trubetskoy())
pal_lv0 <- c("gold", "turquoise", "deeppink", "navy")
p <- plotSpatialFeature(sfe[, sfe$roi],
  features = "Level1", colGeometryName = "cellseg") +
  scale_fill_manual(values = pal_lv1) +
  plotSpatialFeature(sfe[, sfe$roi],
    features = "Level0", colGeometryName = "cellseg") +
  scale_fill_manual(values = pal_lv0) +
  plot_layout(nrow = 1) &
  theme(legend.key.size = ggplot2::unit(.5, "lines"))
fig(p, "17.5_annotation_roi", w = 12, h = 6)

p <- plotCoords(sfe,
  annotate = "Level1", point_size = 0.1) +
  scale_color_manual(values = pal_lv1) +
  plotCoords(sfe,
    annotate = "Level0", point_size = 0.1) +
  scale_color_manual(values = pal_lv0) +
  plot_layout() &
  theme(legend.key.size = ggplot2::unit(0, "pt"))
fig(p, "17.5_annotation_full", w = 12, h = 6)

# ---- 17.6 Downstream ----------------------------------------------------
# wrangling
df <- data.frame(colData(sfe))
df <- filter(df, !is.na(Level0))
fd <- pivot_longer(df, ends_with("um2"))
fd$name <- factor(fd$name, labels = c("nuc.", "cell"))
ws <- c(nlevels(df$Level0), nlevels(df$Level1))
# plotting
p <- ggplot(fd, aes(reorder(Level0, value), value, col = name)) +
  ggplot(fd, aes(reorder(Level1, value), value, col = name)) +
  plot_layout(nrow = 1, widths = ws, guides = "collect") &
  geom_boxplot(key_glyph = "point", outlier.size = 0) &
  labs(col = "area (um2)") &
  theme_bw() & theme(
    axis.title = element_blank(),
    legend.key.size = ggplot2::unit(0, "pt"),
    axis.text.x = element_text(angle = 45, hjust = 1))
fig(p, "17.6_areas_by_celltype", w = 12, h = 5)

message("OSTA chapter reproduced in ", format(Sys.time() - t_start))

# ---- Part 2: iseeSpatial additions (not in the chapter) -----------------
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))

# TODO(verify): assay name written by normalizeRnaCounts.se(); the app needs "logcounts".
if (!"logcounts" %in% assayNames(sfe)) stop("No 'logcounts' assay: check normalizeRnaCounts.se() output name")
logcounts(sfe_unfiltered) <- log1p(counts(sfe_unfiltered))

# gene menus ordered by mean expression
metadata(sfe)$gene_ranks <- rownames(sfe)[order(Matrix::rowMeans(logcounts(sfe)), decreasing = TRUE)]
metadata(sfe_unfiltered)$gene_ranks <- metadata(sfe)$gene_ranks

# PCA for spatial domains
sfe <- scater::runPCA(sfe, ncomponents = 20, ntop = 2000)

sfe <- setPlatform(sfe, "visiumhd")
sfe <- buildSpatialGraph(sfe, k = k)
sfe <- computeNiches(sfe, label_col = "Level1", k = k, n_niches = 8)
sfe <- computeSpatialDomains(sfe, k_spatial = k, lambda = 0.2)
print(table(sfe$niche, useNA = "ifany"))
print(table(sfe$spatial_domain, sfe$Level0, useNA = "ifany"))

saveRDS(sfe, out)
saveRDS(setPlatform(sfe_unfiltered, "visiumhd"), sub("\\.RDS$", "_unfiltered.RDS", out))
message("Saved ", out, " (analysis app) and ", sub("\\.RDS$", "_unfiltered.RDS", out),
        " (QC app). Total time ", format(Sys.time() - t_start))
