# Visium HD vignette: the code of vignettes/articles/visiumhd_colorectal.Rmd, run
# locally (the dataset is 885 MB, too large for the Bioconductor builders),
# with the static figures saved to vignettes/figures/.
# Run from the package root:  Rscript data-raw/visiumhd_vignette.R
# The app screenshots were taken from the objects saved at the end.

library(spiSEE)
library(SpaceTrooper)
library(SingleR)
library(ggplot2)
set.seed(1)
figdir <- "vignettes/articles/figures"
fig <- function(p, name, w = 7, h = 5) {
    ggsave(file.path(figdir, paste0(name, ".jpg")), p, width = w, height = h,
           dpi = 72, quality = 45)
}

# ---- Reading the data
wd <- getwd()
zip <- OSTA.data::OSTA.data_load("VisiumHD_HumanColon_Oliveira")
setwd(wd)    # OSTA.data_load() changes the working directory
dir <- tempfile()
unzip(zip, exdir = dir)
seg <- file.path(dir, "segmented_outputs")
spe <- loadVisiumHD(seg)

# ---- The whole section
xy <- spatialCoords(spe)
markers <- c("CEACAM5", "COL1A1", "DES", "IGKC")
df <- do.call(rbind, lapply(markers, function(g) {
    v <- as.numeric(logcounts(spe)[g, ])
    data.frame(x = xy[, 1], y = xy[, 2], gene = g, logcounts = v)[v > 0, ]
}))
box <- data.frame(xmin = 47000, xmax = 54000, ymin = 9500, ymax = 16500)
p <- ggplot(df, aes(x, y, colour = logcounts)) +
    geom_point(size = 0.01) +
    geom_rect(data = box, inherit.aes = FALSE, colour = "red", fill = NA,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
    facet_wrap(~gene, nrow = 1) +
    scale_colour_viridis_c() +
    coord_fixed() +
    theme_void()
fig(p, "visiumhd_section", w = 12, h = 3.6)

# ---- The region
in_box <- xy[, 1] > 47000 & xy[, 1] < 54000 & xy[, 2] > 9500 & xy[, 2] < 16500
roi <- spe[, in_box]
counts(roi) <- as(counts(roi), "dgCMatrix")
roi <- scuttle::logNormCounts(roi)

polys <- sf::st_read(file.path(seg, "cell_segmentations.geojson"),
                     quiet = TRUE)
id <- as.integer(sub("cellid_0*([0-9]+)-1", "\\1", colnames(roi)))
# Same layout as SpaceTrooper's polygons (an sf table, geometry "global"),
# which plotPolygons() expects.
roi$polygons <- sf::st_sf(
    cell_id = id, global = sf::st_geometry(polys)[match(id, polys$cell_id)])
um_per_px <- jsonlite::fromJSON(
    file.path(seg, "spatial", "scalefactors_json.json"))$microns_per_pixel
roi$um2 <- as.numeric(sf::st_area(sf::st_set_crs(roi$polygons, NA))) *
    um_per_px^2

# ---- Quality control
roi <- scuttle::addPerCellQCMetrics(roi,
    subsets = list(mito = grepl("^MT-", rownames(roi))))
roi <- computeSpatialOutlier(roi, computeBy = "sum", method = "mc")
roi <- computeSpatialOutlier(roi, computeBy = "um2", method = "mc")
roi$high_mito <- scuttle::isOutlier(roi$subsets_mito_percent,
                                    type = "higher")
saveRDS(roi, "visiumhd_region_unfiltered.rds")    # QC app

cells <- roi[, roi$sum_outlier_mc == "NO" & roi$um2_outlier_mc == "NO" &
                !roi$high_mito]

# ---- Annotation
wd <- getwd()
ref_zip <- OSTA.data::OSTA.data_load("Chromium_HumanColon_Oliveira")
setwd(wd)
ref_dir <- tempfile()
unzip(ref_zip, exdir = ref_dir)
ref <- DropletUtils::read10xCounts(
    file.path(ref_dir, "filtered_feature_bc_matrix.h5"), col.names = TRUE)
md <- read.csv(file.path(ref_dir, "cell_metadata.csv"), row.names = 1)
colData(ref)[names(md)] <- md[colnames(ref), ]
rownames(ref) <- make.unique(rowData(ref)$Symbol)
ref <- ref[, ref$QCFilter == "Keep"]
idx <- split(seq_len(ncol(ref)), ref$Level2)
idx <- lapply(idx, function(i) i[sample.int(length(i), min(length(i), 500))])
ref <- ref[, unlist(idx)]
counts(ref) <- as(counts(ref), "dgCMatrix")
ref <- scuttle::logNormCounts(ref)
res <- SingleR(test = cells, ref = ref, labels = ref$Level2,
               aggr.ref = TRUE)
cells$cell_type <- factor(res$pruned.labels)
cells <- cells[, !is.na(cells$cell_type)]

top <- names(sort(table(cells$cell_type), decreasing = TRUE))[1:12]
cells$cell_type_top <- factor(ifelse(cells$cell_type %in% top,
                                     as.character(cells$cell_type), "other"),
                              levels = c(top, "other"))
pal <- c(setNames(unname(pals::glasbey(12)), top), other = "grey85")
p <- plotPolygons(cells, colourBy = "cell_type_top", palette = pal,
                  sampleId = NULL) +
    theme(legend.title = element_blank())
fig(p, "visiumhd_celltypes", w = 9, h = 7)

# ---- Niches and domains
cells <- buildSpatialGraph(cells, k = 15)
cells <- computeNiches(cells, label_col = "cell_type", k = 15, n_niches = 6)
p <- plotPolygons(cells, colourBy = "niche",
                  palette = unname(pals::glasbey(6)), sampleId = NULL) +
    theme(legend.position = "bottom", legend.direction = "vertical",
          legend.title = element_blank())
fig(p, "visiumhd_niches", w = 8, h = 8)

x <- logcounts(cells)
m <- Matrix::rowMeans(x)
v <- Matrix::rowMeans(x^2) - m^2
hvg <- names(sort(v, decreasing = TRUE))[1:2000]
pca <- prcomp(t(as.matrix(x[hvg, ])), rank. = 20)
reducedDim(cells, "PCA") <- pca$x
cells <- computeSpatialDomains(cells, k_spatial = 15, lambda = 0.5,
                               resolution = 0.5)
xy <- spatialCoords(cells)
df <- data.frame(x = unname(xy[, 1]), y = unname(xy[, 2]),
                 cluster = c(as.character(cells$expr_cluster),
                             as.character(cells$spatial_domain)),
                 type = rep(c("expression clusters", "spatial domains"),
                            each = ncol(cells)))
p <- ggplot(df, aes(x, y, colour = cluster)) +
    geom_point(size = 0.2) +
    facet_wrap(~type) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none")
fig(p, "visiumhd_domains", w = 9, h = 4.5)

saveRDS(cells, "visiumhd_region.rds")             # analysis app
