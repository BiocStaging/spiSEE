#' @import SpatialFeatureExperiment
NULL

# ============================================================
# Spatial statistics on a selected region (merged CosMx + Xenium versions).
#
# The neighbour graph is k-NN on cell centroids, rebuilt on the selection so
# that every cell has k neighbours inside the region. k is a parameter: the
# right scale depends on the structures being looked at.
#
# Why k-NN and not polygon adjacency (poly2nb): on CosMx poly2nb took ~10 s
# for ~5,700 cells and left many cells without neighbours; centroids also
# exist on every platform, polygons do not.
#
# `as(se, "SpatialFeatureExperiment")` needs SpatialFeatureExperiment's
# coercion method registered, hence the @import above: relying on load
# order broke a real install (works under load_all, fails after
# R CMD INSTALL). Do not remove it.
# ============================================================

.toSFEWithKNNNeighbors <- function(se, k = 6L) {
  sfe <- as(se, "SpatialFeatureExperiment")
  SpatialFeatureExperiment::colGraph(sfe, "knn") <-
    SpatialFeatureExperiment::findSpatialNeighbors(
      sfe, type = "spatialCoords", method = "knearneigh", k = as.integer(k),
      style = "W", zero.policy = TRUE)
  sfe
}

.joincountDF <- function(sfe, group_col) {
  jc <- spdep::joincount.multi(factor(colData(sfe)[[group_col]]),
                               SpatialFeatureExperiment::colGraph(sfe, "knn"), zero.policy = TRUE)
  df <- as.data.frame(jc)
  colnames(df)[colnames(df) == "z-value"] <- "z_value"
  df$pair <- rownames(jc)
  df
}

.STATS_PLOTS <- c(
  "Spatial neighbour graph (k-NN)"  = "graph_plot",
  "Local Moran's I"                 = "localmoran_plot",
  "Moran scatter plot"              = "moran_plot",
  "Local cluster significance"      = "locclust_plot",
  "Join-count statistics"           = "joincount_plot",
  "Cell-type co-occurrence heatmap" = "joincount_heatmap_plot",
  "L-function (Lest)"               = "lest_plot",
  "L-cross function (Lcross)"       = "lcross_plot")

.computeSpatialStat <- function(se, plot_type, genes = character(0),
                                group_col = NA_character_, k = 6L) {
  has_groups <- !is.na(group_col) && group_col %in% colnames(colData(se)) &&
    !all(is.na(colData(se)[[group_col]]))
  if (has_groups) se <- se[, !is.na(colData(se)[[group_col]])]
  if (ncol(se) < max(5, k + 1))
    return(.notAvailablePlot(sprintf("Too few cells (%d) in this selection for k = %d", ncol(se), k)))
  need_groups <- c("lest_plot", "lcross_plot", "locclust_plot", "joincount_plot", "joincount_heatmap_plot")
  if (plot_type %in% need_groups && !has_groups)
    return(.notAvailablePlot("No cell-type column chosen / available for this statistic"))

  # L-functions work on point patterns; no neighbour graph involved.
  if (plot_type %in% c("lest_plot", "lcross_plot")) {
    counts <- sort(table(droplevels(factor(colData(se)[[group_col]]))), decreasing = TRUE)
    top <- names(counts[counts >= 10])
    n_needed <- if (plot_type == "lcross_plot") 2 else 1
    if (length(top) < n_needed) return(.notAvailablePlot("Not enough cells per cell type in this selection"))
    res <- spatialFDA::calcMetricPerFov(
      spe = se, fun = if (plot_type == "lest_plot") "Lest" else "Lcross",
      subsetby = "sample_id", marks = group_col, selection = top[seq_len(n_needed)],
      rSeq = seq(0, 500, length.out = 100), by = "sample_id", ncores = 1)
    return(spatialFDA::plotMetricPerFov(res, theo = TRUE, correction = "iso", x = "r", imageId = "sample_id"))
  }

  sfe <- .toSFEWithKNNNeighbors(se, k)
  genes <- intersect(genes, rownames(sfe))
  if (plot_type %in% c("localmoran_plot", "moran_plot", "locclust_plot") && length(genes) == 0)
    return(.notAvailablePlot("Choose at least one gene present in the data"))

  switch(plot_type,
    graph_plot = Voyager::plotColGraph(sfe, colGraphName = "knn", colGeometryName = 1) + ggplot2::theme_void(),

    localmoran_plot = {
      sfe <- Voyager::runUnivariate(sfe, features = genes, type = "localmoran", zero.policy = TRUE,
                                    colGraphName = "knn", colGeometryName = 1)
      Voyager::plotLocalResult(sfe, name = "localmoran", features = genes, ncol = 2,
                               colGeometryName = 1, divergent = TRUE, diverge_center = 0)
    },

    moran_plot = {
      sfe <- Voyager::runUnivariate(sfe, feature = genes[1], zero.policy = TRUE, type = "moran.plot",
                                    colGraphName = "knn", colGeometryName = 1)
      Voyager::moranPlot(sfe, feature = genes[1], graphName = "knn")
    },

    locclust_plot = {
      sfe <- Voyager::runUnivariate(sfe, features = genes[1], type = "localmoran", zero.policy = TRUE,
                                    colGraphName = "knn", colGeometryName = 1)
      # localResults()/localResults<-() are SpatialFeatureExperiment's own
      # generic (Voyager imports but does not re-export it -- Voyager::
      # localResults looks fine at a glance since it resolves via
      # inheritance, but errors at call time: "not an exported object").
      lm <- SpatialFeatureExperiment::localResults(sfe)[["localmoran"]][[1]]
      lm$locClust <- with(lm, ifelse(`-log10p_adj` > -log10(0.05), as.character(mean), "non-significant"))
      SpatialFeatureExperiment::localResults(sfe)[["localmoran"]][[1]] <- lm
      Voyager::plotLocalResult(sfe, features = genes[1], name = "localmoran",
                               attribute = "locClust", colGeometryName = 1)
    },

    joincount_plot = {
      df <- .joincountDF(sfe, group_col)
      df <- df[order(-abs(df$z_value)), ]
      df$pair <- factor(df$pair, levels = rev(df$pair))
      ggplot2::ggplot(df, ggplot2::aes(x = pair, y = z_value, fill = z_value > 0)) +
        ggplot2::geom_col() + ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = c(`TRUE` = "#F7A11A", `FALSE` = "#4ECDC4"), guide = "none") +
        ggplot2::labs(x = NULL, y = "z-value", title = sprintf("Join-count statistics (%s, k = %d)", group_col, k)) +
        ggplot2::theme_minimal()
    },

    joincount_heatmap_plot = {
      df <- .joincountDF(sfe, group_col)
      lv <- levels(droplevels(factor(colData(sfe)[[group_col]])))
      z <- matrix(NA_real_, length(lv), length(lv), dimnames = list(lv, lv))
      parts <- strsplit(df$pair, ":")
      for (i in seq_along(parts)) {
        if (length(parts[[i]]) != 2) next          # skip the "Jtot" total row
        z[parts[[i]][1], parts[[i]][2]] <- z[parts[[i]][2], parts[[i]][1]] <- df$z_value[i]
      }
      zl <- max(abs(z), na.rm = TRUE)
      if (!is.finite(zl) || zl == 0) return(.notAvailablePlot("Not enough cell-type pairs in this selection"))
      pheatmap::pheatmap(z, cluster_rows = FALSE, cluster_cols = FALSE,
        color = grDevices::colorRampPalette(c("#4ECDC4", "white", "#F7A11A"))(101),
        breaks = seq(-zl, zl, length.out = 102), display_numbers = TRUE, number_format = "%.1f",
        na_col = "grey90", main = sprintf("Cell-type co-occurrence (join-count z, k = %d)", k), silent = TRUE)
    },

    .notAvailablePlot(paste("Unknown plot type:", plot_type)))
}

# Moran's I scan over the most expressed genes (or metadata(se)$gene_ranks).
.computeSVGRanking <- function(se, top_n = 200L, k = 6L) {
  empty <- data.frame(gene = character(0), moran_i = numeric(0))
  if (ncol(se) < max(5, k + 1)) return(empty)
  cand <- S4Vectors::metadata(se)$gene_ranks
  if (is.null(cand)) {
    a <- if ("counts" %in% SummarizedExperiment::assayNames(se)) "counts" else "logcounts"
    cand <- rownames(se)[order(Matrix::rowSums(SummarizedExperiment::assay(se, a)), decreasing = TRUE)]
  }
  cand <- intersect(utils::head(cand, top_n), rownames(se))
  if (length(cand) == 0) return(empty)
  sfe <- .toSFEWithKNNNeighbors(se, k)
  m <- Voyager::calculateMoransI(sfe, features = cand, colGraphName = "knn",
                                 colGeometryName = 1, zero.policy = TRUE)
  tbl <- data.frame(gene = rownames(m), moran_i = m$moran, row.names = NULL)
  tbl <- tbl[!is.na(tbl$moran_i), , drop = FALSE]
  tbl[order(-tbl$moran_i), , drop = FALSE]
}
