#' iseeSpatial: interactive exploration of imaging- and sequencing-based spatial transcriptomics
#'
#' iSEE panels for spatial data at cellular resolution (CosMx, Xenium,
#' Visium HD). Select a region on an overview map and every panel follows:
#' spatial plots, QC histograms, spatial statistics, spatially variable
#' genes, cellular niches and spatial domains.
#'
#' Start with [iseeSpatial()]; the panel classes can also be mixed with any
#' standard iSEE panel.
#'
#' @import methods
#' @import iSEE
#' @import SpatialExperiment
#' @importFrom SummarizedExperiment colData colData<- assay
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom ggplot2 .data
#' @importFrom DT dataTableOutput renderDataTable datatable formatRound
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c("x", "y", "value", "niche", "prop", "celltype", "fov",
    "name", "pair", "z_value", "xmin", "xmax", "ymin", "ymax", "label"))
