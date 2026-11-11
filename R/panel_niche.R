#' @import iSEE
#' @importFrom shiny plotOutput selectInput numericInput sliderInput actionButton helpText observeEvent checkboxInput
NULL

utils::globalVariables(c("x", "y", "value", "niche", "prop", "celltype"))

# ============================================================
# SpatialNichePanel - cellular niches and spatial domains, for any platform.
#
# Follows the map selection through the `Selection` slot (SpatialPanelBase),
# on any platform.
#
# In the app nothing is computed until Compute is pressed:
#  - Before Compute: the full-data result is shown on the whole tissue
#    (stored by computeNiches()/computeSpatialDomains() offline, see
#    inst/scripts/precompute_niches.R, else computed once on all cells).
#  - Compute: runs on the map selection at that moment with the panel's
#    parameters (k, number of niches, lambda). Later selections or parameter
#    changes wait for the next Compute; switching the plot type reuses the
#    result.
# ============================================================

.NICHE_PLOTS <- c(
  "Niches on tissue"                       = "niche_map",
  "Niche x cell-type enrichment (heatmap)" = "niche_heatmap",
  "Niche composition (bar)"                = "niche_bar",
  "Neighbourhood entropy on tissue"        = "entropy_map",
  "Spatial domains on tissue"              = "domain_map",
  "Expression clusters on tissue"          = "cluster_map",
  "Domains vs expression clusters"         = "domain_vs_cluster"
)

.tissuePlot <- function(se, values, title, discrete = TRUE, use_polygons = FALSE) {
  has_poly <- "polygons" %in% colnames(colData(se))
  if (use_polygons && !has_poly)
    return(.notAvailablePlot("No cell polygons loaded for this dataset: untick 'Draw cell polygons'"))
  if (use_polygons) {
    # Real cell-shape polygons instead of a centroid scatter; same geom_sf
    # mechanism as SpatialPlotPanel's UsePolygons (see panel_plot.R), but
    # hand-rolled here rather than via SpaceTrooper::plotPolygons() since
    # this panel's "on tissue" plots are already a shared ggplot helper
    # used by several PlotTypes, not just one colourBy column.
    # colData(se)$polygons is an sf data frame; data.frame() would split it
    # into polygons.* columns, so build an sf object from its geometry.
    df <- sf::st_sf(value = values, geometry = sf::st_geometry(colData(se)[["polygons"]]))
    p <- ggplot2::ggplot(df, ggplot2::aes(fill = value)) +
      ggplot2::geom_sf(colour = NA) +
      ggplot2::theme_void() + ggplot2::labs(title = title, fill = NULL)
    if (discrete) {
      lv <- levels(factor(values))
      return(p + ggplot2::scale_fill_manual(values = stats::setNames(.discretePalette(length(lv)), lv),
                                            na.value = "grey85"))
    }
    return(p + ggplot2::scale_fill_viridis_c(na.value = "grey85"))
  }
  xy <- as.data.frame(SpatialExperiment::spatialCoords(se))
  colnames(xy) <- c("x", "y")
  xy$value <- values
  pt <- if (nrow(xy) > 50000) list(shape = ".") else list(size = 0.6)
  p <- ggplot2::ggplot(xy, ggplot2::aes(x = x, y = y, colour = value)) +
    do.call(ggplot2::geom_point, pt) +
    ggplot2::coord_fixed() + ggplot2::theme_void() +
    ggplot2::labs(title = title, colour = NULL)
  if (discrete) {
    lv <- levels(factor(values))
    p + ggplot2::scale_colour_manual(values = stats::setNames(.discretePalette(length(lv)), lv),
                                     na.value = "grey85") +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3, shape = 16)))
  } else {
    p + ggplot2::scale_colour_viridis_c(na.value = "grey85")
  }
}

.enrichmentHeatmap <- function(rows, cols, title) {
  keep <- !is.na(rows) & !is.na(cols)
  enr <- .log2Enrichment(as.character(rows[keep]), as.character(cols[keep]))
  zl <- max(abs(enr), 1e-6)
  pheatmap::pheatmap(enr, cluster_rows = nrow(enr) > 2, cluster_cols = ncol(enr) > 2,
    color = grDevices::colorRampPalette(c("#4ECDC4", "white", "#F7A11A"))(101),
    breaks = seq(-zl, zl, length.out = 102), main = title, silent = TRUE)
}

setMethod("initialize", "SpatialNichePanel", function(.Object, ...) {
  .Object <- callNextMethod(.Object, ...)
  if (length(.Object@PlotType)       == 0) .Object@PlotType       <- "niche_map"
  if (length(.Object@LabelCol)       == 0) .Object@LabelCol       <- NA_character_
  if (length(.Object@K)              == 0) .Object@K              <- 20L
  if (length(.Object@NNiches)        == 0) .Object@NNiches        <- 6L
  if (length(.Object@Lambda)         == 0) .Object@Lambda         <- 0.2
  if (length(.Object@ComputeTrigger) == 0) .Object@ComputeTrigger <- ""
  # NA = resolved by .refineParameters() against the dataset, same
  # pattern as SpatialPlotPanel's UsePolygons (see panel_plot.R).
  if (length(.Object@UsePolygons)    == 0) .Object@UsePolygons    <- NA
  .Object
})

#' @rdname SpatialNichePanel
#' @export
setMethod(".refineParameters", "SpatialNichePanel", function(x, se) {
  x <- callNextMethod()
  if (!is.null(x) && isTRUE(is.na(x@UsePolygons)))
    x@UsePolygons <- "polygons" %in% colnames(colData(se)) && ncol(se) <= .POLYGON_AUTO_MAX
  x
})

#' @rdname SpatialNichePanel
#' @export
setMethod(".fullName",   "SpatialNichePanel", function(x) "Spatial Niches & Domains")
#' @rdname SpatialNichePanel
#' @export
setMethod(".panelColor", "SpatialNichePanel", function(x) "#7B2CBF")
#' @rdname SpatialNichePanel
#' @export
setMethod(".defineOutput", "SpatialNichePanel", function(x)
  shiny::plotOutput(.getEncodedName(x), height = "600px"))

#' @rdname SpatialNichePanel
#' @export
setMethod(".defineInterface", "SpatialNichePanel", function(x, se, select_info) {
  enc <- .getEncodedName(x)
  labs <- setdiff(.labelColumns(se), c("niche", "spatial_domain", "expr_cluster"))
  sel_lab <- if (!is.na(x@LabelCol) && x@LabelCol %in% labs) x@LabelCol else .defaultLabelCol(se)
  stored <- c(niches = "niche" %in% colnames(colData(se)),
              domains = "spatial_domain" %in% colnames(colData(se)))
  status <- sprintf("Precomputed on full data: niches %s, domains %s.",
                    if (stored["niches"]) "yes" else "no",
                    if (stored["domains"]) "yes" else "no")
  list(iSEE::collapseBox(
    id = paste0(enc, "_DataBoxOpen"), title = "Parameters", open = TRUE,
    shiny::selectInput(paste0(enc, "_PlotType"), "Show:", choices = .NICHE_PLOTS,
                       selected = x@PlotType),
    shiny::selectInput(paste0(enc, "_LabelCol"), "Cell-type column:",
                       choices = labs, selected = sel_lab),
    shiny::checkboxInput(paste0(enc, "_UsePolygons"),
                         "Draw cell polygons on tissue plots (slow on large regions)",
                         value = isTRUE(x@UsePolygons)),
    shiny::numericInput(paste0(enc, "_K"), "Neighbours per cell (k):",
                        value = x@K, min = 2, max = 200, step = 1),
    shiny::numericInput(paste0(enc, "_NNiches"), "Number of niches:",
                        value = x@NNiches, min = 2, max = 30, step = 1),
    shiny::sliderInput(paste0(enc, "_Lambda"), "Spatial weight for domains (lambda):",
                       min = 0, max = 0.9, value = x@Lambda, step = 0.05),
    shiny::helpText(status),
    shiny::actionButton(paste0(enc, "_Compute"), "Compute on current selection",
                        class = "btn btn-primary btn-sm", width = "100%"),
    shiny::helpText("Until Compute is pressed, the full-data result is shown on the",
      "whole tissue. Compute runs on the current map selection; a new selection",
      "or new parameters wait for the next Compute.")
  ))
})

# Run (or reuse) the analysis needed by the requested plot.
.nicheData <- function(x, se) {
  se <- .subsetBySelection(se, x@Selection)
  lab <- if (!is.na(x@LabelCol)) x@LabelCol else .defaultLabelCol(se)
  wants_domains <- x@PlotType %in% c("domain_map", "cluster_map", "domain_vs_cluster")
  recompute <- nchar(x@ComputeTrigger) > 0
  if (wants_domains) {
    md <- S4Vectors::metadata(se)$iseeSpatial_domains
    stored_ok <- "spatial_domain" %in% colnames(colData(se)) && !is.null(md) &&
      identical(as.integer(md$k_spatial), x@K) && isTRUE(all.equal(md$lambda, x@Lambda))
    if (recompute || !stored_ok)
      se <- computeSpatialDomains(se, k_spatial = x@K, lambda = x@Lambda)
  } else {
    # The stored (offline, full-data) result is reused only if it was made
    # with the parameters currently shown in the panel.
    md <- S4Vectors::metadata(se)$iseeSpatial_niches
    stored_ok <- "niche" %in% colnames(colData(se)) && !is.null(md) &&
      identical(md$label_col, lab) && identical(as.integer(md$k), x@K) &&
      identical(as.integer(md$n_niches), x@NNiches)
    if (recompute || !stored_ok)
      se <- computeNiches(se, label_col = lab, k = x@K, n_niches = x@NNiches)
  }
  list(se = se, lab = lab)
}

# Checks, analysis (or a cached one) and plot. `cache` holds the results of
# one Compute, per analysis kind, so switching the plot type does not re-run
# it; `note` is added to the plot to say what the result was computed on.
.nicheOutput <- function(x, se, cache = new.env(), note = NULL) {
  if (ncol(.subsetBySelection(se, x@Selection)) < 50)
    return(.notAvailablePlot("Too few cells in this selection (need >= 50)"))
  wants_domains <- x@PlotType %in% c("domain_map", "cluster_map", "domain_vs_cluster")
  if (is.na(x@LabelCol) && is.na(.defaultLabelCol(se)) && !wants_domains)
    return(.notAvailablePlot("No categorical cell-type column found in colData"))
  kind <- if (wants_domains) "domains" else "niches"
  d <- cache[[kind]]
  if (is.null(d)) {
    d <- tryCatch(.nicheData(x, se), error = function(e) e)
    if (inherits(d, "error")) return(.notAvailablePlot(conditionMessage(d)))
    cache[[kind]] <- d
  }
  .nichePlot(x, d, note)
}

#' @rdname SpatialNichePanel
#' @export
setMethod(".generateOutput", "SpatialNichePanel", function(x, se, ...) .nicheOutput(x, se))

.nichePlot <- function(x, d, note = NULL) {
  s <- d$se; cd <- colData(s)
  up <- isTRUE(x@UsePolygons)
  # heatmaps are not ggplots: the note goes under the title
  ht <- function(title) if (is.null(note)) title else paste0(title, "\n", note)
  p <- switch(x@PlotType,
    niche_map     = .tissuePlot(s, cd$niche, sprintf("Niches (k = %d, by %s)", x@K, d$lab), use_polygons = up),
    entropy_map   = .tissuePlot(s, cd$niche_entropy, "Neighbourhood entropy (bits)", discrete = FALSE, use_polygons = up),
    niche_heatmap = .enrichmentHeatmap(cd$niche, cd[[d$lab]], ht("Cell types per niche (log2 obs/exp)")),
    niche_bar     = {
      comp <- SingleCellExperiment::reducedDim(s, "NicheComposition")
      ok <- !is.na(cd$niche)
      m <- stats::aggregate(comp[ok, , drop = FALSE], list(niche = cd$niche[ok]), mean)
      long <- stats::reshape(m, direction = "long", varying = colnames(comp), v.names = "prop",
                             timevar = "celltype", times = colnames(comp), idvar = "niche")
      ggplot2::ggplot(long, ggplot2::aes(x = niche, y = prop, fill = celltype)) +
        ggplot2::geom_col() + ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = stats::setNames(.discretePalette(ncol(comp)), colnames(comp))) +
        ggplot2::labs(x = NULL, y = "Mean neighbour proportion", fill = NULL,
                      title = "Average neighbourhood per niche") +
        ggplot2::theme_minimal()
    },
    domain_map    = .tissuePlot(s, cd$spatial_domain, sprintf("Spatial domains (lambda = %.2f)", x@Lambda), use_polygons = up),
    cluster_map   = .tissuePlot(s, cd$expr_cluster, "Expression-only clusters (lambda = 0)", use_polygons = up),
    domain_vs_cluster = .enrichmentHeatmap(cd$spatial_domain, cd$expr_cluster,
                                           ht("Domains vs expression clusters (log2 obs/exp)")),
    .notAvailablePlot(paste("Unknown plot type:", x@PlotType)))
  if (!is.null(note) && inherits(p, "ggplot")) p <- p + ggplot2::labs(caption = note)
  p
}

#' @rdname SpatialNichePanel
#' @export
setMethod(".renderOutput", "SpatialNichePanel", function(x, se, output, pObjects, rObjects) {
  plot_name <- .getEncodedName(x)
  # Selection and parameters are not in the key: they apply on Compute only.
  key <- .panelKey(plot_name, pObjects, rObjects,
                   function(pn) list(pn@PlotType, pn@UsePolygons, pn@ComputeTrigger))
  last <- new.env()   # this session's last Compute: panel snapshot + results
  output[[plot_name]] <- shiny::renderPlot({
    key()
    x <- pObjects$memory[[plot_name]]
    if (!identical(last$trigger, x@ComputeTrigger)) {
      snap <- x
      if (!nzchar(x@ComputeTrigger)) snap@Selection <- ""
      last$trigger <- x@ComputeTrigger; last$snap <- snap; last$cache <- new.env()
    }
    snap <- last$snap
    snap@PlotType <- x@PlotType; snap@UsePolygons <- x@UsePolygons
    n <- ncol(.subsetBySelection(se, snap@Selection))
    note <- if (!nzchar(snap@ComputeTrigger))
      "All cells (full-data result). Press Compute to run on the map selection."
    else sprintf("Computed on %d cells (map selection at last Compute).", n)
    .drawPlot(.nicheOutput(snap, se, last$cache, note))
  })
})

#' @rdname SpatialNichePanel
#' @export
setMethod(".createObservers", "SpatialNichePanel", function(x, se, input, session, pObjects, rObjects) {
  callNextMethod(); enc <- .getEncodedName(x)
  # PlotType redraws immediately (cheap when results are stored); parameters
  # only take effect on Compute, like the other stats panels.
  shiny::observeEvent(input[[paste0(enc, "_PlotType")]], {
    v <- input[[paste0(enc, "_PlotType")]]
    if (!identical(v, pObjects$memory[[enc]]@PlotType)) {
      pObjects$memory[[enc]]@PlotType <- v; .spatialUpdate(enc, rObjects) }
  }, ignoreInit = TRUE)
  .slotObserver(enc, "UsePolygons", input, pObjects, rObjects, cast = as.logical)
  setters <- list(LabelCol = as.character, K = as.integer, NNiches = as.integer, Lambda = as.numeric)
  for (nm in names(setters)) local({
    field <- nm; cast <- setters[[nm]]
    shiny::observeEvent(input[[paste0(enc, "_", field)]], {
      v <- cast(input[[paste0(enc, "_", field)]])
      if (length(v) == 1 && !is.na(v)) methods::slot(pObjects$memory[[enc]], field) <- v
    }, ignoreInit = TRUE)
  })
  shiny::observeEvent(input[[paste0(enc, "_Compute")]], {
    pObjects$memory[[enc]]@ComputeTrigger <- as.character(Sys.time())
    .spatialUpdate(enc, rObjects)
  }, ignoreInit = TRUE)
})

#' Spatial niches and domains panel
#'
#' iSEE panel showing cellular niches (from the cell-type composition of
#' each cell's spatial neighbours) and spatial domains (clustering on
#' neighbourhood-augmented expression), on any platform. It follows the
#' selection of the [SpatialMapPanel()]. Results precomputed offline with
#' [computeNiches()] / [computeSpatialDomains()] are shown instantly when the
#' panel's parameters match them; *Compute* re-runs on the selection.
#'
#' @param PlotType One of the values of the "Show" menu, e.g. `"niche_map"`.
#' @param LabelCol `colData` column with cell-type labels (auto-detected if
#'   `NA`).
#' @param K,NNiches,Lambda Neighbours per cell, number of niches, spatial
#'   weight for domains.
#' @param UsePolygons Draw real cell-shape polygons instead of a centroid
#'   scatter on the "on tissue" plot types (needs a `polygons` column, as
#'   created by SpaceTrooper). `NA` (the default) follows the dataset: on
#'   at load time if a `polygons` column exists, off otherwise.
#' @param ... Other slots, e.g. `Selection`.
#' @return A `SpatialNichePanel` object.
#' @examples
#' SpatialNichePanel(K = 30, NNiches = 8)
#' @export
SpatialNichePanel <- function(PlotType = "niche_map", LabelCol = NA_character_,
                              K = 20L, NNiches = 6L, Lambda = 0.2, UsePolygons = NA, ...)
  new("SpatialNichePanel", PlotType = PlotType, LabelCol = LabelCol,
      K = as.integer(K), NNiches = as.integer(NNiches), Lambda = Lambda,
      UsePolygons = UsePolygons, ...)
