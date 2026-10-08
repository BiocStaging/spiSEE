#' @import iSEE
#' @importFrom shiny plotOutput selectInput numericInput sliderInput
#'   actionButton helpText observeEvent
NULL

utils::globalVariables(c("x", "y", "value", "niche", "prop", "celltype"))

# ============================================================
# SpatialNichePanel - cellular niches and spatial domains, for any platform.
#
# Follows the map selection through the `Selection` slot (SpatialPanelBase),
# on any platform.
#
# Two modes, chosen automatically:
#  - Stored: if computeNiches()/computeSpatialDomains() were run offline
#    (inst/scripts/precompute_niches.R), the stored colData columns are shown
#    instantly for the current selection.
#  - Recompute: pressing Compute re-runs the analysis on the current
#    selection with the parameters in the panel (k, number of niches, lambda).
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

.tissuePlot <- function(se, values, title, discrete = TRUE) {
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
  .Object
})

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
    shiny::numericInput(paste0(enc, "_K"), "Neighbours per cell (k):",
                        value = x@K, min = 2, max = 200, step = 1),
    shiny::numericInput(paste0(enc, "_NNiches"), "Number of niches:",
                        value = x@NNiches, min = 2, max = 30, step = 1),
    shiny::sliderInput(paste0(enc, "_Lambda"), "Spatial weight for domains (lambda):",
                       min = 0, max = 0.9, value = x@Lambda, step = 0.05),
    shiny::helpText(status),
    shiny::actionButton(paste0(enc, "_Compute"), "Compute on current selection",
                        class = "btn btn-primary btn-sm", width = "100%"),
    shiny::helpText("Without Compute, the precomputed full-data result is shown",
      "for the selected cells (instant). Compute re-runs on the selection only.")
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

#' @rdname SpatialNichePanel
#' @export
setMethod(".generateOutput", "SpatialNichePanel", function(x, se, ...) {
  if (ncol(.subsetBySelection(se, x@Selection)) < 50)
    return(.notAvailablePlot("Too few cells in this selection (need >= 50)"))
  if (is.na(x@LabelCol) && is.na(.defaultLabelCol(se)) && !x@PlotType %in% c("domain_map", "cluster_map", "domain_vs_cluster"))
    return(.notAvailablePlot("No categorical cell-type column found in colData"))
  d <- tryCatch(.nicheData(x, se), error = function(e) e)
  if (inherits(d, "error")) return(.notAvailablePlot(conditionMessage(d)))
  s <- d$se; cd <- colData(s)
  switch(x@PlotType,
    niche_map     = .tissuePlot(s, cd$niche, sprintf("Niches (k = %d, by %s)", x@K, d$lab)),
    entropy_map   = .tissuePlot(s, cd$niche_entropy, "Neighbourhood entropy (bits)", discrete = FALSE),
    niche_heatmap = .enrichmentHeatmap(cd$niche, cd[[d$lab]], "Cell types per niche (log2 obs/exp)"),
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
    domain_map    = .tissuePlot(s, cd$spatial_domain, sprintf("Spatial domains (lambda = %.2f)", x@Lambda)),
    cluster_map   = .tissuePlot(s, cd$expr_cluster, "Expression-only clusters (lambda = 0)"),
    domain_vs_cluster = .enrichmentHeatmap(cd$spatial_domain, cd$expr_cluster,
                                           "Domains vs expression clusters (log2 obs/exp)"),
    .notAvailablePlot(paste("Unknown plot type:", x@PlotType)))
})

#' @rdname SpatialNichePanel
#' @export
setMethod(".renderOutput", "SpatialNichePanel", function(x, se, output, pObjects, rObjects) {
  plot_name <- .getEncodedName(x)
  # Same renderCachedPlot pattern as the stats panels: rObjects$modified is
  # shared, so cache on this panel's own state to avoid re-rendering on
  # unrelated panel changes.
  output[[plot_name]] <- shiny::renderCachedPlot({
    force(rObjects$modified[[plot_name]])
    p <- .generateOutput(pObjects$memory[[plot_name]], se)
    if (inherits(p, "pheatmap") || inherits(p, "gtable")) {
      grid::grid.newpage(); grid::grid.draw(if (inherits(p, "pheatmap")) p$gtable else p)
      invisible(NULL)
    } else p
  }, cacheKeyExpr = {
    force(rObjects$modified[[plot_name]])
    pn <- pObjects$memory[[plot_name]]
    list(plot_name, pn@Selection, pn@PlotType, pn@LabelCol,
         pn@K, pn@NNiches, pn@Lambda, pn@ComputeTrigger)
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
      pObjects$memory[[enc]]@PlotType <- v; .requestUpdate(enc, rObjects) }
  }, ignoreInit = TRUE)
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
    .requestUpdate(enc, rObjects)
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
#' @param ... Other slots, e.g. `Selection`.
#' @return A `SpatialNichePanel` object.
#' @examples
#' SpatialNichePanel(K = 30, NNiches = 8)
#' @export
SpatialNichePanel <- function(PlotType = "niche_map", LabelCol = NA_character_,
                              K = 20L, NNiches = 6L, Lambda = 0.2, ...)
  new("SpatialNichePanel", PlotType = PlotType, LabelCol = LabelCol,
      K = as.integer(K), NNiches = as.integer(NNiches), Lambda = Lambda, ...)
