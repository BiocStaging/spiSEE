# ============================================================
# Shared helpers: plot fallbacks, cached rendering, input observers.
# ============================================================

.notAvailablePlot <- function(msg) {
  ggplot2::ggplot() + ggplot2::labs(title = msg) + ggplot2::theme_void()
}

.FORCE_CATEGORICAL <- c("fov", "SingleR_label", "low_qcscore", "ex")

.forceCategorical <- function(se, colorby) {
  if (colorby %in% .FORCE_CATEGORICAL && colorby %in% colnames(colData(se)))
    colData(se)[[colorby]] <- factor(colData(se)[[colorby]])
  se
}

.plotFOVOverview <- function(se, point_size = 0) {
  coords <- as.data.frame(SpatialExperiment::spatialCoords(se))
  coords$fov <- as.integer(as.character(colData(se)[["fov"]]))
  xcol <- colnames(coords)[1]; ycol <- colnames(coords)[2]
  ggplot2::ggplot(coords, ggplot2::aes(x = .data[[xcol]], y = .data[[ycol]], color = fov)) +
    ggplot2::geom_point(size = point_size, alpha = 0.6) +
    ggplot2::scale_color_gradientn(colors = pals::jet(), name = "FOV") +
    ggplot2::coord_fixed() + ggplot2::theme_void() + ggplot2::ggtitle("FOV overview")
}

.discretePalette <- function(n) {
  if (n <= 32) pals::glasbey(32)[seq_len(n)] else grDevices::hcl.colors(n, "Dynamic")
}

# pheatmap / gtable results must be drawn explicitly inside renderPlot.
.drawPlot <- function(p) {
  if (inherits(p, "pheatmap") || inherits(p, "gtable")) {
    grid::grid.newpage(); grid::grid.draw(if (inherits(p, "pheatmap")) p$gtable else p)
    invisible(NULL)
  } else p
}

# Cached plot rendering shared by every plot panel.
#
# iSEE's rObjects$modified is reassigned as a whole list on every update, so
# ANY panel's .requestUpdate() invalidates every panel's render block. The
# cache key is built from this panel's own slots: an unrelated change is a
# cache hit (no recompute, no re-transmission, no flicker).
.renderCachedSpatialPlot <- function(x, se, output, pObjects, rObjects, key_slots) {
  plot_name <- .getEncodedName(x)
  output[[plot_name]] <- shiny::renderCachedPlot({
    force(rObjects$modified[[plot_name]])
    .drawPlot(.generateOutput(pObjects$memory[[plot_name]], se))
  }, cacheKeyExpr = {
    force(rObjects$modified[[plot_name]])
    pn <- pObjects$memory[[plot_name]]
    c(list(plot_name), lapply(key_slots, function(s) methods::slot(pn, s)))
  })
}

# One observer per simple input <-> slot pair.
#   cast:     coerce the input value to the slot type
#   redraw:   request a re-render (FALSE for parameters that only apply on Compute)
#   multiple: NULL input means "nothing selected" (multi-select) instead of "ignore"
.slotObserver <- function(enc, field, input, pObjects, rObjects,
                          cast = identity, redraw = TRUE, multiple = FALSE) {
  id <- paste0(enc, "_", field)
  shiny::observeEvent(input[[id]], {
    v <- input[[id]]
    if (is.null(v)) {
      if (!multiple) return()
      v <- character(0)
    }
    v <- cast(v)
    if (!multiple && (length(v) != 1 || is.na(v))) return()
    if (!identical(v, methods::slot(pObjects$memory[[enc]], field))) {
      methods::slot(pObjects$memory[[enc]], field) <- v
      if (redraw) .requestUpdate(enc, rObjects)
    }
  }, ignoreInit = TRUE, ignoreNULL = !multiple)
}

# Compute button: stamps a trigger slot, which is part of the cache key.
.computeObserver <- function(enc, input, pObjects, rObjects) {
  shiny::observeEvent(input[[paste0(enc, "_Compute")]], {
    pObjects$memory[[enc]]@ComputeTrigger <- as.character(Sys.time())
    .requestUpdate(enc, rObjects)
  }, ignoreInit = TRUE)
}

.computeButton <- function(enc, label = "Compute") {
  shiny::actionButton(paste0(enc, "_Compute"), label, class = "btn btn-primary btn-sm", width = "100%")
}

.selectionHelp <- function(se) {
  what <- if (.selectionUnit(se) == "fov") "FOVs" else "cells"
  shiny::helpText(sprintf("Region: %s selected on the overview map (press Validate there). Empty = whole dataset.", what))
}

# Gene menus can hold >18,000 entries (Visium HD): fill them server-side.
.serverSelectize <- function(session, enc, field, choices, selected) {
  shiny::updateSelectizeInput(session, paste0(enc, "_", field),
                              choices = choices, selected = selected, server = TRUE)
}
