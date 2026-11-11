# ============================================================
# Shared helpers: plot fallbacks, cached rendering, input observers.
# ============================================================

# Placeholder shown instead of a plot: the message centred in the panel, in
# large grey text, so an empty panel reads as "waiting for you" rather than
# broken (F9). The message is also the plot title, for tests and exports.
.notAvailablePlot <- function(msg) {
  wrapped <- paste(strwrap(msg, width = 45), collapse = "\n")
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = wrapped, size = 5.5,
                      colour = "grey45", lineheight = 1.1) +
    ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
    ggplot2::labs(title = msg) + ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_blank())
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

# Update requests for iseeSpatial panels.
#
# iSEE's .requestUpdate() is not used for them: its update observer
# (.create_child_propagation_observer) then calls .generateOutput() on the
# panel itself, and the panel's render computes it again, so every Validate
# or Clear computed each panel twice, including the ones whose result did
# not change (niches, SVG ranking, statistics). Our panels do not send iSEE
# multi-selections, so they need none of that observer's work. Instead each
# panel has its own counter in rObjects, read only by its .panelKey().
.updateName <- function(enc) paste0(enc, "_iseeSpatialUpdate")

.spatialUpdate <- function(enc, rObjects) {
  nm <- .updateName(enc)
  n <- shiny::isolate(rObjects[[nm]])
  rObjects[[nm]] <- if (is.null(n)) 1L else n + 1L
  invisible(NULL)
}

# Reactive render key of one panel. It changes only when the panel was asked
# to update (our counter, or iSEE's rObjects$modified for changes made by
# iSEE itself) AND the slots in key_fun() differ; reactiveVal() ignores
# identical values, so outputs that read it re-run only for their own
# changes. (rObjects$modified is reassigned as a whole list on every update,
# so outputs must not read it directly.)
.panelKey <- function(enc, pObjects, rObjects, key_fun) {
  rv <- shiny::reactiveVal(key_fun(pObjects$memory[[enc]]))
  shiny::observe({
    rObjects[[.updateName(enc)]]
    rObjects$modified
    rv(key_fun(pObjects$memory[[enc]]))
  })
  rv
}

# Runs a panel's output code. An R error is shown in the panel as its
# message (validate(), grey text) and logged to the console, instead of
# Shiny's generic error output (seen as "[object Object]", F3). Shiny's own
# silent errors (req(), validate()) pass through unchanged. Errors are not
# cached by renderCachedPlot, so the panel retries on its next update.
.safeOutput <- function(expr, panel) {
  tryCatch(expr, error = function(e) {
    if (inherits(e, "shiny.silent.error")) stop(e)
    message("[iseeSpatial] ", panel, ": ", conditionMessage(e))
    shiny::validate(shiny::need(FALSE, paste("Error:", conditionMessage(e))))
  })
}

# Cached plot rendering shared by every plot panel; re-renders only when one
# of key_slots changes (see .panelKey).
.renderCachedSpatialPlot <- function(x, se, output, pObjects, rObjects, key_slots) {
  plot_name <- .getEncodedName(x)
  key <- .panelKey(plot_name, pObjects, rObjects,
                   function(pn) lapply(key_slots, function(s) methods::slot(pn, s)))
  fp <- .fingerprintOf(se)
  ver <- as.character(utils::packageVersion("iseeSpatial"))
  output[[plot_name]] <- shiny::renderCachedPlot({
    key()
    .safeOutput(.drawPlot(.generateOutput(pObjects$memory[[plot_name]], se)), plot_name)
  }, cacheKeyExpr = list(fp, ver, plot_name, key()), cache = .plotCache())
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
      if (redraw) .spatialUpdate(enc, rObjects)
    }
  }, ignoreInit = TRUE, ignoreNULL = !multiple)
}

# Compute button: stamps a trigger slot, which is part of the cache key.
.computeObserver <- function(enc, input, pObjects, rObjects) {
  shiny::observeEvent(input[[paste0(enc, "_Compute")]], {
    pObjects$memory[[enc]]@ComputeTrigger <- as.character(Sys.time())
    .spatialUpdate(enc, rObjects)
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
