# ============================================================
# SpatialMapPanel - overview map where the user builds a selection.
#
# One class for every platform; what is selectable comes from the platform
# (platform.R):
#   - "fov"  (CosMx): FOV tiles drawn as plotly rectangles on a cached image
#                      of the tissue; clicks toggle tiles.
#   - "cell" (Xenium, Visium HD): cell centroids as an invisible WebGL point
#                      trace on a cached image; clicks/box-select pick cells.
#
# Selection is entirely client-side (JS attached with htmlwidgets::onRender):
# clicking, box-select, "All" and highlight changes never touch R. Only
# Validate talks to the server: it writes the pending set into the
# `Selection` slot of every panel that has one (.findSelectablePanels) and
# requests their update. Validate REPLACES the previous selection; Clear
# resets it to the whole dataset.
#
# Both background images are rendered once and cached on disk
# (tools::R_user_dir), keyed by a digest of the data's shape -- rasterising
# ~50k points costs ~1 s, which used to be paid on every click.
# ============================================================

setMethod("initialize", "SpatialMapPanel", function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    if (length(.Object@PendingIDs)       == 0) .Object@PendingIDs       <- ""
    if (length(.Object@PendingRemoveIDs) == 0) .Object@PendingRemoveIDs <- ""
    if (length(.Object@CommittedIDs)     == 0) .Object@CommittedIDs     <- ""
    .Object
})

#' @rdname SpatialMapPanel
#' @export
setMethod(".fullName",   "SpatialMapPanel", function(x) "Spatial Overview Map")
#' @rdname SpatialMapPanel
#' @export
setMethod(".panelColor", "SpatialMapPanel", function(x) "#118AB2")
#' @rdname SpatialMapPanel
#' @export
setMethod(".defineOutput", "SpatialMapPanel", function(x)
    .withBusy(plotly::plotlyOutput(.getEncodedName(x), height = "600px")))

#' @rdname SpatialMapPanel
#' @export
setMethod(".defineInterface", "SpatialMapPanel", function(x, se, select_info) {
    enc  <- .getEncodedName(x)
    what <- if (.selectionUnit(se) == "fov") "FOVs" else "cells"
    list(iSEE::collapseBox(
        id = paste0(enc, "_DataBoxOpen"), title = "Selection", open = TRUE,
        shiny::tags$p(
            shiny::tags$b(sprintf("Click %s", what)),
            " or drag a rectangle to build a selection. ",
            shiny::tags$span(style = "color:#3A86FF; font-weight:bold;",
                             "Blue"), " = pending, ",
            shiny::tags$span(style = "color:#FF6B35; font-weight:bold;",
                             "orange"), " = applied, ",
            shiny::tags$span(style = "color:#E63946; font-weight:bold;",
                             "red"), " = pending removal. ",
            shiny::tags$b("Validate"),
            " applies the pending set to all panels (it replaces the previous selection); ",
            shiny::tags$b("Clear"), " goes back to the whole dataset.",
            style = "font-size:12px; color:#555; margin-bottom:8px;"),
        shiny::fluidRow(
            shiny::column(4,
                          shiny::actionButton(paste0(enc, "_Validate"),
                                              "Validate",
                class = "btn btn-primary btn-sm", width = "100%")),
            shiny::column(4,
                          shiny::actionButton(paste0(enc, "_ClearPending"),
                                              "Clear",
                class = "btn btn-default btn-sm", width = "100%")),
            shiny::column(4,
                          shiny::actionButton(paste0(enc, "_SelectAll"),
                                              paste("All", what),
                class = "btn btn-default btn-sm", width = "100%")))
    ))
})

# ── Background images (built once, cached on disk) ──────────────────────

.mapCacheKey <- function(se, unit) {
    xy <- SpatialExperiment::spatialCoords(se)
    digest::digest(list(unit, dim(se), range(xy[, 1]), range(xy[, 2]),
        utils::head(colnames(se), 50)), algo = "xxhash32")
}

.cachedMapBase <- function(se) {
    unit <- .selectionUnit(se)
    cache_dir <- tools::R_user_dir("iseeSpatial", "cache")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    f <- file.path(cache_dir,
                   sprintf("map_%s_%s.rds", unit, .mapCacheKey(se, unit)))
    if (file.exists(f)) return(readRDS(f))
    base <- if (unit == "fov") .buildTileMapBase(se) else .buildPointMapBase(se)
    base$unit <- unit
    saveRDS(base, f)
    base
}

.savePNG <- function(p, xr, yr) {
    # image aspect must match the data's,
    asp <- diff(yr) / diff(xr)
    # or overlaid shapes drift (letterboxing)
    tmp <- tempfile(fileext = ".png")
    ggplot2::ggsave(tmp, p, width = 7, height = 7 * asp, dpi = 150,
                    bg = "white")
    on.exit(unlink(tmp))
    readBin(tmp, "raw", file.info(tmp)$size)
}

# Tiles: reuse SpaceTrooper's FOV layout. Only plain data frames and PNG
# bytes are cached -- never a ggplot object, whose environments would drag
# the whole dataset into the cache file.
.buildTileMapBase <- function(se) {
    built <- ggplot2::ggplot_build(do.call(SpaceTrooper::plotCellsFovs,
        c(list(se), .scaleBarArgs(SpaceTrooper::plotCellsFovs, se))))
    points <- built$data[[1]]
    rects <- built$data[[2]]
    labs <- built$data[[3]]
    xr <- built$layout$panel_params[[1]]$x.range
    yr <- built$layout$panel_params[[1]]$y.range
    pt <- points[1, ]
    r1 <- rects[1, ]
    l1 <- labs[1, ]
    bg <- ggplot2::ggplot() +
        ggplot2::geom_point(data = points, mapping = ggplot2::aes(x = x, y = y),
            colour = pt$colour, alpha = pt$alpha, size = pt$size,
            shape = pt$shape, stroke = pt$stroke) +
        ggplot2::geom_rect(data = rects,
            mapping = ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin,
                                   ymax = ymax),
            colour = r1$colour, fill = NA, linewidth = r1$linewidth,
            linetype = r1$linetype, alpha = r1$alpha) +
        ggplot2::geom_text(data = labs,
                           mapping = ggplot2::aes(x = x, y = y, label = label),
            colour = l1$colour, size = l1$size, alpha = l1$alpha,
            fontface = l1$fontface,
            family = l1$family, hjust = l1$hjust, vjust = l1$vjust) +
        ggplot2::coord_fixed(xlim = xr, ylim = yr,
                             expand = FALSE) + ggplot2::theme_void()
    list(png_raw = .savePNG(bg, xr, yr), xr = xr, yr = yr, rects = rects,
        ids = as.character(labs$label), px = labs$x, py = labs$y)
}

.buildPointMapBase <- function(se) {
    xy <- SpatialExperiment::spatialCoords(se)
    xr <- range(xy[, 1])
    yr <- range(xy[, 2])
    xr <- xr + c(-1, 1) * diff(xr) * 0.05
    yr <- yr + c(-1, 1) * diff(yr) * 0.05
    df <- data.frame(x = xy[, 1], y = xy[, 2])
    size <- if (nrow(df) > 1e5) 0.3 else if (nrow(df) > 2e4) 0.6 else 1.2
    bg <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_point(size = size, alpha = 0.7, colour = "grey35") +
        ggplot2::coord_fixed(xlim = xr, ylim = yr,
                             expand = FALSE) + ggplot2::theme_void()
    list(png_raw = .savePNG(bg, xr, yr), xr = xr, yr = yr,
        ids = colnames(se), px = unname(xy[, 1]), py = unname(xy[, 2]))
}

# ── Widget ──────────────────────────────────────────────────────────────

# Client-side selection logic shared by both modes. Only `redraw()` differs:
# tiles restyle rectangle outlines, cells restyle a highlight trace.
.MAP_JS <- "
function(el, x, data) {
  var enc = data.enc, ids = data.ids, unit = data.unit;
  var pending = new Set(), pendingRemove = new Set(), committed = new Set();
  var idToXY = {};
  if (unit === 'cell') ids.forEach(function(id, i) { idToXY[id] = [data.px[i], data.py[i]]; });

  function tileStyle(id) {
    if (pending.has(id))       return {color: '#3A86FF', dash: 'dash',  width: 2};
    if (pendingRemove.has(id)) return {color: '#E63946', dash: 'dash',  width: 2};
    if (committed.has(id))     return {color: '#FF6B35', dash: 'solid', width: 1.5};
    return {color: '#000000', dash: 'solid', width: 1};
  }
  function redraw() {
    if (unit === 'fov') {
      var upd = {};
      ids.forEach(function(id, i) {
        var s = tileStyle(id);
        upd['shapes[' + i + '].line.color'] = s.color;
        upd['shapes[' + i + '].line.dash']  = s.dash;
        upd['shapes[' + i + '].line.width'] = s.width;
      });
      Plotly.relayout(el, upd);
    } else {
      var xs = [], ys = [], cols = [], cd = [];
      function push(id, c) { var p = idToXY[id]; xs.push(p[0]); ys.push(p[1]); cols.push(c); cd.push(id); }
      committed.forEach(function(id) { if (!pendingRemove.has(id)) push(id, '#FF6B35'); });
      pendingRemove.forEach(function(id) { push(id, '#E63946'); });
      pending.forEach(function(id) { push(id, '#3A86FF'); });
      Plotly.restyle(el, {x: [xs], y: [ys], 'marker.color': [cols], customdata: [cd]}, [1]);
    }
  }
  function sync() {
    Shiny.setInputValue(enc + '_pendingSync', Array.from(pending), {priority: 'event'});
    Shiny.setInputValue(enc + '_pendingRemoveSync', Array.from(pendingRemove), {priority: 'event'});
  }
  // A click toggles the id in whichever set it belongs to, so a second click undoes the first.
  function toggle(id) {
    if (committed.has(id)) { if (pendingRemove.has(id)) pendingRemove.delete(id); else pendingRemove.add(id); }
    else                   { if (pending.has(id)) pending.delete(id); else pending.add(id); }
    redraw(); sync();
  }
  if (unit === 'cell') redraw();   // clears the seeded placeholder point
  el.on('plotly_click', function(evt) {
    var pt = evt.points && evt.points[0];
    if (!pt || pt.customdata === undefined || pt.customdata === '__placeholder__') return;
    toggle(pt.customdata);
  });
  el.on('plotly_selected', function(evt) {          // box select: always additive
    if (!evt || !evt.points || evt.points.length === 0) return;
    evt.points.forEach(function(pt) {
      if (pt.customdata !== undefined && pt.customdata !== '__placeholder__' && !committed.has(pt.customdata)) pending.add(pt.customdata);
    });
    Plotly.relayout(el, {selections: []});           // drop plotly's own selection rectangle
    redraw(); sync();
  });
  Shiny.addCustomMessageHandler(enc + '_selectAll', function(msg) {
    pendingRemove.clear();
    ids.forEach(function(id) { if (!committed.has(id)) pending.add(id); });
    redraw(); sync();
  });
  // The server sends back what was applied (empty list on Clear): adopt it.
  Shiny.addCustomMessageHandler(enc + '_markCommitted', function(applied) {
    committed = new Set(applied); pending.clear(); pendingRemove.clear(); redraw();
  });
}"

.buildMapWidget <- function(base, enc) {
    img <- list(list(
        source = paste0("data:image/png;base64,",
                        base64enc::base64encode(base$png_raw)),
        xref = "x", yref = "y", x = base$xr[1], y = base$yr[2],
        sizex = diff(base$xr), sizey = diff(base$yr), sizing = "stretch",
        layer = "below"))
    axes <- list(xaxis = list(range = base$xr, visible = FALSE,
                              scaleanchor = "y", constrain = "domain"),
        yaxis = list(range = base$yr, visible = FALSE))
    if (base$unit == "fov") {
        shapes <- lapply(seq_along(base$ids), function(i) list(
            type = "rect", x0 = base$rects$xmin[i], x1 = base$rects$xmax[i],
            y0 = base$rects$ymin[i], y1 = base$rects$ymax[i],
            line = list(color = "black", width = 1, dash = "solid"),
            fillcolor = "rgba(0,0,0,0)"))
        p <- plotly::plot_ly(x = base$px, y = base$py, type = "scatter",
                             mode = "markers",
            marker = list(opacity = 0, size = 40), customdata = base$ids,
            hoverinfo = "text", text = base$ids, source = enc)
        p <- plotly::layout(p, images = img, shapes = shapes,
                            xaxis = axes$xaxis, yaxis = axes$yaxis,
            margin = list(l = 0, r = 0, t = 0, b = 0), dragmode = "select")
    } else {
        p <- plotly::plot_ly(x = base$px, y = base$py, type = "scattergl",
                             mode = "markers",
            marker = list(opacity = 0, size = 16), customdata = base$ids,
            hoverinfo = "skip", source = enc)
        # Trace 1 holds the highlighted cells; seeded with one point, cleared
        # on mount.
        p <- plotly::add_trace(p, data = NULL, inherit = FALSE, x = base$px[1],
                               y = base$py[1],
            type = "scattergl", mode = "markers",
            marker = list(size = 7, color = "#3A86FF"), hoverinfo = "skip",
            showlegend = FALSE, customdata = "__placeholder__")
        p <- plotly::layout(p, images = img, xaxis = axes$xaxis,
                            yaxis = axes$yaxis,
            margin = list(l = 0, r = 0, t = 0, b = 0), dragmode = "select")
    }
    p <- plotly::config(p, displayModeBar = FALSE, scrollZoom = FALSE)
    data <- list(enc = enc, unit = base$unit, ids = as.list(base$ids))
    if (base$unit == "cell") {
        data$px <- as.list(base$px)
        data$py <- as.list(base$py)
    }
    htmlwidgets::onRender(p, .MAP_JS, data = data)
}

# Static fallback used by iSEE's export/snapshot machinery only.
#' @rdname SpatialMapPanel
#' @export
setMethod(".generateOutput", "SpatialMapPanel", function(x, se, ...) {
    base <- .cachedMapBase(se)
    img <- png::readPNG(base$png_raw)
    ggplot2::ggplot() +
        ggplot2::annotation_raster(img, xmin = base$xr[1], xmax = base$xr[2],
            ymin = base$yr[1], ymax = base$yr[2]) +
        ggplot2::coord_fixed(xlim = base$xr, ylim = base$yr,
                             expand = FALSE) + ggplot2::theme_void()
})

#' @rdname SpatialMapPanel
#' @export
setMethod(".renderOutput", "SpatialMapPanel",
          function(x, se, output, pObjects, rObjects) {
    plot_name <- .getEncodedName(x)
    base <- .cachedMapBase(se)
    # Built once; deliberately no dependency on rObjects$modified: all
    # interaction is client-side, so it never needs to re-render.
    output[[plot_name]] <- plotly::renderPlotly(.buildMapWidget(base,
        plot_name))
})

#' @rdname SpatialMapPanel
#' @export
setMethod(".createObservers", "SpatialMapPanel",
          function(x, se, input, session, pObjects, rObjects) {
    callNextMethod()
    enc <- .getEncodedName(x)
    all_ids <- .allRegionIDs(se)

    # Pending state mirrored from the browser: stored only, nothing re-renders.
    shiny::observeEvent(input[[paste0(enc, "_pendingSync")]],
        {
            pObjects$memory[[enc]]@PendingIDs <- paste(input[[paste0(enc,
                "_pendingSync")]], collapse = ",")
        },
        ignoreNULL = FALSE,
        ignoreInit = TRUE)
    shiny::observeEvent(input[[paste0(enc, "_pendingRemoveSync")]],
        {
            pObjects$memory[[enc]]@PendingRemoveIDs <- paste(input[[paste0(enc,
                "_pendingRemoveSync")]], collapse = ",")
        },
        ignoreNULL = FALSE,
        ignoreInit = TRUE)

    # `ids` must be the FULL new committed set, not just this round's pending
    # clicks -- it becomes both the panels' Selection and the browser's
    # `committed` set, so anything left out here is dropped everywhere.
    .apply <- function(ids) {
        sel <- paste(ids, collapse = ",")
        for (target in .findSelectablePanels(pObjects)) {
            pObjects$memory[[target]]@Selection <- sel
            .spatialUpdate(target, rObjects)
        }
        pObjects$memory[[enc]]@CommittedIDs <- sel
        pObjects$memory[[enc]]@PendingIDs <- ""
        pObjects$memory[[enc]]@PendingRemoveIDs <- ""
        session$sendCustomMessage(paste0(enc, "_markCommitted"), as.list(ids))
    }

    shiny::observeEvent(input[[paste0(enc, "_Validate")]],
        {
            # Previously committed tiles aren't re-sent by the browser on every
            # click
            # (only the deltas are, via PendingIDs/PendingRemoveIDs), so the new
            # committed set must be folded in from the panel's own CommittedIDs
            # --
            # otherwise a second Validate round wipes out everything from the
            # first.
            committed <- .splitIDs(pObjects$memory[[enc]]@CommittedIDs)
            pending   <- .splitIDs(pObjects$memory[[enc]]@PendingIDs)
            removed   <- .splitIDs(pObjects$memory[[enc]]@PendingRemoveIDs)
            .apply(setdiff(union(committed, pending), removed))
        },
        ignoreInit = TRUE)

    shiny::observeEvent(input[[paste0(enc, "_ClearPending")]],
    .apply(character(0)), ignoreInit = TRUE)

    shiny::observeEvent(input[[paste0(enc, "_SelectAll")]],
        {
            pObjects$memory[[enc]]@PendingIDs <- paste(all_ids, collapse = ",")
            session$sendCustomMessage(paste0(enc, "_selectAll"), list())
        },
        ignoreInit = TRUE)
})

#' Overview map for building a region selection
#'
#' Shows the whole tissue and lets the user pick FOVs (CosMx) or cells
#' (Xenium, Visium HD) by clicking or dragging a rectangle. Selection happens
#' in the browser; pressing *Validate* applies it to every iseeSpatial panel
#' in the app (all classes extending [SpatialPanelBase-class]).
#'
#' @param ... Slots passed to `new()`, e.g. `PanelWidth`.
#' @return A `SpatialMapPanel` object.
#' @examples
#' SpatialMapPanel()
#' @export
SpatialMapPanel <- function(...) new("SpatialMapPanel", ...)
