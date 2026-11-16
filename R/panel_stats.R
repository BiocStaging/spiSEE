# SpatialStatsPanel and SVGPanel: heavy computations, run on Compute only.

setMethod("initialize", "SpatialStatsPanel", function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    if (length(.Object@PlotType) == 0) .Object@PlotType <- "graph_plot"
    if (length(.Object@Genes) == 0) .Object@Genes <- character(0)
    if (length(.Object@LabelCol) == 0) .Object@LabelCol <- NA_character_
    if (length(.Object@K) == 0) .Object@K <- 6L
    if (length(.Object@Alpha) == 0) .Object@Alpha <- 0.05
    if (length(.Object@ComputeTrigger) == 0) .Object@ComputeTrigger <- ""
    .Object
})

# Number of neighbours: one whole number >= 2 (NA allowed where it means
# "from the stored result").
.badK <- function(k, na_ok = FALSE) {
    length(k) != 1 || (is.na(k) && !na_ok) || (!is.na(k) && k < 2)
}

setValidity("SpatialStatsPanel", function(object) {
    a <- object@Alpha
    msg <- c(
        if (length(a) != 1 || is.na(a) || a <= 0 || a >= 1)
            "'Alpha' must be a single number in (0, 1)",
        if (.badK(object@K)) "'K' must be a single whole number >= 2")
    if (length(msg)) msg else TRUE
})

setValidity("SVGPanel", function(object) {
    n <- object@TopN
    msg <- c(
        if (length(n) != 1 || is.na(n) || n < 1)
            "'TopN' must be a single whole number >= 1",
        if (.badK(object@K)) "'K' must be a single whole number >= 2")
    if (length(msg)) msg else TRUE
})

#' @rdname SpatialStatsPanel
#' @export
setMethod(".fullName",   "SpatialStatsPanel", function(x) "Spatial Statistics")
#' @rdname SpatialStatsPanel
#' @export
setMethod(".panelColor", "SpatialStatsPanel", function(x) "#4ECDC4")
#' @rdname SpatialStatsPanel
#' @export
setMethod(".defineOutput", "SpatialStatsPanel", function(x)
    .withBusy(shiny::plotOutput(.getEncodedName(x), height = "550px")))

.kInput <- function(enc, value, label = "Neighbours per cell (k):")
    shiny::numericInput(paste0(enc, "_K"), label, value = value, min = 2,
                        max = 200, step = 1)

#' @rdname SpatialStatsPanel
#' @export
setMethod(".defineInterface", "SpatialStatsPanel",
            function(x, se, select_info) {
    enc <- .getEncodedName(x)
    labs <- .labelColumns(se)
    lab <- .chosenLabel(x, se, labs)
    list(iSEE::collapseBox(
        id = paste0(enc, "_DataBoxOpen"), title = "Parameters", open = TRUE,
        shiny::selectInput(paste0(enc, "_PlotType"), "Statistic:",
                            choices = .STATS_PLOTS, selected = x@PlotType),
        shiny::selectizeInput(paste0(enc, "_Genes"),
                                "Gene(s) (Moran's I, local clusters):",
            choices = x@Genes, selected = x@Genes, multiple = TRUE),
        shiny::selectInput(paste0(enc, "_LabelCol"),
                            "Cell-type column (join-count, L-functions):",
            choices = labs, selected = lab),
        .kInput(enc, x@K),
        shiny::helpText("k sets the spatial scale: small k = immediate",
                        "neighbours, large k = wider region."),
        shiny::numericInput(paste0(enc, "_Alpha"),
                            "Significance level (local clusters, adjusted p):",
            value = x@Alpha, min = 0.001, max = 0.5, step = 0.01),
        .selectionHelp(se),
        .computeButton(enc),
        shiny::helpText("Nothing recomputes until you press Compute.")))
})

#' @rdname SpatialStatsPanel
#' @export
setMethod(".generateOutput", "SpatialStatsPanel", function(x, se, ...) {
    if (!nzchar(x@ComputeTrigger))
        return(.notAvailablePlot(
            "Choose a statistic and a region, then press Compute."))
    lab <- if (!is.na(x@LabelCol)) x@LabelCol else .defaultLabelCol(se)
    .computeSpatialStat(.subsetBySelection(se, x@Selection), x@PlotType,
        genes = x@Genes, group_col = lab, k = x@K, alpha = x@Alpha)
})

#' @rdname SpatialStatsPanel
#' @export
setMethod(".refineParameters", "SpatialStatsPanel", function(x, se) {
    x <- callNextMethod()
    if (is.null(x)) return(x)
    .keyComputeTrigger(x)
})

#' @rdname SpatialStatsPanel
#' @export
setMethod(".renderOutput", "SpatialStatsPanel",
            function(x, se, output, pObjects, rObjects)
    .renderCachedSpatialPlot(x, se, output, pObjects, rObjects,
                            "ComputeTrigger"))

#' @rdname SpatialStatsPanel
#' @export
setMethod(".createObservers", "SpatialStatsPanel",
            function(x, se, input, session, pObjects, rObjects) {
    callNextMethod()
    enc <- .getEncodedName(x)
    genes <- .geneChoices(se)
    sel <- intersect(x@Genes, genes)
    if (!length(sel)) sel <- genes[1]
    .serverSelectize(session, enc, "Genes", genes, sel)
    # Parameters only take effect on Compute (the cache key is ComputeTrigger).
    .slotObserver(enc, "PlotType", input, pObjects, rObjects,
                    cast = as.character, redraw = FALSE)
    .slotObserver(enc, "Genes", input, pObjects, rObjects, cast = as.character,
                    redraw = FALSE, multiple = TRUE)
    .slotObserver(enc, "LabelCol", input, pObjects, rObjects,
                    cast = as.character, redraw = FALSE)
    .slotObserver(enc, "K", input, pObjects, rObjects, cast = as.integer,
                    redraw = FALSE)
    .slotObserver(enc, "Alpha", input, pObjects, rObjects, redraw = FALSE,
        cast = function(v) {
            v <- as.numeric(v)
            if (length(v) == 1 && !is.na(v) && v > 0 && v < 1) v else NA_real_
        })
    .computeObserver(enc, input, pObjects, rObjects)
})

#' Spatial statistics on the selected region
#'
#' Neighbour graph, global and local Moran's I, local cluster significance,
#' join-count statistics and L-functions, computed on the region selected
#' on the [SpatialMapPanel()] when the user presses *Compute*.
#'
#' The neighbour graph is a k-nearest-neighbour graph on the cell
#' centroids of the region. Moran's I and the local clusters come from
#' Voyager, the join counts from spdep. The L-function (`lest_plot`) is
#' computed for the most frequent cell type and the L-cross function
#' (`lcross_plot`) for the two most frequent ones (types with at least 10
#' cells), with spatialFDA, up to a radius of a quarter of the shorter
#' side of the region.
#'
#' Changing a parameter does not recompute anything: results are computed
#' on *Compute* and cached by what they compute (statistic, parameters,
#' region and dataset), so the same analysis opens at once the next time.
#'
#' @param PlotType One of `"graph_plot"`, `"localmoran_plot"`,
#'   `"moran_plot"`, `"locclust_plot"`, `"joincount_plot"`,
#'   `"joincount_heatmap_plot"`, `"lest_plot"`, `"lcross_plot"`.
#' @param Genes Genes for the Moran-type statistics.
#' @param LabelCol Cell-type column (`NA` = auto-detect).
#' @param K Number of spatial neighbours per cell.
#' @param Alpha Significance level for the local cluster plot: a cell is
#'   called significant when its adjusted p-value (local Moran's I) is
#'   below `Alpha`. Shown on the plot.
#' @param ... Other slots: `Selection`, and `ComputeTrigger`, which is set
#'   when *Compute* is pressed. Give any non-empty value (e.g.
#'   `ComputeTrigger = "launch"`) to compute the panel when the app opens.
#' @return A `SpatialStatsPanel` object.
#' @seealso [SVGPanel()] to rank genes by Moran's I, [SpatialNichePanel()]
#'   for niches and domains.
#' @examples
#' SpatialStatsPanel(PlotType = "joincount_plot", K = 10)
#' @export
SpatialStatsPanel <- function(PlotType = "graph_plot", Genes = character(0),
                                LabelCol = NA_character_, K = 6L, Alpha = 0.05,
                                ...)
    new("SpatialStatsPanel", PlotType = PlotType, Genes = Genes,
        LabelCol = LabelCol,
        K = as.integer(K), Alpha = as.numeric(Alpha), ...)

# SVGPanel

setMethod("initialize", "SVGPanel", function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    if (length(.Object@TopN) == 0) .Object@TopN <- 200L
    if (length(.Object@K) == 0) .Object@K <- 6L
    if (length(.Object@LinkedPanel) == 0) .Object@LinkedPanel <- ""
    if (length(.Object@ComputeTrigger) == 0) .Object@ComputeTrigger <- ""
    .Object
})

# Last computed table per panel, so a row click can find its gene without
# recomputing. Kept on the session's own pObjects (a package-level cache
# would be shared between simultaneous sessions of one server).
.svgTables <- function(pObjects) {
    if (is.null(pObjects$spiSEE_svg_tables))
        pObjects$spiSEE_svg_tables <- new.env(parent = emptyenv())
    pObjects$spiSEE_svg_tables
}

#' @rdname SVGPanel
#' @export
setMethod(".fullName",   "SVGPanel", function(x) "Spatially Variable Genes")
#' @rdname SVGPanel
#' @export
setMethod(".panelColor", "SVGPanel", function(x) "#F7A11A")
#' @rdname SVGPanel
#' @export
setMethod(".defineOutput", "SVGPanel",
            function(x) .withBusy(DT::dataTableOutput(.getEncodedName(x))))

#' @rdname SVGPanel
#' @export
setMethod(".defineInterface", "SVGPanel", function(x, se, select_info) {
    enc <- .getEncodedName(x)
    list(iSEE::collapseBox(
        id = paste0(enc, "_DataBoxOpen"), title = "Parameters", open = TRUE,
        shiny::numericInput(paste0(enc, "_TopN"),
                            "Scan the top N most expressed genes:",
            value = x@TopN, min = 10, max = 5000, step = 10),
        .kInput(enc, x@K),
        .selectionHelp(se),
        .computeButton(enc),
        shiny::helpText("Click a row to show that gene in the linked",
                        "Spatial Plot panel.")))
})

#' @rdname SVGPanel
#' @export
setMethod(".generateOutput", "SVGPanel", function(x, se, ...) {
    if (!nzchar(x@ComputeTrigger))
        return(structure(data.frame(gene = character(0), moran_i = numeric(0)),
            reason = "Press Compute to rank genes in the current selection."))
    .computeSVGRanking(.subsetBySelection(se, x@Selection), top_n = x@TopN,
                        k = x@K)
})

#' @rdname SVGPanel
#' @export
setMethod(".refineParameters", "SVGPanel", function(x, se) {
    x <- callNextMethod()
    if (is.null(x)) return(x)
    .keyComputeTrigger(x)
})

#' @rdname SVGPanel
#' @export
setMethod(".renderOutput", "SVGPanel",
            function(x, se, output, pObjects, rObjects) {
    enc <- .getEncodedName(x)
    # Rebuilt only on Compute: re-sending the table resets its page, order and
    # selected row.
    key <- .panelKey(enc, pObjects, rObjects, function(pn) pn@ComputeTrigger)
    fp <- .fingerprintOf(se)
    ver <- as.character(utils::packageVersion("spiSEE"))
    output[[enc]] <- DT::renderDataTable({
        trig <- key()
        x <- .keyPanel(key)
        tbl <- .safeOutput(if (nzchar(trig))
            .cachedResult(digest::digest(list(fp, ver, "svg", trig),
                                        algo = "xxhash64"),
                            .generateOutput(x, se))
        else .generateOutput(x, se), enc)
        assign(enc, tbl, envir = .svgTables(pObjects))
        why <- attr(tbl, "reason")
        names(tbl) <- c("Gene", "Moran's I")
        DT::formatRound(DT::datatable(tbl, selection = "single",
                                        rownames = FALSE,
            options = list(pageLength = 15, order = list(list(1, "desc")),
                language = list(
                    emptyTable = if (is.null(why)) "No result." else why))),
        "Moran's I", 3)
    })
})

#' @rdname SVGPanel
#' @export
setMethod(".createObservers", "SVGPanel",
            function(x, se, input, session, pObjects, rObjects) {
    callNextMethod()
    enc <- .getEncodedName(x)
    .slotObserver(enc, "TopN", input, pObjects, rObjects, cast = as.integer,
                    redraw = FALSE)
    .slotObserver(enc, "K", input, pObjects, rObjects, cast = as.integer,
                    redraw = FALSE)
    .computeObserver(enc, input, pObjects, rObjects)
    shiny::observeEvent(input[[paste0(enc, "_rows_selected")]],
        {
            idx <- input[[paste0(enc, "_rows_selected")]]
            tbl <- get0(enc, envir = .svgTables(pObjects), inherits = FALSE)
            if (is.null(tbl) || !length(idx) || idx > nrow(tbl)) return()
            target <- pObjects$memory[[enc]]@LinkedPanel
            if (!nzchar(target) || is.null(pObjects$memory[[target]]) ||
                !"ColorBy" %in% methods::slotNames(pObjects$memory[[target]])) {
                message("[spiSEE] SVGPanel: LinkedPanel '", target,
                        "' not found or has no ColorBy.")
                return()
            }
            gene <- tbl$gene[idx]
            pObjects$memory[[target]]@ColorBy <- gene
            shiny::updateSelectizeInput(session, paste0(target, "_ColorBy"),
                                        selected = gene,
                choices = .colorByChoices(se), server = TRUE)
            .spatialUpdate(target, rObjects)
        },
        ignoreInit = TRUE,
        ignoreNULL = TRUE)
})

#' Ranking of spatially variable genes (Moran's I)
#'
#' iSEE panel with a table of genes ranked by Moran's I on the region
#' selected on the [SpatialMapPanel()], computed on *Compute* (Voyager, on
#' a k-nearest-neighbour graph of the cell centroids). Clicking a row shows
#' the gene in a linked [SpatialPlotPanel()].
#'
#' @param TopN Number of most expressed genes to scan (or the first `TopN`
#'   of `metadata(se)$gene_ranks`).
#' @param K Number of spatial neighbours per cell.
#' @param LinkedPanel Name of a [SpatialPlotPanel()] (e.g.
#'   `"SpatialPlotPanel2"`) that shows the gene clicked in the table.
#' @param ... Other slots: `Selection`, and `ComputeTrigger` (see
#'   [SpatialStatsPanel()]).
#' @return A `SVGPanel` object.
#' @seealso [SpatialStatsPanel()] for Moran's I of chosen genes.
#' @examples
#' SVGPanel(TopN = 100, LinkedPanel = "SpatialPlotPanel1")
#' @export
SVGPanel <- function(TopN = 200L, K = 6L, LinkedPanel = "", ...)
    new("SVGPanel", TopN = as.integer(TopN), K = as.integer(K),
        LinkedPanel = LinkedPanel, ...)
