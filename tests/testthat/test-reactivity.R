test_that(".panelKey reacts only to its own panel and to real slot changes", {
    shiny::testServer(function(input, output, session) {
        pObjects <- new.env()
        pObjects$memory <- list(A = SpatialPlotPanel(ColorBy = "a"),
                                B = SpatialPlotPanel(ColorBy = "b"))
        rObjects <- shiny::reactiveValues(modified = list())
        key <- .panelKey("A", pObjects, rObjects, function(pn) pn@ColorBy)
        seen <- new.env()
        seen$runs <- 0
        shiny::observe({
            seen$key <- key()
            seen$runs <- seen$runs + 1
        })
    }, {
        session$flushReact()
        expect_equal(seen$runs, 1)
        # another panel updated, then iSEE's reset
        rObjects$modified <- list(B = character(0))
        session$flushReact()
        rObjects$modified <- list()
        session$flushReact()
        expect_equal(seen$runs, 1)
        # this panel marked but its key slots unchanged
        rObjects$modified <- list(A = character(0))
        session$flushReact()
        rObjects$modified <- list()
        session$flushReact()
        expect_equal(seen$runs, 1)
        # this panel's slot really changed
        pObjects$memory$A@ColorBy <- "c"
        rObjects$modified <- list(A = character(0))
        session$flushReact()
        expect_equal(seen$runs, 2)
        expect_identical(seen$key, "c")
        # our own update signal: re-runs on a real change, and leaves iSEE's
        # rObjects$modified alone (iSEE would recompute the panel itself)
        pObjects$memory$A@ColorBy <- "d"
        before <- rObjects$modified
        .spatialUpdate("A", rObjects)
        session$flushReact()
        expect_equal(seen$runs, 3)
        expect_identical(seen$key, "d")
        expect_identical(rObjects$modified, before)
        .spatialUpdate("A", rObjects)
        session$flushReact()
        expect_equal(seen$runs, 3)
    })
})

test_that("an empty ColorBy from a loading menu is ignored", {
    expect_identical(.castColorBy(""), NA_character_)
    expect_identical(.castColorBy("Malat1"), "Malat1")
})

test_that("panel errors become a readable validation message", {
    expect_identical(.safeOutput(42, "P"), 42)
    expect_message(
        err <- tryCatch(.safeOutput(stop("boom"), "P"), error = function(e) e),
        "P: boom")
    expect_s3_class(err, "shiny.silent.error")
    expect_match(conditionMessage(err), "Error: boom")
    # req() / validate() from inside a panel pass through untouched
    err2 <- tryCatch(.safeOutput(shiny::req(FALSE), "P"), error = function(e) e)
    expect_s3_class(err2, "shiny.silent.error")
})

test_that("map selection: browser state is stored, Validate sets every panel", {
    shiny::testServer(function(input, output, session) {
        pObjects <- new.env()
        pObjects$memory <- list(SpatialMapPanel1 = SpatialMapPanel(),
                                SpatialPlotPanel1 = SpatialPlotPanel(),
                                SpatialQCPanel1 = SpatialQCPanel())
        rObjects <- shiny::reactiveValues(modified = list())
        .mapPendingObservers("SpatialMapPanel1", input, pObjects)
    }, {
        session$flushReact()   # observers use ignoreInit = TRUE
        session$setInputs(SpatialMapPanel1_pendingSync = c("3", "5"),
                          SpatialMapPanel1_pendingRemoveSync = "7")
        map <- pObjects$memory$SpatialMapPanel1
        expect_identical(map@PendingIDs, "3,5")
        expect_identical(map@PendingRemoveIDs, "7")

        .applyMapSelection("SpatialMapPanel1", c("3", "5"), session,
                           pObjects, rObjects)
        expect_identical(pObjects$memory$SpatialPlotPanel1@Selection, "3,5")
        expect_identical(pObjects$memory$SpatialQCPanel1@Selection, "3,5")
        map <- pObjects$memory$SpatialMapPanel1
        expect_identical(map@CommittedIDs, "3,5")
        expect_identical(map@PendingIDs, "")
        expect_identical(map@PendingRemoveIDs, "")
        expect_false(is.null(shiny::isolate(
            rObjects[[.updateName("SpatialPlotPanel1")]])))
    })
})

test_that("outputs are drawn from the panel as it was when its key changed", {
    shiny::testServer(function(input, output, session) {
        pObjects <- new.env()
        pObjects$memory <- list(A = SpatialStatsPanel(
            PlotType = "joincount_plot", ComputeTrigger = "k1"))
        rObjects <- shiny::reactiveValues(modified = list())
        key <- .panelKey("A", pObjects, rObjects,
                            function(pn) pn@ComputeTrigger)
    }, {
        session$flushReact()
        # menu changed, Compute not pressed: the drawn panel stays the same
        pObjects$memory$A@PlotType <- "graph_plot"
        .spatialUpdate("A", rObjects)
        session$flushReact()
        expect_identical(.keyPanel(key)@PlotType, "joincount_plot")
        # Compute: key and drawn panel change together
        pObjects$memory$A@ComputeTrigger <- "k2"
        .spatialUpdate("A", rObjects)
        session$flushReact()
        expect_identical(key(), "k2")
        expect_identical(.keyPanel(key)@PlotType, "graph_plot")
    })
})

test_that("a removed panel's key observer stops quietly", {
    shiny::testServer(function(input, output, session) {
        pObjects <- new.env()
        pObjects$memory <- list(A = SpatialPlotPanel(ColorBy = "a"))
        rObjects <- shiny::reactiveValues(modified = list())
        key <- .panelKey("A", pObjects, rObjects, function(pn) pn@ColorBy)
    }, {
        session$flushReact()
        pObjects$memory$A <- NULL
        rObjects$modified <- list(B = character(0))
        expect_no_error(session$flushReact())
        rObjects$modified <- list()
        expect_no_error(session$flushReact())
    })
})
