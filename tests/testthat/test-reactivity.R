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

test_that("panel errors become a readable validation message (F3)", {
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
