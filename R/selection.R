# Selection: one slot (`Selection`), interpreted by the platform.

.splitIDs <- function(x) {
    if (length(x) == 0 || !nzchar(x)) character(0)
    else strsplit(x, ",", fixed = TRUE)[[1]]
}

# Subset to the current map selection: comma-separated FOV ids (tile
# platforms) or cell ids (cell platforms). Empty = no subsetting.
.subsetBySelection <- function(se, selection) {
    ids <- .splitIDs(selection)
    if (length(ids) == 0) return(se)
    if (.selectionUnit(se) == "fov") {
        se[, as.character(colData(se)[[.platformInfo(se)$region_col]]) %in% ids]
    } else {
        se[, colnames(se) %in% ids]
    }
}

# Every panel with a `Selection` slot receives the map's selection. Found by
# slot rather than a fixed list, so panels added at runtime are included.
.findSelectablePanels <- function(pObjects) {
    has_sel <- function(nm)
        "Selection" %in% methods::slotNames(pObjects$memory[[nm]])
    Filter(has_sel, names(pObjects$memory))
}

# All selectable region ids of the dataset, in display order.
.allRegionIDs <- function(se) {
    if (.selectionUnit(se) == "fov") {
        f <- unique(as.character(colData(se)[[.platformInfo(se)$region_col]]))
        if (all(grepl("^-?[0-9]+(\\.[0-9]+)?$",
                        f))) f[order(as.numeric(f))] else sort(f)
    } else colnames(se)
}

# A panel added while the app runs (Organize panels) comes with an empty
# Selection: give it the selection in force on the map.
.mapSelection <- function(pObjects) {
    is_map <- vapply(pObjects$memory, methods::is, logical(1),
                        "SpatialMapPanel")
    if (!any(is_map)) return(NULL)
    pObjects$memory[[which(is_map)[1]]]@CommittedIDs
}

#' @param x A panel extending `SpatialPanelBase`.
#' @param se The `SpatialExperiment` being explored.
#' @param input,session,pObjects,rObjects Shiny and iSEE state, passed by
#'   iSEE.
#' @return `.createObservers()` is called for its side effects (it gives a
#'   panel added while the app runs the selection of the map) and returns
#'   `NULL`.
#' @rdname SpatialPanelBase-class
#' @export
setMethod(".createObservers", "SpatialPanelBase",
            function(x, se, input, session, pObjects, rObjects) {
    callNextMethod()
    enc <- .getEncodedName(x)
    sel <- .mapSelection(pObjects)
    if (!is.null(sel) && !is.null(pObjects$memory[[enc]]) &&
        !identical(pObjects$memory[[enc]]@Selection, sel))
        pObjects$memory[[enc]]@Selection <- sel
    invisible(NULL)
})
