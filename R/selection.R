# ============================================================
# Selection: one slot (`Selection`), interpreted by the platform.
# ============================================================

.splitIDs <- function(x) if (length(x) == 0 || !nzchar(x)) character(0) else strsplit(x, ",", fixed = TRUE)[[1]]

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
  Filter(function(nm) "Selection" %in% methods::slotNames(pObjects$memory[[nm]]),
         names(pObjects$memory))
}

# All selectable region ids of the dataset, in display order.
.allRegionIDs <- function(se) {
  if (.selectionUnit(se) == "fov") {
    f <- unique(as.character(colData(se)[[.platformInfo(se)$region_col]]))
    num <- suppressWarnings(as.numeric(f))
    if (!anyNA(num)) f[order(num)] else sort(f)
  } else colnames(se)
}
