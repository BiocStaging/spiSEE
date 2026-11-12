# iseeSpatial 0.99.0

* Submission version for Bioconductor.
* One class per panel for all platforms (`SpatialMapPanel`, `SpatialPlotPanel`,
  `FOVZoomPanel`, `SpatialQCPanel`, `FOVStatsPanel`, `SpatialStatsPanel`,
  `SVGPanel`, `SpatialNichePanel`). Platform differences live in one table
  (`setPlatform()` / `getPlatform()`).
* New platform: Visium HD (cell-segmented), `load_visiumhd()` and the OSTA
  workflow script in `inst/scripts/`.
* Cellular niches, neighbourhood entropy and spatial domains
  (`buildSpatialGraph()`, `computeNiches()`, `computeSpatialDomains()`).
* The number of neighbours `k` is a parameter everywhere a neighbour graph is
  built (statistics, SVG ranking, niches).
* One launcher, `iseeSpatial(spe, mode = c("analysis", "qc"))`. **Breaking:**
  `iseeSpatial(spe)` used to open the CosMx QC app; it now opens the analysis
  app (use `mode = "qc"`).
* The old platform-specific `CosMx*`/`Xenium*` panel constructors and
  launchers are removed; use the classes above and `iseeSpatial()`.
* Plot images are cached on disk (`spatialCacheDir()`, 1 GB, least recently
  used removed first) so views already seen reopen instantly, also after a
  restart; `cache = FALSE` keeps them in memory, `clearSpatialCache()` empties
  the folder.
* Niches, spatially variable genes and statistics are computed only when
  *Compute* is pressed.
* Statistics, gene rankings and niches are cached by what they compute
  (panel, parameters, region, dataset): a repeated analysis reopens from the
  cache, also after a restart.
* `plotSpatialCells()`: static maps of the cells as polygons, coloured by
  annotations, QC metrics or genes, with the colours of the app.
* Vignettes: *Get started* (SpaceTrooper's small CosMx example) and a
  complete CosMx workflow on the mouse dentate gyrus (OSTA.data).
