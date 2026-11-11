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
