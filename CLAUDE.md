# iseeSpatial (modular) — handover for Claude Code

Written 8 Oct 2026. This folder is a **new copy** (`iseeSpatial_modular`).
Do NOT touch `../iseeSpatial`, `../iseeSpatial-cosmx` or `../iseeSpatial_test`.
Do not push, deploy to Policlinico, or delete anything without asking Emile.
Emile edits slides himself: never modify .pptx files.

## What changed vs `iseeSpatial` (branch xenium-support, ff6943e) + niches
- 21 platform-specific files → one class per panel (`R/AllClasses.R`):
  SpatialMapPanel, SpatialPlotPanel (= old Spatial + Gene panels), FOVZoomPanel,
  SpatialQCPanel, FOVStatsPanel, SpatialStatsPanel, SVGPanel, SpatialNichePanel.
  All but the map extend virtual `SpatialPanelBase` (slot `Selection`).
- Platform differences in one table: `R/platform.R` (`.PLATFORMS`, `getPlatform()`,
  `setPlatform()`). Units: "fov" (CosMx tiles) or "cell" (Xenium, Visium HD).
- Selection: map writes `Selection` into every panel that has the slot
  (`.findSelectablePanels`), interpreted by `.subsetBySelection()` (selection.R).
- Map JS merged into one script (`.MAP_JS` in panel_map.R) for both modes.
  Logic was checked in Node with mocked Plotly/Shiny (click, second click,
  box select, committed→pending removal, select all, placeholder ignored).
- `k` is a panel parameter in SpatialStatsPanel, SVGPanel, SpatialNichePanel.
- Niche panel reuses stored full-data niches only if label/k/n_niches match the panel.
- Gene menus are server-side selectize (Visium HD has ~18k genes).
- Visium HD: `load_visiumhd()` + `inst/scripts/visiumhd_osta_workflow.R`
  (OSTA chapter verbatim + gene ranks, PCA, graph, niches, domains).
- MERSCOPE removed (Emile's decision).
- One launcher `iseeSpatial(spe, mode = "analysis"|"qc")`; old CosMx*/Xenium*
  constructors + launchers are wrappers in `R/deprecated.R`.
- NAMESPACE was written by hand to match the roxygen tags: run `devtools::document()`
  and check the diff; delete stale files in `man/` if any.

## What was NOT run (no iSEE/SpatialExperiment in the build sandbox)
Only parsing, the JS mock test and the pure-matrix niche tests were run.
Everything below needs a real R session.

## Checklist (in order; fix in this folder only)
0. `git init && git add -A && git commit -m 'modular baseline'` so every fix is a reviewable commit.
1. `BiocManager::install(c("BiocNeighbors","igraph","VisiumIO","BiocStyle"))` if missing.
2. `devtools::document()` → fix roxygen errors; compare NAMESPACE to the hand-written one.
3. `devtools::test()` → all pass (real-data tests run if `../cos_isee_ready.RDS` exists).
4. Render without the app, for each platform, every panel's `.generateOutput()`
   (e.g. `iSEE:::.generateOutput(SpatialPlotPanel(ColorBy="SingleR_label", Selection="72,73,74"), spe)`)
   and LOOK at the images. "No error" is not a pass.
5. Launch: `Rscript run_app.R analysis ../cos_isee_niches.RDS --dev`, then qc on
   `data-raw/cos_unfiltered.RDS`, then Xenium (`xenium-analysis`, `xenium-qc`).
   Check by real mouse clicks: map click/box/All/Validate/Clear reach ALL panels,
   panels added via "Organize panels" also receive the selection, SVG row click
   updates SpatialPlotPanel2, Compute buttons, niche k change + Compute.
6. Visium HD (never run yet):
   `Rscript inst/scripts/visiumhd_osta_workflow.R ../visiumhd_isee_ready.RDS 20 osta_visiumhd_figures 2>&1 | tee visiumhd_osta.log`
   (on the server if memory is short). Part 1 = every chunk of OSTA ch. 17 (figures in
   `osta_visiumhd_figures/`): compare each figure and printed table with the book page;
   any difference = report, don't "fix" the chapter code. Then launch
   `Rscript run_app.R analysis ../visiumhd_isee_ready.RDS --platform=visiumhd --dev` and
   `Rscript run_app.R qc ../visiumhd_isee_ready_unfiltered.RDS --platform=visiumhd --dev`.
7. `R CMD build .` + `R CMD check --no-manual iseeSpatial_0.99.0.tar.gz`, then
   `BiocCheck::BiocCheck("iseeSpatial_0.99.0.tar.gz", "new-package" = TRUE)`.
   Target: 0 errors, 0 warnings; justify any note.
8. Log each day in `../LAB_NOTEBOOK.txt` (existing style).

## Things to verify (flagged, not confirmed)
- Every iSEE method has `@rdname <Class>`, so it lands on the class page; check `R CMD check` reports no undocumented S4 methods.
- `scrapper::normalizeRnaCounts.se()` writes an assay named `logcounts` (script stops if not).
- Visium HD `segmented_outputs` contains `cell_segmentations.geojson` (used by run_app.R auto-detect).
- `SpaceTrooper::plotCentroids/plotPolygons/plotMetricHist` on Visium HD objects
  (never tried; Xenium calls passed no `sampleId`, CosMx calls passed `sampleId = NULL` — kept as is).
- `SpaceTrooper::plotMetricHist` with Visium HD QC columns (`sum`, `um2`...).
- Point map with >100k cells (Visium HD): widget size / browser speed. If slow,
  subsample the clickable trace or switch to box-select only.
- `renderCachedPlot` keys now include the panel name (avoids cache sharing between panels).
- Example in `?iseeSpatial` builds an app on a toy object during R CMD check: confirm it does not error.
- Deployed apps on Policlinico call the old launchers; they still work through
  `R/deprecated.R`, but `iseeSpatial(spe)` now opens the ANALYSIS app (was QC).

## Open TODOs before Bioconductor submission
- Authors@R: add co-authors (Dario Righelli? Davide Risso?) and roles — Emile decides.
- URL / BugReports once the GitHub repo is public.
- Small real example data for the vignette (STexampleData / ExperimentHub?) instead of simulation.
- GitHub Actions (biocthis) for R CMD check + BiocCheck.
- Native iSEE multi-selection (tissue ↔ UMAP brushing) — next feature.
