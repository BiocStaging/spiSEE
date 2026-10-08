# iseeSpatial

Interactive exploration of cell-resolution spatial transcriptomics (CosMx,
Xenium, Visium HD) with [iSEE](https://bioconductor.org/packages/iSEE).
Select a region on the overview map; every panel follows: spatial plots, QC
histograms, spatial statistics, spatially variable genes, cellular niches and
spatial domains.

```r
library(iseeSpatial)
spe <- load_visiumhd("segmented_outputs")    # or load_xenium(), load_cosmx()
iseeSpatial(spe, mode = "analysis")          # or mode = "qc"
```

See the vignette for niches/domains and `inst/scripts/` for full workflows
(CosMx preparation, Xenium and Visium HD OSTA workflows, niche precomputation).
