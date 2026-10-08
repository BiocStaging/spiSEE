#!/usr/bin/env Rscript
# One-command launcher (run from the package root; the run_*_app.command
# files do this for you).
#
#   Rscript run_app.R [qc|analysis] [path] [--dev] [--platform=cosmx|xenium|visiumhd]
#
#   path: an .RDS SpatialExperiment, a Xenium output directory, or a Visium HD
#         segmented_outputs directory. Defaults per mode below.
#   Old mode names still work: xenium-qc, xenium-analysis.
#   --dev: load from source (devtools::load_all) instead of the installed package.
options(shiny.launch.browser = TRUE)

args <- commandArgs(trailingOnly = TRUE)
dev_mode <- "--dev" %in% args
plat_arg <- grep("^--platform=", args, value = TRUE)
platform <- if (length(plat_arg)) sub("^--platform=", "", plat_arg[1]) else NULL
args <- args[!grepl("^--", args)]

aliases <- c(`xenium-qc` = "qc", `xenium-analysis` = "analysis")
mode <- if (length(args) && args[1] %in% c("qc", "analysis", names(aliases))) args[1] else "analysis"
if (mode %in% names(aliases)) { if (is.null(platform)) platform <- "xenium"; mode <- aliases[[mode]] }
rest <- if (length(args) && args[1] %in% c("qc", "analysis", names(aliases))) args[-1] else args

default_path <- switch(paste(mode, if (is.null(platform)) "cosmx" else platform),
  `qc cosmx`          = "data-raw/cos_unfiltered.RDS",
  `analysis cosmx`    = "../cos_isee_ready.RDS",
  `qc xenium`         = "data-raw/xenium_unfiltered.RDS",
  `analysis xenium`   = "data-raw/xenium_isee_ready.RDS",
  `qc visiumhd`       = "../visiumhd_isee_ready_unfiltered.RDS",
  `analysis visiumhd` = "../visiumhd_isee_ready.RDS",
  "../cos_isee_ready.RDS")
path <- if (length(rest)) rest[1] else default_path
if (!file.exists(path) && identical(platform, "xenium") && dir.exists("data-raw/xenium_test_data"))
  path <- "data-raw/xenium_test_data"

if (dev_mode || !requireNamespace("iseeSpatial", quietly = TRUE)) {
  suppressMessages(devtools::load_all(".", quiet = TRUE))
} else suppressMessages(library(iseeSpatial))

spe <- if (dir.exists(path)) {
  if (file.exists(file.path(path, "cell_segmentations.geojson")) || grepl("segmented_outputs", path)) load_visiumhd(path)
  else load_xenium(path)
} else readRDS(path)
if (!is.null(platform)) spe <- setPlatform(spe, platform)
cat(sprintf("Launching %s app (%s): %s, %d cells, %d genes\n",
            mode, getPlatform(spe), path, ncol(spe), nrow(spe)))
shiny::runApp(iseeSpatial(spe, mode = mode))
