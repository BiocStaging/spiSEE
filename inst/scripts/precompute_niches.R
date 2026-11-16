#!/usr/bin/env Rscript
# Offline step: build the spatial graph once, compute niches and spatial
# domains on the FULL dataset, save an app-ready RDS. The niche panel then
# shows these at once as long as its parameters match; Compute re-runs on
# the current selection with the panel's parameters.
#
#   Rscript precompute_niches.R in.rds out.rds label_col [k] [n_niches] [lambda]
#
# The script's path: system.file("scripts", "precompute_niches.R",
#                                package = "spiSEE")
# k = neighbours per cell. There is no universal value: a small k describes
# the immediate contact zone, a large k a wider region. Choose it from the
# structures studied, and give it in figures.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
    stop("usage: precompute_niches.R in.rds out.rds label_col ",
        "[k] [n_niches] [lambda]")
suppressPackageStartupMessages(library(spiSEE))
spe <- readRDS(args[1])
label_col <- args[3]
k <- if (length(args) >= 4) as.integer(args[4]) else 20L
n_niches <- if (length(args) >= 5) as.integer(args[5]) else 6L
lambda <- if (length(args) >= 6) as.numeric(args[6]) else 0.2
message(sprintf("cells: %d | labels: %s | k = %d | niches = %d | lambda = %.2f",
    ncol(spe), label_col, k, n_niches, lambda))
timed <- function(what, expr) {
    t0 <- Sys.time()
    out <- force(expr)
    message(what, ": ", format(Sys.time() - t0))
    out
}
spe <- timed("graph", buildSpatialGraph(spe, k = k))
spe <- timed("niches", computeNiches(spe, label_col, k = k,
    n_niches = n_niches))
spe <- timed("domains", computeSpatialDomains(spe, k_spatial = k,
    lambda = lambda))
print(table(spe$niche, useNA = "ifany"))
print(table(spe$spatial_domain))
saveRDS(spe, args[2])
