#!/usr/bin/env Rscript
# Offline step: build the spatial graph once, compute niches and spatial
# domains on the FULL dataset, save an app-ready RDS. The niche panel then
# shows these instantly as long as its k / number of niches match; Compute
# re-runs on the current selection with the panel's parameters.
#
#   Rscript inst/scripts/precompute_niches.R in.RDS out.RDS [label_col] [k] [n_niches] [lambda]
#
# k = neighbours per cell. No universal value: small k describes the
# immediate contact zone, large k a wider region. Choose it from the
# structures you look at (and say which k in figures).
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2)
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
spe <- readRDS(args[1])
label_col <- if (length(args) >= 3 && nzchar(args[3])) args[3] else iseeSpatial:::.defaultLabelCol(spe)
k        <- if (length(args) >= 4) as.integer(args[4]) else 20L
n_niches <- if (length(args) >= 5) as.integer(args[5]) else 8L
lambda   <- if (length(args) >= 6) as.numeric(args[6]) else 0.2
message(sprintf("cells: %d | labels: %s | k = %d | niches = %d | lambda = %.2f",
                ncol(spe), label_col, k, n_niches, lambda))
t0 <- Sys.time(); spe <- buildSpatialGraph(spe, k = k);                     message("graph:   ", format(Sys.time() - t0))
t0 <- Sys.time(); spe <- computeNiches(spe, label_col, k = k, n_niches = n_niches); message("niches:  ", format(Sys.time() - t0))
t0 <- Sys.time(); spe <- computeSpatialDomains(spe, k_spatial = k, lambda = lambda); message("domains: ", format(Sys.time() - t0))
print(table(spe$niche, useNA = "ifany")); print(table(spe$spatial_domain))
saveRDS(spe, args[2])
