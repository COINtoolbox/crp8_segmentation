#!/usr/bin/env Rscript

script_path <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])
script_dir <- dirname(normalizePath(script_path))
source(file.path(script_dir, "figure_03_segmentation_panels_common.R"))
run_figure_03_segmentation_variant("psfmatched")
