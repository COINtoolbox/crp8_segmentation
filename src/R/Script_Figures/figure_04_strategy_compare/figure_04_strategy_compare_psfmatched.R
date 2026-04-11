#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tools)
})

base_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation"
figure_id <- "figure_04_strategy_compare"
tag <- "sagui10"
ncomp <- 10L

voronoi_script <- file.path(base_dir, "src/R/pipelines/compare_sagui_voronoi_demo.R")
pixedfit_script <- file.path(base_dir, "src/R/pipelines/compare_sagui_pixedfit_r_demo.R")

stage_root <- file.path(base_dir, "results/paper_repro_stage", figure_id)
final_dir <- file.path(base_dir, "results/figures/paper_repro", figure_id, paste0(tag, "_psfmatched"))
dir.create(stage_root, recursive = TRUE, showWarnings = FALSE)
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

run_compare <- function(script, env_vars) {
  status <- system2("Rscript", script, env = env_vars)
  if (!identical(status, 0L)) {
    stop("Run failed for script: ", script)
  }
}

build_run_tag <- function(tag, ncomp, suffix) {
  sprintf("%s_n%02d_%s", tag, ncomp, suffix)
}

copy_required <- function(from, to) {
  if (!file.exists(from)) {
    stop("Expected file not found: ", from)
  }
  ok <- file.copy(from, to, overwrite = TRUE)
  if (!ok) {
    stop("Failed to copy ", from, " -> ", to)
  }
}

run_compare(
  voronoi_script,
  c(
    paste0("SAGUI_COMPARE_TAG=", tag),
    paste0("SAGUI_COMPARE_NCOMP=", ncomp),
    "SAGUI_COMPARE_SUFFIX=figure04_psfmatched_sagui",
    "SAGUI_COMPARE_CUBE_MODE=psfmatched",
    paste0("SAGUI_COMPARE_OUTPUT_DIR=", stage_root)
  )
)

run_compare(
  voronoi_script,
  c(
    paste0("SAGUI_COMPARE_TAG=", tag),
    paste0("SAGUI_COMPARE_NCOMP=", ncomp),
    "SAGUI_COMPARE_SUFFIX=figure04_psfmatched_vornomask",
    "SAGUI_COMPARE_CUBE_MODE=psfmatched",
    "SAGUI_COMPARE_VORONOI_MASK=none",
    paste0("SAGUI_COMPARE_OUTPUT_DIR=", stage_root)
  )
)

run_compare(
  pixedfit_script,
  c(
    paste0("SAGUI_COMPARE_TAG=", tag),
    paste0("SAGUI_COMPARE_NCOMP=", ncomp),
    "SAGUI_COMPARE_SUFFIX=figure04_psfmatched_pixedfitr",
    "SAGUI_COMPARE_CUBE_MODE=psfmatched",
    paste0("SAGUI_COMPARE_OUTPUT_DIR=", stage_root)
  )
)

sagui_run <- build_run_tag(tag, ncomp, "figure04_psfmatched_sagui")
voronoi_run <- build_run_tag(tag, ncomp, "figure04_psfmatched_vornomask")
pix_run <- build_run_tag(tag, ncomp, "figure04_psfmatched_pixedfitr")

sagui_src <- file.path(stage_root, "results/figures/sagui_compare", sagui_run, "regions_plot_sagui.png")
voronoi_src <- file.path(stage_root, "results/figures/sagui_compare", voronoi_run, "regions_plot_voronoi.png")
pix_src <- file.path(stage_root, "results/figures/sagui_compare", pix_run, "regions_plot_pixedfit_r.png")

sagui_out <- file.path(final_dir, "regions_plot_sagui.png")
voronoi_out <- file.path(final_dir, "regions_plot_voronoi_nomask.png")
pix_out <- file.path(final_dir, "regions_plot_pixedfit_r.png")

copy_required(sagui_src, sagui_out)
copy_required(voronoi_src, voronoi_out)
copy_required(pix_src, pix_out)

manifest <- data.frame(
  method = c("voronoi_nomask", "pixedfit_r", "sagui"),
  tag = tag,
  ncomp = ncomp,
  cube_mode = "psfmatched",
  run_tag = c(voronoi_run, pix_run, sagui_run),
  source_png = c(voronoi_src, pix_src, sagui_src),
  final_png = c(voronoi_out, pix_out, sagui_out),
  stringsAsFactors = FALSE
)

utils::write.csv(manifest, file.path(final_dir, "manifest.csv"), row.names = FALSE)

message("Finished ", figure_id)
message("Outputs: ", final_dir)
