#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tools)
})

base_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation"
figure_id <- "figure_05_benchmark_grid"
tag <- Sys.getenv("SAGUI_FIG5_TAG", unset = "sagui8")
groups_env <- Sys.getenv("SAGUI_FIG5_GROUPS", unset = "10,20,30,40")
groups <- as.integer(strsplit(groups_env, ",", fixed = TRUE)[[1]])
groups <- groups[is.finite(groups) & groups > 1]
if (!length(groups)) {
  stop("SAGUI_FIG5_GROUPS must contain at least one integer greater than 1.")
}

voronoi_script <- file.path(base_dir, "src/R/pipelines/compare_sagui_voronoi_demo.R")
pixedfit_script <- file.path(base_dir, "src/R/pipelines/compare_sagui_pixedfit_r_demo.R")

stage_root <- file.path(base_dir, "results/paper_repro_stage", figure_id, tag)
final_dir <- file.path(base_dir, "results/figures/paper_repro", figure_id, paste0(tag, "_psfmatched_3x4"))
overleaf_dir <- file.path(base_dir, "results/figures/overleaf_upload", paste0(tag, "_benchmark_3x4"))
dir.create(stage_root, recursive = TRUE, showWarnings = FALSE)
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(overleaf_dir, recursive = TRUE, showWarnings = FALSE)

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

manifest_rows <- vector("list", length(groups) * 3L)
row_idx <- 1L

for (ncomp in groups) {
  run_compare(
    voronoi_script,
    c(
      paste0("SAGUI_COMPARE_TAG=", tag),
      paste0("SAGUI_COMPARE_NCOMP=", ncomp),
      "SAGUI_COMPARE_SUFFIX=figure05_psfmatched_sagui",
      "SAGUI_COMPARE_CUBE_MODE=psfmatched",
      paste0("SAGUI_COMPARE_OUTPUT_DIR=", stage_root)
    )
  )

  run_compare(
    voronoi_script,
    c(
      paste0("SAGUI_COMPARE_TAG=", tag),
      paste0("SAGUI_COMPARE_NCOMP=", ncomp),
      "SAGUI_COMPARE_SUFFIX=figure05_psfmatched_vornomask",
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
      "SAGUI_COMPARE_SUFFIX=figure05_psfmatched_pixedfitr",
      "SAGUI_COMPARE_CUBE_MODE=psfmatched",
      paste0("SAGUI_COMPARE_OUTPUT_DIR=", stage_root)
    )
  )

  sagui_run <- build_run_tag(tag, ncomp, "figure05_psfmatched_sagui")
  voronoi_run <- build_run_tag(tag, ncomp, "figure05_psfmatched_vornomask")
  pix_run <- build_run_tag(tag, ncomp, "figure05_psfmatched_pixedfitr")

  specs <- list(
    list(
      method = "sagui",
      run_tag = sagui_run,
      source = "regions_plot_sagui.png",
      target = sprintf("sagui_n%02d.png", ncomp),
      row = 1L,
      overleaf_subdir = sprintf("%s_n%02d_paperexact", tag, ncomp)
    ),
    list(
      method = "voronoi_nomask",
      run_tag = voronoi_run,
      source = "regions_plot_voronoi.png",
      target = sprintf("voronoi_n%02d.png", ncomp),
      row = 2L,
      overleaf_subdir = sprintf("%s_n%02d_paperexact_vornomask", tag, ncomp)
    ),
    list(
      method = "pixedfit_r",
      run_tag = pix_run,
      source = "regions_plot_pixedfit_r.png",
      target = sprintf("pixedfit_r_n%02d.png", ncomp),
      row = 3L,
      overleaf_subdir = sprintf("%s_n%02d_pixedfitr", tag, ncomp)
    )
  )

  for (spec in specs) {
    src <- file.path(stage_root, "results/figures/sagui_compare", spec$run_tag, spec$source)
    dst <- file.path(final_dir, spec$target)
    overleaf_subdir <- file.path(overleaf_dir, spec$overleaf_subdir)
    overleaf_dst <- file.path(overleaf_subdir, spec$source)
    dir.create(overleaf_subdir, recursive = TRUE, showWarnings = FALSE)
    copy_required(src, dst)
    copy_required(src, overleaf_dst)

    manifest_rows[[row_idx]] <- data.frame(
      method = spec$method,
      tag = tag,
      ncomp = ncomp,
      cube_mode = "psfmatched",
      row_index = spec$row,
      col_index = match(ncomp, groups),
      run_tag = spec$run_tag,
      source_png = src,
      final_png = dst,
      overleaf_png = overleaf_dst,
      stringsAsFactors = FALSE
    )
    row_idx <- row_idx + 1L
  }
}

manifest <- do.call(rbind, manifest_rows)
utils::write.csv(manifest, file.path(final_dir, "manifest.csv"), row.names = FALSE)
utils::write.csv(manifest, file.path(overleaf_dir, "manifest.csv"), row.names = FALSE)

message("Finished ", figure_id, " for ", tag)
message("Outputs: ", final_dir)
message("Overleaf bundle: ", overleaf_dir)
