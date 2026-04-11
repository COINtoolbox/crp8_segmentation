#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(FITSio)
  library(dplyr)
})

base_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation"
figure_id <- "figure_04_strategy_compare"
tag <- "sagui10"
ncomp <- 8L

voronoi_script <- file.path(base_dir, "src/R/pipelines/compare_sagui_voronoi_demo.R")
pixedfit_script <- file.path(base_dir, "src/R/pipelines/compare_sagui_pixedfit_r_demo.R")
region_photometry_path <- file.path(base_dir, "src/R/utils/RegionPhotometry.R")

stage_root <- file.path(base_dir, "results/paper_repro_stage", paste0(figure_id, "_n08"))
figure_out_dir <- file.path(base_dir, "results/figures/paper_repro", figure_id, paste0(tag, "_psfmatched_n08"))
seg_out_dir <- file.path(base_dir, "results/segmentation/paper_repro", figure_id, paste0(tag, "_psfmatched_n08"))
flux_out_dir <- file.path(base_dir, "results/flux_per_region/paper_repro", figure_id, paste0(tag, "_psfmatched_n08"))

dir.create(stage_root, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(seg_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(flux_out_dir, recursive = TRUE, showWarnings = FALSE)

source(region_photometry_path)

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
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(from, to, overwrite = TRUE)
  if (!ok) {
    stop("Failed to copy ", from, " -> ", to)
  }
}

resolve_raw_cube <- function(tag) {
  dashed <- sub("^sagui", "sagui-", gsub("_", "-", tag))
  candidates <- c(
    file.path(base_dir, "data/raw", paste0("datacube_", tag, ".fits")),
    file.path(base_dir, "data/raw", paste0("datacube_", dashed, ".fits"))
  )
  path <- candidates[file.exists(candidates)][1]
  if (is.na(path)) {
    stop("Raw cube not found for ", tag)
  }
  path
}

write_region_photometry_csv <- function(fits_obj, cluster_map, path, error_fallback = c("mad_sky", "flux_over_sqrt_n")) {
  error_fallback <- match.arg(error_fallback)
  sed <- RegionPhotometry(
    fits_obj,
    cluster_map,
    error_fallback = error_fallback
  )

  nircam_filters <- c(
    "F090W", "F115W", "F150W", "F182M", "F200W",
    "F210M", "F277W", "F335M", "F356W", "F410M",
    "F430M", "F444W", "F460M", "F480M"
  )

  raw_names <- names(sed$flux_wide)
  band_cols <- raw_names[grepl("^\\s*[0-9]+\\s*$", raw_names)]
  band_cols <- band_cols[!grepl("_err$|_n_eff$", band_cols)]
  err_cols <- paste0(band_cols, "_err")
  neff_cols <- paste0(band_cols, "_n_eff")
  filters <- nircam_filters[seq_along(band_cols)]

  flux_ordered_jy <- sed$flux_wide |>
    dplyr::select(
      dplyr::any_of(c("region", "n_pix")),
      dplyr::all_of(band_cols),
      dplyr::all_of(err_cols),
      dplyr::all_of(neff_cols)
    ) |>
    dplyr::mutate(
      dplyr::across(dplyr::all_of(band_cols), ~ .x * 1e-8),
      dplyr::across(dplyr::all_of(err_cols), ~ .x * 1e-8)
    ) |>
    dplyr::rename_with(~ filters, .cols = dplyr::all_of(band_cols)) |>
    dplyr::rename_with(~ paste0(filters, "_err"), .cols = dplyr::all_of(err_cols)) |>
    dplyr::rename_with(~ paste0(filters, "_n_eff"), .cols = dplyr::all_of(neff_cols))

  utils::write.csv(flux_ordered_jy, path, row.names = FALSE)
  invisible(flux_ordered_jy)
}

run_compare(
  voronoi_script,
  c(
    paste0("SAGUI_COMPARE_TAG=", tag),
    paste0("SAGUI_COMPARE_NCOMP=", ncomp),
    "SAGUI_COMPARE_SUFFIX=figure04_psfmatched_n08_sagui",
    "SAGUI_COMPARE_CUBE_MODE=psfmatched",
    paste0("SAGUI_COMPARE_OUTPUT_DIR=", stage_root)
  )
)

run_compare(
  voronoi_script,
  c(
    paste0("SAGUI_COMPARE_TAG=", tag),
    paste0("SAGUI_COMPARE_NCOMP=", ncomp),
    "SAGUI_COMPARE_SUFFIX=figure04_psfmatched_n08_vornomask",
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
    "SAGUI_COMPARE_SUFFIX=figure04_psfmatched_n08_pixedfitr",
    "SAGUI_COMPARE_CUBE_MODE=psfmatched",
    paste0("SAGUI_COMPARE_OUTPUT_DIR=", stage_root)
  )
)

sagui_run <- build_run_tag(tag, ncomp, "figure04_psfmatched_n08_sagui")
voronoi_run <- build_run_tag(tag, ncomp, "figure04_psfmatched_n08_vornomask")
pix_run <- build_run_tag(tag, ncomp, "figure04_psfmatched_n08_pixedfitr")

sagui_fig_dir <- file.path(stage_root, "results/figures/sagui_compare", sagui_run)
voronoi_fig_dir <- file.path(stage_root, "results/figures/sagui_compare", voronoi_run)
pix_fig_dir <- file.path(stage_root, "results/figures/sagui_compare", pix_run)

sagui_seg_dir <- file.path(stage_root, "results/segmentation/sagui_compare", sagui_run)
voronoi_seg_dir <- file.path(stage_root, "results/segmentation/sagui_compare", voronoi_run)
pix_seg_dir <- file.path(stage_root, "results/segmentation/sagui_compare", pix_run)

copy_required(file.path(sagui_fig_dir, "regions_plot_sagui.png"), file.path(figure_out_dir, "regions_plot_sagui.png"))
copy_required(file.path(voronoi_fig_dir, "regions_plot_voronoi.png"), file.path(figure_out_dir, "regions_plot_voronoi_nomask.png"))
copy_required(file.path(pix_fig_dir, "regions_plot_pixedfit_r.png"), file.path(figure_out_dir, "regions_plot_pixedfit_r.png"))

copy_required(file.path(sagui_seg_dir, "region_map_sagui.fits"), file.path(seg_out_dir, "region_map_sagui.fits"))
copy_required(file.path(voronoi_seg_dir, "region_map_voronoi.fits"), file.path(seg_out_dir, "region_map_voronoi_nomask.fits"))
copy_required(file.path(pix_seg_dir, "region_map_pixedfit_r.fits"), file.path(seg_out_dir, "region_map_pixedfit_r.fits"))

copy_required(file.path(sagui_seg_dir, "region_metrics_sagui.csv"), file.path(seg_out_dir, "region_metrics_sagui.csv"))
copy_required(file.path(voronoi_seg_dir, "region_metrics_voronoi.csv"), file.path(seg_out_dir, "region_metrics_voronoi_nomask.csv"))
copy_required(file.path(pix_seg_dir, "region_metrics_pixedfit_r.csv"), file.path(seg_out_dir, "region_metrics_pixedfit_r.csv"))
copy_required(file.path(sagui_seg_dir, "comparison_summary.csv"), file.path(seg_out_dir, "comparison_summary_sagui_vs_voronoi.csv"))
copy_required(file.path(pix_seg_dir, "comparison_summary.csv"), file.path(seg_out_dir, "comparison_summary_sagui_vs_pixedfit_r.csv"))

raw_cube_path <- resolve_raw_cube(tag)
raw_fits <- FITSio::readFITS(raw_cube_path)

region_maps <- list(
  sagui = FITSio::readFITS(file.path(seg_out_dir, "region_map_sagui.fits"))$imDat,
  voronoi_nomask = FITSio::readFITS(file.path(seg_out_dir, "region_map_voronoi_nomask.fits"))$imDat,
  pixedfit_r = FITSio::readFITS(file.path(seg_out_dir, "region_map_pixedfit_r.fits"))$imDat
)

for (method_name in names(region_maps)) {
  map <- region_maps[[method_name]]
  write_region_photometry_csv(
    raw_fits,
    map,
    file.path(flux_out_dir, sprintf("SED_flux_wide_%s_%s_mad_sky.csv", tag, method_name)),
    error_fallback = "mad_sky"
  )
  write_region_photometry_csv(
    raw_fits,
    map,
    file.path(flux_out_dir, sprintf("SED_flux_wide_%s_%s_flux_over_sqrt_n.csv", tag, method_name)),
    error_fallback = "flux_over_sqrt_n"
  )
}

sn_summary <- utils::read.csv("/tmp/sagui10_sn_sweep_flux_over_sqrt_n_summary.csv", stringsAsFactors = FALSE)
sn_mad_summary <- utils::read.csv("/tmp/sagui10_sn_sweep_mad_sky_summary.csv", stringsAsFactors = FALSE)
utils::write.csv(sn_summary, file.path(seg_out_dir, "sn_sweep_flux_over_sqrt_n_summary.csv"), row.names = FALSE)
utils::write.csv(sn_mad_summary, file.path(seg_out_dir, "sn_sweep_mad_sky_summary.csv"), row.names = FALSE)

manifest <- data.frame(
  tag = tag,
  ncomp = ncomp,
  figure_dir = figure_out_dir,
  segmentation_dir = seg_out_dir,
  flux_dir = flux_out_dir,
  raw_cube = raw_cube_path,
  stage_root = stage_root,
  stringsAsFactors = FALSE
)
utils::write.csv(manifest, file.path(figure_out_dir, "manifest.csv"), row.names = FALSE)

message("Finished ", figure_id, " N=", ncomp)
message("Figures: ", figure_out_dir)
message("Segmentation FITS: ", seg_out_dir)
message("Flux tables: ", flux_out_dir)
