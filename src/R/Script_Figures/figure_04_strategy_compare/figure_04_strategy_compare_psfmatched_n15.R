#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(FITSio)
  library(dplyr)
})

base_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation"
workspace_dir <- "/Users/rd23aag/Documents/GitHub/sagui"
figure_id <- "figure_04_strategy_compare"
tag <- Sys.getenv("SAGUI_FIG4_TAG", unset = "sagui10")
ncomp <- suppressWarnings(as.integer(Sys.getenv("SAGUI_FIG4_NCOMP", unset = "15")))
if (!is.finite(ncomp) || ncomp <= 1L) {
  stop("SAGUI_FIG4_NCOMP must be an integer greater than 1.")
}

voronoi_script <- file.path(workspace_dir, "experiments/paper_repro/patches/compare_sagui_voronoi_demo.R")
pixedfit_script <- file.path(workspace_dir, "experiments/paper_repro/patches/compare_sagui_pixedfit_r_demo.R")
region_photometry_path <- file.path(base_dir, "src/R/utils/RegionPhotometry.R")

slug <- sprintf("%s_psfmatched_n%02d", tag, ncomp)
stage_root <- file.path(base_dir, "results/paper_repro_stage", paste0(figure_id, sprintf("_n%02d", ncomp)))
figure_out_dir <- file.path(base_dir, "results/figures/paper_repro", figure_id, slug)
seg_out_dir <- file.path(base_dir, "results/segmentation/paper_repro", figure_id, slug)
flux_out_dir <- file.path(base_dir, "results/flux_per_region/paper_repro", figure_id, slug)

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

extract_header_filters <- function(fits_obj) {
  parse_fits_scalar <- function(line) {
    quoted <- regmatches(line, regexpr("'[^']*'", line))
    if (length(quoted) == 1L && nzchar(quoted)) {
      return(trimws(gsub("^'|'$", "", quoted)))
    }
    rhs <- sub("^[^=]+=\\s*", "", line)
    rhs <- sub("\\s*/.*$", "", rhs)
    trimws(rhs)
  }

  hdr_lines <- fits_obj$header
  filters <- character(0)
  if (!is.null(hdr_lines)) {
    filter_lines <- hdr_lines[grepl("^FILTER[0-9]+\\s*=", hdr_lines)]
    if (length(filter_lines)) {
      filters <- vapply(filter_lines, parse_fits_scalar, character(1))
      filters <- filters[nzchar(filters)]
    }
  }

  if (!length(filters) && !is.null(fits_obj$hdr)) {
    hdr <- fits_obj$hdr
    keys <- hdr[seq(1L, length(hdr) - 1L, by = 2L)]
    vals <- hdr[seq(2L, length(hdr), by = 2L)]
    filters <- trimws(vals[grepl("^FILTER[0-9]+$", keys)])
    filters <- filters[nzchar(filters)]
  }

  if (!length(filters)) {
    stop("No FILTERn keywords found in FITS header.")
  }
  unname(filters)
}

write_region_photometry_csv <- function(fits_obj, cluster_map, path, error_fallback = c("mad_sky", "flux_over_sqrt_n")) {
  error_fallback <- match.arg(error_fallback)
  sed <- RegionPhotometry(
    fits_obj,
    cluster_map,
    error_fallback = error_fallback
  )

  raw_names <- names(sed$flux_wide)
  band_cols <- raw_names[grepl("^\\s*[0-9]+\\s*$", raw_names)]
  band_cols <- band_cols[!grepl("_err$|_n_eff$", band_cols)]
  err_cols <- paste0(band_cols, "_err")
  neff_cols <- paste0(band_cols, "_n_eff")
  filters <- extract_header_filters(fits_obj)
  if (length(filters) != length(band_cols)) {
    stop(
      "Header filter count (", length(filters),
      ") does not match photometry band count (", length(band_cols), ")."
    )
  }

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
    sprintf("SAGUI_COMPARE_SUFFIX=figure04_psfmatched_n%02d_sagui", ncomp),
    "SAGUI_COMPARE_CUBE_MODE=psfmatched",
    paste0("SAGUI_COMPARE_OUTPUT_DIR=", stage_root)
  )
)

run_compare(
  voronoi_script,
  c(
    paste0("SAGUI_COMPARE_TAG=", tag),
    paste0("SAGUI_COMPARE_NCOMP=", ncomp),
    sprintf("SAGUI_COMPARE_SUFFIX=figure04_psfmatched_n%02d_vornomask", ncomp),
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
    sprintf("SAGUI_COMPARE_SUFFIX=figure04_psfmatched_n%02d_pixedfitr", ncomp),
    "SAGUI_COMPARE_CUBE_MODE=psfmatched",
    paste0("SAGUI_COMPARE_OUTPUT_DIR=", stage_root)
  )
)

sagui_run <- build_run_tag(tag, ncomp, sprintf("figure04_psfmatched_n%02d_sagui", ncomp))
voronoi_run <- build_run_tag(tag, ncomp, sprintf("figure04_psfmatched_n%02d_vornomask", ncomp))
pix_run <- build_run_tag(tag, ncomp, sprintf("figure04_psfmatched_n%02d_pixedfitr", ncomp))

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
