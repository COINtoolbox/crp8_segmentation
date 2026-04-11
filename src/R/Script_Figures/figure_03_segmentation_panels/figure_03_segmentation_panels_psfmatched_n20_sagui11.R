#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(FITSio)
  library(guara)
  library(capivara)
  library(ragg)
  library(dplyr)
  library(readr)
})

base_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation"
source(file.path(base_dir, "src/R/Script_Figures/figure_03_segmentation_panels/figure_03_segmentation_panels_common.R"))
source(file.path(base_dir, "src/R/utils/RegionPhotometry.R"))

tag <- "sagui11"
mode <- "psfmatched"
ncomp <- 20L
mask_mode <- Sys.getenv("SAGUI11_MASK_MODE", unset = "copula_bridge")
figure_id <- "figure_03_segmentation_panels"
slug <- Sys.getenv("SAGUI11_OUTPUT_SLUG", unset = "psfmatched_n20")

figure_out_dir <- file.path(base_dir, "results/figures/paper_repro", figure_id, slug)
png_dir <- file.path(figure_out_dir, "png")
seg_out_dir <- file.path(base_dir, "results/segmentation/paper_repro", figure_id, slug)
flux_out_dir <- file.path(base_dir, "results/flux_per_region/paper_repro", figure_id, slug)
dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(seg_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(flux_out_dir, recursive = TRUE, showWarnings = FALSE)

parse_fits_scalar <- function(line) {
  quoted <- regmatches(line, regexpr("'[^']*'", line))
  if (length(quoted) == 1L && nzchar(quoted)) {
    return(trimws(gsub("^'|'$", "", quoted)))
  }
  rhs <- sub("^[^=]+=\\s*", "", line)
  rhs <- sub("\\s*/.*$", "", rhs)
  trimws(rhs)
}

extract_header_filters <- function(fits_obj) {
  filter_lines <- fits_obj$header[grepl("^FILTER[0-9]+\\s*=", fits_obj$header)]
  vals <- vapply(filter_lines, parse_fits_scalar, character(1))
  unname(vals[nzchar(vals)])
}

write_region_photometry_csv_dynamic <- function(fits_obj, cluster_map, path, error_fallback = c("mad_sky", "flux_over_sqrt_n")) {
  error_fallback <- match.arg(error_fallback)
  filters <- extract_header_filters(fits_obj)
  sed <- RegionPhotometry(
    fits_obj,
    cluster_map,
    band_values = filters,
    error_fallback = error_fallback
  )

  err_cols <- paste0(filters, "_err")
  neff_cols <- paste0(filters, "_n_eff")

  flux_ordered_jy <- sed$flux_wide |>
    dplyr::select(
      dplyr::any_of(c("region", "n_pix")),
      dplyr::all_of(filters),
      dplyr::any_of(err_cols),
      dplyr::any_of(neff_cols)
    ) |>
    dplyr::mutate(
      dplyr::across(dplyr::all_of(filters), ~ .x * 1e-8),
      dplyr::across(dplyr::any_of(err_cols), ~ .x * 1e-8)
    )

  readr::write_csv(flux_ordered_jy, path)
  invisible(path)
}

keep_large_components <- function(mask, min_area = 120L) {
  labs <- EBImage::bwlabel(EBImage::Image(mask * 1))
  lab_mat <- image_to_matrix(labs)
  ids <- sort(unique(as.integer(lab_mat[lab_mat > 0])))
  if (!length(ids)) {
    return(mask)
  }

  keep_ids <- ids[vapply(
    ids,
    function(id) sum(lab_mat == id, na.rm = TRUE) >= min_area,
    logical(1)
  )]

  matrix(lab_mat %in% keep_ids, nrow = nrow(lab_mat), ncol = ncol(lab_mat))
}

fill_enclosed_holes <- function(mask) {
  image_to_matrix(EBImage::fillHull(EBImage::Image(mask * 1))) > 0
}

build_mask_for_mode <- function(cube, mask_mode) {
  if (identical(mask_mode, "starlet")) {
    info <- build_starlet_mask(cube)
    info$variant_note <- "starlet_mask_on_flux_cube"
    return(info)
  }

  if (identical(mask_mode, "copula")) {
    info <- build_copula_mask(cube)
    info$variant_note <- "copula_default_mask_on_flux_cube"
    return(info)
  }

  if (identical(mask_mode, "copula_relaxed")) {
    info <- build_copula_mask(cube, keep_scales = 2:5, denoise_k = 1L)
    info$variant_note <- "copula_relaxed_mask_on_flux_cube"
    return(info)
  }

  if (identical(mask_mode, "copula_bridge")) {
    allowed_info <- build_copula_mask(cube, keep_scales = 2:5, denoise_k = 1L)
    allowed_mask <- close_open_mask(allowed_info$mask, close_size = 5L, open_size = 1L)
    seed_mask <- keep_large_components(allowed_info$mask, min_area = 120L)
    bridge_mask <- grow_mask_constrained(seed_mask, allowed_mask, max_iter = 24L)
    bridge_mask <- fill_enclosed_holes(bridge_mask)
    allowed_info$mask_seed <- seed_mask
    allowed_info$mask_allowed <- allowed_mask
    allowed_info$mask <- bridge_mask
    allowed_info$variant_note <- "copula_relaxed_bridge_holefill_mask_on_flux_cube"
    return(allowed_info)
  }

  stop("Unsupported mask mode: ", mask_mode)
}

cube_path <- resolve_input_cube(tag, mode)
raw_cube_path <- resolve_raw_cube(tag)
X <- FITSio::readFITS(cube_path)
X_raw <- FITSio::readFITS(raw_cube_path)
cube_seg <- X$imDat

mask_info <- build_mask_for_mode(cube_seg, mask_mode)
seg <- build_segmentation(cube_seg, mask_info$mask, N = ncomp)

png_path <- file.path(png_dir, paste0(tag, ".png"))
fits_path <- file.path(seg_out_dir, paste0(tag, ".fits"))
save_segmentation_png(seg, png_path, N = ncomp)
write_cluster_fits(seg, fits_path)

if (all(dim(seg$cluster_map) == dim(X_raw$imDat)[1:2])) {
  cluster_map_for_photometry <- seg$cluster_map
  match_score <- NA_real_
  match_row <- NA_integer_
  match_col <- NA_integer_
  template_width <- ncol(cluster_map_for_photometry)
  template_height <- nrow(cluster_map_for_photometry)
} else {
  remap <- remap_psf_labels_to_raw(seg$cluster_map, X_raw, X)
  cluster_map_for_photometry <- remap$labels_raw
  match_score <- remap$match_score
  match_row <- remap$match_row
  match_col <- remap$match_col
  template_width <- remap$template_width
  template_height <- remap$template_height
}

mad_path <- file.path(flux_out_dir, paste0("SED_flux_wide_", tag, "_mad_sky.csv"))
sn_path <- file.path(flux_out_dir, paste0("SED_flux_wide_", tag, "_flux_over_sqrt_n.csv"))

write_region_photometry_csv_dynamic(X_raw, cluster_map_for_photometry, mad_path, error_fallback = "mad_sky")
write_region_photometry_csv_dynamic(X_raw, cluster_map_for_photometry, sn_path, error_fallback = "flux_over_sqrt_n")

manifest_path <- file.path(figure_out_dir, "manifest_sagui11.csv")
manifest_row <- data.frame(
  tag = tag,
  mode = mode,
  mask_mode = mask_mode,
  variant_note = if (!is.null(mask_info$variant_note)) mask_info$variant_note else NA_character_,
  ncomp = ncomp,
  source_cube = cube_path,
  photometry_cube = raw_cube_path,
  png_path = png_path,
  fits_path = fits_path,
  mask_pixels = sum(mask_info$mask, na.rm = TRUE),
  mask_seed_pixels = if (!is.null(mask_info$mask_seed)) sum(mask_info$mask_seed, na.rm = TRUE) else NA_integer_,
  mask_allowed_pixels = if (!is.null(mask_info$mask_allowed)) sum(mask_info$mask_allowed, na.rm = TRUE) else NA_integer_,
  labeled_pixels = sum(is.finite(seg$cluster_map) & seg$cluster_map > 0),
  registration_score = match_score,
  registration_row = match_row,
  registration_col = match_col,
  registration_template_width = template_width,
  registration_template_height = template_height,
  stringsAsFactors = FALSE
)
readr::write_csv(manifest_row, manifest_path)

message("Created Sagui-11 PSF-matched N=20 outputs")
message("Mask mode: ", mask_mode)
message("PNG: ", png_path)
message("Segmentation FITS: ", fits_path)
message("Flux CSVs: ", mad_path, " | ", sn_path)
