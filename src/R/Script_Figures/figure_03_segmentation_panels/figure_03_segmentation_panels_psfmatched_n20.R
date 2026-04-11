#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(FITSio)
  library(guara)
  library(capivara)
  library(ragg)
  library(dplyr)
})

base_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation"
source(file.path(base_dir, "src/R/Script_Figures/figure_03_segmentation_panels/figure_03_segmentation_panels_common.R"))
source(file.path(base_dir, "src/R/utils/RegionPhotometry.R"))

figure_id <- "figure_03_segmentation_panels"
mode <- "psfmatched"
ncomp <- 20L
slug <- "psfmatched_n20"

figure_out_dir <- file.path(base_dir, "results/figures/paper_repro", figure_id, slug)
png_dir <- file.path(figure_out_dir, "png")
seg_out_dir <- file.path(base_dir, "results/segmentation/paper_repro", figure_id, slug)
flux_out_dir <- file.path(base_dir, "results/flux_per_region/paper_repro", figure_id, slug)

dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(seg_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(flux_out_dir, recursive = TRUE, showWarnings = FALSE)

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

write_region_photometry_csv_custom <- function(fits_obj, cluster_map, path, error_fallback = c("mad_sky", "flux_over_sqrt_n")) {
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

manifest <- vector("list", length(panel_tags))

for (i in seq_along(panel_tags)) {
  tag <- panel_tags[i]
  cube_path <- resolve_input_cube(tag, mode)
  X <- FITSio::readFITS(cube_path)
  cube <- X$imDat
  raw_cube_path <- resolve_raw_cube(tag)
  X_raw <- FITSio::readFITS(raw_cube_path)

  mask_info <- build_starlet_mask(cube)
  seg <- build_segmentation(cube, mask_info$mask, N = ncomp)

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

  write_region_photometry_csv_custom(
    X_raw,
    cluster_map_for_photometry,
    file.path(flux_out_dir, paste0("SED_flux_wide_", tag, "_mad_sky.csv")),
    error_fallback = "mad_sky"
  )
  write_region_photometry_csv_custom(
    X_raw,
    cluster_map_for_photometry,
    file.path(flux_out_dir, paste0("SED_flux_wide_", tag, "_flux_over_sqrt_n.csv")),
    error_fallback = "flux_over_sqrt_n"
  )

  manifest[[i]] <- data.frame(
    tag = tag,
    mode = mode,
    ncomp = ncomp,
    source_cube = cube_path,
    photometry_cube = raw_cube_path,
    png_path = png_path,
    fits_path = fits_path,
    mask_pixels = sum(mask_info$mask, na.rm = TRUE),
    labeled_pixels = sum(is.finite(seg$cluster_map) & seg$cluster_map > 0),
    registration_score = match_score,
    registration_row = match_row,
    registration_col = match_col,
    registration_template_width = template_width,
    registration_template_height = template_height,
    stringsAsFactors = FALSE
  )
}

manifest_df <- do.call(rbind, manifest)
utils::write.csv(manifest_df, file.path(figure_out_dir, "manifest.csv"), row.names = FALSE)
message("Finished ", figure_id, " variant: ", slug)
message("Figures: ", png_dir)
message("Segmentation FITS: ", seg_out_dir)
message("Flux tables: ", flux_out_dir)
