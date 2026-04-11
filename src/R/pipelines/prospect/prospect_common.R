#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ProSpect)
  library(ProSpectData)
})

prospect_base_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation"

resolve_region_tag <- function(tag) {
  lookup <- c(
    sagui1_2_3 = "sagui1_2_3",
    sagui4 = "sagui4",
    sagui5_6 = "sagui5_6",
    sagui7 = "sagui7",
    sagui8 = "sagui8",
    sagui9 = "sagui9",
    sagui10 = "sagui10",
    sagui11 = "sagui11"
  )
  if (!tag %in% names(lookup)) {
    stop("Unknown tag: ", tag)
  }
  lookup[[tag]]
}

lookup_redshift <- function(tag, fallback = NULL) {
  z_lookup <- c(
    sagui1 = 0.7649,
    sagui2 = 1.0385,
    sagui3 = 0.9965,
    sagui4 = 0.3664,
    sagui5 = 0.3583,
    sagui6 = 0.3583,
    sagui7 = 1.0878,
    sagui8 = 0.6222,
    sagui9 = 0.6651,
    sagui10 = 0.4148
  )

  if (tag %in% names(z_lookup)) {
    return(unname(z_lookup[[tag]]))
  }

  if (identical(tag, "sagui5_6")) {
    return(0.3583)
  }

  if (identical(tag, "sagui1_2_3")) {
    if (!is.null(fallback)) return(fallback)
    stop("No single redshift is defined for tag '", tag, "'. Use an explicit override.")
  }

  if (!is.null(fallback)) return(fallback)
  stop("No redshift lookup available for tag: ", tag)
}

raw_cube_name <- function(tag) {
  dashed <- sub("^sagui", "sagui-", gsub("_", "-", tag))
  file.path(prospect_base_dir, "data/raw", paste0("datacube_", dashed, ".fits"))
}

psf_cube_name <- function(tag) {
  dashed <- sub("^sagui", "sagui-", gsub("_", "-", tag))
  file.path(prospect_base_dir, "data/PSF_matched", paste0("datacube_", dashed, "_psfmatched.fits"))
}

segmentation_fits_path <- function(tag) {
  file.path(
    prospect_base_dir,
    "results/segmentation/paper_repro/figure_03_segmentation_panels/psfmatched",
    paste0(tag, ".fits")
  )
}

region_csv_path <- function(tag) {
  file.path(
    prospect_base_dir,
    "results/flux_per_region/paper_repro/figure_03_segmentation_panels/psfmatched",
    paste0("SED_flux_wide_", tag, ".csv")
  )
}

current_results_csv_path <- function(tag) {
  candidates <- c(
    file.path(prospect_base_dir, "results/sed_fitting", paste0(tag, "_sed_results_new.csv")),
    file.path(prospect_base_dir, "results/sed_fitting", paste0(tag, "_sed_results.csv")),
    file.path(prospect_base_dir, "results/sed_fitting", paste0(sub("^sagui", "sagui-", tag), "_sed_results_new.csv")),
    file.path(prospect_base_dir, "results/sed_fitting", paste0(sub("^sagui", "sagui-", tag), "_sed_results.csv")),
    file.path(prospect_base_dir, "results/sed_fitting", paste0(gsub("_", "-", sub("^sagui", "sagui-", tag)), "_sed_results_new.csv")),
    file.path(prospect_base_dir, "results/sed_fitting", paste0(gsub("_", "-", sub("^sagui", "sagui-", tag)), "_sed_results.csv"))
  )
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) return(NA_character_)
  hits[[1]]
}

known_filters <- c(
  "F090W", "F115W", "F150W", "F182M", "F200W", "F210M",
  "F277W", "F335M", "F356W", "F410M", "F430M", "F444W",
  "F460M", "F480M"
)

infer_filters_from_table <- function(df) {
  known_filters[known_filters %in% names(df)]
}

prospect_filter_symbols <- function(filters) {
  paste0("filt_", filters, "_JWST")
}

make_output_dirs <- function(tag, photometry_mode = c("raw", "psfmatched")) {
  photometry_mode <- match.arg(photometry_mode)
  list(
    fitting = file.path(prospect_base_dir, "results/sed_fitting/prospect", tag, photometry_mode),
    maps = file.path(prospect_base_dir, "results/sed_maps/prospect", tag, photometry_mode)
  )
}

effective_region_csv_path <- function(tag, photometry_mode = c("raw", "psfmatched")) {
  photometry_mode <- match.arg(photometry_mode)
  if (identical(photometry_mode, "raw")) {
    region_csv_path(tag)
  } else {
    file.path(
      prospect_base_dir,
      "results/flux_per_region/paper_repro/figure_03_segmentation_panels/psfmatched",
      paste0("SED_flux_wide_", tag, "_psfmatched.csv")
    )
  }
}

read_region_table <- function(tag) {
  path <- region_csv_path(tag)
  if (!file.exists(path)) {
    stop("Region CSV not found: ", path)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

read_region_table_mode <- function(tag, photometry_mode = c("raw", "psfmatched")) {
  photometry_mode <- match.arg(photometry_mode)
  path <- effective_region_csv_path(tag, photometry_mode)
  if (!file.exists(path)) {
    if (identical(photometry_mode, "psfmatched")) {
      stop("PSF-matched region CSV not found: ", path, ". Generate it first or use photometry_mode='raw'.")
    }
    stop("Region CSV not found: ", path)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

read_region_table_override <- function(tag,
                                       photometry_mode = c("raw", "psfmatched"),
                                       region_csv_override = NULL) {
  photometry_mode <- match.arg(photometry_mode)
  if (!is.null(region_csv_override) && nzchar(region_csv_override)) {
    if (!file.exists(region_csv_override)) {
      stop("Region CSV override not found: ", region_csv_override)
    }
    return(readr::read_csv(region_csv_override, show_col_types = FALSE))
  }
  read_region_table_mode(tag, photometry_mode)
}

prospect_filter_objects <- function(filters) {
  lapply(filters, function(ff) get(paste0("filt_", ff, "_JWST")))
}

prospect_cenwaves <- function(filters) {
  vapply(
    filters,
    function(filter_name) {
      obj <- get(paste0("filt_", filter_name, "_JWST"))
      sum(obj$wave * obj$response, na.rm = TRUE) / sum(obj$response, na.rm = TRUE)
    },
    numeric(1)
  )
}

region_flux_table <- function(tag, region_id, photometry_mode = c("raw", "psfmatched"), systematic_frac = 0.05) {
  photometry_mode <- match.arg(photometry_mode)
  region_df <- read_region_table_mode(tag, photometry_mode)
  filters <- infer_filters_from_table(region_df)
  err_cols <- paste0(filters, "_err")
  row <- region_df |> filter(.data$region == region_id) |> slice(1)
  if (!nrow(row)) {
    stop("Region ", region_id, " not found in ", effective_region_csv_path(tag, photometry_mode))
  }

  flux <- as.numeric(row[1, filters])
  fluxerr <- as.numeric(row[1, err_cols])
  fluxerr_eff <- sqrt(pmax(fluxerr, 0)^2 + (systematic_frac * pmax(flux, 0))^2)

  tibble::tibble(
    filter = filters,
    cenwave = prospect_cenwaves(filters),
    flux = flux,
    fluxerr = pmax(fluxerr_eff, 1e-30)
  )
}

simple_current_start <- function(tag, region_id) {
  current_path <- current_results_csv_path(tag)
  if (!is.character(current_path) || is.na(current_path) || !file.exists(current_path)) {
    return(c(log10(0.10), 2.5, 1.0, 0.5))
  }
  cur <- readr::read_csv(current_path, show_col_types = FALSE) |>
    filter(.data$region == region_id) |>
    slice(1)
  if (!nrow(cur)) {
    return(c(log10(0.10), 2.5, 1.0, 0.5))
  }
  c(
    log10(pmax(cur$sfr_100myr, 1e-4)),
    pmin(pmax(cur$tage, 0.2), 8.5),
    pmin(pmax(cur$tau, 0.2), 5.0),
    pmin(pmax(cur$dust2 * 0.9, 0.0), 2.5)
  )
}

write_manifest <- function(tag, photometry_mode = c("raw", "psfmatched")) {
  photometry_mode <- match.arg(photometry_mode)
  out_dirs <- make_output_dirs(tag, photometry_mode)
  dir.create(out_dirs$fitting, recursive = TRUE, showWarnings = FALSE)
  df <- read_region_table(tag)
  filters <- infer_filters_from_table(df)
  filters_str <- paste(filters, collapse = ",")
  prospect_filters_str <- paste(prospect_filter_symbols(filters), collapse = ",")

  manifest <- tibble::tibble(
    tag = tag,
    photometry_mode = photometry_mode,
    z_used = lookup_redshift(tag, fallback = NA_real_),
    segmentation_fits = segmentation_fits_path(tag),
    photometry_cube = if (photometry_mode == "raw") raw_cube_name(tag) else psf_cube_name(tag),
    region_csv = region_csv_path(tag),
    n_regions = nrow(df),
    n_filters = length(filters),
    filters = filters_str,
    prospect_filters = prospect_filters_str
  )

  readr::write_csv(manifest, file.path(out_dirs$fitting, "manifest.csv"))
  invisible(manifest)
}
