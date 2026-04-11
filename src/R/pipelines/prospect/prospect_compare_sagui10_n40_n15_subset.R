#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(FITSio)
  library(dplyr)
  library(purrr)
  library(readr)
  library(ggplot2)
})

base_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation"

source(file.path(base_dir, "src/R/pipelines/prospect/prospect_learning_helpers.R"))
source(file.path(base_dir, "src/R/utils/RegionPhotometry.R"))

tag <- "sagui10"
systematic_frac <- 0.05
fixed_Z <- 0.02
free_Z <- FALSE
error_fallback <- "flux_over_sqrt_n"
subset_k <- 10L

out_dir <- file.path(
  base_dir,
  "results/sed_fitting/prospect_compare",
  "sagui10_n40_vs_n15_flux_over_sqrt_n_subset10"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

write_region_photometry_csv_custom <- function(fits_obj, cluster_map, path, error_fallback = c("flux_over_sqrt_n", "mad_sky")) {
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

  readr::write_csv(flux_ordered_jy, path)
  invisible(path)
}

select_representative_regions <- function(region_df, k = 10L) {
  region_sizes <- region_df |>
    distinct(region, n_pix) |>
    arrange(desc(n_pix), region)
  n <- nrow(region_sizes)
  if (n <= k) {
    return(region_sizes)
  }
  idx <- unique(round(seq(1, n, length.out = k)))
  region_sizes[idx, , drop = FALSE]
}

fit_region_safe <- function(run_label, ncomp, region_csv, region_id) {
  message("Fitting ", run_label, " region ", region_id)
  tryCatch(
    fit_simple_prospect_region(
      tag = tag,
      region_id = region_id,
      photometry_mode = "raw",
      systematic_frac = systematic_frac,
      fixed_Z = fixed_Z,
      free_Z = free_Z,
      region_csv_override = region_csv
    )$summary |>
      mutate(run_label = run_label, ncomp = ncomp),
    error = function(e) {
      tibble(
        tag = tag,
        region = region_id,
        n_pix = NA_integer_,
        photometry_mode = "raw",
        z_used = lookup_redshift(tag, fallback = NA_real_),
        Z_mode = if (isTRUE(free_Z)) "free" else "fixed",
        Z_fixed = if (isTRUE(free_Z)) NA_real_ else fixed_Z,
        Z_fit = NA_real_,
        logzsol = NA_real_,
        systematic_frac = systematic_frac,
        convergence = NA_integer_,
        objective = NA_real_,
        rms_sigma = NA_real_,
        logMformed = NA_real_,
        age_mw_gyr = NA_real_,
        sfr_recent = NA_real_,
        mSFR = NA_real_,
        mpeak = NA_real_,
        mperiod = NA_real_,
        tau_screen = NA_real_,
        run_label = run_label,
        ncomp = ncomp,
        fit_error = conditionMessage(e)
      )
    }
  )
}

raw_cube_path <- resolve_raw_cube(tag)
raw_fits <- FITSio::readFITS(raw_cube_path)

run_defs <- tibble::tribble(
  ~run_label, ~ncomp, ~seg_path,
  "n40", 40L, file.path(base_dir, "results/segmentation/paper_repro/figure_03_segmentation_panels/psfmatched", "sagui10.fits"),
  "n15", 15L, file.path(base_dir, "results/segmentation/paper_repro/figure_03_segmentation_panels/psfmatched_n15", "sagui10.fits")
) |>
  mutate(region_csv = file.path(out_dir, paste0("SED_flux_wide_", tag, "_", run_label, "_", error_fallback, ".csv")))

for (i in seq_len(nrow(run_defs))) {
  seg_map <- FITSio::readFITS(run_defs$seg_path[i])$imDat
  write_region_photometry_csv_custom(
    raw_fits,
    seg_map,
    run_defs$region_csv[i],
    error_fallback = error_fallback
  )
}

selected_regions <- purrr::map_dfr(seq_len(nrow(run_defs)), function(i) {
  region_df <- read_region_table_override(tag, "raw", run_defs$region_csv[i])
  sel <- select_representative_regions(region_df, k = subset_k)
  sel$run_label <- run_defs$run_label[i]
  sel$ncomp <- run_defs$ncomp[i]
  sel$region_csv <- run_defs$region_csv[i]
  sel
})

results <- purrr::pmap_dfr(
  list(selected_regions$run_label, selected_regions$ncomp, selected_regions$region_csv, selected_regions$region),
  fit_region_safe
)

results <- results |>
  left_join(selected_regions |> select(run_label, region, sampled_n_pix = n_pix), by = c("run_label", "region"))

summary_tbl <- results |>
  group_by(run_label, ncomp) |>
  summarise(
    n_regions = n(),
    fit_success_frac = mean(is.finite(rms_sigma), na.rm = TRUE),
    convergence_frac = mean(convergence == 0, na.rm = TRUE),
    median_n_pix = median(sampled_n_pix, na.rm = TRUE),
    median_rms_sigma = median(rms_sigma, na.rm = TRUE),
    mean_rms_sigma = mean(rms_sigma, na.rm = TRUE),
    p25_rms_sigma = quantile(rms_sigma, 0.25, na.rm = TRUE),
    p75_rms_sigma = quantile(rms_sigma, 0.75, na.rm = TRUE),
    frac_rms_sigma_lt_2 = mean(rms_sigma < 2, na.rm = TRUE),
    frac_rms_sigma_lt_3 = mean(rms_sigma < 3, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(selected_regions, file.path(out_dir, "selected_regions.csv"))
readr::write_csv(results, file.path(out_dir, "prospect_compare_regions_subset.csv"))
readr::write_csv(summary_tbl, file.path(out_dir, "prospect_compare_summary_subset.csv"))
readr::write_csv(run_defs, file.path(out_dir, "region_csv_manifest.csv"))

p_box <- ggplot(results, aes(x = run_label, y = rms_sigma, fill = run_label)) +
  geom_boxplot(width = 0.58, alpha = 0.70, outlier.shape = NA) +
  geom_jitter(width = 0.08, alpha = 0.55, size = 1.8, colour = "#1D1D1D") +
  scale_fill_manual(values = c(n40 = "#8DA0CB", n15 = "#FC8D62")) +
  labs(
    x = NULL,
    y = "Region fit RMS (sigma units)",
    title = "ProSpect subset comparison: sagui10",
    subtitle = "Representative regions across the n_pix distribution"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "prospect_compare_rms_boxplot_subset.png"), p_box, width = 6.5, height = 4.8, dpi = 220)

print(summary_tbl)
