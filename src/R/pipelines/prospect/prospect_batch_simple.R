#!/usr/bin/env Rscript

# Batch version of the simple ProSpect learning setup.
# It uses the same conservative model as prospect_learn_one_region.R, but fits
# every region in the selected segmentation table and writes one summary CSV.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
})

source("/Users/rd23aag/Documents/GitHub/crp8_segmentation/src/R/pipelines/prospect/prospect_learning_helpers.R")

tag <- Sys.getenv("SAGUI_PROSPECT_TAG", unset = "sagui10")
photometry_mode <- Sys.getenv("SAGUI_PROSPECT_PHOTOMETRY_MODE", unset = "raw")
systematic_frac <- as.numeric(Sys.getenv("SAGUI_PROSPECT_SYSTEMATIC_FRAC", unset = "0.05"))
fixed_Z <- as.numeric(Sys.getenv("SAGUI_PROSPECT_FIXED_Z_METAL", unset = "0.02"))
fixed_z_override <- Sys.getenv("SAGUI_PROSPECT_FIXED_REDSHIFT", unset = "")
free_Z <- tolower(Sys.getenv("SAGUI_PROSPECT_FREE_Z", unset = "0")) %in% c("1", "true", "yes")
logzsol_lo <- as.numeric(Sys.getenv("SAGUI_PROSPECT_LOGZSOL_LO", unset = "-2.0"))
logzsol_hi <- as.numeric(Sys.getenv("SAGUI_PROSPECT_LOGZSOL_HI", unset = "0.3"))
region_csv_override <- Sys.getenv("SAGUI_PROSPECT_REGION_CSV", unset = "")
region_limit <- as.integer(Sys.getenv("SAGUI_PROSPECT_REGION_LIMIT", unset = "0"))
run_suffix <- Sys.getenv("SAGUI_PROSPECT_RUN_SUFFIX", unset = "")

region_df <- read_region_table_override(tag, photometry_mode, if (nzchar(region_csv_override)) region_csv_override else NULL)
regions <- sort(unique(region_df$region))
if (is.finite(region_limit) && region_limit > 0) {
  regions <- head(regions, region_limit)
}

out_dir <- file.path(
  prospect_base_dir,
  "results/sed_fitting/prospect_learning",
  tag,
  photometry_mode,
  if (nzchar(run_suffix)) paste0("batch_simple_", run_suffix) else "batch_simple"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fit_region_safe <- function(region_id) {
  message("Fitting region ", region_id)
  tryCatch(
    fit_simple_prospect_region(
      tag = tag,
      region_id = region_id,
      photometry_mode = photometry_mode,
      systematic_frac = systematic_frac,
      fixed_z = if (nzchar(fixed_z_override)) as.numeric(fixed_z_override) else NULL,
      fixed_Z = fixed_Z,
      free_Z = free_Z,
      logzsol_bounds = c(logzsol_lo, logzsol_hi),
      region_csv_override = if (nzchar(region_csv_override)) region_csv_override else NULL
    )$summary |>
      dplyr::mutate(Av = 1.086 * .data$tau_screen),
    error = function(e) {
      tibble::tibble(
        tag = tag,
        region = region_id,
        n_pix = region_df$n_pix[region_df$region == region_id][[1]],
        photometry_mode = photometry_mode,
        z_used = if (nzchar(fixed_z_override)) as.numeric(fixed_z_override) else lookup_redshift(tag, fallback = NA_real_),
        Z_mode = if (isTRUE(free_Z)) "free" else "fixed",
        Z_fixed = fixed_Z,
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
        Av = NA_real_,
        fit_error = conditionMessage(e)
      )
    }
  )
}

results <- purrr::map_dfr(regions, fit_region_safe)

summary_path <- file.path(out_dir, paste0(tag, "_batch_simple_summary.csv"))
readr::write_csv(results, summary_path)

cat("Wrote:\n", summary_path, "\n")
