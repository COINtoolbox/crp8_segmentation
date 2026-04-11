#!/usr/bin/env Rscript

# A small, readable ProSpect example for learning.
# It fits one segmented region, saves a diagnostic plot, and writes a compact
# summary table so you can inspect what the code is actually doing.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

source("/Users/rd23aag/Documents/GitHub/crp8_segmentation/src/R/pipelines/prospect/prospect_learning_helpers.R")

tag <- Sys.getenv("SAGUI_PROSPECT_TAG", unset = "sagui10")
region_id <- as.integer(Sys.getenv("SAGUI_PROSPECT_REGION", unset = "1"))
photometry_mode <- Sys.getenv("SAGUI_PROSPECT_PHOTOMETRY_MODE", unset = "raw")
systematic_frac <- as.numeric(Sys.getenv("SAGUI_PROSPECT_SYSTEMATIC_FRAC", unset = "0.05"))
fixed_Z <- as.numeric(Sys.getenv("SAGUI_PROSPECT_FIXED_Z_METAL", unset = "0.02"))
fixed_z_override <- Sys.getenv("SAGUI_PROSPECT_FIXED_REDSHIFT", unset = "")
free_Z <- tolower(Sys.getenv("SAGUI_PROSPECT_FREE_Z", unset = "0")) %in% c("1", "true", "yes")
logzsol_lo <- as.numeric(Sys.getenv("SAGUI_PROSPECT_LOGZSOL_LO", unset = "-2.0"))
logzsol_hi <- as.numeric(Sys.getenv("SAGUI_PROSPECT_LOGZSOL_HI", unset = "0.3"))
region_csv_override <- Sys.getenv("SAGUI_PROSPECT_REGION_CSV", unset = "")

out_dir <- file.path(
  prospect_base_dir,
  "results/sed_fitting/prospect_learning",
  tag,
  photometry_mode,
  paste0("region_", sprintf("%02d", region_id))
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fit_result <- fit_simple_prospect_region(
  tag = tag,
  region_id = region_id,
  photometry_mode = photometry_mode,
  systematic_frac = systematic_frac,
  fixed_z = if (nzchar(fixed_z_override)) as.numeric(fixed_z_override) else NULL,
  fixed_Z = fixed_Z,
  free_Z = free_Z,
  logzsol_bounds = c(logzsol_lo, logzsol_hi),
  region_csv_override = if (nzchar(region_csv_override)) region_csv_override else NULL
)

summary_tbl <- fit_result$summary |>
  dplyr::mutate(Av = 1.086 * .data$tau_screen)
plot_obj <- plot_simple_prospect_region(fit_result)

png_path <- file.path(out_dir, paste0(tag, "_region_", sprintf("%02d", region_id), "_simple_fit.png"))
pdf_path <- file.path(out_dir, paste0(tag, "_region_", sprintf("%02d", region_id), "_simple_fit.pdf"))
summary_path <- file.path(out_dir, paste0(tag, "_region_", sprintf("%02d", region_id), "_simple_summary.csv"))
photom_path <- file.path(out_dir, paste0(tag, "_region_", sprintf("%02d", region_id), "_model_photometry.csv"))

obs_df <- tibble::as_tibble(fit_result$data$flux) |>
  rename(obs_flux = "flux", obs_fluxerr = "fluxerr")
model_df <- tibble::tibble(
  filter = obs_df$filter,
  cenwave = obs_df$cenwave,
  model_flux = fit_result$checked$SEDout$Photom
)

photom_df <- left_join(obs_df, model_df, by = c("filter", "cenwave"))

readr::write_csv(summary_tbl, summary_path)
readr::write_csv(photom_df, photom_path)

grDevices::png(png_path, width = 1800, height = 1250, res = 220, bg = "white")
print(plot_obj)
grDevices::dev.off()

grDevices::pdf(pdf_path, width = 10.5, height = 7.3, bg = "white")
print(plot_obj)
grDevices::dev.off()

cat("Wrote:\n", summary_path, "\n", photom_path, "\n", png_path, "\n", pdf_path, "\n")
