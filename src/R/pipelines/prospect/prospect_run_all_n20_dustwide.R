#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(glue)
})

base_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation"
batch_script <- file.path(base_dir, "src/R/pipelines/prospect/prospect_batch_simple.R")
paint_script <- file.path(base_dir, "src/R/pipelines/prospect/prospect_paint_maps.R")

run_suffix <- Sys.getenv("SAGUI_PROSPECT_RUN_SUFFIX", unset = "n20_freez_sys10_dustwide_madsky")
style_suffix <- Sys.getenv("SAGUI_PROSPECT_STYLE_SUFFIX", unset = "legacy_cosmic_n20_freez_sys10_dustwide_madsky")
systematic_frac <- Sys.getenv("SAGUI_PROSPECT_SYSTEMATIC_FRAC", unset = "0.10")
tau_screen_hi <- Sys.getenv("SAGUI_PROSPECT_TAU_SCREEN_HI", unset = "5.5")
do_paint <- tolower(Sys.getenv("SAGUI_PROSPECT_DO_PAINT", unset = "1")) %in% c("1", "true", "yes")

tags <- c("sagui4", "sagui5_6", "sagui7", "sagui8", "sagui9", "sagui10", "sagui11")
fixed_redshift <- c(sagui11 = "1.1")

run_one <- function(tag) {
  region_csv <- file.path(base_dir, "results/flux_per_region", glue("SED_flux_wide_{tag}.csv"))
  segmentation_fits <- file.path(
    base_dir,
    "results/segmentation/paper_repro/figure_03_segmentation_panels/psfmatched_n20",
    glue("{tag}.fits")
  )
  if (!file.exists(region_csv)) {
    stop("Missing region CSV: ", region_csv)
  }
  if (!file.exists(segmentation_fits)) {
    stop("Missing segmentation FITS: ", segmentation_fits)
  }

  env <- c(
    glue("SAGUI_PROSPECT_TAG={tag}"),
    "SAGUI_PROSPECT_PHOTOMETRY_MODE=raw",
    glue("SAGUI_PROSPECT_REGION_CSV={region_csv}"),
    glue("SAGUI_PROSPECT_SYSTEMATIC_FRAC={systematic_frac}"),
    "SAGUI_PROSPECT_FREE_Z=1",
    glue("SAGUI_PROSPECT_TAU_SCREEN_HI={tau_screen_hi}"),
    "SAGUI_PROSPECT_LOGZSOL_LO=-2.0",
    "SAGUI_PROSPECT_LOGZSOL_HI=0.3",
    glue("SAGUI_PROSPECT_RUN_SUFFIX={run_suffix}")
  )
  if (tag %in% names(fixed_redshift)) {
    env <- c(env, glue("SAGUI_PROSPECT_FIXED_REDSHIFT={fixed_redshift[[tag]]}"))
  }

  message("Running batch fit for ", tag)
  status <- system2("Rscript", batch_script, env = env)
  if (!identical(status, 0L)) {
    stop("Batch fit failed for ", tag)
  }

  if (do_paint) {
    results_csv <- file.path(
      base_dir,
      "results/sed_fitting/prospect_learning",
      tag,
      "raw",
      glue("batch_simple_{run_suffix}"),
      glue("{tag}_batch_simple_summary.csv")
    )
    paint_env <- c(
      glue("SAGUI_PROSPECT_TAG={tag}"),
      "SAGUI_PROSPECT_PHOTOMETRY_MODE=raw",
      glue("SAGUI_PROSPECT_RESULTS_CSV={results_csv}"),
      glue("SAGUI_PROSPECT_SEGMENTATION_FITS={segmentation_fits}"),
      glue("SAGUI_PROSPECT_STYLE_SUFFIX={style_suffix}"),
      "SAGUI_PROSPECT_PALETTE=legacy_cosmic_dusk"
    )
    message("Painting maps for ", tag)
    status <- system2("Rscript", paint_script, env = paint_env)
    if (!identical(status, 0L)) {
      stop("Map painting failed for ", tag)
    }
  }
}

for (tag in tags) {
  run_one(tag)
}

message("Completed ProSpect N=20 dust-wide run for: ", paste(tags, collapse = ", "))
