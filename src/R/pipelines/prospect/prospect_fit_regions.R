#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ProSpect)
  library(ProSpectData)
  library(readr)
  library(dplyr)
  library(purrr)
})

source("/Users/rd23aag/Documents/GitHub/crp8_segmentation/src/R/pipelines/prospect/prospect_common.R")

tag <- Sys.getenv("SAGUI_PROSPECT_TAG", unset = "sagui10")
photometry_mode <- Sys.getenv("SAGUI_PROSPECT_PHOTOMETRY_MODE", unset = "raw")
photometry_mode <- match.arg(photometry_mode, c("raw", "psfmatched"))
current_results_override <- Sys.getenv("SAGUI_PROSPECT_CURRENT_RESULTS_CSV", unset = "")
out_path <- file.path(
  prospect_base_dir,
  "results/sed_fitting/prospect",
  tag,
  photometry_mode,
  paste0(tag, "_prospect_results.csv")
)

input_df <- readr::read_csv(
  file.path(prospect_base_dir, "results/sed_fitting/prospect", tag, photometry_mode, paste0(tag, "_prospect_input.csv")),
  show_col_types = FALSE
)
default_current_path <- file.path(prospect_base_dir, "results/sed_fitting", paste0(tag, "_sed_results_new.csv"))
current_df <- if (nzchar(current_results_override) && file.exists(current_results_override)) {
  readr::read_csv(current_results_override, show_col_types = FALSE)
} else if (file.exists(default_current_path)) {
  readr::read_csv(default_current_path, show_col_types = FALSE)
} else {
  tibble::tibble()
}

filters <- infer_filters_from_table(input_df)
cenwaves <- vapply(
  filters,
  function(filter_name) {
    obj <- get(paste0("filt_", filter_name, "_JWST"))
    sum(obj$wave * obj$response, na.rm = TRUE) / sum(obj$response, na.rm = TRUE)
  },
  numeric(1)
)
filtout <- lapply(filters, function(ff) get(paste0("filt_", ff, "_JWST")))

make_data <- function(region_id) {
  row <- input_df |> filter(.data$region == region_id)
  flux_df <- data.frame(
    filter = filters,
    cenwave = cenwaves,
    flux = as.numeric(row[1, filters]),
    fluxerr = pmax(as.numeric(row[1, paste0(filters, "_err")]), 1e-12)
  )
  list(
    flux = flux_df,
    SFH = SFHfunc,
    speclib = BC03lr,
    Dale = FALSE,
    filtout = filtout,
    fit = "optim",
    like = "st",
    parm.names = c("mSFR", "mpeak", "mperiod", "mskew", "tau_birth", "tau_screen", "Z"),
    logged = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
    intervals = list(
      lo = c(-4, 0.05, 0.10, -3.0, 0.00, 0.00, log10(1e-4)),
      hi = c(4, 10.0, 6.00, 3.0, 6.00, 6.00, log10(0.05))
    ),
    mon.names = c("masstot", "SFRburst"),
    verbose = FALSE,
    arglist = list(
      z = 0.42,
      ref = "Planck18",
      massfunc = massfunc_snorm
    )
  )
}

make_current_start <- function(region_id) {
  if (!nrow(current_df) || !"region" %in% names(current_df)) {
    return(fallback_start)
  }
  cur <- current_df |> filter(.data$region == region_id) |> slice(1)
  if (!nrow(cur)) {
    return(fallback_start)
  }
  c(
    log10(pmax(cur$sfr_100myr, 1e-4)),
    pmin(pmax(cur$tage, 0.1), 8.5),
    pmin(pmax(cur$tau, 0.15), 5.0),
    0,
    pmin(pmax(cur$dust2 * 1.2, 0.0), 5.5),
    pmin(pmax(cur$dust2 * 0.7, 0.0), 5.5),
    log10(pmax(0.02 * 10^cur$logzsol, 1e-4))
  )
}

fallback_start <- c(log10(0.10), 2.5, 1.0, 0.0, 1.0, 0.5, log10(0.010))

age_mass_weighted_gyr <- function(stars) {
  sum(stars$agevec * stars$massvec, na.rm = TRUE) / sum(stars$massvec, na.rm = TRUE) / 1e9
}

recent_sfr_from_sfh <- function(stars, window_yr = 1e8) {
  idx <- stars$agevec <= window_yr
  if (!any(idx, na.rm = TRUE)) return(stars$SFRburst)
  stats::weighted.mean(stars$SFR[idx], w = pmax(stars$massvec[idx], 0), na.rm = TRUE)
}

fit_once <- function(data_obj, start) {
  objective_fn <- function(par) ProSpectSEDlike(par, data_obj)
  tryCatch(
    nlminb(
      start = start,
      objective = objective_fn,
      lower = data_obj$intervals$lo,
      upper = data_obj$intervals$hi,
      control = list(iter.max = 180, eval.max = 450)
    ),
    error = function(e) NULL
  )
}

fit_region <- function(region_id) {
  message("Fitting region ", region_id)
  data_obj <- make_data(region_id)
  starts <- list(make_current_start(region_id), fallback_start)
  fits <- purrr::map(starts, ~ fit_once(data_obj, .x))
  fits <- purrr::compact(fits)
  if (!length(fits)) stop("All fits failed for region ", region_id)

  score <- function(x) if (is.finite(x$objective)) x$objective else Inf
  best <- fits[[which.min(vapply(fits, score, numeric(1)))]]

  data_obj$fit <- "check"
  checked <- ProSpectSEDlike(best$par, data_obj)
  stars <- checked$SEDout$Stars
  z_linear <- 10^best$par[7]
  logm_formed <- log10(pmax(stars$masstot, .Machine$double.eps))
  sfr_recent <- recent_sfr_from_sfh(stars)

  tibble::tibble(
    region = region_id,
    n_pix = input_df$n_pix[input_df$region == region_id][[1]],
    convergence = best$convergence,
    objective = best$objective,
    logMformed = logm_formed,
    logMstar = logm_formed,
    logzsol = log10(pmax(z_linear, 1e-8) / 0.02),
    Av = 1.086 * best$par[6],
    age_mw = age_mass_weighted_gyr(stars),
    sfr = sfr_recent,
    ssfr = sfr_recent / pmax(stars$masstot, .Machine$double.eps),
    mSFR = 10^best$par[1],
    mpeak = best$par[2],
    mperiod = best$par[3],
    mskew = best$par[4],
    tau_birth = best$par[5],
    tau_screen = best$par[6],
    Z = z_linear,
    z_used = 0.42
  )
}

regions <- sort(unique(input_df$region))
results <- vector("list", length(regions))

for (i in seq_along(regions)) {
  results[[i]] <- fit_region(regions[[i]])
  partial <- dplyr::bind_rows(results[seq_len(i)])
  readr::write_csv(partial, out_path)
}

final <- dplyr::bind_rows(results)
readr::write_csv(final, out_path)
cat("Wrote:", out_path, "\n")
