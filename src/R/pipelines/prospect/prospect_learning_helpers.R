#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ProSpect)
  library(ProSpectData)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

source("/Users/rd23aag/Documents/GitHub/crp8_segmentation/src/R/pipelines/prospect/prospect_common.R")

throughput_nircam_path <- "/Users/rd23aag/Documents/GitHub/sagui/inst/extdata/throughput_nircam.csv"

pretty_angstrom_ticks <- function(x) {
  out <- rep("", length(x))
  tick_map <- c(
    "10000" = "1×10^4",
    "20000" = "2×10^4",
    "30000" = "3×10^4",
    "40000" = "4×10^4"
  )
  key <- as.character(round(x))
  out[key %in% names(tick_map)] <- unname(tick_map[key[key %in% names(tick_map)]])
  out
}

make_simple_prospect_data <- function(tag,
                                      region_id,
                                      photometry_mode = c("raw", "psfmatched"),
                                      systematic_frac = 0.05,
                                      fixed_z = NULL,
                                      fixed_Z = 0.02,
                                      free_Z = FALSE,
                                      logzsol_bounds = c(-2.0, 0.3),
                                      region_csv_override = NULL) {
  photometry_mode <- match.arg(photometry_mode)
  region_df <- read_region_table_override(tag, photometry_mode, region_csv_override)
  filters <- infer_filters_from_table(region_df)
  err_cols <- paste0(filters, "_err")
  row <- region_df |> filter(.data$region == region_id) |> slice(1)
  if (!nrow(row)) {
    stop("Region ", region_id, " not found in region table.")
  }

  flux <- as.numeric(row[1, filters])
  fluxerr <- as.numeric(row[1, err_cols])
  fluxerr_eff <- sqrt(pmax(fluxerr, 0)^2 + (systematic_frac * pmax(flux, 0))^2)

  flux_tbl <- tibble::tibble(
    filter = filters,
    cenwave = prospect_cenwaves(filters),
    flux = flux,
    fluxerr = pmax(fluxerr_eff, 1e-30)
  )
  filtout <- prospect_filter_objects(filters)
  z_use <- if (is.null(fixed_z)) lookup_redshift(tag) else fixed_z
  tau_screen_lo <- as.numeric(Sys.getenv("SAGUI_PROSPECT_TAU_SCREEN_LO", unset = "0.00"))
  tau_screen_hi <- as.numeric(Sys.getenv("SAGUI_PROSPECT_TAU_SCREEN_HI", unset = "5.50"))
  mpeak_lo <- as.numeric(Sys.getenv("SAGUI_PROSPECT_MPEAK_LO", unset = "0.30"))
  mpeak_hi <- as.numeric(Sys.getenv("SAGUI_PROSPECT_MPEAK_HI", unset = "8.50"))
  mperiod_lo <- as.numeric(Sys.getenv("SAGUI_PROSPECT_MPERIOD_LO", unset = "0.20"))
  mperiod_hi <- as.numeric(Sys.getenv("SAGUI_PROSPECT_MPERIOD_HI", unset = "5.00"))

  parm_names <- c("mSFR", "mpeak", "mperiod", "tau_screen")
  logged <- c(TRUE, FALSE, FALSE, FALSE)
  lo <- c(-3.0, mpeak_lo, mperiod_lo, tau_screen_lo)
  hi <- c(3.0, mpeak_hi, mperiod_hi, tau_screen_hi)
  arglist <- list(
    z = z_use,
    ref = "Planck18",
    massfunc = massfunc_snorm,
    tau_birth = 1.0
  )

  if (isTRUE(free_Z)) {
    parm_names <- c(parm_names, "Z")
    logged <- c(logged, TRUE)
    lo <- c(lo, log10(0.02 * 10^logzsol_bounds[1]))
    hi <- c(hi, log10(0.02 * 10^logzsol_bounds[2]))
  } else {
    arglist$Z <- fixed_Z
  }

  data_obj <- list(
    flux = as.data.frame(flux_tbl),
    SFH = SFHfunc,
    speclib = BC03lr,
    Dale = FALSE,
    filtout = filtout,
    fit = "optim",
    like = "st",
    parm.names = parm_names,
    logged = logged,
    intervals = list(lo = lo, hi = hi),
    mon.names = c("masstot", "SFRburst"),
    verbose = FALSE,
    arglist = arglist
  )

  attr(data_obj, "tag") <- tag
  attr(data_obj, "region_id") <- region_id
  attr(data_obj, "photometry_mode") <- photometry_mode
  attr(data_obj, "filters") <- filters
  attr(data_obj, "fixed_z") <- z_use
  attr(data_obj, "fixed_Z") <- fixed_Z
  attr(data_obj, "free_Z") <- isTRUE(free_Z)
  attr(data_obj, "region_csv_override") <- region_csv_override
  data_obj
}

default_simple_starts <- function(tag, region_id, free_Z = FALSE) {
  current_base <- simple_current_start(tag, region_id)
  current_path <- current_results_csv_path(tag)
  z_start <- log10(0.02)
  tau_hi <- as.numeric(Sys.getenv("SAGUI_PROSPECT_TAU_SCREEN_HI", unset = "5.50"))
  dust_mid <- min(1.5, tau_hi)
  dust_high <- min(3.0, tau_hi)
  dust_max <- min(4.5, tau_hi)
  if (is.character(current_path) && !is.na(current_path) && file.exists(current_path)) {
    cur <- readr::read_csv(current_path, show_col_types = FALSE) |>
      filter(.data$region == region_id) |>
      slice(1)
    if (nrow(cur) && "logzsol" %in% names(cur)) {
      z_start <- log10(pmin(pmax(0.02 * 10^cur$logzsol, 0.02 * 10^-2), 0.02 * 10^0.3))
    }
  }

  starts <- list(
    current_base,
    c(log10(0.10), 2.5, 1.0, 0.5),
    c(log10(0.30), 1.5, 0.8, 0.2),
    c(log10(0.10), 2.5, 1.0, dust_mid),
    c(log10(0.30), 1.5, 0.8, dust_high),
    c(log10(0.70), 1.0, 0.5, dust_max)
  )

  if (isTRUE(free_Z)) {
    starts <- lapply(starts, function(x) c(x, z_start))
  }

  starts
}

fit_simple_prospect_region <- function(tag,
                                       region_id,
                                       photometry_mode = c("raw", "psfmatched"),
                                       systematic_frac = 0.05,
                                       fixed_z = NULL,
                                       fixed_Z = 0.02,
                                       free_Z = FALSE,
                                       logzsol_bounds = c(-2.0, 0.3),
                                       region_csv_override = NULL) {
  photometry_mode <- match.arg(photometry_mode)
  data_obj <- make_simple_prospect_data(
    tag = tag,
    region_id = region_id,
    photometry_mode = photometry_mode,
    systematic_frac = systematic_frac,
    fixed_z = fixed_z,
    fixed_Z = fixed_Z,
    free_Z = free_Z,
    logzsol_bounds = logzsol_bounds,
    region_csv_override = region_csv_override
  )

  objective_fn <- function(par) {
    out <- suppressWarnings(
      tryCatch(
        ProSpectSEDlike(par, data_obj),
        error = function(e) Inf
      )
    )
    if (!is.finite(out)) Inf else out
  }

  starts <- default_simple_starts(tag, region_id, free_Z = free_Z)
  fits <- lapply(starts, function(start) {
    tryCatch(
      nlminb(
        start = start,
        objective = objective_fn,
        lower = data_obj$intervals$lo,
        upper = data_obj$intervals$hi,
        control = list(iter.max = 400, eval.max = 1200)
      ),
      error = function(e) NULL
    )
  })
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) {
    stop("All simple ProSpect fits failed for region ", region_id)
  }

  scores <- vapply(fits, function(x) if (is.finite(x$objective)) x$objective else Inf, numeric(1))
  best <- fits[[which.min(scores)]]

  data_obj$fit <- "check"
  checked <- ProSpectSEDlike(best$par, data_obj)
  stars <- checked$SEDout$Stars
  z_linear <- if (isTRUE(free_Z)) 10^best$par[5] else fixed_Z

  obs_flux <- data_obj$flux$flux
  model_flux <- checked$SEDout$Photom
  rms_sigma <- sqrt(mean(((obs_flux - model_flux) / pmax(data_obj$flux$fluxerr, 1e-30))^2, na.rm = TRUE))

  region_tbl <- read_region_table_override(tag, photometry_mode, region_csv_override)
  n_pix <- region_tbl$n_pix[region_tbl$region == region_id][[1]]

  summary <- tibble::tibble(
    tag = tag,
    region = region_id,
    n_pix = n_pix,
    photometry_mode = photometry_mode,
    z_used = attr(data_obj, "fixed_z"),
    Z_mode = if (isTRUE(free_Z)) "free" else "fixed",
    Z_fixed = if (isTRUE(free_Z)) NA_real_ else fixed_Z,
    Z_fit = z_linear,
    logzsol = log10(pmax(z_linear, 1e-12) / 0.02),
    systematic_frac = systematic_frac,
    convergence = best$convergence,
    objective = best$objective,
    rms_sigma = rms_sigma,
    logMformed = log10(pmax(stars$masstot, .Machine$double.eps)),
    age_mw_gyr = sum(stars$agevec * stars$massvec, na.rm = TRUE) / sum(stars$massvec, na.rm = TRUE) / 1e9,
    sfr_recent = {
      idx <- stars$agevec <= 1e8
      if (any(idx, na.rm = TRUE)) {
        stats::weighted.mean(stars$SFR[idx], w = pmax(stars$massvec[idx], 0), na.rm = TRUE)
      } else {
        stars$SFRburst
      }
    },
    mSFR = 10^best$par[1],
    mpeak = best$par[2],
    mperiod = best$par[3],
    tau_screen = best$par[4]
  )

  list(
    data = data_obj,
    fit = best,
    checked = checked,
    summary = summary
  )
}

plot_simple_prospect_region <- function(fit_result,
                                        title = NULL,
                                        subtitle = NULL) {
  data_obj <- fit_result$data
  checked <- fit_result$checked
  summary <- fit_result$summary
  filters <- attr(data_obj, "filters")

  throughput_df <- readr::read_csv(throughput_nircam_path, show_col_types = FALSE) |>
    filter(.data$filter %in% filters) |>
    mutate(
      wavelength = .data$wavelength * 1e4,
      filter = factor(.data$filter, levels = filters)
    )

  filter_cols <- setNames(viridisLite::turbo(length(filters)), filters)

  model_curve <- tibble::tibble(
    wave_A = checked$SEDout$FinalFlux[, "wave"],
    flux = checked$SEDout$FinalFlux[, "flux"]
  )

  model_photom <- tibble::tibble(
    filter = filters,
    cenwave = data_obj$flux$cenwave,
    flux = checked$SEDout$Photom
  )

  obs_df <- tibble::as_tibble(data_obj$flux)

  x_lo <- min(obs_df$cenwave) * 0.90
  x_hi <- max(obs_df$cenwave) * 1.05
  model_curve_visible <- model_curve |>
    filter(.data$wave_A >= x_lo, .data$wave_A <= x_hi)
  y_pos <- c(
    obs_df$flux[is.finite(obs_df$flux) & obs_df$flux > 0],
    model_photom$flux[is.finite(model_photom$flux) & model_photom$flux > 0],
    model_curve_visible$flux[is.finite(model_curve_visible$flux) & model_curve_visible$flux > 0 & model_curve_visible$flux >= 0.05 * max(obs_df$flux, na.rm = TRUE)]
  )
  y_min <- min(y_pos, na.rm = TRUE)
  y_max <- max(y_pos, na.rm = TRUE)
  trans_base <- y_min * 0.48
  trans_top <- y_min * 0.90
  y_lo <- trans_base * 0.82
  y_hi <- y_max * 1.20

  thr_scaled <- throughput_df |>
    group_by(.data$filter) |>
    mutate(
      throughput_norm = .data$throughput / max(.data$throughput, na.rm = TRUE),
      y = trans_base * (trans_top / trans_base)^pmax(.data$throughput_norm, 0)
    ) |>
    ungroup() |>
    filter(.data$y >= y_lo, .data$y <= y_hi)

  model_curve_visible <- model_curve_visible |>
    filter(.data$flux >= y_lo, .data$flux <= y_hi)

  if (is.null(title)) {
    title <- paste0("Simple ProSpect fit: ", summary$tag, " region ", summary$region)
  }
  if (is.null(subtitle)) {
    subtitle <- paste0(
      "z = ", sprintf("%.4f", summary$z_used),
      " | Z ", if (identical(summary$Z_mode, "free")) {
        paste0("fit = ", sprintf("%.4f", summary$Z_fit), " (logZ/Zsun = ", sprintf("%.2f", summary$logzsol), ")")
      } else {
        paste0("fixed = ", sprintf("%.4f", summary$Z_fixed))
      },
      " | RMS = ", sprintf("%.2f", summary$rms_sigma)
    )
  }

  ggplot() +
    geom_line(
      data = thr_scaled,
      aes(x = .data$wavelength, y = .data$y, colour = .data$filter, group = .data$filter),
      linewidth = 1.4,
      alpha = 0.95,
      show.legend = FALSE
    ) +
    geom_line(
      data = model_curve_visible,
      aes(x = .data$wave_A, y = .data$flux, colour = "Fitted spectrum"),
      linewidth = 1.15,
      alpha = 0.95,
      show.legend = TRUE
    ) +
    geom_errorbar(
      data = obs_df,
      aes(x = .data$cenwave, ymin = pmax(.data$flux - .data$fluxerr, trans_base * 0.98), ymax = .data$flux + .data$fluxerr),
      width = 0,
      linewidth = 0.42,
      colour = "#2D2D2D"
    ) +
    geom_point(
      data = obs_df,
      aes(x = .data$cenwave, y = .data$flux, colour = "Observed photometry"),
      shape = 21,
      fill = "white",
      stroke = 0.85,
      size = 2.8,
      show.legend = TRUE
    ) +
    geom_point(
      data = model_photom,
      aes(x = .data$cenwave, y = .data$flux),
      shape = 16,
      colour = "#7F63E3",
      size = 1.5,
      alpha = 0.9
    ) +
    scale_colour_manual(
      values = c(
        setNames(unname(filter_cols), names(filter_cols)),
        "Fitted spectrum" = "#B197FC",
        "Observed photometry" = "#2D2D2D"
      ),
      breaks = c("Fitted spectrum", "Observed photometry")
    ) +
    scale_x_continuous(
      breaks = c(10000, 20000, 30000, 40000),
      labels = pretty_angstrom_ticks,
      limits = c(x_lo, x_hi),
      expand = c(0, 0)
    ) +
    scale_y_log10(
      limits = c(y_lo, y_hi),
      expand = c(0, 0)
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Wavelength [Å]",
      y = "Flux density [Jy]",
      colour = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = c(0.15, 0.90),
      legend.background = element_rect(fill = scales::alpha("white", 0.9), colour = "#D0D0D0"),
      axis.line = element_line(colour = "#3A3A3A", linewidth = 0.45),
      panel.border = element_rect(colour = "#4A4A4A", fill = NA, linewidth = 0.35),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 11)
    )
}
