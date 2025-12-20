
plot_sed_with_filters <- function(
    reg_long,                   # data.frame/tibble: region, band, flux, flux_err
    region_id,                  # integer region ID
    throughput_tbl,             # data.frame: wavelength, throughput, filter
    nircam_lambda_um,           # named numeric: F090W = 0.902, ...
    scale_mode = c("normalized","auto"),
    band_height_frac = 0.15,
    errors_are_fractional = FALSE,
    show_legend = TRUE
){
  scale_mode <- match.arg(scale_mode)

  ## --------------------------------------------------------
  ## 0) Infer filter order from nircam_lambda_um + throughput
  ## --------------------------------------------------------
  filters_all  <- names(nircam_lambda_um)
  filters_thr  <- unique(as.character(throughput_tbl$filter))
  filters_used <- filters_all[filters_all %in% filters_thr]

  if (length(filters_used) == 0L) {
    stop("No overlap between names(nircam_lambda_um) and throughput_tbl$filter.")
  }

  ## --------------------------------------------------------
  ## 1) Subset SED to the requested region (base R)
  ## --------------------------------------------------------
  reg_long <- as.data.frame(reg_long)

  sed <- reg_long[reg_long$region == region_id, , drop = FALSE]
  if (nrow(sed) == 0L) {
    stop("No rows found for region == ", region_id)
  }

  # band is "0","1",... -> integer
  sed$band_idx <- as.integer(as.character(sed$band))

  # sanity: warn if any band_idx out of range
  bad_idx <- is.na(sed$band_idx) | sed$band_idx < 0 | sed$band_idx >= length(filters_used)
  if (any(bad_idx)) {
    warning("Some band indices are out of range for filters_used; dropping those rows.")
    sed <- sed[!bad_idx, , drop = FALSE]
  }

  sed$filt      <- filters_used[sed$band_idx + 1L]
  sed$lambda_um <- nircam_lambda_um[sed$filt]

  sed <- sed[!is.na(sed$lambda_um), , drop = FALSE]
  if (nrow(sed) == 0L) {
    stop("No valid lambda_um after mapping band indices to filters.")
  }

  ## --------------------------------------------------------
  ## 2) Throughput subset (base R)
  ## --------------------------------------------------------
  thr <- as.data.frame(throughput_tbl)
  thr$filter <- as.character(thr$filter)
  thr        <- thr[thr$filter %in% filters_used, , drop = FALSE]

  if (!all(c("wavelength", "throughput") %in% names(thr))) {
    stop("throughput_tbl must have columns 'wavelength' and 'throughput'. Got: ",
         paste(names(thr), collapse = ", "))
  }

  ## --------------------------------------------------------
  ## 3) Scaling
  ## --------------------------------------------------------
  if (scale_mode == "normalized") {

    K <- max(sed$flux, na.rm = TRUE)

    sed$y <- sed$flux / K
    if (errors_are_fractional) {
      sed$ymin <- (sed$flux * (1 - sed$flux_err)) / K
      sed$ymax <- (sed$flux * (1 + sed$flux_err)) / K
    } else {
      sed$ymin <- (sed$flux - sed$flux_err) / K
      sed$ymax <- (sed$flux + sed$flux_err) / K
    }
    sed$ymin <- pmax(sed$ymin, 0)

    thr$y <- thr$throughput * band_height_frac

    ylab    <- "Normalized Flux"
    ylimits <- c(-0.02, 1.05)

  } else {

    sed_scale <- stats::quantile(sed$flux, 0.95, na.rm = TRUE)

    thr$y <- thr$throughput * sed_scale * band_height_frac

    sed$y <- sed$flux
    if (errors_are_fractional) {
      sed$ymin <- pmax(sed$flux * (1 - sed$flux_err), 0)
      sed$ymax <- sed$flux * (1 + sed$flux_err)
    } else {
      sed$ymin <- pmax(sed$flux - sed$flux_err, 0)
      sed$ymax <- sed$flux + sed$flux_err
    }

    ylab    <- "Flux"
    ylimits <- NULL
  }

  ## --------------------------------------------------------
  ## 4) Plot (ggplot only)
  ## --------------------------------------------------------
  thr$filter <- factor(thr$filter, levels = filters_used)

  p <- ggplot2::ggplot() +
    ggplot2::geom_area(
      data  = thr,
      ggplot2::aes(x = wavelength, y = y, fill = filter, group = filter),
      position = "identity",
      alpha    = 0.35,
      color    = NA,
      linewidth = 0
    ) +
    ggplot2::geom_errorbar(
      data  = sed,
      ggplot2::aes(x = lambda_um, ymin = ymin, ymax = ymax),
      width = 0.06,
      color = "black"
    ) +
    ggplot2::geom_point(
      data  = sed,
      ggplot2::aes(x = lambda_um, y = y),
      size  = 2.6,
      shape = 21,
      fill  = "orange2",
      color = "black"
    ) +
    ggplot2::labs(
      x     = "Wavelength (µm)",
      y     = ylab,
      title = paste("Region", region_id, "- JWST/NIRCam SED")
    ) +
    ggplot2::scale_fill_viridis_d(option = "turbo", name = "Filter") +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position = if (show_legend) "right" else "none",
      panel.grid      = ggplot2::element_blank(),
      axis.line       = ggplot2::element_line(color = "black")
    )

  if (!is.null(ylimits)) {
    p <- p + ggplot2::coord_cartesian(ylim = ylimits)
  }

  p
}


## ------------------------------------------------------------
## Region photometry based on SED labels (from hclust)
## ------------------------------------------------------------

SED <- RegionPhotometry(
  X,
  seg_cap$cluster_map,
  error_fallback = "poisson"
)

num_cols  <- as.character(0:9)        # adjust to nbands if needed
err_cols  <- paste0(num_cols, "_err")
neff_cols <- paste0(num_cols, "_n_eff")

# Rename only the numeric bands 0..9 here
flux_wide_named <- SED$flux_wide %>%
  rename_with(
    ~ nircam_filters[match(.x, num_cols)],
    .cols = all_of(num_cols)
  )

flux_wide_named_jy <- SED$flux_wide %>%
  # 1) flux [nJy] -> Jy
  mutate(across(all_of(num_cols), ~ .x * 1e-8)) %>%
  # 2) errors [nJy] -> Jy
  mutate(across(all_of(err_cols), ~ .x * 1e-8)) %>%
  # 3) rename band columns
  rename_with(
    ~ nircam_filters[match(.x, num_cols)],
    .cols = all_of(num_cols)
  ) %>%
  # 4) rename error columns: 0_err -> F090W_err
  rename_with(
    ~ paste0(nircam_filters, "_err")[match(.x, err_cols)],
    .cols = all_of(err_cols)
  ) %>%
  # 5) rename n_eff columns: 0_n_eff -> F090W_n_eff
  rename_with(
    ~ paste0(nircam_filters, "_n_eff")[match(.x, neff_cols)],
    .cols = all_of(neff_cols)
  )

readr::write_csv(flux_wide_named_jy, "SED_flux_wide_Jy_reg8.csv")









filters  <- readr::read_csv("throughput_nircam.csv", show_col_types = FALSE)

p1 <- plot_sed_with_filters(
  reg_long         = SED$flux_long,
  region       = 22,
  throughput_tbl   = throughput_all,
  nircam_lambda_um = nircam_lambda_um,
  scale_mode       = "normalized",
  band_height_frac = 0.99,
  errors_are_fractional = FALSE,
  show_legend      = TRUE
)
p1





