#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(FITSio)
  library(readr)
  library(sagui)
})

source("/Users/rd23aag/Documents/GitHub/crp8_segmentation/src/R/pipelines/prospect/prospect_common.R")

tag <- Sys.getenv("SAGUI_PROSPECT_TAG", unset = "sagui10")
photometry_mode <- Sys.getenv("SAGUI_PROSPECT_PHOTOMETRY_MODE", unset = "raw")
photometry_mode <- match.arg(photometry_mode, c("raw", "psfmatched"))
results_override <- Sys.getenv("SAGUI_PROSPECT_RESULTS_CSV", unset = "")
segmentation_override <- Sys.getenv("SAGUI_PROSPECT_SEGMENTATION_FITS", unset = "")
photometry_cube_override <- Sys.getenv("SAGUI_PROSPECT_PHOTOMETRY_CUBE", unset = "")
palette_mode <- Sys.getenv("SAGUI_PROSPECT_PALETTE", unset = "legacy_cosmic_dusk")
style_suffix <- Sys.getenv("SAGUI_PROSPECT_STYLE_SUFFIX", unset = palette_mode)

use_smooth <- TRUE
smooth_lambda <- 10
smooth_adjacency <- "queen"
na_color <- "black"
crop_padding_px <- 4L
default_contrast_quantiles <- c(0.02, 0.98)
show_isophote_contours <- tolower(Sys.getenv("SAGUI_PROSPECT_SHOW_CONTOURS", unset = "1")) %in% c("1", "true", "yes")
contour_bins <- 8L
contour_alpha <- 0.42
contour_linewidth <- 0.22
property_contrast_quantiles <- list(
  logM = c(0.05, 0.95),
  logZ = c(0.10, 0.90),
  Av = c(0.02, 0.995),
  age_mw = c(0.05, 0.95),
  sfr = c(0.03, 0.97),
  sigma_sfr = c(0.03, 0.97),
  ssfr = c(0.03, 0.97)
)

out_dirs <- make_output_dirs(tag, photometry_mode)
if (nzchar(style_suffix)) {
  out_dirs$maps <- paste0(out_dirs$maps, "_", style_suffix)
}
figure_out_dir <- file.path(
  prospect_base_dir,
  "results/figures/sed_maps/prospect",
  tag,
  if (nzchar(style_suffix)) paste0(photometry_mode, "_", style_suffix) else photometry_mode
)
dir.create(out_dirs$maps, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_out_dir, recursive = TRUE, showWarnings = FALSE)

fit_results_path <- if (nzchar(results_override)) results_override else {
  file.path(out_dirs$fitting, paste0(tag, "_prospect_results.csv"))
}
segmentation_path <- if (nzchar(segmentation_override)) {
  segmentation_override
} else {
  segmentation_fits_path(tag)
}

if (!file.exists(segmentation_path)) {
  stop("Segmentation FITS not found: ", segmentation_path)
}
if (!file.exists(fit_results_path)) {
  stop("Fit results CSV not found: ", fit_results_path)
}

property_specs <- list(
  logM = list(
    label = expression(log[10](M/M[odot])),
    units = "dex",
    scale = "log10(M/Msun)"
  ),
  logZ = list(
    label = expression(log[10](Z/Z[odot])),
    units = "dex",
    scale = "log10(Z/Zsun)"
  ),
  Av = list(
    label = expression(A[V]~"[mag]"),
    units = "mag",
    scale = "Av"
  ),
  age_mw = list(
    label = expression(phantom(.)*"<t>"[M]~"[Gyr]"),
    units = "Gyr",
    scale = "mass-weighted age"
  ),
  sfr = list(
    label = expression(SFR~"["*M[odot]~yr^{-1}*"]"),
    units = "Msun/yr",
    scale = "SFR"
  ),
  sigma_sfr = list(
    label = expression(Sigma[SFR]~"["*M[odot]~yr^{-1}~kpc^{-2}*"]"),
    units = "Msun/yr/kpc^2",
    scale = "Sigma_SFR"
  ),
  ssfr = list(
    label = expression(sSFR~"[yr"^{-1}*"]"),
    units = "yr^-1",
    scale = "sSFR"
  )
)

palette_stops <- function(mode = c("legacy_cosmic_dusk", "paper_cool", "paper_sea", "paper_softwarm")) {
  mode <- match.arg(mode)
  switch(
    mode,
    legacy_cosmic_dusk = c(
      "#0A1026",
      "#182A6B",
      "#2F5F9A",
      "#4FA7A6",
      "#D8E1D3",
      "#F1C76A",
      "#A96A2A"
    ),
    paper_cool = c(
      "#07111F",
      "#112A5A",
      "#1F5AA6",
      "#49A8CC",
      "#D8F0EC",
      "#F5F1E7"
    ),
    paper_sea = c(
      "#06101C",
      "#19345F",
      "#2E6FA8",
      "#59B8BF",
      "#CDEBE0",
      "#F3EEE0"
    ),
    paper_softwarm = c(
      "#07111F",
      "#163163",
      "#2D6AA9",
      "#78C8C0",
      "#F0E7CF",
      "#E6C98E"
    )
  )
}

palette_cosmic_dusk <- function(n = 256, mode = palette_mode) {
  stops <- palette_stops(mode)
  grDevices::colorRampPalette(stops, space = "Lab")(n)
}

mat_to_df <- function(Z, flip_y = FALSE, transpose = FALSE) {
  if (transpose) Z <- t(Z)
  nx <- nrow(Z)
  ny <- ncol(Z)
  df <- expand.grid(x = seq_len(nx), y = seq_len(ny))
  df$value <- as.vector(Z)
  if (flip_y) df$y <- ny - df$y + 1
  df
}

crop_df_to_data <- function(df, pad = 4L) {
  df_ok <- df[is.finite(df$value), ]
  if (!nrow(df_ok)) return(df)

  x_all <- sort(unique(df$x))
  y_all <- sort(unique(df$y))
  x_rng <- range(df_ok$x)
  y_rng <- range(df_ok$y)

  xmin <- max(min(x_all), x_rng[1] - pad)
  xmax <- min(max(x_all), x_rng[2] + pad)
  ymin <- max(min(y_all), y_rng[1] - pad)
  ymax <- min(max(y_all), y_rng[2] + pad)

  df[df$x >= xmin & df$x <= xmax & df$y >= ymin & df$y <= ymax, ]
}

resolve_contrast_quantiles <- function(property_name, values) {
  probs <- property_contrast_quantiles[[property_name]]
  if (is.null(probs)) return(default_contrast_quantiles)
  vals_ok <- values[is.finite(values)]
  n_unique <- length(unique(signif(vals_ok, 4)))
  if (identical(property_name, "logZ") && n_unique <= 15) return(c(0, 1))
  probs
}

make_scale_spec <- function(vals_ok, property_name = NULL) {
  probs <- resolve_contrast_quantiles(property_name, vals_ok)
  lims <- as.numeric(stats::quantile(vals_ok, probs = probs, na.rm = TRUE))
  if (!all(is.finite(lims)) || diff(lims) == 0) {
    lims <- range(vals_ok, finite = TRUE)
  }
  if (!all(is.finite(lims))) {
    lims <- c(0, 1)
  }
  if (diff(lims) == 0) {
    centre <- lims[[1]]
    pad <- max(abs(centre) * 0.05, 1e-6)
    lims <- c(centre - pad, centre + pad)
  }

  inner_probs <- c(0.00, 0.03, 0.08, 0.16, 0.28, 0.45, 0.62, 0.78, 0.90, 0.97, 1.00)
  qs <- stats::quantile(vals_ok, probs = inner_probs, na.rm = TRUE, names = FALSE)
  qs <- pmin(pmax(qs, lims[1]), lims[2])
  qs <- unique(qs)
  if (length(qs) < 3L) {
    qs <- seq(lims[1], lims[2], length.out = 7L)
  }
  if (diff(range(qs, finite = TRUE)) == 0) {
    qs <- seq(lims[1], lims[2], length.out = 7L)
  }

  n_breaks <- if (identical(property_name, "sigma_sfr") || identical(property_name, "ssfr")) 3 else 4
  raw_breaks <- pretty(lims, n = n_breaks)
  raw_breaks <- raw_breaks[raw_breaks >= lims[1] & raw_breaks <= lims[2]]
  if (length(raw_breaks) < 3L) {
    raw_breaks <- seq(lims[1], lims[2], length.out = n_breaks)
  }

  list(
    lims = lims,
    colours = palette_cosmic_dusk(length(qs), mode = palette_mode),
    values = scales::rescale(qs, from = range(qs, finite = TRUE)),
    breaks = unique(raw_breaks)
  )
}

format_break_labels <- function(breaks, property_name) {
  breaks <- as.numeric(breaks)
  ok <- is.finite(breaks)
  out <- rep("", length(breaks))
  if (!any(ok)) return(out)

  vals <- breaks[ok]
  uniq <- sort(unique(vals))
  min_step <- if (length(uniq) >= 2L) min(diff(uniq)) else NA_real_

  if (identical(property_name, "sigma_sfr") || identical(property_name, "ssfr")) {
    for (digits in 2:5) {
      lab <- scales::label_scientific(digits = digits)(vals)
      if (length(unique(lab)) == length(lab)) {
        out[ok] <- lab
        return(out)
      }
    }
    out[ok] <- scales::label_scientific(digits = 5)(vals)
    return(out)
  }

  if (!is.finite(min_step) || min_step <= 0) {
    acc_candidates <- c(0.1, 0.01, 0.001, 0.0001)
  } else {
    base_acc <- 10^floor(log10(min_step))
    acc_candidates <- unique(c(base_acc, base_acc / 2, base_acc / 5, base_acc / 10, 0.1, 0.01, 0.001, 0.0001))
    acc_candidates <- acc_candidates[acc_candidates > 0]
    acc_candidates <- sort(acc_candidates, decreasing = TRUE)
  }

  if (identical(property_name, "logZ")) {
    acc_candidates <- sort(unique(c(0.01, 0.005, acc_candidates)), decreasing = TRUE)
  }

  for (acc in acc_candidates) {
    lab <- scales::label_number(accuracy = acc, trim = TRUE)(vals)
    if (length(unique(lab)) == length(lab)) {
      out[ok] <- lab
      return(out)
    }
  }

  out[ok] <- scales::label_number(accuracy = min(acc_candidates), trim = TRUE)(vals)
  out
}

crop_matrix_to_df <- function(mat, df_ref) {
  if (!nrow(df_ref)) return(mat_to_df(mat))
  xr <- range(df_ref$x)
  yr <- range(df_ref$y)
  df <- mat_to_df(mat, flip_y = FALSE, transpose = FALSE)
  df[df$x >= xr[1] & df$x <= xr[2] & df$y >= yr[1] & df$y <= yr[2], ]
}

plot_univariate <- function(df_pix,
                            property_name = NULL,
                            contour_df = NULL,
                            show_legend = TRUE,
                            transparent_bg = TRUE,
                            legend_n = 4,
                            legend_position = "bottom") {
  vals_ok <- df_pix$value[is.finite(df_pix$value)]
  if (!length(vals_ok)) stop("No finite values available for plotting.")

  scale_spec <- make_scale_spec(vals_ok, property_name)
  bg_fill <- if (transparent_bg) "transparent" else na_color
  na_val <- if (transparent_bg) NA else na_color

  p <- ggplot(df_pix, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    scale_fill_gradientn(
      colours = scale_spec$colours,
      values = scale_spec$values,
      limits = scale_spec$lims,
      breaks = scale_spec$breaks,
      labels = function(brks) format_break_labels(brks, property_name),
      oob = scales::squish,
      na.value = na_val,
      name = NULL,
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        label.position = "bottom",
        ticks = TRUE,
        barheight = grid::unit(9, "mm"),
        barwidth = grid::unit(140, "mm")
      )
    ) +
    theme_void() +
    theme(
      legend.position = if (show_legend) legend_position else "none",
      plot.margin = margin(0, 0, 0, 0),
      panel.background = element_rect(fill = bg_fill, colour = NA),
      plot.background = element_rect(fill = bg_fill, colour = NA),
      legend.background = element_rect(fill = bg_fill, colour = NA),
      legend.key = element_rect(fill = bg_fill, colour = NA),
      legend.box.spacing = grid::unit(1.5, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.text = element_text(size = 19)
    )

  if (!is.null(contour_df)) {
    p <- p +
      geom_contour(
        data = contour_df,
        aes(x = x, y = y, z = value),
        inherit.aes = FALSE,
        bins = contour_bins,
        colour = scales::alpha("white", contour_alpha),
        linewidth = contour_linewidth
      )
  }

  p
}

region_to_pixels <- function(region_map, region_df) {
  idmat <- region_map
  idmat[idmat == 0L] <- NA_integer_

  vals <- region_df$value
  names(vals) <- as.character(region_df$region_id)

  out <- matrix(NA_real_, nrow(idmat), ncol(idmat))
  ok <- is.finite(idmat)
  out[ok] <- vals[as.character(idmat[ok])]
  out
}

map_region_field_to_pixels <- function(region_map, reg_field, property_name) {
  vals_ok <- reg_field$value[is.finite(reg_field$value)]
  if (isTRUE(use_smooth) && length(unique(signif(vals_ok, 4))) > 2) {
    res <- sagui::smooth_region_field_laplacian(
      region_id_mat = region_map,
      region_values = reg_field,
      adjacency = smooth_adjacency,
      lambda = smooth_lambda
    )
    return(res$interpolated_matrix)
  }
  region_to_pixels(region_map, reg_field)
}

pick_first_existing <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (!length(hit)) return(NULL)
  hit[[1]]
}

get_hdr_value <- function(fits_obj, key, default = NA_character_) {
  idx <- match(key, fits_obj$hdr)
  if (is.na(idx) || idx >= length(fits_obj$hdr)) {
    return(default)
  }
  fits_obj$hdr[idx + 1L]
}

pixel_scale_arcsec <- function(fits_obj) {
  cdelt1 <- as.numeric(get_hdr_value(fits_obj, "CDELT1", "1"))
  cdelt2 <- as.numeric(get_hdr_value(fits_obj, "CDELT2", "1"))
  pc11 <- as.numeric(get_hdr_value(fits_obj, "PC1_1", "1"))
  pc22 <- as.numeric(get_hdr_value(fits_obj, "PC2_2", "1"))
  c(
    x = abs(cdelt1 * pc11) * 3600,
    y = abs(cdelt2 * pc22) * 3600
  )
}

angular_diameter_distance_mpc <- function(z, H0 = 70, omega_m = 0.3) {
  if (!is.finite(z) || z <= 0) return(NA_real_)
  omega_l <- 1 - omega_m
  ez_inv <- function(zz) 1 / sqrt(omega_m * (1 + zz)^3 + omega_l)
  dc_mpc <- (299792.458 / H0) * stats::integrate(ez_inv, lower = 0, upper = z, subdivisions = 2000L)$value
  dc_mpc / (1 + z)
}

kpc_per_arcsec_at_z <- function(z) {
  da_mpc <- angular_diameter_distance_mpc(z)
  if (!is.finite(da_mpc)) return(NA_real_)
  da_mpc * 1000 / 206265
}

infer_results_redshift <- function(df, tag) {
  if ("z_used" %in% names(df)) {
    z_vals <- unique(round(as.numeric(df$z_used[is.finite(df$z_used)]), 8))
    if (length(z_vals)) return(z_vals[[1]])
  }
  lookup_redshift(tag, fallback = NA_real_)
}

region_area_table <- function(region_map, cube_fits_obj, z_used) {
  ids <- sort(unique(as.integer(region_map[is.finite(region_map) & region_map > 0])))
  if (!length(ids) || !is.finite(z_used)) {
    return(tibble::tibble(region_id = integer(), n_pix = integer(), area_kpc2 = numeric()))
  }

  pix_scale <- pixel_scale_arcsec(cube_fits_obj)
  kpc_per_arcsec <- kpc_per_arcsec_at_z(z_used)
  pixel_area_kpc2 <- (pix_scale["x"] * kpc_per_arcsec) * (pix_scale["y"] * kpc_per_arcsec)

  tibble::tibble(
    region_id = ids,
    n_pix = vapply(ids, function(id) sum(region_map == id, na.rm = TRUE), integer(1)),
    area_kpc2 = n_pix * pixel_area_kpc2
  )
}

build_property_table <- function(df) {
  region_col <- pick_first_existing(df, c("region", "region_id"))
  if (is.null(region_col)) stop("Results CSV must contain a region identifier column.")

  out <- tibble::tibble(region_id = as.integer(df[[region_col]]))

  mass_col <- pick_first_existing(df, c("logMstar", "logmass", "logM", "logMformed"))
  if (!is.null(mass_col)) {
    out$logM <- as.numeric(df[[mass_col]])
  } else if ("mass" %in% names(df)) {
    out$logM <- log10(pmax(as.numeric(df[["mass"]]), .Machine$double.eps))
  }

  logz_col <- pick_first_existing(df, c("logzsol", "logZ", "metallicity_log"))
  if (!is.null(logz_col)) out$logZ <- as.numeric(df[[logz_col]])

  av_col <- pick_first_existing(df, c("Av", "A_V"))
  if (!is.null(av_col)) {
    out$Av <- as.numeric(df[[av_col]])
  } else if ("dust2" %in% names(df)) {
    out$Av <- 1.086 * as.numeric(df[["dust2"]])
  }

  age_col <- pick_first_existing(df, c("age_mw", "age_mass_weighted", "age_mw_gyr", "tage"))
  if (!is.null(age_col)) out$age_mw <- as.numeric(df[[age_col]])

  sfr_col <- pick_first_existing(df, c("sfr", "sfr_recent", "sfr_100myr", "sfr_10myr"))
  if (!is.null(sfr_col)) out$sfr <- as.numeric(df[[sfr_col]])

  ssfr_col <- pick_first_existing(df, c("ssfr", "ssfr_100myr", "ssfr_10myr"))
  if (!is.null(ssfr_col)) {
    out$ssfr <- as.numeric(df[[ssfr_col]])
  } else if ("sfr" %in% names(out) && "logM" %in% names(out)) {
    out$ssfr <- out$sfr / pmax(10^out$logM, .Machine$double.eps)
  }

  out
}

results_df <- readr::read_csv(fit_results_path, show_col_types = FALSE)
property_df <- build_property_table(results_df)

seg <- FITSio::readFITS(segmentation_path)
region_map <- seg$imDat
region_map[region_map == 0L] <- NA_integer_

cube_path <- NULL
if (nzchar(photometry_cube_override)) {
  cube_path <- photometry_cube_override
} else {
  manifest_path <- file.path(out_dirs$fitting, "manifest.csv")
  if (file.exists(manifest_path)) {
    manifest <- readr::read_csv(manifest_path, show_col_types = FALSE)
    if ("photometry_cube" %in% names(manifest) && nrow(manifest)) {
      cube_path <- manifest$photometry_cube[[1]]
    }
  }
}
if (is.null(cube_path) || !nzchar(cube_path)) {
  cube_path <- if (identical(photometry_mode, "raw")) raw_cube_name(tag) else psf_cube_name(tag)
}
contour_df <- NULL
if (show_isophote_contours && file.exists(cube_path)) {
  cube_fits <- FITSio::readFITS(cube_path)
  cube <- cube_fits$imDat
  collapsed <- tryCatch(
    guara::collapse_white_light(cube, kclip = 2),
    error = function(e) apply(cube, c(1, 2), function(v) sum(pmax(v, 0), na.rm = TRUE))
  )
  collapsed <- asinh(pmax(collapsed, 0))
}

z_used_map <- infer_results_redshift(results_df, tag)
if (exists("cube_fits", inherits = FALSE) && "sfr" %in% names(property_df)) {
  area_df <- region_area_table(region_map, cube_fits, z_used_map)
  if (nrow(area_df)) {
    property_df <- property_df |>
      dplyr::left_join(area_df, by = "region_id") |>
      dplyr::mutate(
        sigma_sfr = ifelse(is.finite(.data$area_kpc2) & .data$area_kpc2 > 0,
                           .data$sfr / .data$area_kpc2,
                           NA_real_)
      ) |>
      dplyr::select(-tidyselect::any_of(c("n_pix", "area_kpc2")))
  }
}

property_names <- setdiff(names(property_df), "region_id")
if (!length(property_names)) {
  stop("No plottable properties found in results CSV.")
}

manifest_rows <- vector("list", length(property_names))

for (ii in seq_along(property_names)) {
  property_name <- property_names[[ii]]
  reg_field <- property_df |>
    transmute(region_id = .data$region_id, value = .data[[property_name]])

  pix_mat <- map_region_field_to_pixels(region_map, reg_field, property_name)
  pix_df <- mat_to_df(pix_mat, flip_y = FALSE, transpose = FALSE)
  pix_df <- crop_df_to_data(pix_df, pad = crop_padding_px)
  if (show_isophote_contours && exists("collapsed", inherits = FALSE)) {
    contour_df <- crop_matrix_to_df(collapsed, pix_df)
  }

  png_path <- file.path(figure_out_dir, paste0(property_name, "_smooth.png"))
  fits_path <- file.path(out_dirs$maps, paste0(property_name, "_smooth.fits"))
  csv_path <- file.path(out_dirs$maps, paste0(property_name, "_regions.csv"))

  p <- plot_univariate(
    df_pix = pix_df,
    property_name = property_name,
    contour_df = contour_df,
    show_legend = TRUE,
    transparent_bg = TRUE,
    legend_n = 4,
    legend_position = "bottom"
  )

  grDevices::png(png_path, width = 1800, height = 1800, res = 220, bg = "transparent")
  print(p)
  grDevices::dev.off()

  FITSio::writeFITSim(pix_mat, file = fits_path)
  readr::write_csv(reg_field, csv_path)

  manifest_rows[[ii]] <- tibble::tibble(
    tag = tag,
    photometry_mode = photometry_mode,
    property = property_name,
    property_scale = if (!is.null(property_specs[[property_name]])) property_specs[[property_name]]$scale else property_name,
    property_units = if (!is.null(property_specs[[property_name]])) property_specs[[property_name]]$units else NA_character_,
    source_results = fit_results_path,
    segmentation_fits = segmentation_path,
    output_png = png_path,
    output_fits = fits_path,
    output_regions_csv = csv_path
  )
}

manifest <- dplyr::bind_rows(manifest_rows)
readr::write_csv(manifest, file.path(out_dirs$maps, paste0(tag, "_paint_manifest.csv")))
readr::write_csv(manifest, file.path(figure_out_dir, "manifest.csv"))
readr::write_csv(
  tibble::tibble(
    property = property_names,
    scale = vapply(property_names, function(x) property_specs[[x]]$scale, character(1)),
    units = vapply(property_names, function(x) property_specs[[x]]$units, character(1))
  ),
  file.path(out_dirs$maps, paste0(tag, "_property_units.csv"))
)

message("Painted ", length(property_names), " property maps for ", tag)
message("Figure outputs: ", figure_out_dir)
message("Map products: ", out_dirs$maps)
