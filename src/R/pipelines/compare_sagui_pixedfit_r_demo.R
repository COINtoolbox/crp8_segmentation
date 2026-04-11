suppressPackageStartupMessages({
  library(FITSio)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ragg)
  library(EBImage)
})

# ------------------------------------------------------------------
# User config
# ------------------------------------------------------------------
base_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation"
sagui_pkg_dir <- "/Users/rd23aag/Documents/GitHub/sagui"
output_root <- Sys.getenv("SAGUI_COMPARE_OUTPUT_DIR", unset = base_dir)
run_suffix <- Sys.getenv("SAGUI_COMPARE_SUFFIX", unset = "pixedfitr")

tag <- Sys.getenv("SAGUI_COMPARE_TAG", unset = "sagui10")
ncomp <- suppressWarnings(as.integer(Sys.getenv("SAGUI_COMPARE_NCOMP", unset = "12")))
if (!is.finite(ncomp) || ncomp <= 1) {
  stop("SAGUI_COMPARE_NCOMP must be an integer greater than 1.")
}

collapse_kclip <- 2
cluster_pretransform <- "none"
starlet_J <- 5L
starlet_scales <- 2:5
starlet_denoise_k <- 2L
cube_mode <- Sys.getenv("SAGUI_COMPARE_CUBE_MODE", unset = "raw")
cube_path_override <- Sys.getenv("SAGUI_COMPARE_CUBE_PATH", unset = "")

pixlike_blur_sigma <- 2
pixlike_otsu_factor <- 0.18
pixlike_min_area <- 30L
pixlike_dmin_bin <- 4
pixlike_del_r <- 2
pixlike_redc_chi2_limit <- 4

# JWST/NIRCam band order used in the sagui photometric cubes
jwst_filters <- c(
  "F090W", "F115W", "F150W", "F182M", "F200W",
  "F210M", "F277W", "F335M", "F356W", "F410M",
  "F430M", "F444W", "F460M", "F480M"
)

jwst_lambda_um <- c(
  F090W = 0.902, F115W = 1.154, F150W = 1.501, F182M = 1.842,
  F200W = 1.989, F210M = 2.099, F277W = 2.770, F335M = 3.365,
  F356W = 3.563, F410M = 4.082, F430M = 4.286, F444W = 4.421,
  F460M = 4.620, F480M = 4.828
)

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
resolve_compare_cube_path <- function(base_dir, tag, cube_mode = c("raw", "psfmatched"), override = "") {
  cube_mode <- match.arg(cube_mode)
  if (nzchar(override)) {
    return(override)
  }

  dashed <- sub("^sagui", "sagui-", gsub("_", "-", tag))
  candidates <- if (identical(cube_mode, "psfmatched")) {
    c(
      file.path(base_dir, "data/PSF_matched", paste0("datacube_", dashed, "_psfmatched.fits"))
    )
  } else {
    c(
      file.path(base_dir, "data/raw", paste0("datacube_", tag, ".fits")),
      file.path(base_dir, "data/raw", paste0("datacube_", dashed, ".fits"))
    )
  }

  path <- candidates[file.exists(candidates)][1]
  if (is.na(path)) {
    stop("Cube not found for tag=", tag, ", cube_mode=", cube_mode, ". Tried: ", paste(candidates, collapse = ", "))
  }
  path
}

cube_path <- resolve_compare_cube_path(base_dir, tag, cube_mode = cube_mode, override = cube_path_override)
run_tag <- if (nzchar(run_suffix)) {
  sprintf("%s_n%02d_%s", tag, ncomp, gsub("[^A-Za-z0-9]+", "_", tolower(run_suffix)))
} else {
  sprintf("%s_n%02d", tag, ncomp)
}
figure_dir <- file.path(output_root, "results/figures/sagui_compare", run_tag)
data_dir <- file.path(output_root, "results/segmentation/sagui_compare", run_tag)
flux_table_path <- file.path(base_dir, "results/flux_per_region", paste0("SED_flux_wide_", tag, ".csv"))

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

palette_van_gogh_div <- function(n = 256) {
  stops <- c(
    "#0A1026",
    "#182A6B",
    "#2F5F9A",
    "#4FA7A6",
    "#D8E1D3",
    "#F1C76A",
    "#A96A2A"
  )
  grDevices::colorRampPalette(stops, space = "Lab")(n)
}

sanitize_name <- function(x) {
  tolower(gsub("[^A-Za-z0-9]+", "_", x))
}

save_plot_png <- function(plot, path, width = 7, height = 7, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ragg::agg_png(path, width = width, height = height, units = "in", res = dpi, background = "white")
  print(plot)
  grDevices::dev.off()
}

safe_band_axis <- function(axDat, nbands) {
  vals <- tryCatch(FITSio::axVec(3, axDat), error = function(e) NULL)
  if (!is.null(vals) && length(vals) == nbands) {
    return(vals)
  }
  seq_len(nbands)
}

infer_wavelength_axis <- function(band_axis, nbands) {
  finite_axis <- band_axis[is.finite(band_axis)]
  if (length(finite_axis) == nbands && max(finite_axis, na.rm = TRUE) > 100) {
    return(as.numeric(band_axis))
  }
  if (length(finite_axis) == nbands && max(finite_axis, na.rm = TRUE) > 0 && max(finite_axis, na.rm = TRUE) < 20) {
    return(as.numeric(band_axis) * 1e4)
  }
  seq_len(nbands)
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

  unname(filters[nzchar(filters)])
}

load_sagui_package <- function(path) {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(path, export_all = FALSE, helpers = FALSE, attach_testthat = FALSE, quiet = TRUE)
    return(invisible(TRUE))
  }
  if (!requireNamespace("sagui", quietly = TRUE)) {
    stop("Need either pkgload or an installed sagui package to run this demo.")
  }
  invisible(TRUE)
}

row_scale_matrix <- function(x) {
  out <- t(apply(x, 1, function(v) {
    xf <- v[is.finite(v)]
    if (!length(xf)) {
      return(rep(0, length(v)))
    }
    center <- stats::median(xf)
    scale <- stats::mad(xf, center = center, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(scale) || scale <= 0) {
      scale <- stats::sd(xf, na.rm = TRUE)
    }
    if (!is.finite(scale) || scale <= 0) {
      scale <- 1
    }
    sv <- (v - center) / scale
    sv[!is.finite(sv)] <- 0
    sv
  }))
  storage.mode(out) <- "double"
  out
}

compute_cluster_snr <- function(cluster_map, signal_map, noise_map) {
  ids <- sort(unique(stats::na.omit(as.integer(cluster_map))))
  values <- vapply(ids, function(id) {
    idx <- which(cluster_map == id)
    s <- sum(signal_map[idx], na.rm = TRUE)
    n <- sqrt(sum(noise_map[idx]^2, na.rm = TRUE))
    if (!is.finite(n) || n <= 0) {
      return(NA_real_)
    }
    s / n
  }, numeric(1))
  stats::setNames(values, as.character(ids))
}

compute_segmentation_diagnostics <- function(cluster_map,
                                             cube,
                                             wavelength_axis,
                                             signal_map,
                                             noise_map,
                                             method_label) {
  spectra_all <- sagui::cube_to_matrix(cube)
  flat_map <- as.vector(cluster_map)
  valid <- !is.na(flat_map)
  labels <- as.integer(flat_map[valid])
  raw_spectra <- spectra_all[valid, , drop = FALSE]
  scaled_spectra <- row_scale_matrix(raw_spectra)

  pixel_values <- rep(NA_real_, sum(valid))
  region_ids <- sort(unique(labels))
  cluster_snr <- compute_cluster_snr(cluster_map, signal_map, noise_map)

  per_region <- vector("list", length(region_ids))
  names(per_region) <- as.character(region_ids)

  for (i in seq_along(region_ids)) {
    region_id <- region_ids[i]
    idx <- which(labels == region_id)
    region_spectra <- scaled_spectra[idx, , drop = FALSE]

    if (nrow(region_spectra) == 1L) {
      pixel_var <- 0
      avg_corr <- NA_real_
    } else {
      template <- apply(region_spectra, 2, stats::median, na.rm = TRUE)
      diffs <- sweep(region_spectra, 2, template, FUN = "-")
      pixel_var <- rowMeans(diffs^2, na.rm = TRUE)
      corr_mat <- suppressWarnings(stats::cor(t(region_spectra), use = "pairwise.complete.obs"))
      avg_corr <- mean(corr_mat[upper.tri(corr_mat)], na.rm = TRUE)
    }

    pixel_values[idx] <- pixel_var
    per_region[[i]] <- data.frame(
      method = method_label,
      region = region_id,
      n_pixels = length(idx),
      avg_corr = avg_corr,
      mean_pixel_intravariance = mean(pixel_var, na.rm = TRUE),
      median_pixel_intravariance = stats::median(pixel_var, na.rm = TRUE),
      cluster_snr = unname(cluster_snr[as.character(region_id)]),
      stringsAsFactors = FALSE
    )
  }

  intravariance_map <- matrix(NA_real_, nrow = nrow(cluster_map), ncol = ncol(cluster_map))
  intravariance_map[valid] <- pixel_values

  summary_row <- data.frame(
    method = method_label,
    n_regions = length(region_ids),
    valid_pixels = sum(valid),
    mean_pixels_per_region = mean(vapply(per_region, `[[`, numeric(1), "n_pixels")),
    median_pixels_per_region = stats::median(vapply(per_region, `[[`, numeric(1), "n_pixels")),
    mean_cluster_snr = mean(vapply(per_region, `[[`, numeric(1), "cluster_snr"), na.rm = TRUE),
    mean_avg_corr = mean(vapply(per_region, `[[`, numeric(1), "avg_corr"), na.rm = TRUE),
    median_avg_corr = stats::median(vapply(per_region, `[[`, numeric(1), "avg_corr"), na.rm = TRUE),
    mean_pixel_intravariance = mean(pixel_values, na.rm = TRUE),
    median_pixel_intravariance = stats::median(pixel_values, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  sed_long <- as.data.frame(raw_spectra)
  names(sed_long) <- paste0("band_", seq_len(ncol(raw_spectra)))
  sed_long$sed_id <- seq_len(nrow(sed_long))
  sed_long$cluster <- factor(labels, levels = region_ids)

  sed_long <- tidyr::pivot_longer(
    sed_long,
    cols = starts_with("band_"),
    names_to = "band_id",
    values_to = "flux"
  ) |>
    mutate(
      band_index = as.integer(sub("^band_", "", band_id)),
      wavelength = wavelength_axis[band_index]
    )

  sed_median <- sed_long |>
    group_by(cluster, wavelength) |>
    summarise(median_flux = stats::median(flux, na.rm = TRUE), .groups = "drop")

  cluster_summary <- bind_rows(per_region) |>
    mutate(
      cluster = factor(region, levels = region_ids),
      label = paste0(
        "n: ", n_pixels,
        "\nAvgCorr: ", sprintf("%.2f", avg_corr),
        "\nMeanVar: ", sprintf("%.2f", mean_pixel_intravariance)
      )
    )

  list(
    per_region = bind_rows(per_region),
    summary = summary_row,
    intravariance_map = intravariance_map,
    sed_long = sed_long,
    sed_median = sed_median,
    cluster_summary = cluster_summary
  )
}

plot_scalar_map <- function(mat,
                            palette = palette_van_gogh_div(256),
                            na_color = "white",
                            quantiles = c(0.02, 0.98)) {
  df <- expand.grid(
    row = seq_len(nrow(mat)),
    col = seq_len(ncol(mat))
  )
  df$value <- as.vector(mat)
  vals <- df$value[is.finite(df$value)]
  limits <- stats::quantile(vals, probs = quantiles, na.rm = TRUE)
  if (!all(is.finite(limits)) || limits[1] == limits[2]) {
    center <- vals[1]
    limits <- c(center - 0.5, center + 0.5)
  }

  ggplot(df, aes(x = row, y = col, fill = value)) +
    geom_raster() +
    coord_fixed() +
    scale_fill_gradientn(
      colours = palette,
      na.value = na_color,
      limits = limits,
      oob = scales::squish
    ) +
    theme_void() +
    labs(fill = NULL) +
    theme(
      legend.position = "bottom",
      legend.key.width = grid::unit(1.2, "cm"),
      legend.text = element_text(size = 9)
    )
}

plot_cluster_sed_shape_variance <- function(diagnostics,
                                            palette = "magma") {
  region_levels <- levels(diagnostics$cluster_summary$cluster)
  n_clusters <- length(region_levels)
  colors <- if (length(palette) == 1L) {
    viridis::viridis(n_clusters, option = palette)
  } else {
    grDevices::colorRampPalette(palette)(n_clusters)
  }
  names(colors) <- region_levels

  ggplot(diagnostics$sed_long, aes(x = wavelength, y = flux, group = sed_id)) +
    geom_line(color = "gray25", linewidth = 0.22, alpha = 0.12) +
    geom_line(
      data = diagnostics$sed_median,
      aes(x = wavelength, y = median_flux),
      inherit.aes = FALSE,
      color = "#C8862F",
      linewidth = 1.0
    ) +
    geom_point(
      data = diagnostics$sed_median,
      aes(x = wavelength, y = median_flux),
      inherit.aes = FALSE,
      size = 1.4,
      shape = 21,
      stroke = 0.25,
      fill = "#FFD77A",
      color = "black"
    ) +
    facet_wrap(~ cluster, scales = "free_y") +
    geom_text(
      data = diagnostics$cluster_summary,
      aes(x = Inf, y = Inf, label = label),
      hjust = 1.05,
      vjust = 1.1,
      inherit.aes = FALSE,
      size = 2.7,
      lineheight = 1.0
    ) +
    scale_x_continuous(
      name = expression(Wavelength ~ (ring(A))),
      breaks = scales::pretty_breaks(n = 5)
    ) +
    scale_y_continuous(
      name = "Flux",
      labels = scales::label_scientific(digits = 2)
    ) +
    theme_classic() +
    theme(
      strip.background = element_rect(fill = "gray92", color = "gray30"),
      strip.text = element_text(face = "bold", size = 9),
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 10, color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
      panel.spacing = grid::unit(0.45, "lines")
    )
}

plot_summary_metric <- function(summary_df) {
  long_df <- summary_df |>
    select(method, mean_cluster_snr, mean_avg_corr, mean_pixel_intravariance) |>
    tidyr::pivot_longer(
      cols = -method,
      names_to = "metric",
      values_to = "value"
    ) |>
    mutate(metric = factor(metric,
                           levels = c("mean_cluster_snr", "mean_avg_corr", "mean_pixel_intravariance"),
                           labels = c("Mean cluster S/N", "Mean intra-cluster corr.", "Mean pixel intravariance")))

  ggplot(long_df, aes(x = method, y = value, fill = method)) +
    geom_col(width = 0.72, color = NA) +
    facet_wrap(~ metric, scales = "free_y") +
    scale_fill_manual(values = palette_van_gogh_div(max(3, nrow(summary_df)))) +
    theme_bw() +
    labs(x = NULL, y = NULL) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "gray92", color = "gray30"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 15, hjust = 1)
    )
}

build_pixlike_mask <- function(collapse_map,
                               blur_sigma = 2,
                               otsu_factor = 0.18,
                               min_area = 30L) {
  img <- collapse_map
  img[!is.finite(img)] <- 0
  img <- img - min(img, na.rm = TRUE)
  peak <- max(img, na.rm = TRUE)
  if (!is.finite(peak) || peak <= 0) {
    stop("Collapsed image is empty; cannot build pixlike mask.")
  }
  img <- img / peak
  smoothed <- if (blur_sigma > 0) EBImage::gblur(img, sigma = blur_sigma) else img
  threshold <- min(1, EBImage::otsu(smoothed) * otsu_factor)

  mask <- smoothed > threshold
  mask <- EBImage::opening(mask, EBImage::makeBrush(3, shape = "disc"))
  mask <- EBImage::closing(mask, EBImage::makeBrush(5, shape = "disc"))

  label_map <- EBImage::bwlabel(mask)
  center_y <- ceiling(nrow(mask) / 2)
  center_x <- ceiling(ncol(mask) / 2)
  label_id <- label_map[center_y, center_x]

  if (!is.finite(label_id) || label_id <= 0) {
    brightest <- which.max(as.vector(smoothed))
    rc <- arrayInd(brightest, .dim = dim(smoothed))
    label_id <- label_map[rc[1], rc[2]]
  }

  if (!is.finite(label_id) || label_id <= 0) {
    return(mask * FALSE)
  }

  keep <- label_map == label_id
  if (sum(keep) < min_area) {
    return(mask * FALSE)
  }

  keep
}

estimate_band_errors <- function(cube, gal_region) {
  nbands <- dim(cube)[3]
  errs <- numeric(nbands)

  for (bb in seq_len(nbands)) {
    img <- cube[, , bb]
    background <- img[!gal_region & is.finite(img)]
    if (length(background) < 50) {
      background <- as.vector(img[is.finite(img)])
    }
    center <- stats::median(background, na.rm = TRUE)
    sig <- stats::mad(background, center = center, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(sig) || sig <= 0) {
      sig <- stats::sd(background, na.rm = TRUE)
    }
    if (!is.finite(sig) || sig <= 0) {
      sig <- max(abs(center), 1e-6)
    }
    errs[bb] <- sig
  }

  errs
}

make_error_cube <- function(cube, band_errors) {
  arr <- array(0, dim = dim(cube))
  for (bb in seq_len(dim(cube)[3])) {
    arr[, , bb] <- band_errors[bb]
  }
  arr
}

extract_spectra <- function(cube, coords) {
  nbands <- dim(cube)[3]
  out <- matrix(NA_real_, nrow = nrow(coords), ncol = nbands)
  for (bb in seq_len(nbands)) {
    out[, bb] <- cube[, , bb][coords]
  }
  out
}

redchi2_two_seds <- function(sed1_f, sed1_ferr, sed2_f, sed2_ferr) {
  denom <- sed1_ferr^2 + sed2_ferr^2
  denom[!is.finite(denom) | denom <= 0] <- NA_real_
  top <- sum(sed2_f * sed1_f / denom, na.rm = TRUE)
  bottom <- sum(sed1_f^2 / denom, na.rm = TRUE)
  if (!is.finite(bottom) || bottom <= 0) {
    return(Inf)
  }
  norm <- top / bottom
  sum((sed2_f - norm * sed1_f)^2 / denom, na.rm = TRUE) / length(sed1_f)
}

prepare_flux_correction <- function(cube, gal_region) {
  flux_corr <- cube
  nbands <- dim(cube)[3]
  for (bb in seq_len(nbands)) {
    img <- cube[, , bb]
    pos_vals <- img[gal_region & is.finite(img) & img > 0]
    if (length(pos_vals)) {
      lowest <- min(pos_vals, na.rm = TRUE)
      idx <- which(gal_region & is.finite(img) & img < 0)
      flux_corr[, , bb][idx] <- lowest
    }
  }
  flux_corr
}

estimate_systematic_factor <- function(flux_corr,
                                       err_cube,
                                       gal_region,
                                       ref_band) {
  valid_coords <- which(gal_region, arr.ind = TRUE)
  ref_vals <- flux_corr[, , ref_band][gal_region]
  rc <- valid_coords[which.max(ref_vals), , drop = FALSE]
  yc <- rc[1, 1]
  xc <- rc[1, 2]

  yseq <- seq(max(1, yc - 2), min(nrow(gal_region), yc + 2))
  xseq <- seq(max(1, xc - 2), min(ncol(gal_region), xc + 2))
  neighbors <- expand.grid(row = yseq, col = xseq)
  neighbors <- neighbors[!(neighbors$row == yc & neighbors$col == xc), , drop = FALSE]
  neighbors <- neighbors[gal_region[cbind(neighbors$row, neighbors$col)], , drop = FALSE]

  if (!nrow(neighbors)) {
    return(0.03)
  }

  sed1_f <- flux_corr[yc, xc, ]
  factor <- 0.01
  repeat {
    sed1_ferr <- sqrt(err_cube[yc, xc, ]^2 + (factor * sed1_f)^2)
    chi2_vals <- apply(neighbors, 1, function(coord) {
      sed2_f <- flux_corr[coord[1], coord[2], ]
      sed2_ferr <- sqrt(err_cube[coord[1], coord[2], ]^2 + (factor * sed2_f)^2)
      redchi2_two_seds(sed1_f = sed1_f, sed1_ferr = sed1_ferr, sed2_f = sed2_f, sed2_ferr = sed2_ferr)
    })
    chi2_med <- stats::median(chi2_vals[is.finite(chi2_vals)], na.rm = TRUE)
    if (!is.finite(chi2_med) || chi2_med <= 2 || factor >= 0.50) {
      break
    }
    factor <- factor + 0.01
  }
  factor
}

pix_coords_in_radius <- function(center_row, center_col, radius, dims, lower = 0) {
  yseq <- seq(max(1, floor(center_row - radius - 1)), min(dims[1], ceiling(center_row + radius + 1)))
  xseq <- seq(max(1, floor(center_col - radius - 1)), min(dims[2], ceiling(center_col + radius + 1)))
  grid <- expand.grid(row = yseq, col = xseq)
  dist <- sqrt((grid$row - center_row)^2 + (grid$col - center_col)^2)
  grid[dist > lower & dist <= radius, , drop = FALSE]
}

summed_flux_and_error <- function(cube, err_cube, coords) {
  nbands <- dim(cube)[3]
  flux <- numeric(nbands)
  err2 <- numeric(nbands)
  for (bb in seq_len(nbands)) {
    vals <- cube[, , bb][as.matrix(coords)]
    errs <- err_cube[, , bb][as.matrix(coords)]
    flux[bb] <- sum(vals, na.rm = TRUE)
    err2[bb] <- sum(errs^2, na.rm = TRUE)
  }
  list(flux = flux, err2 = err2)
}

filter_ring_by_shape <- function(ring_coords,
                                 center_coord,
                                 flux_corr,
                                 err_cube_corr,
                                 redc_chi2_limit) {
  if (!nrow(ring_coords)) {
    return(ring_coords)
  }

  sed1_f <- flux_corr[center_coord[1], center_coord[2], ]
  sed1_ferr <- err_cube_corr[center_coord[1], center_coord[2], ]
  keep <- apply(ring_coords, 1, function(coord) {
    sed2_f <- flux_corr[coord[1], coord[2], ]
    sed2_ferr <- err_cube_corr[coord[1], coord[2], ]
    redchi2_two_seds(sed1_f, sed1_ferr, sed2_f, sed2_ferr) <= redc_chi2_limit
  })
  ring_coords[keep, , drop = FALSE]
}

pixlike_pixel_binning <- function(cube,
                                  gal_region,
                                  ref_band = NULL,
                                  Dmin_bin = 4,
                                  SNR = 5,
                                  redc_chi2_limit = 4,
                                  del_r = 2) {
  dims <- dim(cube)
  nbands <- dims[3]

  if (is.null(ref_band)) {
    ref_band <- if (nbands == 1) 1L else floor((nbands - 1) / 2) + 1L
  }
  ref_band <- as.integer(ref_band)

  if (length(SNR) == 1L) {
    SN_threshold <- rep(as.numeric(SNR), nbands)
  } else if (length(SNR) == nbands) {
    SN_threshold <- as.numeric(SNR)
  } else {
    stop("SNR must be length 1 or match the number of bands.")
  }

  band_errors <- estimate_band_errors(cube, gal_region)
  err_cube <- make_error_cube(cube, band_errors)
  flux_corr <- prepare_flux_correction(cube, gal_region)
  systematic_factor <- estimate_systematic_factor(flux_corr, err_cube, gal_region, ref_band)
  err_cube_corr <- sqrt(err_cube^2 + (systematic_factor * flux_corr)^2)

  pixbin_map <- matrix(0L, nrow = dims[1], ncol = dims[2])
  total_npix <- sum(gal_region)
  assigned_npix <- 0L
  count_bin <- 0L

  repeat {
    remaining_mask <- gal_region & pixbin_map == 0L
    if (!any(remaining_mask)) {
      break
    }

    remaining_coords <- which(remaining_mask, arr.ind = TRUE)
    ref_vals <- cube[, , ref_band][remaining_mask]
    center_idx <- which.max(ref_vals)
    center_coord <- remaining_coords[center_idx, , drop = FALSE]
    bin_rad <- 0.5 * Dmin_bin

    candidate <- pix_coords_in_radius(center_coord[1, 1], center_coord[1, 2], bin_rad, dims[1:2])
    candidate <- candidate[remaining_mask[as.matrix(candidate)], , drop = FALSE]
    if (!nrow(candidate)) {
      pixbin_map[center_coord[1, 1], center_coord[1, 2]] <- count_bin + 1L
      count_bin <- count_bin + 1L
      assigned_npix <- assigned_npix + 1L
      next
    }

    totals <- summed_flux_and_error(cube, err_cube, candidate)
    total_snr <- totals$flux / sqrt(totals$err2)

    if (all(total_snr >= SN_threshold, na.rm = FALSE)) {
      count_bin <- count_bin + 1L
      pixbin_map[as.matrix(candidate)] <- count_bin
      assigned_npix <- assigned_npix + nrow(candidate)
      next
    }

    cumulative <- candidate
    cumulative_totals <- totals

    repeat {
      rmin <- bin_rad
      rmax <- rmin + del_r
      ring_coords <- pix_coords_in_radius(center_coord[1, 1], center_coord[1, 2], rmax, dims[1:2], lower = rmin)
      ring_coords <- ring_coords[remaining_mask[as.matrix(ring_coords)], , drop = FALSE]

      if (nrow(ring_coords)) {
        ring_coords <- filter_ring_by_shape(
          ring_coords = ring_coords,
          center_coord = center_coord[1, ],
          flux_corr = flux_corr,
          err_cube_corr = err_cube_corr,
          redc_chi2_limit = redc_chi2_limit
        )
      }

      if (nrow(ring_coords)) {
        cumulative <- unique(rbind(cumulative, ring_coords))
        ring_totals <- summed_flux_and_error(cube, err_cube, ring_coords)
        cumulative_totals$flux <- cumulative_totals$flux + ring_totals$flux
        cumulative_totals$err2 <- cumulative_totals$err2 + ring_totals$err2
      }

      total_snr <- cumulative_totals$flux / sqrt(cumulative_totals$err2)
      if (all(total_snr >= SN_threshold, na.rm = FALSE)) {
        count_bin <- count_bin + 1L
        pixbin_map[as.matrix(cumulative)] <- count_bin
        assigned_npix <- assigned_npix + nrow(cumulative)
        break
      }

      remaining_coords <- which(remaining_mask, arr.ind = TRUE)
      remaining_totals <- summed_flux_and_error(cube, err_cube, remaining_coords)
      remaining_snr <- remaining_totals$flux / sqrt(remaining_totals$err2)

      if (nrow(remaining_coords) == nrow(cumulative) ||
          !all(remaining_snr >= SN_threshold, na.rm = FALSE)) {
        count_bin <- count_bin + 1L
        pixbin_map[as.matrix(remaining_coords)] <- count_bin
        assigned_npix <- assigned_npix + nrow(remaining_coords)
        break
      }

      if (assigned_npix + nrow(cumulative) >= total_npix) {
        count_bin <- count_bin + 1L
        pixbin_map[as.matrix(cumulative)] <- count_bin
        assigned_npix <- assigned_npix + nrow(cumulative)
        break
      }

      bin_rad <- bin_rad + del_r
      if (bin_rad > max(dims[1:2])) {
        count_bin <- count_bin + 1L
        pixbin_map[as.matrix(cumulative)] <- count_bin
        assigned_npix <- assigned_npix + nrow(cumulative)
        break
      }
    }
  }

  pixbin_map[pixbin_map == 0L] <- NA_integer_
  list(
    cluster_map = pixbin_map,
    band_errors = band_errors,
    systematic_factor = systematic_factor,
    nbins = length(unique(stats::na.omit(as.vector(pixbin_map))))
  )
}

find_target_snr_pixlike <- function(cube,
                                    gal_region,
                                    target_ngroups,
                                    Dmin_bin = 4,
                                    redc_chi2_limit = 4,
                                    del_r = 2,
                                    lower = 1,
                                    upper = 40,
                                    tol = 1L,
                                    max_iter = 10L) {
  cache <- new.env(parent = emptyenv())

  run_once <- function(sn_value) {
    key <- sprintf("%.6f", sn_value)
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }
    res <- pixlike_pixel_binning(
      cube = cube,
      gal_region = gal_region,
      Dmin_bin = Dmin_bin,
      SNR = sn_value,
      redc_chi2_limit = redc_chi2_limit,
      del_r = del_r
    )
    assign(key, res, envir = cache)
    res
  }

  lower_res <- run_once(lower)
  upper_res <- run_once(upper)
  while (upper_res$nbins > target_ngroups + tol && upper < 1000) {
    upper <- upper * 1.5
    upper_res <- run_once(upper)
  }

  sn_linear <- seq(lower, upper, length.out = max(12L, max_iter * 3L))
  sn_log <- exp(seq(log(max(lower, 1e-3)), log(max(upper, lower + 1e-3)), length.out = max(12L, max_iter * 3L)))
  candidates <- sort(unique(round(c(lower, upper, sn_linear, sn_log), 6)))

  candidate_results <- lapply(candidates, run_once)
  candidate_bins <- vapply(candidate_results, `[[`, numeric(1), "nbins")
  deltas <- abs(candidate_bins - target_ngroups)

  best_idx <- which(deltas == min(deltas))
  if (length(best_idx) > 1L) {
    best_idx <- best_idx[which.min(abs(candidate_bins[best_idx] - target_ngroups))]
  }
  best_idx <- best_idx[1]

  list(
    target_sn = candidates[best_idx],
    result = candidate_results[[best_idx]],
    iterations = length(candidates),
    search_table = data.frame(
      target_sn = candidates,
      nbins = candidate_bins,
      delta = deltas
    )
  )
}

# ------------------------------------------------------------------
# Run comparison
# ------------------------------------------------------------------
if (!file.exists(cube_path)) {
  stop("Cube not found: ", cube_path)
}
if (!dir.exists(sagui_pkg_dir) && !requireNamespace("sagui", quietly = TRUE)) {
  stop("Local sagui package directory not found and sagui is not installed.")
}

set.seed(42)
load_sagui_package(sagui_pkg_dir)

cube_fits <- FITSio::readFITS(cube_path)
cube <- cube_fits$imDat
band_axis <- safe_band_axis(cube_fits$axDat, dim(cube)[3])
filters_for_tag <- extract_header_filters(cube_fits)
wavelength_axis <- if (length(filters_for_tag) == dim(cube)[3]) {
  unname(jwst_lambda_um[filters_for_tag]) * 1e4
} else {
  infer_wavelength_axis(band_axis, dim(cube)[3])
}
collapse_fn <- function(x) guara::collapse_white_light(x, kclip = collapse_kclip)

collapsed_original <- guara::collapse_white_light(cube, kclip = collapse_kclip)
starlet_dec <- guara::starlet_mask(collapsed_original, J = starlet_J)
starlet_rec <- guara::starlet_reconstruct(
  starlet_dec,
  keep_scales = starlet_scales,
  include_coarse = FALSE,
  denoise_k = starlet_denoise_k,
  mode = "soft"
)
shared_mask <- is.finite(starlet_rec) & (starlet_rec > 0)

masked_cube <- guara::mask_cube(cube, shared_mask, mode = "na")
masked_fits <- cube_fits
masked_fits$imDat <- masked_cube

sagui_seg <- sagui::segment_regions(
  input = masked_fits,
  Ncomp = ncomp,
  cluster_pretransform = cluster_pretransform,
  use_starlet_mask = FALSE,
  collapse_fn = collapse_fn,
  mask_mode = "na"
)
sagui_seg$mask <- shared_mask
sagui_seg$collapsed <- collapsed_original

pixlike_mask <- build_pixlike_mask(
  collapse_map = collapsed_original,
  blur_sigma = pixlike_blur_sigma,
  otsu_factor = pixlike_otsu_factor,
  min_area = pixlike_min_area
)

if (!any(pixlike_mask)) {
  stop("pixlike foreground mask is empty.")
}

pixlike_target <- find_target_snr_pixlike(
  cube = cube,
  gal_region = pixlike_mask,
  target_ngroups = ncomp,
  Dmin_bin = pixlike_dmin_bin,
  redc_chi2_limit = pixlike_redc_chi2_limit,
  del_r = pixlike_del_r
)

pixlike_res <- pixlike_target$result
method_maps <- list(
  sagui = sagui_seg$cluster_map,
  pixedfit_r = pixlike_res$cluster_map
)

diagnostics <- purrr::imap(method_maps, function(map, method_name) {
  compute_segmentation_diagnostics(
    cluster_map = map,
    cube = cube_fits,
    wavelength_axis = wavelength_axis,
    signal_map = collapsed_original,
    noise_map = sqrt(pmax(collapsed_original, 0)),
    method_label = method_name
  )
})

summary_df <- bind_rows(purrr::map(diagnostics, "summary"))
summary_df$target_sn <- NA_real_
summary_df$target_sn[summary_df$method == "pixedfit_r"] <- pixlike_target$target_sn
summary_df$mask_mode <- ifelse(summary_df$method == "sagui", "paperexact_starlet", "pixlike_foreground")
summary_df$systematic_factor <- NA_real_
summary_df$systematic_factor[summary_df$method == "pixedfit_r"] <- pixlike_res$systematic_factor

write.csv(summary_df, file.path(data_dir, "comparison_summary.csv"), row.names = FALSE)

for (method_name in names(method_maps)) {
  method_key <- sanitize_name(method_name)
  map <- method_maps[[method_name]]
  diag <- diagnostics[[method_name]]

  write.csv(
    diag$per_region,
    file.path(data_dir, paste0("region_metrics_", method_key, ".csv")),
    row.names = FALSE
  )

  map_to_write <- map
  map_to_write[!is.finite(map_to_write)] <- 0
  FITSio::writeFITSim(map_to_write, file = file.path(data_dir, paste0("region_map_", method_key, ".fits")))
  FITSio::writeFITSim(diag$intravariance_map, file = file.path(data_dir, paste0("pixel_intravariance_", method_key, ".fits")))

  region_plot <- sagui::plot_region_map(
    list(cluster_map = map),
    palette = palette_van_gogh_div(length(unique(stats::na.omit(as.vector(map))))),
    border_color = "black",
    border_linewidth = 0.9,
    background_color = "white"
  )
  intravariance_plot <- plot_scalar_map(diag$intravariance_map)
  sed_plot <- plot_cluster_sed_shape_variance(diag)

  save_plot_png(region_plot, file.path(figure_dir, paste0("regions_plot_", method_key, ".png")), width = 6, height = 6)
  save_plot_png(intravariance_plot, file.path(figure_dir, paste0("pixel_intravariance_", method_key, ".png")), width = 6, height = 6)
  save_plot_png(sed_plot, file.path(figure_dir, paste0("sed_combo_long_", method_key, ".png")), width = 14, height = 10)
}

sagui_mask_plot <- plot_scalar_map(shared_mask * 1, palette = c("white", "#0A1026"), quantiles = c(0, 1))
pixlike_mask_plot <- plot_scalar_map(pixlike_mask * 1, palette = c("white", "#0A1026"), quantiles = c(0, 1))
collapsed_plot <- plot_scalar_map(collapsed_original)
summary_plot <- plot_summary_metric(summary_df)

save_plot_png(sagui_mask_plot, file.path(figure_dir, "sagui_starlet_mask.png"), width = 6, height = 6)
save_plot_png(pixlike_mask_plot, file.path(figure_dir, "pixedfit_r_foreground_mask.png"), width = 6, height = 6)
save_plot_png(collapsed_plot, file.path(figure_dir, "collapsed_white_light.png"), width = 6, height = 6)
save_plot_png(summary_plot, file.path(figure_dir, "comparison_metrics.png"), width = 10, height = 5)

message("Finished R-translated piXedfit-style comparison for ", run_tag)
message("Figures: ", figure_dir)
message("Data: ", data_dir)
