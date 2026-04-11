suppressPackageStartupMessages({
  library(FITSio)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(reticulate)
  library(ragg)
})

# ------------------------------------------------------------------
# User config
# ------------------------------------------------------------------
base_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation"
capivara_paper_dir <- "/Users/rd23aag/Documents/Capivara_MNRAS_paperI"
sagui_pkg_dir <- "/Users/rd23aag/Documents/GitHub/sagui"
output_root <- Sys.getenv("SAGUI_COMPARE_OUTPUT_DIR", unset = base_dir)
run_suffix <- Sys.getenv("SAGUI_COMPARE_SUFFIX", unset = "papermask")

tag <- Sys.getenv("SAGUI_COMPARE_TAG", unset = "sagui7")
ncomp <- suppressWarnings(as.integer(Sys.getenv("SAGUI_COMPARE_NCOMP", unset = "8")))
if (!is.finite(ncomp) || ncomp <= 1) {
  stop("SAGUI_COMPARE_NCOMP must be an integer greater than 1.")
}

collapse_kclip <- 2
mask_pretransform <- "none"
cluster_pretransform <- "none"
voronoi_mask_mode <- Sys.getenv("SAGUI_COMPARE_VORONOI_MASK", unset = "shared")
starlet_J <- 5L
starlet_scales <- 2:5
starlet_denoise_k <- 2L
cube_mode <- Sys.getenv("SAGUI_COMPARE_CUBE_MODE", unset = "raw")
cube_path_override <- Sys.getenv("SAGUI_COMPARE_CUBE_PATH", unset = "")

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

# Optional external region map (for example, a piXedfit bin map) with the
# same spatial dimensions as the FITS cube. Supported formats: FITS, RDS, CSV.
external_map_path <- NA_character_
external_label <- "piXedfit"

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
voronoi_py_path <- file.path(capivara_paper_dir, "voronoi_2d_binning.py")
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
  if (nbands == length(jwst_filters)) {
    return(unname(jwst_lambda_um[jwst_filters]) * 1e4)
  }

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

load_external_cluster_map <- function(path, dims) {
  if (!is.character(path) || length(path) != 1 || is.na(path) || !nzchar(path)) {
    return(NULL)
  }
  if (!file.exists(path)) {
    stop("External map not found: ", path)
  }

  ext <- tolower(tools::file_ext(path))
  if (ext == "fits") {
    obj <- FITSio::readFITS(path)
    map <- obj$imDat
  } else if (ext == "rds") {
    map <- readRDS(path)
  } else if (ext == "csv") {
    dat <- read.csv(path, stringsAsFactors = FALSE)
    needed <- c("row", "col", "cluster")
    if (!all(needed %in% names(dat))) {
      stop("CSV external map must contain columns: row, col, cluster")
    }
    map <- matrix(NA_real_, nrow = dims[1], ncol = dims[2])
    ok <- dat$row >= 1 & dat$row <= dims[1] & dat$col >= 1 & dat$col <= dims[2]
    map[cbind(dat$row[ok], dat$col[ok])] <- dat$cluster[ok]
  } else {
    stop("Unsupported external map format: ", ext)
  }

  if (!is.matrix(map) || !identical(dim(map), dims)) {
    stop("External map dimensions do not match the cube dimensions.")
  }

  map[!is.finite(map)] <- NA_real_
  map[map <= 0] <- NA_real_
  as.matrix(map)
}

source_voronoi_code <- local({
  loaded <- FALSE
  function(py_path) {
    if (loaded) {
      return(invisible(TRUE))
    }
    mpl_dir <- file.path(tempdir(), "mplconfig")
    dir.create(mpl_dir, recursive = TRUE, showWarnings = FALSE)
    Sys.setenv(MPLCONFIGDIR = mpl_dir)
    reticulate::source_python(py_path, envir = .GlobalEnv)
    loaded <<- TRUE
    invisible(TRUE)
  }
})

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

prepare_voronoi_data <- function(signal_map, support_mask) {
  stopifnot(is.matrix(signal_map), is.matrix(support_mask), identical(dim(signal_map), dim(support_mask)))

  n_rows <- nrow(signal_map)
  n_cols <- ncol(signal_map)

  signal_map[!support_mask] <- 0
  signal_map[!is.finite(signal_map) | signal_map <= 0] <- 0
  noise_map <- sqrt(signal_map)
  noise_map[!is.finite(noise_map) | noise_map <= 0] <- Inf

  valid <- signal_map > 0 & support_mask
  valid_flat <- as.vector(valid)
  coords <- expand.grid(y = seq_len(n_rows), x = seq_len(n_cols))

  list(
    valid_coords = coords[valid_flat, , drop = FALSE],
    signal_vector = as.vector(signal_map)[valid_flat],
    noise_vector = as.vector(noise_map)[valid_flat],
    signal_map = signal_map,
    noise_map = noise_map,
    valid_mask = valid,
    n_rows = n_rows,
    n_cols = n_cols
  )
}

perform_voronoi_binning_local <- function(voronoi_data, target_sn) {
  np <- reticulate::import("numpy")
  result <- voronoi_2d_binning(
    x = np$array(as.numeric(voronoi_data$valid_coords$x), dtype = "float64"),
    y = np$array(as.numeric(voronoi_data$valid_coords$y), dtype = "float64"),
    signal = np$array(as.numeric(voronoi_data$signal_vector), dtype = "float64"),
    noise = np$array(as.numeric(voronoi_data$noise_vector), dtype = "float64"),
    target_sn = as.numeric(target_sn),
    pixelsize = 1,
    plot = FALSE,
    quiet = TRUE
  )

  bin_ids <- as.integer(reticulate::py_to_r(result[[1]]))
  out <- matrix(NA_integer_, nrow = voronoi_data$n_rows, ncol = voronoi_data$n_cols)
  out[cbind(voronoi_data$valid_coords$y, voronoi_data$valid_coords$x)] <- bin_ids + 1L
  out
}

find_target_sn_local <- function(voronoi_data, target_ngroups, tol = 1L, max_iter = 30L) {
  sn_seed <- voronoi_data$signal_vector / voronoi_data$noise_vector
  sn_seed <- sn_seed[is.finite(sn_seed) & sn_seed > 0]
  lower_sn <- 0.5
  upper_sn <- max(20, stats::quantile(sn_seed, 0.99, na.rm = TRUE) * 5)
  best <- list(target_sn = NA_real_, ngroups = Inf, delta = Inf)

  upper_map <- perform_voronoi_binning_local(voronoi_data, upper_sn)
  upper_ngroups <- length(unique(stats::na.omit(as.vector(upper_map))))
  while (upper_ngroups > target_ngroups + tol && upper_sn < 1e5) {
    upper_sn <- upper_sn * 2
    upper_map <- perform_voronoi_binning_local(voronoi_data, upper_sn)
    upper_ngroups <- length(unique(stats::na.omit(as.vector(upper_map))))
  }

  for (iter in seq_len(max_iter)) {
    current_sn <- (lower_sn + upper_sn) / 2
    candidate_map <- perform_voronoi_binning_local(voronoi_data, current_sn)
    current_ngroups <- length(unique(stats::na.omit(as.vector(candidate_map))))
    delta <- abs(current_ngroups - target_ngroups)

    if (delta < best$delta) {
      best <- list(target_sn = current_sn, ngroups = current_ngroups, delta = delta)
    }
    if (delta <= tol) {
      return(list(target_sn = current_sn, ngroups = current_ngroups, iterations = iter))
    }

    if (current_ngroups > target_ngroups) {
      lower_sn <- current_sn
    } else {
      upper_sn <- current_sn
    }
  }

  list(
    target_sn = best$target_sn,
    ngroups = best$ngroups,
    iterations = max_iter
  )
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
                            title = NULL,
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
                                            title = NULL,
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

# ------------------------------------------------------------------
# Run comparison
# ------------------------------------------------------------------
if (!file.exists(cube_path)) {
  stop("Cube not found: ", cube_path)
}
if (!file.exists(voronoi_py_path)) {
  stop("Voronoi Python implementation not found: ", voronoi_py_path)
}
if (!dir.exists(sagui_pkg_dir) && !requireNamespace("sagui", quietly = TRUE)) {
  stop("Local sagui package directory not found and sagui is not installed.")
}

set.seed(42)
load_sagui_package(sagui_pkg_dir)
source_voronoi_code(voronoi_py_path)

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

# Match the published paper workflow for the foreground mask:
# white-light collapse on the original cube, starlet reconstruction on scales 2:5,
# soft thresholding with denoise_k = 2, then keep strictly positive pixels.
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
sagui_seg$starlet <- list(
  decomposition = starlet_dec,
  reconstruction = starlet_rec
)

voronoi_support <- switch(
  tolower(voronoi_mask_mode),
  shared = shared_mask,
  none = matrix(TRUE, nrow = nrow(shared_mask), ncol = ncol(shared_mask)),
  stop("Unsupported SAGUI_COMPARE_VORONOI_MASK: ", voronoi_mask_mode)
)

voronoi_data <- prepare_voronoi_data(collapsed_original, voronoi_support)
target_sn_result <- find_target_sn_local(voronoi_data, target_ngroups = ncomp)
voronoi_map <- perform_voronoi_binning_local(voronoi_data, target_sn_result$target_sn)

external_map <- load_external_cluster_map(external_map_path, dims = dim(shared_mask))

method_maps <- list(
  sagui = sagui_seg$cluster_map,
  voronoi = voronoi_map
)
if (!is.null(external_map)) {
  method_maps[[external_label]] <- external_map
}

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
summary_df$target_sn[summary_df$method == "voronoi"] <- target_sn_result$target_sn
summary_df$voronoi_mask_mode <- ifelse(summary_df$method == "voronoi", voronoi_mask_mode, "shared")

write.csv(summary_df, file.path(data_dir, "comparison_summary.csv"), row.names = FALSE)

for (method_name in names(method_maps)) {
  method_key <- sanitize_name(method_name)
  map <- method_maps[[method_name]]
  diag <- diagnostics[[method_name]]

  write.csv(diag$per_region,
            file.path(data_dir, paste0("region_metrics_", method_key, ".csv")),
            row.names = FALSE)

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

  sed_plot <- plot_cluster_sed_shape_variance(diagnostics = diag)

  save_plot_png(region_plot, file.path(figure_dir, paste0("regions_plot_", method_key, ".png")), width = 6, height = 6)
  save_plot_png(intravariance_plot, file.path(figure_dir, paste0("pixel_intravariance_", method_key, ".png")), width = 6, height = 6)
  save_plot_png(sed_plot, file.path(figure_dir, paste0("sed_combo_long_", method_key, ".png")), width = 14, height = 10)
}

mask_plot <- plot_scalar_map(shared_mask * 1, palette = c("white", "#0A1026"), quantiles = c(0, 1))
voronoi_support_plot <- plot_scalar_map(voronoi_support * 1, palette = c("white", "#0A1026"), quantiles = c(0, 1))
collapsed_plot <- plot_scalar_map(collapsed_original)
summary_plot <- plot_summary_metric(summary_df)

save_plot_png(mask_plot, file.path(figure_dir, "shared_starlet_mask.png"), width = 6, height = 6)
save_plot_png(mask_plot, file.path(figure_dir, "sagui_starlet_mask.png"), width = 6, height = 6)
save_plot_png(voronoi_support_plot, file.path(figure_dir, "voronoi_support.png"), width = 6, height = 6)
save_plot_png(collapsed_plot, file.path(figure_dir, "collapsed_white_light.png"), width = 6, height = 6)
save_plot_png(summary_plot, file.path(figure_dir, "comparison_metrics.png"), width = 10, height = 5)

message("Finished comparison demo for ", run_tag)
message("Figures: ", figure_dir)
message("Data: ", data_dir)
