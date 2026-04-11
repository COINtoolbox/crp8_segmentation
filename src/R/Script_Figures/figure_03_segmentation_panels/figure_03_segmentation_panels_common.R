#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(guara)
  library(capivara)
  library(FITSio)
  library(ragg)
  library(EBImage)
})

set.seed(42)

figure_id <- "figure_03_segmentation_panels"
base_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation"

panel_tags <- c(
  "sagui1_2_3",
  "sagui4",
  "sagui5_6",
  "sagui7",
  "sagui8",
  "sagui9",
  "sagui10"
)

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

variant_slug <- function(mode) {
  switch(
    mode,
    original = "original",
    psfmatched = "psfmatched",
    copula = "copula_gaussian",
    contourlet_psfmatched = "contourlet_psfmatched",
    stop("Unknown mode: ", mode)
  )
}

resolve_input_cube <- function(tag, mode) {
  if (identical(mode, "original") || identical(mode, "copula")) {
    dashed <- sub("^sagui", "sagui-", gsub("_", "-", tag))
    candidates <- c(
      file.path(base_dir, "data/raw", paste0("datacube_", tag, ".fits")),
      file.path(base_dir, "data/raw", paste0("datacube_", dashed, ".fits"))
    )
  } else if (identical(mode, "psfmatched") || identical(mode, "contourlet_psfmatched")) {
    dashed <- sub("^sagui", "sagui-", gsub("_", "-", tag))
    candidates <- c(
      file.path(base_dir, "data/PSF_matched", paste0("datacube_", dashed, "_psfmatched.fits"))
    )
  } else {
    stop("Unknown mode: ", mode)
  }

  path <- candidates[file.exists(candidates)][1]
  if (is.na(path)) {
    stop("Input cube not found for ", tag, ". Tried: ", paste(candidates, collapse = ", "))
  }

  path
}

resolve_raw_cube <- function(tag) {
  dashed <- sub("^sagui", "sagui-", gsub("_", "-", tag))
  candidates <- c(
    file.path(base_dir, "data/raw", paste0("datacube_", tag, ".fits")),
    file.path(base_dir, "data/raw", paste0("datacube_", dashed, ".fits"))
  )
  path <- candidates[file.exists(candidates)][1]
  if (is.na(path)) {
    stop("Raw cube not found for ", tag, ". Tried: ", paste(candidates, collapse = ", "))
  }
  path
}

binary_brush <- function(size = 3L, shape = "disc") {
  EBImage::makeBrush(size, shape = shape, step = FALSE) > 0
}

image_to_matrix <- function(x) {
  arr <- as.array(x)
  if (length(dim(arr)) <= 2L) {
    return(arr)
  }
  arr[, , 1]
}

close_open_mask <- function(mask, close_size = 5L, open_size = 3L) {
  img <- EBImage::Image(mask * 1)
  closed <- EBImage::closing(img, binary_brush(close_size))
  opened <- EBImage::opening(closed, binary_brush(open_size))
  image_to_matrix(opened) > 0
}

oriented_gaussian_kernel <- function(theta_deg, sigma_major, sigma_minor) {
  radius <- ceiling(3 * max(sigma_major, sigma_minor))
  grid <- seq(-radius, radius)
  xy <- expand.grid(x = grid, y = grid)
  theta <- theta_deg * pi / 180
  xr <- xy$x * cos(theta) + xy$y * sin(theta)
  yr <- -xy$x * sin(theta) + xy$y * cos(theta)
  ker <- exp(-0.5 * ((xr / sigma_major)^2 + (yr / sigma_minor)^2))
  ker <- matrix(ker, nrow = length(grid), ncol = length(grid), byrow = FALSE)
  ker <- ker / sum(ker)
  ker - mean(ker)
}

contourlet_like_response <- function(image, angles_deg, scales) {
  nr <- nrow(image)
  nc <- ncol(image)
  best <- matrix(0, nrow = nr, ncol = nc)

  for (scale in scales) {
    sigma_major <- scale[1]
    sigma_minor <- scale[2]
    for (theta in angles_deg) {
      ker <- oriented_gaussian_kernel(theta, sigma_major, sigma_minor)
      resp <- EBImage::filter2(image, ker)
      best <- pmax(best, resp, 0, na.rm = TRUE)
    }
  }

  best
}

robust_threshold <- function(x, k) {
  vals <- x[is.finite(x)]
  med <- stats::median(vals, na.rm = TRUE)
  madv <- stats::mad(vals, center = med, constant = 1, na.rm = TRUE)
  med + k * madv
}

positive_quantile <- function(x, prob, fallback_k = 2) {
  vals <- x[is.finite(x) & x > 0]
  if (length(vals) >= 10) {
    return(as.numeric(stats::quantile(vals, probs = prob, na.rm = TRUE, names = FALSE)))
  }
  robust_threshold(x, fallback_k)
}

grow_mask_constrained <- function(seed, allowed, max_iter = 12L) {
  brush <- binary_brush(3L)
  grown <- seed

  for (iter in seq_len(max_iter)) {
    dil <- EBImage::dilate(EBImage::Image(grown * 1), brush)
    next_mask <- (image_to_matrix(dil) > 0) & allowed
    next_mask <- next_mask | seed
    if (identical(next_mask, grown)) {
      break
    }
    grown <- next_mask
  }

  grown
}

empirical_copula_cube <- function(cube, gaussian = TRUE, eps = 1e-6) {
  dims <- dim(cube)
  mat <- matrix(cube, nrow = dims[1] * dims[2], ncol = dims[3])
  out <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))

  for (j in seq_len(ncol(mat))) {
    col <- mat[, j]
    ok <- is.finite(col)
    n_ok <- sum(ok)
    if (!n_ok) {
      next
    }
    ranks <- rank(col[ok], ties.method = "average")
    u <- (ranks - 0.5) / n_ok
    u <- pmin(pmax(u, eps), 1 - eps)
    out[ok, j] <- if (gaussian) stats::qnorm(u) else u
  }

  array(out, dim = dims, dimnames = dimnames(cube))
}

build_starlet_mask <- function(cube, collapse_kclip = 2, J = 5L, keep_scales = 2:5, denoise_k = 2L) {
  collapsed <- collapse_white_light(cube, kclip = collapse_kclip)
  dec <- guara::starlet_mask(collapsed, J = J)
  rec <- starlet_reconstruct(
    dec,
    keep_scales = keep_scales,
    include_coarse = FALSE,
    denoise_k = denoise_k,
    mode = "soft"
  )

  list(
    collapsed = collapsed,
    mask = is.finite(rec) & rec > 0
  )
}

build_copula_mask <- function(raw_cube, collapse_kclip = 2, J = 5L, keep_scales = 2:5, denoise_k = 2L) {
  transformed_cube <- empirical_copula_cube(raw_cube, gaussian = TRUE)
  collapsed_copula <- collapse_white_light(transformed_cube, kclip = collapse_kclip)
  dec <- guara::starlet_mask(collapsed_copula, J = J)
  rec <- starlet_reconstruct(
    dec,
    keep_scales = keep_scales,
    include_coarse = FALSE,
    denoise_k = denoise_k,
    mode = "soft"
  )

  list(
    collapsed = collapsed_copula,
    mask = is.finite(rec) & rec > 0,
    transformed_cube = transformed_cube
  )
}

build_contourlet_mask <- function(cube, collapse_kclip = 2) {
  collapsed <- collapse_white_light(cube, kclip = collapse_kclip)
  dec <- guara::starlet_mask(collapsed, J = 5L)
  rec <- starlet_reconstruct(
    dec,
    keep_scales = 2:5,
    include_coarse = FALSE,
    denoise_k = 2L,
    mode = "soft"
  )
  seed_mask <- is.finite(rec) & rec > 0

  white_light_prepped <- asinh(pmax(collapsed, 0))
  background_map <- EBImage::gblur(white_light_prepped, sigma = 8)
  contrast_map <- pmax(white_light_prepped - background_map, 0)

  arm_score <- contourlet_like_response(
    image = contrast_map,
    angles_deg = seq(0, 165, by = 15),
    scales = list(
      c(4.5, 1.2),
      c(7.5, 1.8),
      c(11.0, 2.6)
    )
  )

  seed_dilated <- image_to_matrix(
    EBImage::dilate(EBImage::Image(seed_mask * 1), binary_brush(5L))
  ) > 0
  outside_seed <- !seed_dilated & is.finite(arm_score)

  score_low <- positive_quantile(arm_score[outside_seed], 0.92, fallback_k = 1.6)
  score_high <- positive_quantile(arm_score[outside_seed], 0.98, fallback_k = 2.5)
  contrast_low <- positive_quantile(contrast_map[outside_seed], 0.80, fallback_k = 1.2)
  contrast_high <- positive_quantile(contrast_map[outside_seed], 0.93, fallback_k = 1.8)

  strong_arm <- (arm_score > score_high) & (contrast_map > contrast_high)
  allowed_growth <- seed_mask | ((arm_score > score_low) & (contrast_map > contrast_low))
  allowed_growth <- close_open_mask(allowed_growth, close_size = 3L, open_size = 3L)

  grown <- grow_mask_constrained(
    seed = seed_mask | strong_arm,
    allowed = allowed_growth,
    max_iter = 8L
  )
  grown <- close_open_mask(grown, close_size = 3L, open_size = 3L)

  list(
    collapsed = collapsed,
    mask = grown,
    seed_mask = seed_mask,
    added_pixels = grown & !seed_mask
  )
}

build_segmentation <- function(cube, mask, N = 40L) {
  cube_na <- guara::mask_cube(cube, mask, mode = "na")
  seg <- capivara::segment(list(imDat = cube_na), N = N)
  seg$mask <- mask
  seg
}

write_cluster_fits <- function(seg, path) {
  out <- seg$cluster_map
  out[!is.finite(out)] <- 0L
  FITSio::writeFITSim(out, file = path, header = seg$header, axDat = seg$axDat)
}

save_segmentation_png <- function(seg, path, N = 40L) {
  plot_obj <- plot_cluster_voronoi_style(
    seg,
    palette = palette_van_gogh_div(N),
    border_color = "black",
    border_linewidth = 1,
    background_color = "transparent"
  )

  ragg::agg_png(path, width = 2000, height = 2000, res = 300, background = "transparent")
  print(plot_obj)
  grDevices::dev.off()
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

resize_matrix <- function(x, width, height) {
  image_to_matrix(EBImage::resize(EBImage::Image(x), w = width, h = height))
}

match_template_ncc <- function(search, template) {
  search <- as.matrix(search)
  template <- as.matrix(template)
  sh <- dim(search)
  th <- dim(template)
  nr <- sh[1] - th[1] + 1L
  nc <- sh[2] - th[2] + 1L
  if (nr < 1L || nc < 1L) {
    stop("Template is larger than search image.")
  }

  template0 <- template - mean(template, na.rm = TRUE)
  template_norm <- sqrt(sum(template0^2, na.rm = TRUE))
  best_score <- -Inf
  best_row <- 1L
  best_col <- 1L

  for (row in seq_len(nr)) {
    for (col in seq_len(nc)) {
      patch <- search[row:(row + th[1] - 1L), col:(col + th[2] - 1L)]
      patch0 <- patch - mean(patch, na.rm = TRUE)
      patch_norm <- sqrt(sum(patch0^2, na.rm = TRUE))
      score <- if (template_norm == 0 || patch_norm == 0) {
        -Inf
      } else {
        sum(template0 * patch0, na.rm = TRUE) / (template_norm * patch_norm)
      }

      if (is.finite(score) && score > best_score) {
        best_score <- score
        best_row <- row
        best_col <- col
      }
    }
  }

  list(
    score = best_score,
    row = best_row,
    col = best_col,
    template_dim = dim(template),
    search_dim = dim(search)
  )
}

remap_psf_labels_to_raw <- function(labels_psf, raw_fits_obj, psf_fits_obj, collapse_kclip = 2) {
  raw_cube <- raw_fits_obj$imDat
  psf_cube <- psf_fits_obj$imDat

  raw_white <- asinh(pmax(collapse_white_light(raw_cube, kclip = collapse_kclip), 0))
  psf_white <- asinh(pmax(collapse_white_light(psf_cube, kclip = collapse_kclip), 0))

  raw_scale <- pixel_scale_arcsec(raw_fits_obj)
  psf_scale <- pixel_scale_arcsec(psf_fits_obj)
  scale_x <- raw_scale["x"] / psf_scale["x"]
  scale_y <- raw_scale["y"] / psf_scale["y"]

  template_width <- max(1L, as.integer(round(ncol(raw_white) * scale_x)))
  template_height <- max(1L, as.integer(round(nrow(raw_white) * scale_y)))
  raw_template <- resize_matrix(raw_white, width = template_width, height = template_height)
  match_info <- match_template_ncc(psf_white, raw_template)

  raw_rows <- seq_len(nrow(raw_white))
  raw_cols <- seq_len(ncol(raw_white))
  template_rows <- if (length(raw_rows) == 1L) {
    rep(1, 1L)
  } else {
    1 + (raw_rows - 1) * (template_height - 1) / (length(raw_rows) - 1)
  }
  template_cols <- if (length(raw_cols) == 1L) {
    rep(1, 1L)
  } else {
    1 + (raw_cols - 1) * (template_width - 1) / (length(raw_cols) - 1)
  }

  psf_row_idx <- round(match_info$row + template_rows - 1)
  psf_col_idx <- round(match_info$col + template_cols - 1)
  psf_row_idx <- pmin(pmax(psf_row_idx, 1L), nrow(labels_psf))
  psf_col_idx <- pmin(pmax(psf_col_idx, 1L), ncol(labels_psf))

  labels_raw <- matrix(0L, nrow = nrow(raw_white), ncol = ncol(raw_white))
  for (i in seq_along(psf_row_idx)) {
    labels_raw[i, ] <- labels_psf[psf_row_idx[i], psf_col_idx]
  }

  list(
    labels_raw = labels_raw,
    match_score = match_info$score,
    match_row = match_info$row,
    match_col = match_info$col,
    template_width = template_width,
    template_height = template_height,
    scale_x = scale_x,
    scale_y = scale_y
  )
}

write_region_photometry_csv <- function(fits_obj, cluster_map, path) {
  SED <- RegionPhotometry(
    fits_obj,
    cluster_map,
    error_fallback = "poisson"
  )

  band_cols <- names(SED$flux_wide)[grepl("^[0-9]+$", names(SED$flux_wide))]
  band_cols <- as.character(sort(as.integer(band_cols)))
  err_cols <- paste0(band_cols, "_err")
  neff_cols <- paste0(band_cols, "_n_eff")

  err_cols <- err_cols[err_cols %in% names(SED$flux_wide)]
  neff_cols <- neff_cols[neff_cols %in% names(SED$flux_wide)]
  hdr <- fits_obj$header
  if (is.null(hdr)) {
    hdr <- fits_obj$hdr
  }
  if (is.null(hdr)) {
    stop("FITS object does not contain header cards.")
  }
  keys <- hdr[seq(1L, length(hdr) - 1L, by = 2L)]
  vals <- hdr[seq(2L, length(hdr), by = 2L)]
  filters <- trimws(vals[grepl("^FILTER[0-9]+$", keys)])
  filters <- filters[nzchar(filters)]
  if (length(filters) != length(band_cols)) {
    stop(
      "Header filter count (", length(filters),
      ") does not match photometry band count (", length(band_cols), ")."
    )
  }

  flux_ordered_jy <- SED$flux_wide |>
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

  utils::write.csv(flux_ordered_jy, path, row.names = FALSE)
  invisible(flux_ordered_jy)
}

run_figure_03_segmentation_variant <- function(mode = c("original", "psfmatched", "copula", "contourlet_psfmatched")) {
  mode <- match.arg(mode)
  slug <- variant_slug(mode)

  figure_out_dir <- file.path(base_dir, "results/figures/paper_repro", figure_id, slug)
  png_dir <- file.path(figure_out_dir, "png")
  seg_out_dir <- file.path(base_dir, "results/segmentation/paper_repro", figure_id, slug)
  flux_out_dir <- file.path(base_dir, "results/flux_per_region/paper_repro", figure_id, slug)
  dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(seg_out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(flux_out_dir, recursive = TRUE, showWarnings = FALSE)

  manifest <- vector("list", length(panel_tags))

  for (i in seq_along(panel_tags)) {
    tag <- panel_tags[i]
    cube_path <- resolve_input_cube(tag, mode)
    X <- FITSio::readFITS(cube_path)
    cube <- X$imDat
    raw_cube_path <- resolve_raw_cube(tag)
    X_raw <- FITSio::readFITS(raw_cube_path)

    if (identical(mode, "original") || identical(mode, "psfmatched")) {
      mask_info <- build_starlet_mask(cube)
      seg <- build_segmentation(cube, mask_info$mask, N = 40L)
      variant_note <- "starlet_mask"
    } else if (identical(mode, "copula")) {
      mask_info <- build_copula_mask(cube)
      seg <- build_segmentation(cube, mask_info$mask, N = 40L)
      variant_note <- "copula_gaussian_mask_on_raw_cube"
    } else if (identical(mode, "contourlet_psfmatched")) {
      mask_info <- build_contourlet_mask(cube)
      seg <- build_segmentation(cube, mask_info$mask, N = 40L)
      variant_note <- "contourlet_like_mask"
    } else {
      stop("Unsupported mode: ", mode)
    }

    png_path <- file.path(png_dir, paste0(tag, ".png"))
    fits_path <- file.path(seg_out_dir, paste0(tag, ".fits"))
    csv_path <- file.path(flux_out_dir, paste0("SED_flux_wide_", tag, ".csv"))
    save_segmentation_png(seg, png_path, N = 40L)
    write_cluster_fits(seg, fits_path)

    if (all(dim(seg$cluster_map) == dim(X_raw$imDat)[1:2])) {
      cluster_map_for_photometry <- seg$cluster_map
      match_score <- NA_real_
      match_row <- NA_integer_
      match_col <- NA_integer_
      template_width <- ncol(cluster_map_for_photometry)
      template_height <- nrow(cluster_map_for_photometry)
    } else if (identical(mode, "psfmatched") || identical(mode, "contourlet_psfmatched")) {
      remap <- remap_psf_labels_to_raw(seg$cluster_map, X_raw, X)
      cluster_map_for_photometry <- remap$labels_raw
      match_score <- remap$match_score
      match_row <- remap$match_row
      match_col <- remap$match_col
      template_width <- remap$template_width
      template_height <- remap$template_height
    } else {
      cluster_map_for_photometry <- seg$cluster_map
      match_score <- NA_real_
      match_row <- NA_integer_
      match_col <- NA_integer_
      template_width <- ncol(cluster_map_for_photometry)
      template_height <- nrow(cluster_map_for_photometry)
    }

    write_region_photometry_csv(X_raw, cluster_map_for_photometry, csv_path)

    manifest[[i]] <- data.frame(
      tag = tag,
      mode = mode,
      variant_note = variant_note,
      source_cube = cube_path,
      photometry_cube = raw_cube_path,
      png_path = png_path,
      fits_path = fits_path,
      csv_path = csv_path,
      mask_pixels = sum(mask_info$mask, na.rm = TRUE),
      labeled_pixels = sum(is.finite(seg$cluster_map) & seg$cluster_map > 0),
      registration_score = match_score,
      registration_row = match_row,
      registration_col = match_col,
      registration_template_width = template_width,
      registration_template_height = template_height,
      stringsAsFactors = FALSE
    )
  }

  manifest_df <- do.call(rbind, manifest)
  utils::write.csv(manifest_df, file.path(figure_out_dir, "manifest.csv"), row.names = FALSE)
  message("Finished ", figure_id, " variant: ", mode)
  message("Figures: ", png_dir)
  message("Segmentation FITS: ", seg_out_dir)
  invisible(manifest_df)
}
