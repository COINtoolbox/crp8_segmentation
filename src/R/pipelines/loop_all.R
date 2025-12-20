suppressPackageStartupMessages({
  library(guara)
  library(capivara)
  library(FITSio)
  library(viridis)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

set.seed(42)

# -------------------------------------------------------------------
# Fixed output folders (your requested destinations)
# -------------------------------------------------------------------
out_csv_dir     <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation/results/flux_per_region"
out_seg_dir     <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation/results/figures/sagui_seg"
out_starlet_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation/results/figures/starlet"

# Optional safety: keep or delete these three lines
dir.create(out_csv_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(out_seg_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(out_starlet_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------
palette_van_gogh <- function(n) {
  base <- c("#263C8B", "#547FFF", "#405CFF",
            "#FFFAA3", "#FFDE38", "#BFA524")
  grDevices::colorRampPalette(base)(n)
}

mat_to_df <- function(mat, label) {
  m <- mat
  m[m <= 0] <- NA_real_
  df <- as.data.frame(as.table(m))
  names(df) <- c("y","x","value")
  df$y <- as.integer(df$y)
  df$x <- as.integer(df$x)
  df$panel <- label
  df
}

normalize_panel <- function(v) {
  s <- asinh_stretch(v)
  rng <- range(s, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[2] == rng[1]) {
    return(rep(NA_real_, length(v)))
  }
  (s - rng[1]) / (rng[2] - rng[1])
}

clean_prefix <- function(fits_path) {
  # datacube_sagui1_2_3.fits -> sagui1_2_3
  prefix <- tools::file_path_sans_ext(basename(fits_path))
  prefix <- sub("^datacube_", "", prefix)
  prefix
}

# -------------------------------------------------------------------
# Main per-file processor
# -------------------------------------------------------------------
process_one_cube <- function(fits_path,
                             J = 5,
                             keep_scales_for_mask = 2:5,
                             denoise_k_mask = 2,
                             denoise_k_scales = 1,
                             kclip_collapse = 2,
                             N = 40,
                             background_color = "#d9bea1") {

  prefix <- clean_prefix(fits_path)
  message("Processing: ", prefix)

  # 1) Read cube
  X    <- FITSio::readFITS(fits_path)
  cube <- X$imDat
  stopifnot(length(dim(cube)) == 3L)

  # 2) Collapse (white-light)
  collapsed <- collapse_white_light(cube, kclip = kclip_collapse)

  # 3) Starlet decomposition + mosaic
  dec <- guara::starlet_mask(collapsed, J = J)

  rec_all <- starlet_reconstruct(
    dec,
    keep_scales    = keep_scales_for_mask,
    include_coarse = FALSE,
    denoise_k      = denoise_k_mask,
    mode           = "soft"
  )

  df_list <- list()
  df_list[[1]] <- mat_to_df(collapsed, "Original")

  for (j in seq_len(J)) {
    rec_j <- starlet_reconstruct(
      dec,
      keep_scales    = j,
      include_coarse = FALSE,
      denoise_k      = denoise_k_scales,
      mode           = "soft"
    )
    df_list[[length(df_list) + 1]] <- mat_to_df(rec_j, paste0("Scale j = ", j))
  }

  df_all <- do.call(rbind, df_list) %>%
    group_by(panel) %>%
    mutate(value_norm = normalize_panel(value)) %>%
    ungroup()

  pdf(file.path(out_starlet_dir, paste0(prefix, "_starlet_scales.pdf")),
      height = 4, width = 9)
  ggplot(df_all, aes(x = y, y = x, fill = value_norm)) +
    geom_raster() +
    coord_fixed() +
    scale_fill_gradientn(
      colours  = palette_van_gogh(256),
      na.value = "black",
      limits   = c(0, 1)
    ) +
    facet_wrap(~ panel, ncol = 3) +
    theme_bw() +
    theme(
      strip.background  = element_rect(fill = "gray80", colour = NA),
      panel.grid        = element_blank(),
      panel.background  = element_rect(fill = "black", colour = NA),
      axis.title        = element_blank(),
      axis.text         = element_blank(),
      axis.ticks        = element_blank(),
      legend.position   = "none",
      strip.text        = element_text(size = 10, face = "bold")
    )
  dev.off()

  # 4) Capivara segmentation masked by starlet reconstruction
  mask_rec <- is.finite(rec_all) & (rec_all > 0)
  cube_na  <- guara::mask_cube(cube, mask_rec, mode = "na")
  seg_cap  <- capivara::segment(list(imDat = cube_na), N = N)

  pdf(file.path(out_seg_dir, paste0(prefix, "_segmentation.pdf")),
      height = 6.5, width = 5.5)
  plot_cluster_voronoi_style(
    seg_cap,
    palette = palette_van_gogh(N),
    border_color = "black",
    border_linewidth = 1,
    background_color = background_color
  )
  dev.off()

  # 5) Region photometry -> CSV
  if (!exists("nircam_filters")) {
    stop("Object `nircam_filters` not found. Define it before calling process_one_cube().")
  }

  SED <- RegionPhotometry(
    X,
    seg_cap$cluster_map,
    error_fallback = "poisson"
  )

  # Keep your current convention; change if nbands differs
  num_cols  <- as.character(0:9)
  err_cols  <- paste0(num_cols, "_err")
  neff_cols <- paste0(num_cols, "_n_eff")

  flux_wide_named_jy <- SED$flux_wide %>%
    mutate(across(all_of(num_cols), ~ .x * 1e-8)) %>%   # nJy -> Jy
    mutate(across(all_of(err_cols), ~ .x * 1e-8)) %>%   # nJy -> Jy
    rename_with(~ nircam_filters[match(.x, num_cols)],
                .cols = all_of(num_cols)) %>%
    rename_with(~ paste0(nircam_filters, "_err")[match(.x, err_cols)],
                .cols = all_of(err_cols)) %>%
    rename_with(~ paste0(nircam_filters, "_n_eff")[match(.x, neff_cols)],
                .cols = all_of(neff_cols))

  out_csv <- file.path(out_csv_dir, paste0(prefix, "_SED_flux_wide_Jy.csv"))
  readr::write_csv(flux_wide_named_jy, out_csv)

  invisible(list(
    prefix     = prefix,
    starlet_pdf= file.path(out_starlet_dir, paste0(prefix, "_starlet_scales.pdf")),
    seg_pdf    = file.path(out_seg_dir, paste0(prefix, "_segmentation.pdf")),
    sed_csv    = out_csv
  ))
}

# -------------------------------------------------------------------
# Loop over all FITS
# -------------------------------------------------------------------
data_dir  <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation/data"
fits_list <- sort(list.files(data_dir, pattern = "^datacube_sagui.*\\.fits$", full.names = TRUE))

results <- lapply(fits_list, function(fp) {
  tryCatch(
    process_one_cube(fp, J = 5, N = 40),
    error = function(e) {
      message("FAILED on ", basename(fp), " : ", conditionMessage(e))
      NULL
    }
  )
})

results <- Filter(Negate(is.null), results)
