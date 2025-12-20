suppressPackageStartupMessages({
  library(guara)
  library(capivara)
  library(FITSio)
  library(viridis)
  library(dbscan)
  library(dplyr)
  library(ggplot2)
})
palette_van_gogh <- function(n) {
  base <- c("#263C8B", "#547FFF","#405CFF",
            "#FFFAA3", "#FFDE38", "#BFA524"
  )
  grDevices::colorRampPalette(base)(n)
}
require(guara)
set.seed(42)
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

## ============================================================
## 1. Load datacube & collapse white light
## ============================================================



fits_path <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation/data/datacube_sagui5_6.fits"

X    <- FITSio::readFITS(fits_path)
cube <- X$imDat    # [nx, ny, nlam]
nx   <- dim(cube)[1]
ny   <- dim(cube)[2]
nlam <- dim(cube)[3]


collapsed <- collapse_white_light(cube, kclip = 2)


image(asinh_stretch(collapsed), col = viridis::turbo(256))


## ============================================================
## 2. Starlet decomposition & reconstruction
## ============================================================
J   <- 5
dec <- guara::starlet_mask(collapsed, J = J)

# Reconstruction: example with scales 2–6, no coarse, soft denoise
rec_all <- starlet_reconstruct(
  dec,
  keep_scales    = 2:5,
  include_coarse = FALSE,
  denoise_k      = 2,
  mode           = "soft"
)
J <- length(dec$w)  # or set explicitly, e.g. J <- 7
df_list <- list()

## Original image
df_list[[1]] <- mat_to_df(collapsed, "Original")

## For each scale j, reconstruct using only that scale
for (j in 1:5) {
  rec_j <- starlet_reconstruct(
    dec,
    keep_scales    = j,        # <- only this scale
    include_coarse = FALSE,
    denoise_k      = 1,
    mode           = "soft"
  )

  df_list[[length(df_list) + 1]] <-
    mat_to_df(rec_j, paste0("Scale j = ", j))
}

df_all <- do.call(rbind, df_list)

normalize_panel <- function(v) {
  s <-  asinh_stretch(v)
  rng <- range(s, na.rm = TRUE)
  (s - rng[1]) / (rng[2] - rng[1])
}

df_all <- df_all %>%
  dplyr::group_by(panel) %>%
  dplyr::mutate(value_norm = normalize_panel(value)) %>%
  dplyr::ungroup()



pdf("/Users/rd23aag/Documents/GitHub/crp8_segmentation/results/figures/starlet/starlet_2d_sagui5_6.pdf",height = 3.5,width = 7.5)
ggplot(df_all, aes(x = y, y = x, fill = value_norm)) +
  geom_raster() +
  coord_fixed() +
  scale_fill_gradientn(
    colours  = palette_van_gogh(256),
    na.value = "black",
    limits   = c(0, 1)    # very important!
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


## ============================================================
## 4. Capivara segmentation examples
## ============================================================
mask_rec <- is.finite(rec_all) & (rec_all > 0)
## (b) Segment original cube but masked by starlet reconstruction
cube_na  <- guara::mask_cube(cube, mask_rec, mode = "na")
N <- 40
seg_cap  <- capivara::segment(list(imDat = cube_na), N = N)
cluster_map_vis <- seg_cap$cluster_map
cluster_map_vis[!is.finite(cluster_map_vis)] <- 0L   # background/masked -> 0

FITSio::writeFITSim(
  cluster_map_vis,
  file   = "/Users/rd23aag/Documents/GitHub/crp8_segmentation/results/segmentation/starlet_capivara/sagui5_6.fits",
  header = seg_cap$header,
  axDat  = seg_cap$axDat
)




pdf("/Users/rd23aag/Documents/GitHub/crp8_segmentation/results/figures/sagui_seg/sagui_seg5_6.pdf",height = 3.5,width = 6)
plot_cluster_voronoi_style(
  seg_cap,
  palette = palette_van_gogh(N),
  border_color = "black",
  border_linewidth = 1,
  background_color = "black"
)
dev.off()



SED <- RegionPhotometry(
  X,
  seg_cap$cluster_map,
  error_fallback = "poisson"
)


num_cols  <- as.character(0:9)        # adjust to nbands if needed
err_cols  <- paste0(num_cols, "_err")
neff_cols <- paste0(num_cols, "_n_eff")
filters <- nircam_filters[1:10]
stopifnot(length(filters) == length(num_cols))

flux_ordered_jy <- SED$flux_wide %>%
  # convert units first
  dplyr::mutate(
    dplyr::across(dplyr::all_of(num_cols),  ~ .x * 1e-8),
    dplyr::across(dplyr::all_of(err_cols),  ~ .x * 1e-8)
  ) %>%
  # rename columns to filter names
  dplyr::rename_with(~ filters, .cols = dplyr::all_of(num_cols)) %>%
  dplyr::rename_with(~ paste0(filters, "_err"),  .cols = dplyr::all_of(err_cols)) %>%
  dplyr::rename_with(~ paste0(filters, "_n_eff"),.cols = dplyr::all_of(neff_cols)) %>%
  # NOW enforce wavelength-ordered column layout
  dplyr::relocate(
    dplyr::any_of(c("region_id", "region", "cluster", "cluster_id")),  # keep whichever exists
    dplyr::any_of("n_pix"),
    dplyr::all_of(filters),
    dplyr::all_of(paste0(filters, "_err")),
    dplyr::all_of(paste0(filters, "_n_eff"))
  )

readr::write_csv(flux_ordered_jy,
                 "/Users/rd23aag/Documents/GitHub/crp8_segmentation/results/flux_per_region/SED_flux_wide_sagui_5_6.csv")














fits_dir   <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation/data"
fits_paths <- list.files(fits_dir, pattern = "\\.fits$", full.names = TRUE)

seg_df_all <- do.call(rbind, lapply(fits_paths, function(fits_path) {

  X    <- FITSio::readFITS(fits_path)
  cube <- X$imDat

  collapsed <- collapse_white_light(cube, kclip = 1)

  dec <- guara::starlet_mask(collapsed, J = 5)
  rec_all <- starlet_reconstruct(dec, keep_scales = 2:5,
                                 include_coarse = FALSE,
                                 denoise_k = 2, mode = "soft")

  mask_rec <- is.finite(rec_all) & (rec_all > 0)

  cube_na <- guara::mask_cube(cube, mask_rec, mode = "na")
  seg_cap <- capivara::segment(list(imDat = cube_na), N = 40)

  lab <- seg_cap$cluster_map
  dataset_id <- tools::file_path_sans_ext(basename(fits_path))

  df <- as.data.frame(as.table(lab))
  names(df) <- c("y","x","cluster")
  df$y <- as.integer(df$y)
  df$x <- as.integer(df$x)
  df$dataset <- dataset_id

  df
}))

N_global <- max(seg_df_all$cluster, na.rm = TRUE)

ggplot(seg_df_all, aes(x = y, y = x, fill = factor(cluster))) +
  geom_raster() +
#  coord_fixed() +
  scale_fill_manual(values = palette_van_gogh(max(8, N_global)),
                    na.value = "black") +
  facet_wrap(~ dataset, ncol = 3,scales="free") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray80", colour = NA),
    panel.grid       = element_blank(),
    panel.background = element_rect(fill = "black", colour = NA),
    axis.title       = element_blank(),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    legend.position  = "none",
    strip.text       = element_text(size = 10, face = "bold")
  )













sed <- read.csv("sagui-5-6_sed_results.csv")

# SFR(t) for exponential declining SFH
sfr_exp <- function(t, tau) {
  exp(-t / tau)
}
mw_age_exp <- function(tage, tau, n_grid = 1000) {
  # time grid from 0 (start of SF) to tage (now)
  t <- seq(0, tage, length.out = n_grid)
  dt <- t[2] - t[1]

  sfr <- sfr_exp(t, tau)

  # star *formation* time = t
  # stellar age now = tage - t
  num <- sum((tage - t) * sfr) * dt       # ∫ (tage - t) SFR(t) dt
  den <- sum(sfr) * dt                    # ∫ SFR(t) dt

  num / den  # returns in same units as tage (Gyr)
}

current_sfr <- function(mass, tage, tau) {
  # tau and tage in Gyr
  # mass in Msun (living mass or total mass formed)

  # Normalization constant
  A <- mass / (tau * (1 - exp(-tage / tau)))

  # SFR now:
  sfr_now_gyr = A * exp(-tage / tau)  # Msun / Gyr

  # convert to Msun / yr
  sfr_now_yr = sfr_now_gyr / 1e9

  return(sfr_now_yr)
}




sed2 <- sed %>%
  mutate(
    t_mw_gyr = purrr::map2_dbl(tage, tau, mw_age_exp),
    sfr_now = current_sfr(mass, tage, tau)
  )

# columns: mass, logzsol, dust2, tage, tau
# row i == region i


# Example: mass per region
Age <- sed2$t_mw_gyr
names(Age) <- as.character(seq_len(nrow(pars)))  # region 1..N

p_mass <- plot_property_map(
  cluster_map = seg_cap$cluster_map,
  values      = Age,
  palette     = "magma",
  value_label = "Age"
)

p_mass

sfr_log10 <- (sed2$sfr_now)
#sfr_log10[!is.finite(sfr_log10)] <- NA          # just in case

names(sfr_log10) <- as.character(seq_len(nrow(sed2)))

p_sfr <- plot_property_map(
  cluster_map = seg_cap$cluster_map,
  values      = asinh_stretch(sfr_log10),
  palette     = "plasma",
  na_color = "black",
  value_label = expression("SFR [" * M["\u2609"] * "/yr]")
)

p_sfr

writeFITSim(seg_cap$cluster_map,"sagui.fits")


## ------------------------------------------------------------
## NIRCam filter definitions & helpers
## ------------------------------------------------------------

# Central wavelengths (µm)
nircam_lambda_um <- c(
  F090W = 0.902, F115W = 1.154, F150W = 1.501, F182M = 1.842,
  F200W = 1.989, F210M = 2.099, F277W = 2.770, F335M = 3.365,
  F356W = 3.563, F410M = 4.082, F430M = 4.286, F444W = 4.421,
  F460M = 4.620, F480M = 4.828
)

nircam_filters <- names(nircam_lambda_um)

build_filter_spec <- function(
    filters,
    lambda,
    throughput_input = NULL,
    exclude = NULL,
    assume_ordered = TRUE,
    band_index = NULL
){
  if (is.null(exclude)) exclude <- character(0)
  filters_use <- setdiff(filters, exclude)

  if (is.null(band_index)) {
    band_index <- seq_along(filters_use) - 1L   # 0-based
  }

  tibble::tibble(
    filter     = filters_use,
    lambda_um  = as.numeric(lambda[filters_use]),
    band_index = band_index
  )
}

band_lookup <- tibble::tibble(
  band      = as.character(seq_along(nircam_filters) - 1L),  # "0","1",...
  band_idx  = seq_along(nircam_filters) - 1L,
  filter    = nircam_filters,
  lambda_um = nircam_lambda_um[nircam_filters]
)

spec_nircam <- build_filter_spec(
  filters          = nircam_filters,
  lambda           = nircam_lambda_um,
  throughput_input = "/Users/rd23aag/Downloads/nircam_throughputs/mean_throughputs",
  exclude          = c("F150W2"),     # unused here but kept for API symmetry
  assume_ordered   = TRUE,
  band_index       = seq_along(nircam_filters) - 1L
)

make_throughput_table <- function(
    input,                                    # dir path OR character vector of files/URLs
    filters,                                  # e.g., c("F090W","F115W",...)
    nircam_lambda_um,                         # named numeric: central wavelength by filter
    exclude = c("F150W2"),            # odd filenames to ignore
    output = "throughput_nircam.csv",         # output file
    format = c("csv","rds"),                  # write CSV or RDS
    preview = FALSE, preview_path = "throughput_preview.png"
){
  stopifnot(requireNamespace("dplyr", quietly=TRUE))
  stopifnot(requireNamespace("purrr", quietly=TRUE))
  stopifnot(requireNamespace("stringr", quietly=TRUE))
  stopifnot(requireNamespace("ggplot2", quietly=TRUE))
  stopifnot(requireNamespace("readr", quietly=TRUE))

  library(dplyr); library(purrr); library(stringr); library(ggplot2); library(readr)

  format <- match.arg(format)

  # Collect files
  files <- if (length(input) == 1 && dir.exists(input)) {
    list.files(input, pattern = "_mean_system_throughput\\.txt$", full.names = TRUE)
  } else {
    as.character(input)
  }
  files <- files[str_detect(basename(files), "_mean_system_throughput\\.txt$")]
  if (!length(files)) stop("No throughput files found.")

  # Keep files that START with desired filters & not excluded
  keep_pat    <- paste0("^(", paste(filters, collapse="|"), ")_")
  exclude_pat <- paste0("^(", paste(exclude,  collapse="|"), ")_")
  base <- basename(files)
  files <- files[str_detect(base, keep_pat) & !str_detect(base, exclude_pat)]
  if (!length(files)) stop("No files left after filtering.")

  # Reader
  read_one <- function(f){
    df <- suppressMessages(readr::read_table(f, show_col_types = FALSE))
    if (ncol(df) < 2) stop(paste("Bad file:", f))
    nm <- stringr::str_extract(basename(f), "^F[0-9]{3,4}[A-Z]+")
    if (is.na(nm)) nm <- stringr::str_extract(basename(f), "^[^_]+")
    tibble::tibble(
      wavelength = df[[1]],
      throughput = df[[2]],
      filter     = nm,
      lambda_c   = unname(nircam_lambda_um[nm])
    )
  }

  throughput_all <- purrr::map_dfr(files, read_one) %>%
    filter(!filter %in% exclude, filter %in% filters) %>%
    arrange(filter, wavelength) %>%
    mutate(filter = factor(filter, levels = names(sort(nircam_lambda_um))))

  # Write
  if (format == "csv") readr::write_csv(throughput_all, output) else saveRDS(throughput_all, output)

  if (preview) {
    p <- ggplot(throughput_all, aes(wavelength, throughput, color = filter)) +
      geom_line(linewidth=0.6) + theme_minimal(base_size=12) +
      labs(x="Wavelength (μm)", y="System Throughput", color="Filter") +
      scale_color_viridis_d(option="turbo")
    ggsave(preview_path, p, width=8, height=4.6, dpi=150)
  }

  message(sprintf("Throughput table saved to %s  (%d rows, %d filters).",
                  output, nrow(throughput_all), dplyr::n_distinct(throughput_all$filter)))
  invisible(throughput_all)
}

throughput_all <- make_throughput_table(
  input            = "/Users/rd23aag/Downloads/nircam_throughputs/mean_throughputs/",
  filters          = nircam_filters,
  nircam_lambda_um = nircam_lambda_um,
  exclude          = c("F150W2", "150W2"),
  output           = "throughput_nircam.csv",
  format           = "csv",
  preview          = TRUE
)

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

throughput_all <- readr::read_csv("throughput_nircam.csv", show_col_types = FALSE)

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






















## ============================================================
## 7. Region photometry + NIRCam filter spec + SED plot
## ============================================================


build_filter_spec <- function(
    filters,
    lambda,
    throughput_input = NULL,
    exclude = NULL,
    assume_ordered = TRUE,
    band_index = NULL
){
  filters_use <- setdiff(filters, exclude %||% character(0))

  if (is.null(band_index)) {
    band_index <- seq_along(filters_use) - 1L
  }

  tibble::tibble(
    filter     = filters_use,
    lambda_um  = as.numeric(lambda[filters_use]),
    band_index = band_index
  )
}


# NIRCam filter info
filters <- c(
  "F090W","F115W","F150W","F182M","F200W",
  "F210M","F277W","F335M","F356W","F410M",
  "F430M","F444W","F460M","F480M"
)

nircam_lambda_um <- c(
  F090W = 0.902, F115W = 1.154, F150W = 1.501, F182M = 1.842,
  F200W = 1.989, F210M = 2.099, F277W = 2.770, F335M = 3.365,
  F356W = 3.563, F410M = 4.082, F430M = 4.286, F444W = 4.421,
  F460M = 4.620, F480M = 4.828
)

band_lookup <- tibble::tibble(
  band      = as.character(0:13),     # matches your flux_long$band
  band_idx  = 0:13,
  filter    = filters,
  lambda_um = nircam_lambda_um[filters]
)
band_lookup

spec_nircam <- build_filter_spec(
  filters          = filters,
  lambda           = nircam_lambda_um,
  throughput_input = "/Users/rd23aag/Downloads/nircam_throughputs/mean_throughputs",
  exclude          = c("F150W2"),
  assume_ordered   = TRUE,
  band_index       = 0:(length(filters) - 1)
)


throughput_all <- make_throughput_table(
  input  = "/Users/rd23aag/Downloads/nircam_throughputs/mean_throughputs/",
  filters = filters,
  nircam_lambda_um = nircam_lambda_um,
  exclude = c("F150W2","150W2"),
  output  = "throughput_nircam.csv",
  format  = "csv",
  preview = TRUE
)

plot_sed_with_filters <- function(
    reg_long,                   # tibble with: region, band, flux, flux_err
    region_id,                  # integer region ID
    throughput_tbl,             # output of make_throughput_table()
    filters,                    # vector of JWST filter names in correct order
    nircam_lambda_um,           # named numeric vector: F090W = 0.902, ...
    scale_mode = c("normalized","auto"),
    band_height_frac = 0.15,
    errors_are_fractional = FALSE,
    show_legend = TRUE
){
  scale_mode <- match.arg(scale_mode)

  ## ----------- SED for a single region -----------
  sed <- reg_long |>
    dplyr::filter(.data$region == region_id) |>
    dplyr::mutate(
      band_idx  = as.integer(as.character(.data$band)),
      filt      = filters[band_idx + 1L],          # map 0→F090W, 1→F115W,...
      lambda_um = nircam_lambda_um[filt]
    ) |>
    dplyr::filter(!is.na(.data$lambda_um)) |>
    dplyr::arrange(.data$lambda_um)

  stopifnot(nrow(sed) > 0)

  ## ----------- Throughput -----------
  if (!"filter" %in% names(throughput_tbl)) {
    stop("throughput_tbl must have a column named 'filter'. Got: ",
         paste(names(throughput_tbl), collapse = ", "))
  }

  thr <- throughput_tbl |>
    dplyr::mutate(filter = as.character(.data$filter)) |>
    dplyr::filter(.data$filter %in% filters)

  ## ----------- Scaling -----------
  if (scale_mode == "normalized") {

    K <- max(sed$flux, na.rm = TRUE)

    sed <- sed |>
      dplyr::mutate(
        y    = .data$flux / K,
        ymin = if (errors_are_fractional) {
          (.data$flux * (1 - .data$flux_err)) / K
        } else {
          (.data$flux - .data$flux_err) / K
        },
        ymax = if (errors_are_fractional) {
          (.data$flux * (1 + .data$flux_err)) / K
        } else {
          (.data$flux + .data$flux_err) / K
        }
      )

    sed$ymin <- pmax(sed$ymin, 0)

    thr <- thr |>
      dplyr::mutate(y = .data$throughput * band_height_frac)

    ylab    <- "Normalized Flux"
    ylimits <- c(-0.02, 1.05)

  } else {

    sed_scale <- stats::quantile(sed$flux, 0.95, na.rm = TRUE)

    thr <- thr |>
      dplyr::mutate(y = .data$throughput * sed_scale * band_height_frac)

    sed <- sed |>
      dplyr::mutate(
        y    = .data$flux,
        ymin = if (errors_are_fractional) {
          pmax(.data$flux * (1 - .data$flux_err), 0)
        } else {
          pmax(.data$flux - .data$flux_err, 0)
        },
        ymax = if (errors_are_fractional) {
          .data$flux * (1 + .data$flux_err)
        } else {
          .data$flux + .data$flux_err
        }
      )

    ylab    <- "Flux"
    ylimits <- NULL
  }

  ## ----------- Plot -----------
  p <- ggplot2::ggplot() +
    ggplot2::geom_area(
      data  = thr,
      ggplot2::aes(x = wavelength, y = y, fill = filter, group = filter),
      position = "identity",
      alpha    = 0.35,
      color    = NA,
      size     = 0
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
    ggplot2::scale_color_viridis_d(option = "turbo", name = "Filter") +
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


# Region photometry based on SED labels (from hclust)
SED <- RegionPhotometry(
  X,
  seg_cap$cluster_map,
  error_fallback = "poisson"
)


num_cols <- as.character(0:9)   # the numeric-band columns

flux_wide_named <- SED$flux_wide %>%
  rename_with(
    ~ filters[match(.x, num_cols)],   # "0" -> "F090W", etc.
    .cols = all_of(num_cols)
  )

library(dplyr)

# bands numéricas


# erros e n_eff
err_cols  <- paste0(num_cols, "_err")
neff_cols <- paste0(num_cols, "_n_eff")

flux_wide_named_jy <- SED$flux_wide %>%

  # 1) converter fluxos das bandas para Jy
  mutate(across(all_of(num_cols), ~ .x * 1e-8)) %>%   # 10 nJy -> Jy

  # 2) converter erros correspondentes para Jy (se existirem)
  mutate(across(all_of(err_cols), ~ .x * 1e-8)) %>%

  # 3) renomear colunas das bandas: 0->F090W, 1->F115W, etc
  rename_with(
    ~ filters[match(.x, num_cols)],
    .cols = all_of(num_cols)
  ) %>%

  # 4) renomear colunas de erro: 0_err -> F090W_err
  rename_with(
    ~ paste0(filters, "_err")[match(.x, err_cols)],
    .cols = all_of(err_cols)
  ) %>%

  # 5) renomear colunas n_eff: 0_n_eff -> F090W_n_eff
  rename_with(
    ~ paste0(filters, "_n_eff")[match(.x, neff_cols)],
    .cols = all_of(neff_cols)
  )


readr::write_csv(flux_wide_named_jy, "SED_flux_wide_Jy_reg8.csv")




p1 <- plot_sed_with_filters(
  reg_long         = SED$flux_long,
  region_id        = 10,
  throughput_tbl   = throughput_all,
  filters          = nircam_filters,
  nircam_lambda_um = nircam_lambda_um,
  scale_mode       = "normalized",
  band_height_frac = 0.95,
  errors_are_fractional = FALSE,
  show_legend      = FALSE
)

p1



## ============================================================
## 5. HDBSCAN on bright pixels of rec_all (cluster arms/knots)
## ============================================================

df <- which(is.finite(rec_all), arr.ind = TRUE)
df <- data.frame(
  y     = df[, 1],
  x     = df[, 2],
  value = rec_all[df]
)

# Keep only bright pixels (top 2%)
thr <- quantile(df$value, 0.98, na.rm = TRUE)
pts <- df[df$value > thr, ]

# HDBSCAN on (x, y)
pts_scaled    <- scale(pts[, c("x", "y")])
hdb           <- hdbscan(pts_scaled, minPts = 5)
pts$cluster   <- hdb$cluster

# reconstruct cluster-labelled matrix
cluster_mat <- matrix(0L, nrow = nx, ncol = ny)
sel <- pts$cluster > 0
if (any(sel)) {
  cluster_mat[cbind(pts$y[sel], pts$x[sel])] <- pts$cluster[sel]
}

df_lab <- as.data.frame(as.table(cluster_mat))
names(df_lab) <- c("y", "x", "cluster")
df_lab$y       <- as.integer(df_lab$y)
df_lab$x       <- as.integer(df_lab$x)
df_lab$cluster <- as.integer(df_lab$cluster)
df_lab$cluster_f <- factor(ifelse(df_lab$cluster == 0, NA, df_lab$cluster))

ggplot(df_lab, aes(x = x, y = y, fill = cluster_f)) +
  geom_raster(na.rm = FALSE) +
  scale_y_reverse() +
  coord_equal() +
  scale_fill_viridis_d(na.value = "gray") +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle("HDBSCAN clusters on bright pixels")


## ============================================================
## 6. Hierarchical clustering in SED space (guara_hclust_sample_assign)
## ============================================================

res_hc <- guara_hclust_sample_assign(
  cube,
  mask    = mask_rec,
  N       = 10,
  m_limit = 8000,
  seed    = 1
)

lab_img <- res_hc$labels_img
lab_img[lab_img == 0] <- NA

par(mfrow = c(1, 2), mar = c(2, 2, 2, 1))
image(asinh_stretch(collapsed), col = pal, main = "Collapsed")
image(log(lab_img), col = viridis::plasma(10), main = "hclust SED labels")





## ============================================================
## 8. Alternative SED clustering via guara_sed_cluster (k-means / hclust)
## ============================================================

res_sed <- guara_sed_cluster(
  cube,
  mask       = mask_rec,  # hard reconstruction / mask
  N          = 50,        # number of clusters
  method     = "kmeans",  # or "hclust" for small datasets
  max_hclust = 100000,
  seed       = 1
)

lab_img_plot <- res_sed$labels_img
lab_img_plot[is.na(lab_img_plot) | lab_img_plot == 0] <- NA

par(mfrow = c(1, 1), mar = c(2, 2, 2, 1))
image(lab_img_plot,
      col      = viridis::magma(100),
      useRaster = TRUE,
      asp      = 1,
      main     = "SED-based clustering (guara_sed_cluster)")


## ============================================================
## 9. Example: other cube (reg1) & quick Capivara segmentation
## ============================================================

Xm <- FITSio::readFITS("/Users/rd23aag/Documents/GitHub/crp8_segmentation/data/raw/datacube_reg1.fits")
manga_seg <- capivara::segment(Xm, N = 12)

par(mfrow = c(1, 1), mar = c(2, 2, 2, 1))
image(manga_seg$cluster_map,
      col = viridis::plasma(12),
      asp = 1,
      main = "Capivara segmentation (reg1 cube)")
