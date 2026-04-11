# ============================================================
#  batch_sagui_sed_maps.R
#  Purpose:
#    For each sagui tag found in data/raw/*.fits:
#      - read segmentation FITS (region map)
#      - read SED CSV (per-region table)
#      - compute derived properties
#      - map region values -> pixels (optionally smoothed)
#      - save per-property FITS + PNG maps
# ============================================================

suppressPackageStartupMessages({
  library(FITSio)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(classInt)
})

# -----------------------------
# 0) Global config
# -----------------------------
base_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation"

# smoothing + plotting options
use_smooth <- TRUE
smooth_lambda <- 10          # try 0.2, 0.5, 1, 2, 5
smooth_adjacency <- "queen"   # or "rook" depending on your function


use_fisher <- FALSE
fisher_n <- 40
na_color <- "black"

# -----------------------------
# 1) Tag discovery
# -----------------------------
get_sagui_tag <- function(path) {
  fname <- basename(path)
  m <- regexpr("sagui[0-9]+(_[0-9]+){0,2}", fname)
  if (m == -1) stop("Could not extract sagui index from: ", fname)
  regmatches(fname, m)
}

fits_dir <- file.path(base_dir, "data/raw")
fits_files <- list.files(fits_dir, pattern = "\\.fits$", full.names = TRUE)
sagui_tags <- sort(unique(vapply(fits_files, get_sagui_tag, character(1))))
print(sagui_tags)

# -----------------------------
# 2) Path helpers (FIXED for your naming)
#    tag: "sagui10" -> SED: "sagui-10_sed_results.csv"
#    tag: "sagui5_6" -> SED: "sagui-5_6_sed_results.csv"
# -----------------------------
sed_path_from_tag <- function(tag, base_dir) {
  sed_dir <- file.path(base_dir, "results/sed_fitting")
  idx <- sub("^sagui", "", tag)
  file.path(sed_dir, paste0("sagui-", idx, "_sed_results.csv"))
}

paths_from_tag <- function(tag, base_dir) {
  list(
    seg_fits = file.path(base_dir, "results/segmentation/starlet_capivara", paste0(tag, ".fits")),
    sed_csv  = sed_path_from_tag(tag, base_dir),
    out_dir  = file.path(base_dir, "results/figures/sed_maps", tag)
  )
}

# -----------------------------
# 3) Palette
# -----------------------------
palette_cosmic_dusk <- function(n = 256) {
  stops <- c(
    "#0A1026",  # near-black indigo
    "#182A6B",  # deep ultramarine
    "#2F5F9A",  # starry blue
    "#4FA7A6",  # calm cyan
    "#D8E1D3",  # luminous haze
    "#F1C76A",  # warm starlight
    "#A96A2A"   # deep amber
  )
  grDevices::colorRampPalette(stops, space = "Lab")(n)
}
palette_wes_astro <- function(n = 256) {
  stops <- c(
    "#2F3E46",  # charcoal blue
    "#52796F",  # desaturated teal
    "#84A98C",  # soft green
    "#CAD2C5",  # cream
    "#E6B89C",  # peach
    "#C97C5D"   # muted terracotta
  )
  grDevices::colorRampPalette(stops, space = "Lab")(n)
}

# -----------------------------
# 4) Core helpers
# -----------------------------
stop_if_missing <- function(df, cols, where = "") {
  miss <- setdiff(cols, names(df))
  if (length(miss)) stop(where, " missing columns: ", paste(miss, collapse = ", "))
}

mat_to_df <- function(Z, flip_y = FALSE, transpose = FALSE) {
  if (transpose) Z <- t(Z)
  nx <- nrow(Z); ny <- ncol(Z)
  df <- expand.grid(x = seq_len(nx), y = seq_len(ny))
  df$value <- as.vector(Z)
  if (flip_y) df$y <- ny - df$y + 1
  df
}

crop_df_to_data <- function(df) {
  df_ok <- df[is.finite(df$value), ]
  if (!nrow(df_ok)) return(df)
  
  x_rng <- range(df_ok$x)
  y_rng <- range(df_ok$y)
  
  df[df$x >= x_rng[1] & df$x <= x_rng[2] &
       df$y >= y_rng[1] & df$y <= y_rng[2], ]
}


region_to_pixels <- function(region_map, region_df) {
  stop_if_missing(region_df, c("region_id", "value"), "region_to_pixels():")
  idmat <- region_map
  idmat[idmat == 0L] <- NA_integer_
  
  v <- region_df$value
  names(v) <- as.character(region_df$region_id)
  
  out <- matrix(NA_real_, nrow(idmat), ncol(idmat))
  ok <- is.finite(idmat)
  out[ok] <- v[as.character(idmat[ok])]
  out
}

map_region_field_to_pixels <- function(region_map, reg_field) {
  if (isTRUE(use_smooth)) {
    res <- smooth_region_field_laplacian(
      region_id_mat = region_map,
      region_values = reg_field,
      adjacency     = smooth_adjacency,
      lambda        = smooth_lambda
    )
    return(res$interpolated_matrix)
  }
  region_to_pixels(region_map, reg_field)
}


make_legend_breaks <- function(vals_ok, n = 4) {
  rng <- range(vals_ok, finite = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) return(rng)
  seq(rng[1], rng[2], length.out = n)
}

plot_univariate <- function(df_pix, legend_title, title_text,
                            show_legend = FALSE,
                            transparent_bg = TRUE,
                            legend_n = 4,
                            legend_position = "bottom") {
  
  vals_ok <- df_pix$value[is.finite(df_pix$value)]
  if (!length(vals_ok)) stop("No finite values for plot: ", title_text)
  
  brks_leg <- make_legend_breaks(vals_ok, n = legend_n)
  
  bg_fill <- if (transparent_bg) "transparent" else na_color
  na_val  <- if (transparent_bg) NA else na_color
  
  ggplot(df_pix, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    scale_fill_gradientn(
      colours  = palette_cosmic_dusk(256),
      breaks   = brks_leg,
      labels   = scales::label_scientific(digits = 2),
      oob      = scales::squish,
      na.value = na_val,
      name     = NULL,
      guide    = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        label.position = "bottom",
        ticks = TRUE,
        barheight = grid::unit(4, "mm"),
        barwidth  = grid::unit(80, "mm")
      )
    ) +
    theme_void() +
    theme(
      legend.position = if (show_legend) legend_position else "none",
      plot.margin = margin(0, 0, 0, 0),
      
      panel.background  = element_rect(fill = bg_fill, colour = NA),
      plot.background   = element_rect(fill = bg_fill, colour = NA),
      legend.background = element_rect(fill = bg_fill, colour = NA),
      legend.key        = element_rect(fill = bg_fill, colour = NA),
      
      legend.box.spacing = grid::unit(1.5, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.text  = element_text(size = 12),
      legend.title = element_text(size = 15)
    )
}




# -----------------------------
# 5) Derived properties
# -----------------------------
safe_log10 <- function(x) ifelse(is.finite(x) & x > 0, log10(x), NA_real_)

mw_age_exp <- function(tage, tau, n_grid = 1000) {
  if (!is.finite(tage) || !is.finite(tau) || tage <= 0 || tau <= 0) return(NA_real_)
  t <- seq(0, tage, length.out = n_grid)
  dt <- t[2] - t[1]
  sfr <- exp(-t / tau)
  num <- sum((tage - t) * sfr) * dt
  den <- sum(sfr) * dt
  if (den == 0) return(NA_real_)
  num / den
}

current_sfr <- function(mass, tage, tau) {
  if (!is.finite(mass) || !is.finite(tage) || !is.finite(tau) || mass <= 0 || tage <= 0 || tau <= 0) return(NA_real_)
  A <- mass / (tau * (1 - exp(-tage / tau)))   # Msun/Gyr
  (A * exp(-tage / tau)) / 1e9                 # Msun/yr
}

# -----------------------------
# 6) Property registry
# -----------------------------
props <- list(
  logZ = list(
    col   = "logzsol",
    label = expression(log[10](Z/Z[odot]))
  ),
  logM = list(
    col   = "logM",
    label = expression(log[10](M/M[odot]))
  ),
  
  # Dust: concise + standard in papers is A_V (if dust2 ~ A_V in your model)
  Av = list(
    col   = "dust2",
    label = expression(A[V]~"[mag]")
  ),
  
  # Age: show mass-weighted as the headline age map
  age_mw = list(
    col   = "t_mw_gyr",
    label = expression(phantom(.)*"<t>"[M]~"[Gyr]")
  ),
  
  # SFH timescale
  tau = list(
    col   = "tau",
    label = expression(tau~"[Gyr]")
  ),
  
  # SFR and "current" SFR (use 0 rather than now)
  sfr = list(
    col   = "sfr",
    label = expression(SFR~"["*M[odot]*"/yr"*"]")
  ),
  sfr0 = list(
    col   = "sfr_now",
    label = expression(SFR[0]~"["*M[odot]*"/yr"*"]")
  ),
  
  # sSFR
  ssfr = list(
    col   = "ssfr",
    label = expression(sSFR~"[yr"^{-1}*"]")
  )
)



# -----------------------------
# 7) Process one sagui tag
# -----------------------------
process_one_sagui <- function(tag) {
  paths <- paths_from_tag(tag, base_dir)
  
  message("\n==============================")
  message("Processing: ", tag)
  message("==============================")
  
  if (!file.exists(paths$seg_fits)) {
    message("Skipping ", tag, " (missing seg FITS): ", paths$seg_fits)
    return(invisible(FALSE))
  }
  if (!file.exists(paths$sed_csv)) {
    message("Skipping ", tag, " (missing SED csv): ", paths$sed_csv)
    return(invisible(FALSE))
  }
  
  dir.create(paths$out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load segmentation
  seg <- FITSio::readFITS(paths$seg_fits)
  region_map <- seg$imDat
  region_map[region_map == 0L] <- NA_integer_
  
  # Load SED
  sed <- read.csv(paths$sed_csv, stringsAsFactors = FALSE)
  
  # Your SED has `region` (not region_id): standardize here
  stop_if_missing(sed, c("region", "mass", "logzsol", "dust2", "tage", "tau", "sfr"), where = paste0(tag, ": SED"))
  sed <- sed %>% mutate(region_id = as.integer(region))
  
  # Derived
  sed <- sed %>%
    mutate(
      logM       = safe_log10(mass),
      logSFR     = safe_log10(sfr),
      t_mw_gyr   = mapply(mw_age_exp, tage, tau),
      sfr_now    = mapply(current_sfr, mass, tage, tau),
      logSFR_now = safe_log10(sfr_now),
      ssfr       = ifelse(is.finite(sfr) & is.finite(mass) & mass > 0, sfr / mass, NA_real_),
      logSSFR    = safe_log10(ssfr)
    )
  
  # Inner loop: properties
  for (nm in names(props)) {
    col <- props[[nm]]$col
    lab <- props[[nm]]$label
    
    if (!col %in% names(sed)) {
      message("  - ", nm, " (missing column: ", col, ")")
      next
    }
    
    reg_field <- sed %>% transmute(region_id, value = as.numeric(.data[[col]]))
    stop_if_missing(reg_field, c("region_id", "value"), where = paste0(tag, " reg_field(", nm, "):"))
    
    Z <- map_region_field_to_pixels(region_map, reg_field)
    
    suffix <- if (use_smooth) "_smooth" else ""
    fits_out <- file.path(paths$out_dir, paste0(nm, suffix, ".fits"))
    png_out  <- file.path(paths$out_dir, paste0(nm, suffix, ".png"))
    
    FITSio::writeFITSim(Z, fits_out)
    
    df_pix <- mat_to_df(Z, flip_y = FALSE, transpose = FALSE)
    p <- plot_univariate(df_pix, legend_title = lab,
                         title_text = paste0(tag, " — ", nm),
                         show_legend = TRUE,
                         transparent_bg = TRUE)
    
    print(p)
    ggsave(png_out, p, width = 6, height = 6, dpi = 300, bg = "transparent")
    
    
    message("  ✓ wrote ", nm)
  }
  
  invisible(TRUE)
}

# -----------------------------
# 8) Outer loop over all saguis
# -----------------------------
for (tag in sagui_tags) {
  process_one_sagui(tag)
}
