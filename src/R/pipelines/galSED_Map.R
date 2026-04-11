# ======================================================================
#  map_regions.R
#  Purpose:
#    - Read a segmentation FITS (region id per pixel)
#    - Read per-region SED table (one row per region)
#    - Map any region-level variable to pixel grid (optionally smoothed)
#    - Plot univariate maps (continuous or Fisher-binned gradient)
#    - Plot bivariate maps (two variables combined into a 2D color key)
#
#  Notes:
#    - Expects "region_id == row_number()" unless you provide an explicit
#      region_id column in the CSV.
#    - If you want smoothing, you must have smooth_region_field_laplacian()
#      available (e.g. source your function before running, or define it here).
# ======================================================================

palette_van_gogh_seq_continuous <- function(n) {
  stops <- c(
    "#0B1D39",  # deep night blue
    "#123B73",
    "#2E6FA3",
    "#4FA7A1",
    "#E8E2D2",
    "#F2D36B",
    "#D2A43A"
  )
  
  ramp <- grDevices::colorRamp(stops, space = "Lab")
  
  x <- seq(0, 1, length.out = n)
  rgb <- ramp(x)
  
  grDevices::rgb(rgb[,1], rgb[,2], rgb[,3], maxColorValue = 255)
}


# -----------------------------
#  0) Dependencies
# -----------------------------
pkgs <- c(
  "FITSio", "dplyr", "tidyr", "purrr", "ggplot2",
  "scales", "classInt", "grid"
)
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install)

library(dplyr)
library(ggplot2)

# Optional bivariate helper packages (nice but not mandatory)
# install.packages(c("biscale", "cowplot"))
has_biscale <- requireNamespace("biscale", quietly = TRUE)
has_cowplot <- requireNamespace("cowplot", quietly = TRUE)

bi_class_safe <- function(x, y, dim = 3) {
  # Returns a character vector like "1-3", "2-2", ... with NA where missing.
  stopifnot(length(x) == length(y))
  ok <- is.finite(x) & is.finite(y)
  
  out <- rep(NA_character_, length(x))
  if (!any(ok)) return(out)
  
  xok <- x[ok]
  yok <- y[ok]
  
  # If variable is (nearly) constant, assign everyone to the middle bin.
  x_unique <- length(unique(xok))
  y_unique <- length(unique(yok))
  
  # Reduce dim if there aren't enough unique values
  dx <- min(dim, x_unique)
  dy <- min(dim, y_unique)
  
  if (dx < 2) {
    xb <- rep(ceiling(dim / 2), length(xok))
  } else {
    # Use rank to handle ties; then cut into dx groups
    rx <- rank(xok, ties.method = "average")
    xb <- dplyr::ntile(rx, dx)
    # If dx < dim, stretch bins to 1..dim for consistent legend positions
    xb <- ceiling(xb * dim / dx)
  }
  
  if (dy < 2) {
    yb <- rep(ceiling(dim / 2), length(yok))
  } else {
    ry <- rank(yok, ties.method = "average")
    yb <- dplyr::ntile(ry, dy)
    yb <- ceiling(yb * dim / dy)
  }
  
  out[ok] <- paste0(xb, "-", yb)
  out
}


# -----------------------------
#  1) Configuration
# -----------------------------
cfg <- list(
  fits_path = "/Users/rd23aag/Documents/GitHub/crp8_segmentation/results/segmentation/starlet_capivara/sagui9.fits",
  sed_path  = "/Users/rd23aag/Documents/GitHub/crp8_segmentation/results/sed_fitting/sagui-9_sed_results.csv",
  
  # Region map handling
  zero_is_na = TRUE,
  
  # Which property to plot (univariate)
  property1 = "logzsol",     # example: logzsol
  property2 = "sfr",     # for bivariate (must exist in SED table)
  property_label1 = "Log Z",
  property_label2 = "SFR",
  
  # Region smoothing (requires smooth_region_field_laplacian)
  use_smooth = TRUE,
  smooth_adjacency = "queen",
  smooth_lambda = 500,
  
  # Univariate binning (Fisher breaks -> nonlinear gradient mapping)
  use_fisher_binning = TRUE,
  fisher_n = 30,
  
  # Plot appearance
  flip_y = FALSE,
  transpose = FALSE,
  na_color = "black",
  
  # Bivariate options
  do_bivariate = TRUE,
  biv_n = 6,                 # 3x3 or 4x4 are common
  biv_style = "BlueOr",       # biscale styles: "DkBlue","DkCyan","GrPink",...
  biv_show_legend = TRUE
)

# -----------------------------
#  2) Utility helpers
# -----------------------------
stop_if_missing <- function(df, cols) {
  missing <- setdiff(cols, names(df))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
}

# Basic palette placeholder:
# If you already have palette_van_gogh_div(), keep using it.
# Otherwise this provides a robust default (viridis-like without extra deps).
palette_default <- function(n) {
  grDevices::colorRampPalette(c("#0B1D39", "#2E6FA3", "#E8E2D2", "#F2D36B", "#6B3F1D"))(n)
}

mat_to_df <- function(Z, panel = "value", flip_y = TRUE, transpose = FALSE) {
  if (transpose) Z <- t(Z)
  
  nx <- nrow(Z)
  ny <- ncol(Z)
  df <- expand.grid(x = seq_len(nx), y = seq_len(ny))
  df$value <- as.vector(Z)         # column-major
  df$panel <- panel
  
  if (flip_y) df$y <- ny - df$y + 1
  df
}

region_to_pixels <- function(region_map, region_df) {
  stopifnot(all(c("region_id", "value") %in% names(region_df)))
  
  idmat <- region_map
  idmat[idmat == 0L] <- NA_integer_
  
  v <- region_df$value
  names(v) <- as.character(region_df$region_id)
  
  out <- matrix(NA_real_, nrow(idmat), ncol(idmat))
  ok <- is.finite(idmat)
  out[ok] <- v[as.character(idmat[ok])]
  out
}

read_region_map <- function(fits_path, zero_is_na = TRUE) {
  seg <- FITSio::readFITS(fits_path)
  region_map <- seg$imDat
  if (zero_is_na) region_map[region_map == 0L] <- NA_integer_
  region_map
}

read_sed <- function(csv_path) {
  read.csv(csv_path, stringsAsFactors = FALSE)
}

attach_region_id <- function(sed_df) {
  if (!("region_id" %in% names(sed_df))) {
    sed_df %>% mutate(region_id = row_number())
  } else {
    sed_df %>% mutate(region_id = as.integer(region_id))
  }
}

make_region_field <- function(sed_df, property) {
  stop_if_missing(sed_df, c("region_id", property))
  sed_df %>%
    transmute(region_id = as.integer(region_id),
              value = as.numeric(.data[[property]]))
}

# Optional smoothing wrapper.
# Requires: smooth_region_field_laplacian(region_id_mat, region_values, ...)
map_region_field_to_pixels <- function(region_map, reg_field, cfg) {
  if (isTRUE(cfg$use_smooth)) {
    if (!exists("smooth_region_field_laplacian", mode = "function")) {
      stop(
        "cfg$use_smooth=TRUE but smooth_region_field_laplacian() is not available.\n",
        "Define it or source it before running, or set cfg$use_smooth=FALSE."
      )
    }
    res <- smooth_region_field_laplacian(
      region_id_mat = region_map,
      region_values = reg_field,
      adjacency = cfg$smooth_adjacency,
      lambda = cfg$smooth_lambda
    )
    res$interpolated_matrix
  } else {
    region_to_pixels(region_map, reg_field)
  }
}

# Fisher breaks for nonlinear mapping in scale_fill_gradientn(values=...)
fisher_breaks <- function(x, n = 20) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(c(min(x, na.rm = TRUE), max(x, na.rm = TRUE)))
  classInt::classIntervals(x, n = n, style = "fisher")$brks
}

# -----------------------------
#  3) Plotters (univariate)
# -----------------------------
plot_univariate <- function(df, cfg, palette_fun = palette_van_gogh_seq_continuous) {
  vals <- df$value
  vals_ok <- vals[is.finite(vals)]
  
  if (length(vals_ok) == 0) stop("No finite values to plot.")
  
  if (isTRUE(cfg$use_fisher_binning)) {
    brks <- fisher_breaks(vals_ok, n = cfg$fisher_n)
    ggplot(df, aes(x = x, y = y, fill = value)) +
      geom_tile() +
      coord_fixed() +
      scale_fill_gradientn(
        colours = palette_fun(cfg$fisher_n),
        values  = scales::rescale(brks),
        na.value = cfg$na_color,
        name = cfg$property_label1
      ) +
      theme_void() +
      theme(
        legend.position = "right",
        panel.background = element_rect(fill = cfg$na_color, colour = NA)
      )
  } else {
    ggplot(df, aes(x = x, y = y, fill = value)) +
      geom_raster() +
      coord_fixed() +
      scale_fill_gradientn(
        colours = palette_fun(100),
        na.value = cfg$na_color,
        name = cfg$property_label1
      ) +
      theme_void() +
      theme(
        legend.position = "right",
        panel.background = element_rect(fill = cfg$na_color, colour = NA)
      )
  }
}

# -----------------------------
#  4) Bivariate mapping
#     Option A (recommended): biscale + cowplot (if installed)
#     Option B (fallback): simple manual 2D binning -> custom palette
# -----------------------------
bivariate_fallback <- function(df, xvar, yvar, n = 10,
                               xlab = "X", ylab = "Y",
                               na_color = "black") {
  # Make n-quantile bins on x and y and encode as "i-j"
  x <- df[[xvar]]
  y <- df[[yvar]]
  ok <- is.finite(x) & is.finite(y)
  
  xb <- rep(NA_integer_, nrow(df))
  yb <- rep(NA_integer_, nrow(df))
  xb[ok] <- as.integer(cut(x[ok], breaks = quantile(x[ok], probs = seq(0, 1, length.out = n + 1), na.rm = TRUE),
                           include.lowest = TRUE, labels = FALSE))
  yb[ok] <- as.integer(cut(y[ok], breaks = quantile(y[ok], probs = seq(0, 1, length.out = n + 1), na.rm = TRUE),
                           include.lowest = TRUE, labels = FALSE))
  
  df$bin <- ifelse(ok, paste0(xb, "-", yb), NA_character_)
  
  # Simple n x n palette: dark->light across x, and cool->warm across y
  # (Replace this with your own bivariate palette if you want full control.)
  base_x <- grDevices::colorRampPalette(c("#0B1D39", "#F2D36B"))(n) # x gradient
  base_y <- grDevices::colorRampPalette(c("#2E6FA3", "#6B3F1D"))(n) # y gradient
  
  pal <- setNames(character(n * n), character(n * n))
  k <- 1
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      # blend base_x[i] with base_y[j]
      pal[k] <- grDevices::rgb(
        (col2rgb(base_x[i])[1,1] + col2rgb(base_y[j])[1,1]) / 2,
        (col2rgb(base_x[i])[2,1] + col2rgb(base_y[j])[2,1]) / 2,
        (col2rgb(base_x[i])[3,1] + col2rgb(base_y[j])[3,1]) / 2,
        maxColorValue = 255
      )
      names(pal)[k] <- paste0(i, "-", j)
      k <- k + 1
    }
  }
  
  p <- ggplot(df, aes(x = x, y = y, fill = bin)) +
    geom_raster() +
    coord_fixed() +
    scale_fill_manual(values = pal, na.value = na_color, drop = FALSE) +
    theme_void() +
    theme(
      legend.position = "right",
      panel.background = element_rect(fill = na_color, colour = NA)
    ) +
    labs(fill = paste0(xlab, " × ", ylab))
  
  list(plot = p, palette = pal)
}

bivariate_palette_vangogh <- function(n = 6) {
  # More orthogonal axes (easy to distinguish)
  base_x <- grDevices::colorRampPalette(c("orange", "blue"))(n)  # blue -> pale yellow
  base_y <- grDevices::colorRampPalette(c("white", "black"))(n)  # purple -> green
  
  pal <- setNames(character(n * n), character(n * n))
  k <- 1
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      pal[k] <- grDevices::rgb(
        (col2rgb(base_x[i])[1,1] + col2rgb(base_y[j])[1,1]) / 2,
        (col2rgb(base_x[i])[2,1] + col2rgb(base_y[j])[2,1]) / 2,
        (col2rgb(base_x[i])[3,1] + col2rgb(base_y[j])[3,1]) / 2,
        maxColorValue = 255
      )
      names(pal)[k] <- paste0(i, "-", j)
      k <- k + 1
    }
  }
  pal
}

bivariate_legend_square <- function(pal, n, xlab, ylab, text_size = 20) {
  leg_df <- expand.grid(x = 1:n, y = 1:n)
  leg_df$bi_class <- paste0(leg_df$x, "-", leg_df$y)
  
  ggplot(leg_df, aes(x = x, y = y, fill = bi_class)) +
    geom_tile() +
    scale_fill_manual(values = pal, drop = FALSE) +
    coord_fixed() +
    coord_cartesian(xlim = c(-0.8, n + 0.2), ylim = c(-0.8, n + 0.2), clip = "off") +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0),
      panel.border = element_rect(colour = "white", fill = NA, linewidth = 0.4)
    ) +
    annotate("segment", x = 1, y = -0.35, xend = n, yend = -0.35,
             linewidth = 0.4, colour = "white",
             arrow = grid::arrow(length = unit(0.10, "in"))) +
    annotate("segment", x = -0.35, y = 1, xend = -0.35, yend = n,
             linewidth = 0.4, colour = "white",
             arrow = grid::arrow(length = unit(0.10, "in"))) +
    annotate("text", x = (n + 1) / 2, y = -0.70, label = xlab,
             colour = "white", size = text_size / ggplot2::.pt) +
    annotate("text", x = -0.70, y = (n + 1) / 2, label = ylab,
             colour = "white", size = text_size / ggplot2::.pt, angle = 90)
}



plot_bivariate <- function(df, cfg) {
  stop_if_missing(df, c("x", "y", "value1", "value2"))
  if (!has_cowplot) stop("Package 'cowplot' is required to inset the legend.")
  
  d2c <- df %>%
    mutate(bi_class = bi_class_safe(value1, value2, dim = cfg$biv_n))
  
  pal <- bivariate_palette_vangogh(cfg$biv_n)
  
  p_map <- ggplot(d2c, aes(x = x, y = y, fill = bi_class)) +
    geom_raster() +
    coord_fixed() +
    scale_fill_manual(values = pal, na.value = cfg$na_color, drop = FALSE) +
    theme_void() +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = cfg$na_color, colour = NA),
      plot.background  = element_rect(fill = cfg$na_color, colour = NA)
    )
  
  p_leg <- bivariate_legend_square(
    pal = pal, n = cfg$biv_n,
    xlab = cfg$property_label1,
    ylab = cfg$property_label2,
    text_size = 6
  )
  
  composed <- cowplot::ggdraw() +
    cowplot::draw_plot(p_map, 0, 0, 1, 1) +
    cowplot::draw_plot(p_leg, 0.72, 0.06, 0.24, 0.24)  # smaller, readable
  
  list(map = p_map, legend = p_leg, palette = pal, composed = composed)
}




# -----------------------------
#  5) Main workflow
# -----------------------------
region_map <- read_region_map(cfg$fits_path, zero_is_na = cfg$zero_is_na)

sed <- read_sed(cfg$sed_path) %>%
  attach_region_id()
sed <- sed %>%
  mutate(
    logmass = ifelse(is.finite(mass) & mass > 0, log10(mass), NA_real_)
  )


# --- Univariate: property1 ---
reg1 <- make_region_field(sed, cfg$property1)
Z1 <- map_region_field_to_pixels(region_map, reg1, cfg)
df1 <- mat_to_df(Z1, flip_y = cfg$flip_y, transpose = cfg$transpose)

p1 <- plot_univariate(df1, cfg, palette_fun = palette_van_gogh_seq_continuous)
print(p1)



reg2 <- make_region_field(sed, cfg$property2)
Z2 <- map_region_field_to_pixels(region_map, reg2, cfg)
df2 <- mat_to_df(Z2, flip_y = cfg$flip_y, transpose = cfg$transpose)

p2 <- plot_univariate(df2, cfg, palette_fun = palette_van_gogh_seq_continuous)
print(p2)




# --- Bivariate: property1 vs property2 ---

  reg2 <- make_region_field(sed, cfg$property2)
  
  Z2 <- map_region_field_to_pixels(region_map, reg2, cfg)
  df2 <- mat_to_df(Z2, flip_y = cfg$flip_y, transpose = cfg$transpose)
  
  # Merge into one df with two value columns
  df_biv <- df1 %>%
    rename(value1 = value) %>%
    left_join(df2 %>% rename(value2 = value), by = c("x", "y", "panel"))
  
  biv <- plot_bivariate(df_biv,cfg)
  
  
  print(biv$composed)


# ======================================================================
# End of script
# ======================================================================
