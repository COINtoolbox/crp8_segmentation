suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(FITSio)
  library(sagui)
})

# ============================================================
# batch_surface_sfr_maps.R
# Purpose:
#   - For each sagui tag with segmentation + NEW SED results,
#   - use the surface-SFR quantities already present in the
#     table (`sigma_sfr_10myr`, `sigma_sfr_100myr`),
#   - map them back to the segmentation grid (optionally smoothed),
#   - write a clean, dedicated figure set without touching the
#     existing SFR outputs.
# ============================================================

base_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation"

use_smooth <- TRUE
smooth_lambda <- 10
smooth_adjacency <- "queen"
na_color <- "black"
crop_padding_px <- 4L
contrast_quantiles <- c(0.02, 0.98)
surface_sfr_column <- "sigma_sfr_100myr"
property_contrast_quantiles <- list(
  logZ = c(0.10, 0.90),
  Av = c(0.02, 0.995),
  age_mw = c(0.05, 0.95),
  sfr = c(0.03, 0.97)
)

get_sagui_tag <- function(path) {
  fname <- basename(path)
  m <- regexpr("sagui[0-9]+(_[0-9]+){0,2}", fname)
  if (m == -1) stop("Could not extract sagui index from: ", fname)
  regmatches(fname, m)
}

sed_path_from_tag <- function(tag, base_dir) {
  idx <- sub("^sagui", "", tag)
  idx_hyphen <- gsub("_", "-", idx)
  candidates <- unique(c(
    file.path(base_dir, "results/sed_fitting", paste0("sagui-", idx, "_sed_results_new.csv")),
    file.path(base_dir, "results/sed_fitting", paste0("sagui-", idx_hyphen, "_sed_results_new.csv"))
  ))
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) {
    return(existing[[1]])
  }
  candidates[[1]]
}

paths_from_tag <- function(tag, base_dir) {
  list(
    seg_fits = file.path(base_dir, "results/segmentation/starlet_capivara", paste0(tag, ".fits")),
    sed_csv  = sed_path_from_tag(tag, base_dir),
    out_dir  = file.path(base_dir, "results/figures/sed_maps_surface_sfr", tag)
  )
}

fits_dir <- file.path(base_dir, "data/raw")
fits_files <- list.files(fits_dir, pattern = "\\.fits$", full.names = TRUE)
sagui_tags <- sort(unique(vapply(fits_files, get_sagui_tag, character(1))))

stop_if_missing <- function(df, cols, where = "") {
  miss <- setdiff(cols, names(df))
  if (length(miss)) stop(where, " missing columns: ", paste(miss, collapse = ", "))
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

make_legend_breaks <- function(vals_ok, n = 4) {
  rng <- range(vals_ok, finite = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rng)
  seq(rng[1], rng[2], length.out = n)
}

palette_cosmic_dusk <- function(n = 256) {
  stops <- c(
    "#061024",
    "#132A67",
    "#2858A9",
    "#2E9BB8",
    "#D8EEE6",
    "#F6D26C",
    "#B86D1E"
  )
  grDevices::colorRampPalette(stops, space = "Lab")(n)
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
  if (is.null(probs)) {
    return(contrast_quantiles)
  }

  vals_ok <- values[is.finite(values)]
  n_unique <- length(unique(signif(vals_ok, 4)))

  if (identical(property_name, "logZ") && n_unique <= 15) {
    return(c(0, 1))
  }

  if (identical(property_name, "sfr") && n_unique <= 15) {
    return(c(0.02, 0.99))
  }

  probs
}

plot_univariate <- function(df_pix,
                            property_name = NULL,
                            show_legend = TRUE,
                            transparent_bg = TRUE,
                            legend_n = 4,
                            legend_position = "bottom") {
  vals_ok <- df_pix$value[is.finite(df_pix$value)]
  if (!length(vals_ok)) stop("No finite values available for plotting.")

  probs <- resolve_contrast_quantiles(property_name, vals_ok)
  lims <- as.numeric(stats::quantile(vals_ok, probs = probs, na.rm = TRUE))
  if (!all(is.finite(lims)) || diff(lims) == 0) {
    lims <- range(vals_ok, finite = TRUE)
  }
  brks_leg <- seq(lims[1], lims[2], length.out = legend_n)
  bg_fill <- if (transparent_bg) "transparent" else na_color
  na_val  <- if (transparent_bg) NA else na_color

  ggplot(df_pix, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    scale_fill_gradientn(
      colours = palette_cosmic_dusk(256),
      limits = lims,
      breaks = brks_leg,
      labels = scales::label_scientific(digits = 2),
      oob = scales::squish,
      na.value = na_val,
      name = NULL,
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        label.position = "bottom",
        ticks = TRUE,
        barheight = grid::unit(5, "mm"),
        barwidth = grid::unit(96, "mm")
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
      legend.text = element_text(size = 17)
    )
}

region_to_pixels <- function(region_map, region_df) {
  stop_if_missing(region_df, c("region_id", "value"), "region_to_pixels():")
  idmat <- region_map
  idmat[idmat == 0L] <- NA_integer_

  vals <- region_df$value
  names(vals) <- as.character(region_df$region_id)

  out <- matrix(NA_real_, nrow(idmat), ncol(idmat))
  ok <- is.finite(idmat)
  out[ok] <- vals[as.character(idmat[ok])]
  out
}

map_region_field_to_pixels <- function(region_map, reg_field) {
  if (isTRUE(use_smooth)) {
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

should_smooth_property <- function(property_name, values) {
  if (!isTRUE(use_smooth)) {
    return(FALSE)
  }

  vals_ok <- values[is.finite(values)]
  n_unique <- length(unique(signif(vals_ok, 4)))

  if (identical(property_name, "logZ") && n_unique <= 15) {
    return(FALSE)
  }

  if (identical(property_name, "sfr") && n_unique <= 15) {
    return(FALSE)
  }

  TRUE
}

process_one_sagui <- function(tag) {
  paths <- paths_from_tag(tag, base_dir)

  message("\n==============================")
  message("Processing surface SFR: ", tag)
  message("==============================")

  needed <- c(paths$seg_fits, paths$sed_csv)
  missing <- needed[!file.exists(needed)]
  if (length(missing)) {
    message("Skipping ", tag, " (missing inputs):")
    for (path in missing) {
      message("  - ", path)
    }
    return(invisible(FALSE))
  }

  dir.create(paths$out_dir, recursive = TRUE, showWarnings = FALSE)

  seg <- FITSio::readFITS(paths$seg_fits)
  region_map <- seg$imDat
  region_map[region_map == 0L] <- NA_integer_

  sed <- read.csv(paths$sed_csv, stringsAsFactors = FALSE)
  stop_if_missing(
    sed,
    c("region", "logzsol", "dust2", "age_mass_weighted", "sigma_sfr_10myr", "sigma_sfr_100myr", "area_kpc2"),
    where = paste0(tag, ": SED")
  )

  sed <- sed %>% mutate(region_id = as.integer(region))

  props <- list(
    logZ = "logzsol",
    Av = "dust2",
    age_mw = "age_mass_weighted",
    sfr = surface_sfr_column
  )

  suffix <- if (use_smooth) "_smooth" else ""
  for (nm in names(props)) {
    col <- props[[nm]]
    reg_field <- sed %>%
      transmute(region_id, value = as.numeric(.data[[col]]))
    do_smooth <- should_smooth_property(nm, reg_field$value)
    Z <- if (do_smooth) {
      map_region_field_to_pixels(region_map, reg_field)
    } else {
      region_to_pixels(region_map, reg_field)
    }

    fits_out <- file.path(paths$out_dir, paste0(nm, suffix, ".fits"))
    png_out <- file.path(paths$out_dir, paste0(nm, suffix, ".png"))

    FITSio::writeFITSim(Z, fits_out)
    df_pix <- mat_to_df(Z, flip_y = FALSE, transpose = FALSE)
    df_pix <- crop_df_to_data(df_pix, pad = crop_padding_px)
    p <- plot_univariate(df_pix, property_name = nm, show_legend = TRUE, transparent_bg = TRUE)
    ggsave(png_out, p, width = 7, height = 7, dpi = 320, bg = "transparent")
  }

  csv_out <- file.path(paths$out_dir, "surface_sfr_regions.csv")
  write.csv(
    sed %>%
      mutate(selected_surface_sfr = .data[[surface_sfr_column]]) %>%
      select(region_id, n_pix, area_kpc2, sigma_sfr_10myr, sigma_sfr_100myr, selected_surface_sfr),
    csv_out,
    row.names = FALSE
  )

  message("  ✓ wrote new-table surface-SFR products for ", tag, " using ", surface_sfr_column)
  invisible(TRUE)
}

for (tag in sagui_tags) {
  process_one_sagui(tag)
}
