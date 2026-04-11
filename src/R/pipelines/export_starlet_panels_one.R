#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(FITSio)
  library(ggplot2)
  library(dplyr)
  library(ragg)
  library(guara)
})

base_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation"
source(file.path(base_dir, "src/R/pipelines/export_starlet_panels_png.R"))

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

mat_to_df <- function(mat, label) {
  m <- mat
  m[m <= 0] <- NA_real_
  df <- as.data.frame(as.table(m))
  names(df) <- c("y", "x", "value")
  df$y <- as.integer(df$y)
  df$x <- as.integer(df$x)
  df$panel <- label
  df
}

tag <- Sys.getenv("SAGUI_TAG", unset = "sagui11")
J <- suppressWarnings(as.integer(Sys.getenv("SAGUI_STARLET_J", unset = "5")))
if (!is.finite(J) || J < 1L) {
  stop("SAGUI_STARLET_J must be a positive integer.")
}
out_tag <- Sys.getenv("SAGUI_STARLET_OUT_TAG", unset = tag)

fits_path <- resolve_raw_cube(tag)
out_dir <- file.path(base_dir, "results/figures/starlet_png", out_tag)

message("Using FITS cube: ", fits_path)

X <- FITSio::readFITS(fits_path)
cube <- X$imDat
collapsed <- collapse_white_light(cube, kclip = 2)
dec <- guara::starlet_mask(collapsed, J = J)

export_starlet_panels_png(
  collapsed = collapsed,
  dec = dec,
  J = J,
  out_dir = out_dir,
  tag = out_tag,
  prefix = "starlet_2d",
  px = 1800,
  res = 300,
  denoise_k = 1,
  mode = "soft",
  keep_include_coarse = FALSE,
  na_transparent = TRUE
)

message("Exported starlet panels for ", tag, " to ", out_dir, " using J=", J)
