# --- minimal deps ---
library(FITSio)
library(dplyr)
library(viridisLite)

# 1) INPUTS ---------------------------------------------------------------
labels <- seg$cluster_map            # matrix [nx, ny] with region ids (1..N; <=0 = sky)
#df     <- flux_wide_out              # data frame with at least: region, n_pix, and SED columns
df_Sed <- read.csv("~/Documents/GitHub/crp8_segmentation/SEDfitting/flux_wide_capivara_reg1_logMlogZ.csv")
df_Sed_wavelet <- read.csv("~/Documents/GitHub/crp8_segmentation/SEDfitting/flux_wide_wavelet_reg1_logMlogZ.csv")


# --- choose which columns to paint (from flux_wide_out) ---
props <- c("logM","logZ","flag")

# --- paint a property back to the label grid (no orientation change) ---
paint_property <- function(labels, df, prop_col, background_na = TRUE) {
  stopifnot(is.matrix(labels), "region" %in% names(df), prop_col %in% names(df))
  lut <- df[[prop_col]]; names(lut) <- as.character(as.integer(df$region))
  lab_vec <- as.vector(labels)
  out_vec <- unname(lut[as.character(lab_vec)])
  if (background_na) out_vec[!is.finite(lab_vec) | lab_vec <= 0] <- NA_real_
  matrix(out_vec, nrow = nrow(labels), ncol = ncol(labels))
}

paint_many <- function(labels, df, props) {
  arr <- array(NA_real_, dim = c(nrow(labels), ncol(labels), length(props)),
               dimnames = list(NULL, NULL, props))
  for (k in seq_along(props)) arr[,,k] <- paint_property(labels, df, props[k])
  arr
}

# build maps aligned with seg$cluster_map
labels <- seg$cluster_map
df     <- df_Sed
maps   <- paint_many(labels, df, props)

# --- ggplot helper that preserves matrix orientation ---
plot_map_gg <- function(mat, quantity = "value", flip_y_axis = FALSE) {
  dfp <- as.data.frame(as.table(mat), responseName = "val")
  names(dfp) <- c("row","col","val")
  dfp$row <- as.integer(dfp$row)
  dfp$col <- as.integer(dfp$col)
  dfp$val <- as.numeric(dfp$val)
  
  p <- ggplot(filter(dfp), aes(x = row, y = col, fill = val)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    scale_fill_viridis_c(na.value = "transparent") +
    labs(x = "x [px]", y = "y [px]", fill = quantity) +
    theme_bw() +
    theme(legend.position = "right")
  
  if (flip_y_axis) p <- p + scale_y_reverse()
  p
}

# example: plot one layer
# legend title = "logZ"
plot_map_gg(maps[, , "logZ"], quantity = "logZ")

# another example with units
plot_map_gg(maps[, , "logM"], quantity = "logM")




labels <- labels_fix
df_wavelet     <- df_Sed_wavelet
maps   <- paint_many(labels, df_wavelet, props)

# --- ggplot helper that preserves matrix orientation ---
plot_map_gg <- function(mat, quantity = "value", flip_y_axis = FALSE) {
  dfp <- as.data.frame(as.table(mat), responseName = "val")
  names(dfp) <- c("row","col","val")
  dfp$row <- as.integer(dfp$row)
  dfp$col <- as.integer(dfp$col)
  dfp$val <- as.numeric(dfp$val)
  
  p <- ggplot(filter(dfp), aes(x = row, y = col, fill = val)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    scale_fill_viridis_c(na.value = "transparent") +
    labs(x = "x [px]", y = "y [px]", fill = quantity) +
    theme_bw() +
    theme(legend.position = "right")
  
  if (flip_y_axis) p <- p + scale_y_reverse()
  p
}

# example: plot one layer
# legend title = "logZ"
plot_map_gg(maps[, , "logZ"], quantity = "logZ")

# another example with units
plot_map_gg(maps[, , "logM"], quantity = "logM")


# extract the slice you want
logM_map <- maps[, , "logM"]
logZ_map <- maps[, , "logZ"]


# write as FITS
writeFITSim(logM_map, file = "../../sed_maps/logM_wavelet.fits")
writeFITSim(logZ_map, file = "../../sed_maps/logZ_wavelet.fits")




plot_map_gg_mask <- function(mat, labels, sed_df, bad_flags = c(-1), title = "", flip_y_axis = FALSE){
  stopifnot(all(dim(mat) == dim(labels)))
  bad_regions <- sed_df$region[sed_df$flag %in% bad_flags]
  
  dfp <- as.data.frame(as.table(mat), responseName = "val")
  names(dfp) <- c("row","col","val")
  dfp$row <- as.integer(dfp$row); dfp$col <- as.integer(dfp$col)
  # add region per pixel and mask only for plotting
  dfp$region <- labels[cbind(dfp$row, dfp$col)]
  dfp$val[dfp$region %in% bad_regions] <- NA_real_
  
  p <- ggplot(dfp, aes(x = row, y = col, fill = val)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    scale_fill_viridis_c(na.value = "transparent", name = title) +
    labs(x = "x [px]", y = "y [px]") +
    theme_bw() + theme(legend.position = "bottom")
  if (flip_y_axis) p <- p + scale_y_reverse()
  p
}

plot_map_gg_mask(maps[, , "logM"], seg$cluster_map, 
                 df_Sed, bad_flags = c(-1), title = "logM")




# extract the slice you want
logM_map <- maps[, , "logM"]
logZ_map <- maps[, , "logZ"]


# write as FITS
writeFITSim(logM_map, file = "../../sed_maps/logM.fits")
writeFITSim(logZ_map, file = "../../sed_maps/logZ.fits")

