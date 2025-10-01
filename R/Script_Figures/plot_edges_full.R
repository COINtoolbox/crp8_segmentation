# --- Packages ---
library(FITSio)
library(ggplot2)
library(ggtext)     # for element_markdown()
library(patchwork)  # for | and / operators
library(ragg)

# --- Inputs (edit paths as needed) ---
files <- sprintf("..//..//data/raw/datacube_reg%d.fits", c(1,2,4,5,6,7,8))

# --- Small helpers ---
read_cube <- function(path) {
  x <- FITSio::readFITS(path)
  x$imDat
}

# unwrap: make_rgb sometimes returns a list
unwrap_rgb <- function(x) {
  if (is.list(x)) {
    if (!is.null(x$rgb)) return(x$rgb)
    return(x[[1]])
  }
  x
}

# Convert RGB (array or cimg) into dataframe for ggplot
rgb_to_df <- function(cube_rgb, rotate_ccw = TRUE) {
  arr <- unwrap_rgb(cube_rgb)
  arr <- as.array(arr)  # works for array or cimg
  if (length(dim(arr)) == 4) {
    arr <- arr[,,,1]    # drop z if present
  }
  if (length(dim(arr)) == 2) {
    arr <- array(rep(arr, 3), dim = c(dim(arr), 3))
  }
  H <- dim(arr)[1]; W <- dim(arr)[2]
  hex <- grDevices::rgb(arr[,,1], arr[,,2], arr[,,3])
  M <- array(hex, dim = c(H, W))
  if (rotate_ccw) M <- t(apply(M, 2, rev))
  df <- expand.grid(y = seq_len(nrow(M)), x = seq_len(ncol(M)))
  df$hex    <- as.character(M[cbind(df$y, df$x)])
  df$x_plot <- df$x; df$y_plot <- df$y
  df
}

# Plain RGB (no masking). NA pixels -> black
plot_rgb_plain <- function(cube_rgb, rotate_ccw = TRUE) {
  df <- rgb_to_df(cube_rgb, rotate_ccw = rotate_ccw)
  df$hex[is.na(df$hex)] <- "#000000"
  ggplot(df, aes(x = -x_plot, y = y_plot, fill = hex)) +
    geom_raster() +
    scale_fill_identity() +
    coord_equal() + theme_void()
}

# Overlay contour from label matrix
# Overlay contour from label matrix directly on top of the RGB image
add_contour <- function(cube_rgb, L, rotate_ccw = TRUE,
                        col = "#ffcc00", size = 0.7, alpha = 0.95) {
  # RGB base
  df_rgb <- rgb_to_df(cube_rgb, rotate_ccw = rotate_ccw)
  df_rgb$hex[is.na(df_rgb$hex)] <- "#000000"
  
  # Contour data
  M <- if (rotate_ccw) t(apply(L, 2, rev)) else L
  dfc <- expand.grid(y = seq_len(nrow(M)), x = seq_len(ncol(M)))
  v <- as.numeric(M[cbind(dfc$y, dfc$x)] > 0L)
  v[is.na(v)] <- 0
  dfc$val <- v
  
  ggplot() +
    geom_raster(data = df_rgb, aes(x = -x_plot, y = y_plot, fill = hex)) +
    scale_fill_identity() +
    geom_contour(
      data = dfc,
      inherit.aes = FALSE,
      aes(x = x, y = y, z = val),
      breaks = 0.5,
      colour = col,
      linewidth = size,
      alpha = alpha
    ) +
    coord_equal() +
    theme_void()
}

# --- Core pipeline for one cube ---
process_cube <- function(cube,
                         smart_bg_q = 0.2,
                         edge_q = 0.75,
                         drop_frame = TRUE,
                         close_k = 3L, close_iters = 1L,
                         conn = 8L, min_size = 25L,
                         r = 7, g = 4, b = 2,
                         pansharpen = 0.5, guide_band = 2,
                         upscale = 1L, unsharp_sigma = 1.1, unsharp_amount = 0.7,
                         sat = 0.9, gamma = 1.0,
                         rotate_ccw = TRUE) {
  
  P <- smart_sum(cube, method = "localivar", bg_q = smart_bg_q)
  
  L_edge_groups <- segment_sobel(
    P, edge_q = edge_q, drop_frame = drop_frame,
    close_k = close_k, close_iters = close_iters,
    conn = conn, min_size = min_size
  )
  
  cube_rgb <- make_rgb(
    cube, r = r, g = g, b = b,
    pansharpen = pansharpen, guide_band = guide_band,
    upscale = upscale,
    unsharp_sigma = unsharp_sigma, unsharp_amount = unsharp_amount,
    sat = sat, gamma = gamma
  )
  
  cube_rgb <- unwrap_rgb(cube_rgb)
  
  # single plot = rgb + contours
  p_rgb_contours <- add_contour(cube_rgb, L_edge_groups, rotate_ccw = rotate_ccw)
  
  # mask view (optional)
  p_mask <- plot_mask_bw(L_edge_groups, rotate_ccw = rotate_ccw)
  
  list(P = P,
       L = L_edge_groups,
       rgb = cube_rgb,
       p_rgb_contours = p_rgb_contours,
       p_mask = p_mask)
}


# Black/white mask plot; NA -> black
plot_mask_bw <- function(L, rotate_ccw = TRUE) {
  m <- if (rotate_ccw) t(apply(L, 2, rev)) else L
  df <- expand.grid(y = seq_len(nrow(m)), x = seq_len(ncol(m)))
  v <- as.integer(m[cbind(df$y, df$x)]) > 0L
  v[is.na(v)] <- NA
  df$val <- v
  ggplot(df, aes(x = x, y = y, fill = val)) +
    geom_raster() +
    scale_fill_manual(
      values = c("FALSE" = "#000000", "TRUE" = "#FFFFFF"),
      guide = "none",
      na.value = "#000000"
    ) +
    coord_equal() + theme_void()
}

# --- Core pipeline for one cube ---
process_cube <- function(cube,
                         smart_bg_q = 0.2,
                         edge_q = 0.75,
                         drop_frame = TRUE,
                         close_k = 3L, close_iters = 1L,
                         conn = 8L, min_size = 25L,
                         r = 7, g = 4, b = 2,
                         pansharpen = 0.5, guide_band = 2,
                         upscale = 1L, unsharp_sigma = 1.1, unsharp_amount = 0.7,
                         sat = 0.9, gamma = 1.0,
                         rotate_ccw = TRUE) {
  
  P <- smart_sum(cube, method = "localivar", bg_q = smart_bg_q)
  
  L_edge_groups <- segment_sobel(
    P, edge_q = edge_q, drop_frame = drop_frame,
    close_k = close_k, close_iters = close_iters,
    conn = conn, min_size = min_size
  )
  
  cube_rgb <- make_rgb(
    cube, r = r, g = g, b = b,
    pansharpen = pansharpen, guide_band = guide_band,
    upscale = upscale,
    unsharp_sigma = unsharp_sigma, unsharp_amount = unsharp_amount,
    sat = sat, gamma = gamma
  )
  
  cube_rgb <- unwrap_rgb(cube_rgb)
  
  # single plot = rgb + contours
  p_rgb_contours <- add_contour(cube_rgb, L_edge_groups, rotate_ccw = rotate_ccw)
  
  # mask view (optional)
  p_mask <- plot_mask_bw(L_edge_groups, rotate_ccw = rotate_ccw)
  
  list(P = P,
       L = L_edge_groups,
       rgb = cube_rgb,
       p_rgb_contours = p_rgb_contours,
       p_mask = p_mask)
}



# --- Run pipeline over all cubes ---
cubes <- lapply(files, read_cube)
res_list <- lapply(cubes, process_cube)

rgb_plots  <- lapply(res_list, `[[`, "p_rgb_contours")
mask_plots <- lapply(res_list, `[[`, "p_mask")

# --- Layout mosaics ---
layout_rgb  <- (rgb_plots[[1]] | rgb_plots[[4]] | rgb_plots[[7]]) /
  (rgb_plots[[2]] | rgb_plots[[3]] | rgb_plots[[5]] | rgb_plots[[6]]) &
  theme(panel.spacing = grid::unit(0, "pt"),
        plot.margin   = grid::unit(c(0,0,0,0), "pt"))

layout_mask <- (mask_plots[[1]] | mask_plots[[4]] | mask_plots[[7]]) /
  (mask_plots[[2]] | mask_plots[[3]] | mask_plots[[5]] | mask_plots[[6]]) &
  theme(panel.spacing = grid::unit(0, "pt"),
        plot.margin   = grid::unit(c(0,0,0,0), "pt"))

title_block <- plot_annotation(
  title = "**<span style='color:black;'>JWST/NIRCam</span>**"
) & theme(
  plot.title = element_markdown(size = 48, face = "bold", hjust = 0.5,
                                margin = margin(b = 6)),
  plot.title.position = "plot"
)

rgb_mosaic  <- layout_rgb  + title_block
mask_mosaic <- layout_mask + title_block

# --- Save PNGs ---
ragg::agg_png("mosaic_rgb_contours.png", width = 8000, height = 4500,
              units = "px", background = "black")
print(rgb_mosaic)
dev.off()

ragg::agg_png("mosaic_masks_bw.png", width = 8000, height = 4500,
              units = "px", background = "black")
print(mask_mosaic)
dev.off()
