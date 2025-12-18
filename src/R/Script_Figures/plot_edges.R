# --- Packages ---
library(FITSio)
library(ggplot2)
library(ggtext)
library(patchwork)
library(ragg)

# If you use purrr/dplyr habits:
# library(purrr); library(dplyr)

# --- Inputs (edit paths as needed) ---
files <- sprintf("..//data/raw/datacube_reg%d.fits", c(1,2,4,5,6,7,8))

# --- Small helpers ---
read_cube <- function(path) {
  x <- FITSio::readFITS(path)
  x$imDat
}

plot_rgb_masked <- function(cube_rgb, L, rotate_ccw = TRUE) {
  df <- rgb_mask_to_df(cube_rgb, L, rotate_ccw = rotate_ccw)
  ggplot(df, aes(x = -x_plot, y = y_plot, fill = hex)) +
    geom_raster() +
    scale_fill_identity() +
    coord_equal() +
    theme_void()
}

plot_mask_bw <- function(L, rotate_ccw = TRUE) {
  # Apply the same rotation logic as RGB
  m <- if (rotate_ccw) t(apply(L, 2, rev)) else L
  
  # Build a grid in the transformed coordinates
  df <- expand.grid(y_plot = seq_len(nrow(m)),
                    x_plot = seq_len(ncol(m)))
  df$val <- as.integer(m[cbind(df$y_plot, df$x_plot)]) > 0L
  
  ggplot(df, aes(x = -x_plot, y = y_plot, fill = val)) +  # <- mirror x to match RGB
    geom_raster() +
    scale_fill_manual(values = c("FALSE" = "#000000", "TRUE" = "#FFFFFF"), guide = "none") +
    coord_equal() +
    theme_void()
}


# --- Core per-cube pipeline ---
process_cube <- function(cube,
                         # Energy map:
                         smart_bg_q = 0.2,
                         # Sobel child segmentation:
                         edge_q = 0.75, drop_frame = TRUE,
                         close_k = 3L, close_iters = 1L,
                         conn = 8L, min_size = 25L,
                         # RGB:
                         r = 7, g = 4, b = 2,
                         pansharpen = 0.5, guide_band = 2,
                         upscale = 1L, unsharp_sigma = 1.1, unsharp_amount = 0.7,
                         sat = 0.9, gamma = 1.0,
                         rotate_ccw = TRUE) {
  
  H <- dim(cube)[1]; W <- dim(cube)[2]
  
  # Option A (fast, robust in faint regions): local IVAR energy
  P <- smart_sum(cube, method = "localivar", bg_q = smart_bg_q)
  
  # Child segmentation (Sobel groups → clean → fill)
  L_edge_groups <- segment_sobel(
    P, edge_q = edge_q, drop_frame = drop_frame,
    close_k = close_k, close_iters = close_iters,
    conn = conn, min_size = min_size
  )
  
  # Build RGB
  cube_rgb <- make_rgb(
    cube, r = r, g = g, b = b,
    pansharpen = pansharpen, guide_band = guide_band,
    upscale = upscale,
    unsharp_sigma = unsharp_sigma, unsharp_amount = unsharp_amount,
    sat = sat, gamma = gamma
  )
  
  # Masked RGB plot
  p_rgb_masked <- plot_rgb_masked(cube_rgb, L_edge_groups, rotate_ccw = rotate_ccw)
  
  # Pure mask (B/W) plot
  p_mask <- plot_mask_bw(L_edge_groups, rotate_ccw = rotate_ccw) 
  
  list(P = P,
       L = L_edge_groups,
       rgb = cube_rgb,
       p_rgb_masked = p_rgb_masked,
       p_mask = p_mask)
}

# --- Read all cubes ---
cubes <- lapply(files, read_cube)

# --- Run the pipeline over all cubes ---
res_list <- lapply(cubes, process_cube)

# Collect plots
rgb_plots  <- lapply(res_list, `[[`, "p_rgb_masked")
mask_plots <- lapply(res_list, `[[`, "p_mask")

# --- Arrange mosaics (adjust layout to taste) ---
# Example layout (7 panels): (1 | 4 | 7) / (2 | 3 | 5 | 6)
layout_rgb  <- (rgb_plots[[1]] | rgb_plots[[4]] | rgb_plots[[7]]) /
  (rgb_plots[[2]] | rgb_plots[[3]] | rgb_plots[[5]] | rgb_plots[[6]]) &
  theme(panel.spacing = grid::unit(0, "pt"),
        plot.margin   = grid::unit(c(0,0,0,0), "pt"))

layout_mask <- (mask_plots[[1]] | mask_plots[[4]] | mask_plots[[7]]) /
  (mask_plots[[2]] | mask_plots[[3]] | mask_plots[[5]] | mask_plots[[6]]) &
  theme(panel.spacing = grid::unit(0, "pt"),
        plot.margin   = grid::unit(c(0,0,0,0), "pt"))

# Optional title (kept minimal & black per your preference)
title_block <- plot_annotation(
  title = "**<span style='color:black;'>JWST/NIRCam</span>**"
) & theme(
  plot.title = element_markdown(size = 48, face = "bold", hjust = 0.5,
                                margin = margin(b = 6)),
  plot.title.position = "plot"
)

rgb_mosaic  <- layout_rgb  + title_block
mask_mosaic <- layout_mask + title_block

# --- Save high-res PNGs ---
ragg::agg_png("mosaic_rgb_masked.png", width = 8000, height = 4500, units = "px", background = "black")
print(rgb_mosaic)
dev.off()

ragg::agg_png("mosaic_masks_bw.png", width = 8000, height = 4500, units = "px", 
              background = "black")
print(mask_mosaic)
dev.off()
