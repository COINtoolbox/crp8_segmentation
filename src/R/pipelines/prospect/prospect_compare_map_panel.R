#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(png)
  library(grid)
  library(ggplot2)
  library(patchwork)
})

base_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation/results/figures/sed_maps/prospect/sagui10"
current_dir <- file.path(base_dir, "raw_paper_cool_current")
prospect_dir <- file.path(base_dir, "raw_paper_cool_prospect")
out_dir <- file.path(base_dir, "comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

properties <- c("logM", "age_mw", "sfr", "Av")
labels <- c(
  logM = "log M",
  age_mw = "Age_MW",
  sfr = "SFR",
  Av = "A_V"
)

img_panel <- function(path, title = NULL) {
  img <- png::readPNG(path)
  ggplot() +
    annotation_custom(rasterGrob(img, interpolate = TRUE), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    theme_void() +
    labs(title = title) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.margin = margin(2, 2, 2, 2)
    )
}

top_row <- wrap_plots(lapply(properties, function(prop) {
  img_panel(file.path(current_dir, paste0(prop, "_smooth.png")), labels[[prop]])
}), nrow = 1) +
  plot_annotation(title = "Current result")

bottom_row <- wrap_plots(lapply(properties, function(prop) {
  img_panel(file.path(prospect_dir, paste0(prop, "_smooth.png")), labels[[prop]])
}), nrow = 1) +
  plot_annotation(title = "ProSpect first pass")

panel <- top_row / bottom_row

png_path <- file.path(out_dir, "current_vs_prospect_maps.png")
png(png_path, width = 2600, height = 1800, res = 220, bg = "white")
print(panel)
dev.off()

cat("Wrote:", png_path, "\n")
