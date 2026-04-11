suppressPackageStartupMessages({
  library(png)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(ragg)
})

# ============================================================
# assemble_segmentation_mosaic.R
# Purpose:
#   - Rebuild the segmentation mosaic from the existing sagui
#     segmentation PNGs,
#   - keep the current outputs untouched,
#   - write a dedicated paper-panel version in an organized folder.
# ============================================================

base_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation"
seg_dir <- file.path(base_dir, "results/figures/sagui_seg")
out_dir <- file.path(base_dir, "results/figures/paper_panels/segmentation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

panel_order <- c(
  "sagui1_2_3",
  "sagui4",
  "sagui5_6",
  "sagui7",
  "sagui8",
  "sagui9"
)

pretty_tag <- function(tag) {
  sub("^sagui", "sagui-", gsub("_", "/", tag))
}

make_seg_panel <- function(image_path, tag_label) {
  img <- png::readPNG(image_path)
  grob <- grid::rasterGrob(img, interpolate = FALSE)
  label_df <- data.frame(x = 0.03, y = 0.04, label = tag_label)

  ggplot() +
    annotation_custom(grob, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    geom_label(
      data = label_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 0,
      fill = "#7A7A7A",
      colour = "white",
      fontface = "bold",
      family = "serif",
      size = 8.5,
      label.padding = grid::unit(0.2, "lines"),
      label.r = grid::unit(0.12, "lines"),
      linewidth = 0
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    theme_void() +
    theme(
      panel.border = element_rect(colour = "#7C7C7C", fill = NA, linewidth = 0.7),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(1, 1, 1, 1)
    )
}

plots <- list()
for (tag in panel_order) {
  image_path <- file.path(seg_dir, paste0(tag, ".png"))
  if (!file.exists(image_path)) {
    next
  }
  plots[[tag]] <- make_seg_panel(image_path, pretty_tag(tag))
}

if (!length(plots)) {
  stop("No segmentation PNGs found to assemble.")
}

mosaic <- patchwork::wrap_plots(plots, ncol = 3) &
  theme(
    plot.margin = margin(0, 0, 0, 0),
    panel.spacing = grid::unit(0, "pt")
  )

out_path <- file.path(out_dir, "sagui_segmentation_mosaic.png")
ragg::agg_png(out_path, width = 4200, height = 2800, units = "px", background = "white")
print(mosaic)
dev.off()
message("Wrote segmentation mosaic: ", out_path)
