suppressPackageStartupMessages({
  library(png)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(ragg)
})

# ============================================================
# assemble_sed_surface_sfr_panels.R
# Purpose:
#   - Build paper-style 2x2 SED-property panels for each sagui tag,
#   - replacing SFR with surface SFR while keeping the existing
#     single-property PNG products untouched.
# ============================================================

base_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation"
map_root <- file.path(base_dir, "results/figures/sed_maps_surface_sfr")
out_dir <- file.path(base_dir, "results/figures/paper_panels/sed_properties")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tag_dirs <- list.dirs(map_root, recursive = FALSE, full.names = TRUE)

property_specs <- list(
  list(file = "logZ_smooth.png", label = "log10(Z/Z⊙)"),
  list(file = "sigma_sfr_100myr_smooth.png", label = "ΣSFR [M⊙ yr⁻¹ kpc⁻²]"),
  list(file = "age_mw_smooth.png", label = "⟨t⟩ₘ [Gyr]"),
  list(file = "Av_smooth.png", label = "A_V [mag]")
)

make_png_panel <- function(image_path, label_text) {
  img <- png::readPNG(image_path)
  grob <- grid::rasterGrob(img, interpolate = FALSE)
  label_df <- data.frame(x = 0.03, y = 0.97, label = label_text)

  ggplot() +
    annotation_custom(grob, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    geom_label(
      data = label_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      fill = "#6B6B6B",
      colour = "white",
      fontface = "bold",
      family = "serif",
      size = 8.6,
      label.padding = grid::unit(0.18, "lines"),
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

for (tag_dir in tag_dirs) {
  tag <- basename(tag_dir)
  panel_inputs <- file.path(tag_dir, vapply(property_specs, `[[`, character(1), "file"))
  if (!all(file.exists(panel_inputs))) {
    next
  }

  plots <- Map(
    make_png_panel,
    image_path = panel_inputs,
    label_text = vapply(property_specs, `[[`, character(1), "label")
  )

  panel <- patchwork::wrap_plots(plots, ncol = 2) &
    theme(
      plot.margin = margin(0, 0, 0, 0),
      panel.spacing = grid::unit(0, "pt")
    )

  out_path <- file.path(out_dir, paste0(tag, "_surface_sfr_panel.png"))
  ragg::agg_png(out_path, width = 2400, height = 2360, units = "px", background = "white")
  print(panel)
  dev.off()
  message("Wrote property panel: ", out_path)
}
