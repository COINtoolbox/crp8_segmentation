export_starlet_panels_png <- function(
    collapsed,
    dec,
    J,
    out_dir,
    tag,
    prefix = "starlet_2d",
    px = 1600,
    res = 300,
    denoise_k = 1,
    mode = "soft",
    keep_include_coarse = FALSE,
    na_transparent = TRUE   # TRUE => NA pixels become transparent
) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # helper: panel normalization (same as your logic)
  normalize_panel <- function(v) {
    s <- asinh_stretch(v)
    rng <- range(s, na.rm = TRUE)
    (s - rng[1]) / (rng[2] - rng[1])
  }
  
  # helper: build a ggplot for a single matrix
  plot_one <- function(mat) {
    df <- mat_to_df(mat, "tmp")          # your function already sets <=0 to NA
    df$value_norm <- normalize_panel(df$value)
    
    ggplot(df, aes(x = y, y = x, fill = value_norm)) +
      geom_raster() +
      coord_fixed(expand = FALSE) +
      scale_fill_gradientn(
        colours = palette_van_gogh_div(256),
        na.value = if (isTRUE(na_transparent)) NA else "black",
        oob = scales::squish
      ) +
      theme_void() +
      theme(
        legend.position  = "none",
        plot.margin      = grid::unit(c(0, 0, 0, 0), "pt"),
        panel.background = element_rect(fill = NA, colour = NA),
        plot.background  = element_rect(fill = NA, colour = NA)
      )
  }
  
  # file naming helper
  fpath <- function(suffix) {
    file.path(out_dir, sprintf("%s_%s_%s.png", prefix, tag, suffix))
  }
  
  # 1) Original
  p0 <- plot_one(collapsed)
  ragg::agg_png(fpath("Original"), width = px, height = px, res = res, background = "transparent")
  print(p0)
  dev.off()
  
  # 2) Each scale j as its own PNG
  for (j in seq_len(J)) {
    rec_j <- starlet_reconstruct(
      dec,
      keep_scales    = j,
      include_coarse = keep_include_coarse,
      denoise_k      = denoise_k,
      mode           = mode
    )
    
    pj <- plot_one(rec_j)
    ragg::agg_png(fpath(sprintf("j_%02d", j)), width = px, height = px, res = res, background = "transparent")
    print(pj)
    dev.off()
  }
  
  invisible(TRUE)
}
