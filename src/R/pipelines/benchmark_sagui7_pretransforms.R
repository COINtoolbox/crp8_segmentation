suppressPackageStartupMessages({
  library(FITSio)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(readr)
  library(tidyr)
  library(viridis)
})

repo_root <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation"
sagui_root <- "/Users/rd23aag/Documents/GitHub/sagui"

cfg <- list(
  cube_path = file.path(repo_root, "data/raw/datacube_sagui7.fits"),
  figure_dir = file.path(repo_root, "results/figures/sagui_seg/benchmarks"),
  segmentation_dir = file.path(repo_root, "results/segmentation/sagui_benchmarks/sagui7"),
  Ncomp = 20,
  starlet_J = 5,
  starlet_scales = 1:5,
  kclip = 0.5,
  mask_mode = "na",
  hclust_method = "ward.D2",
  seed = 42
)

transform_specs <- list(
  list(label = "none", method = "none"),
  list(label = "clipped_log1p", method = function(x) log1p(pmax(x, 0))),
  list(label = "asinh", method = "asinh"),
  list(label = "copula_uniform", method = "copula_uniform"),
  list(label = "copula_gaussian", method = "copula_gaussian")
)

dir.create(cfg$figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg$segmentation_dir, showWarnings = FALSE, recursive = TRUE)

load_local_sagui <- function(path) {
  if (requireNamespace("pkgload", quietly = TRUE) && dir.exists(path)) {
    pkgload::load_all(path, quiet = TRUE)
    return(invisible(TRUE))
  }

  library(sagui)
  invisible(FALSE)
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

safe_label <- function(x) {
  gsub("[^a-z0-9_]+", "_", tolower(x))
}

matrix_to_df <- function(mat) {
  df <- as.data.frame(as.table(mat))
  names(df) <- c("y", "x", "value")
  df$y <- as.integer(df$y)
  df$x <- as.integer(df$x)
  df
}

normalize_panel <- function(v) {
  s <- asinh_stretch(v)
  rng <- range(s, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(rep(NA_real_, length(v)))
  }
  (s - rng[1]) / diff(rng)
}

expand_mask <- function(mask, iterations = 1L) {
  out <- mask
  nr <- nrow(mask)
  nc <- ncol(mask)

  for (iter in seq_len(iterations)) {
    expanded <- out
    idx <- which(out, arr.ind = TRUE)
    if (!nrow(idx)) {
      return(out)
    }

    for (k in seq_len(nrow(idx))) {
      i <- idx[k, 1]
      j <- idx[k, 2]
      ii <- max(1, i - 1):min(nr, i + 1)
      jj <- max(1, j - 1):min(nc, j + 1)
      expanded[ii, jj] <- TRUE
    }

    out <- expanded
  }

  out
}

build_relaxed_mask_info <- function(input, cfg, pretransform = "none") {
  cube_for_mask <- pretransform_cube(input$imDat, method = pretransform)
  collapsed <- collapse_white_light(cube_for_mask, kclip = cfg$kclip)
  dec <- starlet_mask(collapsed, J = cfg$starlet_J)

  rec_compact <- starlet_reconstruct(
    dec,
    keep_scales = 1:3,
    include_coarse = FALSE,
    denoise_k = 0,
    mode = "soft"
  )

  rec_diffuse <- starlet_reconstruct(
    dec,
    keep_scales = 2:cfg$starlet_J,
    include_coarse = FALSE,
    denoise_k = 0,
    mode = "soft"
  )

  mask <- (is.finite(rec_compact) & rec_compact > 0) |
    (is.finite(rec_diffuse) & rec_diffuse > 0)

  mask <- expand_mask(mask, iterations = 1L)

  list(
    collapsed = collapsed,
    decomposition = dec,
    reconstruction = list(compact = rec_compact, diffuse = rec_diffuse),
    mask = mask
  )
}

cluster_sizes <- function(cluster_map) {
  ids <- sort(unique(as.integer(cluster_map[is.finite(cluster_map) & cluster_map > 0])))
  if (!length(ids)) {
    return(integer())
  }

  vapply(ids, function(id) sum(cluster_map == id, na.rm = TRUE), integer(1))
}

save_scalar_map_png <- function(mat, out_path, title = NULL, mask = FALSE) {
  df <- matrix_to_df(mat)

  if (isTRUE(mask)) {
    df$value <- factor(ifelse(is.finite(df$value) & df$value > 0, "mask", "background"))
    p <- ggplot(df, aes(x = y, y = x, fill = value)) +
      geom_raster() +
      coord_fixed() +
      scale_fill_manual(values = c(background = "#0A1026", mask = "#F1C76A")) +
      theme_void(base_size = 14) +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.background = element_rect(fill = "white", colour = NA)
      )
  } else {
    df$value_norm <- normalize_panel(df$value)
    p <- ggplot(df, aes(x = y, y = x, fill = value_norm)) +
      geom_raster() +
      coord_fixed() +
      scale_fill_gradientn(
        colours = palette_van_gogh_div(256),
        na.value = "#0A1026",
        limits = c(0, 1)
      ) +
      theme_void(base_size = 14) +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.background = element_rect(fill = "white", colour = NA)
      )
  }

  if (!is.null(title)) {
    p <- p + labs(title = title)
  }

  ggplot2::ggsave(
    filename = out_path,
    plot = p,
    width = 5,
    height = 5,
    dpi = 300,
    bg = "white"
  )
}

summarize_segmentation <- function(seg, transform, elapsed_sec) {
  sizes <- cluster_sizes(seg$cluster_map)

  tibble::tibble(
    transform = transform,
    elapsed_sec = elapsed_sec,
    mask_pixels = sum(seg$mask, na.rm = TRUE),
    labeled_pixels = sum(is.finite(seg$cluster_map) & seg$cluster_map > 0),
    n_regions = length(sizes),
    min_region_pixels = if (length(sizes)) min(sizes) else NA_integer_,
    median_region_pixels = if (length(sizes)) as.numeric(stats::median(sizes)) else NA_real_,
    mean_region_pixels = if (length(sizes)) mean(sizes) else NA_real_,
    max_region_pixels = if (length(sizes)) max(sizes) else NA_integer_,
    mean_cluster_snr = mean(seg$cluster_snr, na.rm = TRUE),
    median_cluster_snr = stats::median(seg$cluster_snr, na.rm = TRUE)
  )
}

boundary_from_labels <- function(cluster_map) {
  m <- cluster_map
  nr <- nrow(m)
  nc <- ncol(m)
  edge <- matrix(FALSE, nrow = nr, ncol = nc)
  valid <- is.finite(m) & m > 0

  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      if (!valid[i, j]) {
        next
      }

      ngh_i <- max(1, i - 1):min(nr, i + 1)
      ngh_j <- max(1, j - 1):min(nc, j + 1)
      neigh <- m[ngh_i, ngh_j]
      if (any(is.finite(neigh) & neigh > 0 & neigh != m[i, j])) {
        edge[i, j] <- TRUE
      }
    }
  }

  edge
}

write_cluster_fits <- function(cluster_map, path) {
  out <- cluster_map
  out[!is.finite(out)] <- 0L
  FITSio::writeFITSim(out, file = path)
}

format_transform_label <- function(label, elapsed_sec, n_regions) {
  pretty_label <- gsub("_", " ", label, fixed = TRUE)
  sprintf("%s\n%.2fs | %d regions", pretty_label, elapsed_sec, n_regions)
}

make_segmentation_plot <- function(seg, title_label, ncomp) {
  plot_region_map(
    seg,
    palette = palette_van_gogh_div(ncomp),
    border_color = "black",
    border_linewidth = 0.9,
    background_color = "transparent"
  ) +
    ggplot2::labs(title = title_label) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 11,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 7)
      ),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )
}

make_overlay_plot <- function(seg, collapsed, title_label) {
  base_df <- matrix_to_df(collapsed)
  base_df$value_norm <- normalize_panel(base_df$value)

  boundary_df <- matrix_to_df(boundary_from_labels(seg$cluster_map))
  boundary_df <- boundary_df[isTRUE(boundary_df$value) | boundary_df$value == 1, , drop = FALSE]

  ggplot(base_df, aes(x = y, y = x)) +
    geom_raster(aes(fill = value_norm)) +
    geom_raster(
      data = boundary_df,
      fill = "black"
    ) +
    coord_fixed() +
    scale_fill_gradientn(
      colours = palette_van_gogh_div(256),
      na.value = "#0A1026",
      limits = c(0, 1)
    ) +
    labs(title = title_label) +
    theme_void(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5, margin = margin(b = 7)),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(6, 6, 6, 6)
    )
}

make_segmentation_panel <- function(plot_list) {
  patchwork::wrap_plots(plotlist = plot_list, ncol = 3) +
    patchwork::plot_annotation(
      title = "sagui7 segmentation benchmark",
      subtitle = "Starlet mask built on transformed data; clustering on original cube",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA)
      )
    )
}

make_overlay_panel <- function(plot_list) {
  patchwork::wrap_plots(plotlist = plot_list, ncol = 3) +
    patchwork::plot_annotation(
      title = "sagui7 boundary overlays",
      subtitle = "Mask on transformed data, boundaries over original white-light image",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA)
      )
    )
}

make_metrics_plot <- function(summary_df) {
  metrics_df <- summary_df |>
    dplyr::select(transform, elapsed_sec, labeled_pixels, median_region_pixels, mean_cluster_snr) |>
    tidyr::pivot_longer(
      cols = -transform,
      names_to = "metric",
      values_to = "value"
    ) |>
    dplyr::mutate(
      metric = factor(
        metric,
        levels = c("elapsed_sec", "labeled_pixels", "median_region_pixels", "mean_cluster_snr"),
        labels = c("Runtime (s)", "Labeled Pixels", "Median Region Size", "Mean Cluster SNR")
      ),
      transform = factor(transform, levels = summary_df$transform)
    )

  metric_palette <- setNames(
    palette_van_gogh_div(nlevels(metrics_df$transform)),
    levels(metrics_df$transform)
  )

  ggplot(metrics_df, aes(x = transform, y = value, fill = transform)) +
    geom_col(show.legend = FALSE, colour = "black", linewidth = 0.2) +
    facet_wrap(~ metric, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = metric_palette) +
    theme_bw(base_size = 13) +
    theme(
      strip.background = element_rect(fill = "#d9bea1", colour = NA),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1),
      axis.title = element_blank()
    ) +
    labs(
      title = "sagui7 transform summary",
      subtitle = "Same starlet mask and target number of regions"
    )
}

save_individual_map <- function(seg, transform, figure_dir, ncomp) {
  p <- plot_region_map(
    seg,
    palette = palette_van_gogh_div(ncomp),
    border_color = "black",
    border_linewidth = 0.8,
    background_color = "transparent"
  )

  out_path <- file.path(
    figure_dir,
    paste0("sagui7_", safe_label(transform), ".png")
  )

  ggplot2::ggsave(
    filename = out_path,
    plot = p,
    width = 5,
    height = 5,
    dpi = 300,
    bg = "transparent"
  )

  out_path
}

save_overlay_map <- function(seg, collapsed, transform, figure_dir) {
  out_path <- file.path(
    figure_dir,
    paste0("sagui7_", safe_label(transform), "_overlay.png")
  )

  p <- make_overlay_plot(
    seg,
    collapsed,
    gsub("_", " ", transform, fixed = TRUE)
  )

  ggplot2::ggsave(
    filename = out_path,
    plot = p,
    width = 5,
    height = 5,
    dpi = 300,
    bg = "white"
  )

  out_path
}

set.seed(cfg$seed)
load_local_sagui(sagui_root)

X <- FITSio::readFITS(cfg$cube_path)
collapsed <- collapse_white_light(X$imDat, kclip = cfg$kclip)

FITSio::writeFITSim(
  collapsed,
  file = file.path(cfg$segmentation_dir, "sagui7_collapsed_white_light.fits")
)
save_scalar_map_png(
  collapsed,
  file.path(cfg$figure_dir, "sagui7_collapsed_white_light.png"),
  title = "sagui7 white-light"
)

summary_rows <- vector("list", length(transform_specs))
seg_plots <- vector("list", length(transform_specs))
overlay_plots <- vector("list", length(transform_specs))
figure_paths <- character(length(transform_specs))
overlay_paths <- character(length(transform_specs))
mask_paths <- character(length(transform_specs))

for (i in seq_along(transform_specs)) {
  transform_label <- transform_specs[[i]]$label
  transform_method <- transform_specs[[i]]$method
  message("Running transform: ", transform_label)

  mask_info <- build_relaxed_mask_info(X, cfg, pretransform = transform_method)
  mask_fits_path <- file.path(
    cfg$segmentation_dir,
    paste0("sagui7_", safe_label(transform_label), "_mask.fits")
  )
  FITSio::writeFITSim(ifelse(mask_info$mask, 1L, 0L), file = mask_fits_path)
  mask_paths[[i]] <- file.path(
    cfg$figure_dir,
    paste0("sagui7_", safe_label(transform_label), "_mask.png")
  )
  save_scalar_map_png(
    ifelse(mask_info$mask, 1L, 0L),
    mask_paths[[i]],
    title = paste("mask:", gsub("_", " ", transform_label, fixed = TRUE)),
    mask = TRUE
  )

  masked_input <- X
  masked_input$imDat <- mask_cube(X$imDat, mask_info$mask, mode = cfg$mask_mode)

  timing <- system.time({
    seg <- segment_regions(
      input = masked_input,
      Ncomp = cfg$Ncomp,
      use_starlet_mask = FALSE,
      collapse_fn = function(cube) collapse_white_light(cube, kclip = cfg$kclip),
      starlet_J = cfg$starlet_J,
      starlet_scales = cfg$starlet_scales,
      mask_mode = cfg$mask_mode,
      hclust_method = cfg$hclust_method,
      cluster_pretransform = "none"
    )
  })

  seg$mask <- mask_info$mask
  seg$collapsed <- mask_info$collapsed
  seg$starlet <- list(
    decomposition = mask_info$decomposition,
    reconstruction = mask_info$reconstruction
  )

  summary_rows[[i]] <- summarize_segmentation(seg, transform_label, timing[["elapsed"]])
  panel_label <- format_transform_label(
    transform_label,
    timing[["elapsed"]],
    summary_rows[[i]]$n_regions
  )
  seg_plots[[i]] <- make_segmentation_plot(seg, panel_label, cfg$Ncomp)
  overlay_plots[[i]] <- make_overlay_plot(seg, mask_info$collapsed, panel_label)

  write_cluster_fits(
    seg$cluster_map,
    file.path(cfg$segmentation_dir, paste0("sagui7_", safe_label(transform_label), ".fits"))
  )
  figure_paths[[i]] <- save_individual_map(seg, transform_label, cfg$figure_dir, cfg$Ncomp)
  overlay_paths[[i]] <- save_overlay_map(seg, mask_info$collapsed, transform_label, cfg$figure_dir)
}

summary_df <- dplyr::bind_rows(summary_rows) |>
  dplyr::arrange(match(transform, vapply(transform_specs, `[[`, character(1), "label")))

readr::write_csv(
  summary_df,
  file.path(cfg$segmentation_dir, "sagui7_pretransform_summary.csv")
)
readr::write_csv(
  tibble::tibble(
    transform = vapply(transform_specs, `[[`, character(1), "label"),
    figure_path = figure_paths,
    overlay_path = overlay_paths,
    mask_path = mask_paths
  ),
  file.path(cfg$segmentation_dir, "sagui7_pretransform_figure_index.csv")
)

seg_panel <- make_segmentation_panel(seg_plots)
metrics_panel <- make_metrics_plot(summary_df)
overlay_panel <- make_overlay_panel(overlay_plots)

ggplot2::ggsave(
  filename = file.path(cfg$figure_dir, "sagui7_pretransform_benchmark_panel.png"),
  plot = seg_panel,
  width = 13,
  height = 8.5,
  dpi = 300,
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(cfg$figure_dir, "sagui7_pretransform_benchmark_metrics.png"),
  plot = metrics_panel,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(cfg$figure_dir, "sagui7_pretransform_overlay_panel.png"),
  plot = overlay_panel,
  width = 13,
  height = 8.5,
  dpi = 300,
  bg = "white"
)
