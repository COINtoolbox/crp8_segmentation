# ======================================================================
# 1) Build and save a tidy table of NIRCam filter throughputs
# ======================================================================
make_throughput_table <- function(
    input,                                    # dir path OR character vector of files/URLs
    filters,                                  # e.g., c("F090W","F115W",...)
    nircam_lambda_um,                         # named numeric: central wavelength by filter
    exclude = c("F150W2"),            # odd filenames to ignore
    output = "throughput_nircam.csv",         # output file
    format = c("csv","rds"),                  # write CSV or RDS
    preview = FALSE, preview_path = "throughput_preview.png"
){
  stopifnot(requireNamespace("dplyr", quietly=TRUE))
  stopifnot(requireNamespace("purrr", quietly=TRUE))
  stopifnot(requireNamespace("stringr", quietly=TRUE))
  stopifnot(requireNamespace("ggplot2", quietly=TRUE))
  stopifnot(requireNamespace("readr", quietly=TRUE))
  
  library(dplyr); library(purrr); library(stringr); library(ggplot2); library(readr)
  
  format <- match.arg(format)
  
  # Collect files
  files <- if (length(input) == 1 && dir.exists(input)) {
    list.files(input, pattern = "_mean_system_throughput\\.txt$", full.names = TRUE)
  } else {
    as.character(input)
  }
  files <- files[str_detect(basename(files), "_mean_system_throughput\\.txt$")]
  if (!length(files)) stop("No throughput files found.")
  
  # Keep files that START with desired filters & not excluded
  keep_pat    <- paste0("^(", paste(filters, collapse="|"), ")_")
  exclude_pat <- paste0("^(", paste(exclude,  collapse="|"), ")_")
  base <- basename(files)
  files <- files[str_detect(base, keep_pat) & !str_detect(base, exclude_pat)]
  if (!length(files)) stop("No files left after filtering.")
  
  # Reader
  read_one <- function(f){
    df <- suppressMessages(readr::read_table(f, show_col_types = FALSE))
    if (ncol(df) < 2) stop(paste("Bad file:", f))
    nm <- stringr::str_extract(basename(f), "^F[0-9]{3,4}[A-Z]+")
    if (is.na(nm)) nm <- stringr::str_extract(basename(f), "^[^_]+")
    tibble::tibble(
      wavelength = df[[1]],
      throughput = df[[2]],
      filter     = nm,
      lambda_c   = unname(nircam_lambda_um[nm])
    )
  }
  
  throughput_all <- purrr::map_dfr(files, read_one) %>%
    filter(!filter %in% exclude, filter %in% filters) %>%
    arrange(filter, wavelength) %>%
    mutate(filter = factor(filter, levels = names(sort(nircam_lambda_um))))
  
  # Write
  if (format == "csv") readr::write_csv(throughput_all, output) else saveRDS(throughput_all, output)
  
  if (preview) {
    p <- ggplot(throughput_all, aes(wavelength, throughput, color = filter)) +
      geom_line(linewidth=0.6) + theme_minimal(base_size=12) +
      labs(x="Wavelength (μm)", y="System Throughput", color="Filter") +
      scale_color_viridis_d(option="turbo")
    ggsave(preview_path, p, width=8, height=4.6, dpi=150)
  }
  
  message(sprintf("Throughput table saved to %s  (%d rows, %d filters).",
                  output, nrow(throughput_all), dplyr::n_distinct(throughput_all$filter)))
  invisible(throughput_all)
}

# ======================================================================
# 2) Plot SED + filters for any region
#    scale_mode = "auto" (raw units; filters occupy small data-driven height)
#               or "normalized" (SED to 1; errors normalized consistently)
# ======================================================================
plot_sed_with_filters <- function(
    reg_long,                   # tibble with: region, band, flux, flux_err
    region_id,                  # integer region ID
    throughput_tbl,             # output of make_throughput_table()
    filters,                    # vector of JWST filter names in correct order
    nircam_lambda_um,           # named numeric vector: F090W = 0.902, ...
    scale_mode = c("normalized","auto"),
    band_height_frac = 0.15,
    errors_are_fractional = FALSE,
    show_legend = TRUE
){
  library(dplyr)
  library(ggplot2)
  library(rlang)
  
  scale_mode <- match.arg(scale_mode)
  
  # ----------- SED for a single region -----------
  sed <- reg_long %>%
    filter(region == region_id) %>%
    mutate(
      band_idx  = as.integer(as.character(band)),     # "0","1",... -> 0,1,...
      filt      = filters[band_idx + 1L],             # 0 -> F090W, 1 -> F115W, ...
      lambda_um = nircam_lambda_um[filt]
    ) %>%
    filter(!is.na(lambda_um)) %>%
    arrange(lambda_um)
  
  stopifnot(nrow(sed) > 0)
  
  # ----------- Throughput table -----------
  if (!"filter" %in% names(throughput_tbl)) {
    stop("throughput_tbl must have a column named 'filter'. ",
         "Got: ", paste(names(throughput_tbl), collapse = ", "))
  }
  
  thr <- throughput_tbl %>%
    mutate(filter = as.character(.data$filter)) %>%   # ensure it's a vector, not factor
    filter(.data$filter %in% filters)
  
  # ----------- Scaling -----------
  
  if (scale_mode == "normalized") {
    K <- max(sed$flux, na.rm = TRUE)
    
    sed <- sed %>%
      mutate(
        y = flux / K,
        ymin = if (errors_are_fractional) {
          (flux * (1 - flux_err)) / K
        } else {
          (flux - flux_err) / K
        },
        ymax = if (errors_are_fractional) {
          (flux * (1 + flux_err)) / K
        } else {
          (flux + flux_err) / K
        }
      )
    
    sed$ymin <- pmax(sed$ymin, 0)
    
    thr <- thr %>%
      mutate(y = throughput * band_height_frac)
    
    ylab    <- "Normalized Flux"
    ylimits <- c(-0.02, 1.05)
    
  } else {  # "auto"
    sed_scale <- quantile(sed$flux, 0.95, na.rm = TRUE)
    
    thr <- thr %>%
      mutate(y = throughput * sed_scale * band_height_frac)
    
    sed <- sed %>%
      mutate(
        y = flux,
        ymin = if (errors_are_fractional) {
          pmax(flux * (1 - flux_err), 0)
        } else {
          pmax(flux - flux_err, 0)
        },
        ymax = if (errors_are_fractional) {
          flux * (1 + flux_err)
        } else {
          flux + flux_err
        }
      )
    
    ylab    <- "Flux"
    ylimits <- NULL
  }
  
  # ----------- Plot -----------
  
  p <- ggplot() +
    geom_area(
      data  = thr,
      aes(x = wavelength, y = y, fill = filter, group = filter),
      position = "identity",
      alpha    = 0.35,
      color    = NA,
      size     = 0
    ) +
    geom_errorbar(
      data  = sed,
      aes(x = lambda_um, ymin = ymin, ymax = ymax),
      width = 0.06,
      color = "black"
    ) +
    geom_point(
      data  = sed,
      aes(x = lambda_um, y = y),
      size  = 2.6,
      shape = 21,
      fill  = "orange2",
      color = "black"
    ) +
    labs(
      x     = "Wavelength (µm)",
      y     = ylab,
      title = paste("Region", region_id, "- JWST/NIRCam SED")
    ) +
    scale_color_viridis_d(option = "turbo", name = "Filter") +
    scale_fill_viridis_d(option = "turbo", name = "Filter") +
    theme_bw(base_size = 14) +
    theme(
      legend.position = if (show_legend) "right" else "none",
      panel.grid      = element_blank(),
      axis.line       = element_line(color = "black")
    )
  
  if (!is.null(ylimits)) {
    p <- p + coord_cartesian(ylim = ylimits)
  }
  
  p
}







filters <- c("F090W","F115W","F150W","F182M","F200W",
             "F210M","F277W","F335M","F356W","F410M",
             "F430M","F444W","F480M")

nircam_lambda_um <- c(
  F090W = 0.902, F115W = 1.154, F150W = 1.501, F182M = 1.842,
  F200W = 1.989, F210M = 2.099, F277W = 2.770, F335M = 3.365,
  F356W = 3.563, F410M = 4.082, F430M = 4.286, F444W = 4.421,
  F480M = 4.828
)

# 1) Create and save the filter table (run once, reuse later)
throughput_all <- make_throughput_table(
  input  = "/Users/rd23aag/Downloads/nircam_throughputs/mean_throughputs/",
  filters = filters,
  nircam_lambda_um = nircam_lambda_um,
  exclude = c("F150W2","150W2"),
  output  = "throughput_nircam.csv",
  format  = "csv",
  preview = TRUE
)


plot_sed_with_filters <- function(
    reg_long,                   # data frame with columns: region, band, flux, flux_err
    region_id,                  # integer/label selecting the region
    throughput_tbl,             # output of make_throughput_table()
    nircam_lambda_um,           # named numeric vector
    filters,                    # filters to include (subset of throughput_tbl$filter)
    scale_mode = c("auto","normalized"),
    band_height_frac = 0.12,    # for "auto": fraction of SED scale used by filters
    errors_are_fractional = FALSE, # TRUE if flux_err is relative (fraction)
    show_legend = FALSE
){
  library(dplyr); library(ggplot2)
  
  scale_mode <- match.arg(scale_mode)
  
  # --- Build SED for the region
  sed <- reg_long %>%
    filter(region == region_id, band %in% filters) %>%
    mutate(lambda_um = nircam_lambda_um[band]) %>%
    arrange(lambda_um)
  
  stopifnot(nrow(sed) > 0)
  
  # --- Prepare filters (sorted, grouped, subset)
  thr <- throughput_tbl %>%
    filter(filter %in% filters) %>%
    arrange(filter, wavelength)
  
  # --- Compute y-scaling
  if (scale_mode == "normalized") {
    K <- max(sed$flux, na.rm = TRUE)
    sed <- sed %>%
      mutate(
        y    = flux / K,
        ymin = if (errors_are_fractional) pmax((flux*(1 - flux_err)) / K, 0) else pmax((flux - flux_err) / K, 0),
        ymax = if (errors_are_fractional)        (flux*(1 + flux_err)) / K   else        (flux + flux_err) / K
      )
    thr <- thr %>% mutate(y = throughput * band_height_frac)  # put filters at bottom
    y_lab <- "Normalized Flux"
    y_lim <- c(-0.02, 1.05)
  } else { # "auto" raw-units scaling
    # robust SED scale: use 95th percentile of flux to avoid outliers
    sed_scale <- stats::quantile(sed$flux, 0.95, na.rm = TRUE)
    thr <- thr %>% mutate(y = throughput * (band_height_frac * sed_scale))
    sed <- sed %>%
      mutate(
        y    = flux,
        ymin = pmax(flux - (if (errors_are_fractional) flux*flux_err else flux_err), 0),
        ymax =        flux + (if (errors_are_fractional) flux*flux_err else flux_err)
      )
    y_lab <- "Flux"
    y_lim <- NULL
  }
  
  # --- Plot
  p <- ggplot() +
    geom_area(data = thr,
              aes(x = wavelength, y = y, fill = filter, group = filter),
              alpha = 0.28, color = NA, position = "identity") +
  #  geom_line(data = thr,
 #             aes(x = wavelength, y = y, color = filter, group = filter),
 #             linewidth = 0.7) +
    geom_errorbar(data = sed,
                  aes(x = lambda_um, ymin = ymin, ymax = ymax),
                  width = 0.05, color = "black") +
    geom_point(data = sed,
               aes(x = lambda_um, y = y),
               size = 2, shape = 21, fill = "orange2", color = "black") +
    labs(x = "Wavelength (μm)", y = y_lab) +
    scale_color_viridis_d(option = "turbo", name = "Filter") +
    scale_fill_viridis_d(option = "turbo",  name = "Filter") +
    theme_bw(base_size = 13) +
    theme(
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(face = "bold"),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black")
    )
  
  if (!is.null(y_lim)) p <- p + coord_cartesian(ylim = y_lim)
  return(p)
}

region_id <- 15
p1 <- plot_sed_with_filters(
  reg_long = reg$flux_long,
  region_id = region_id,
  throughput_tbl = throughput_all,
  nircam_lambda_um = nircam_lambda_um,
  filters = filters,
  scale_mode = "normalized",            # or "normalized"
  band_height_frac = 0.95,        # adjust visual height of bands
  errors_are_fractional = FALSE,  # set TRUE if flux_err is fractional
  show_legend = FALSE
)
p1
