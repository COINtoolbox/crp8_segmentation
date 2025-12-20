## ------------------------------------------------------------
## NIRCam filter definitions & helpers
## ------------------------------------------------------------

# Central wavelengths (µm)
nircam_lambda_um <- c(
  F090W = 0.902, F115W = 1.154, F150W = 1.501, F182M = 1.842,
  F200W = 1.989, F210M = 2.099, F277W = 2.770, F335M = 3.365,
  F356W = 3.563, F410M = 4.082, F430M = 4.286, F444W = 4.421,
  F460M = 4.620, F480M = 4.828
)

nircam_filters <- names(nircam_lambda_um)

build_filter_spec <- function(
    filters,
    lambda,
    throughput_input = NULL,
    exclude = NULL,
    assume_ordered = TRUE,
    band_index = NULL
){
  if (is.null(exclude)) exclude <- character(0)
  filters_use <- setdiff(filters, exclude)

  if (is.null(band_index)) {
    band_index <- seq_along(filters_use) - 1L   # 0-based
  }

  tibble::tibble(
    filter     = filters_use,
    lambda_um  = as.numeric(lambda[filters_use]),
    band_index = band_index
  )
}

band_lookup <- tibble::tibble(
  band      = as.character(seq_along(nircam_filters) - 1L),  # "0","1",...
  band_idx  = seq_along(nircam_filters) - 1L,
  filter    = nircam_filters,
  lambda_um = nircam_lambda_um[nircam_filters]
)

spec_nircam <- build_filter_spec(
  filters          = nircam_filters,
  lambda           = nircam_lambda_um,
  throughput_input = "/Users/rd23aag/Downloads/nircam_throughputs/mean_throughputs",
  exclude          = c("F150W2"),     # unused here but kept for API symmetry
  assume_ordered   = TRUE,
  band_index       = seq_along(nircam_filters) - 1L
)

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

throughput_all <- make_throughput_table(
  input            = "/Users/rd23aag/Downloads/nircam_throughputs/mean_throughputs/",
  filters          = nircam_filters,
  nircam_lambda_um = nircam_lambda_um,
  exclude          = c("F150W2", "150W2"),
  output           = "throughput_nircam.csv",
  format           = "csv",
  preview          = TRUE
)
