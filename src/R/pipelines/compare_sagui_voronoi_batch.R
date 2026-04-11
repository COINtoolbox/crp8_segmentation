suppressPackageStartupMessages({
  library(tools)
})

base_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation"
script_path <- file.path(base_dir, "src/R/pipelines/compare_sagui_voronoi_demo.R")

ncomp <- suppressWarnings(as.integer(Sys.getenv("SAGUI_COMPARE_NCOMP", unset = "8")))
run_suffix <- Sys.getenv("SAGUI_COMPARE_SUFFIX", unset = "papermask")
if (!is.finite(ncomp) || ncomp <= 1) {
  stop("SAGUI_COMPARE_NCOMP must be an integer greater than 1.")
}

fits_dir <- file.path(base_dir, "data/raw")
fits_files <- list.files(
  fits_dir,
  full.names = TRUE,
  pattern = "^datacube_sagui.*[.]fits$"
)

get_tag <- function(path) {
  sub("^datacube_(sagui.*)[.]fits$", "\\1", basename(path))
}

tags <- vapply(fits_files, get_tag, character(1))
summary_rows <- vector("list", length(tags))

for (i in seq_along(tags)) {
  tag <- tags[i]
  message("Running comparison for ", tag, " (N=", ncomp, ")")
  cmd <- sprintf(
    "SAGUI_COMPARE_TAG=%s SAGUI_COMPARE_NCOMP=%d SAGUI_COMPARE_SUFFIX=%s Rscript %s",
    shQuote(tag),
    ncomp,
    shQuote(run_suffix),
    shQuote(script_path)
  )
  status <- system(cmd)
  if (!identical(status, 0L)) {
    warning("Comparison failed for ", tag, " with exit status ", status)
    next
  }

  run_tag <- if (nzchar(run_suffix)) {
    sprintf("%s_n%02d_%s", tag, ncomp, gsub("[^A-Za-z0-9]+", "_", tolower(run_suffix)))
  } else {
    sprintf("%s_n%02d", tag, ncomp)
  }
  summary_path <- file.path(
    base_dir,
    "results/segmentation/sagui_compare",
    run_tag,
    "comparison_summary.csv"
  )
  if (file.exists(summary_path)) {
    dat <- utils::read.csv(summary_path, stringsAsFactors = FALSE)
    dat$tag <- tag
    summary_rows[[i]] <- dat
  }
}

batch_summary <- do.call(rbind, Filter(Negate(is.null), summary_rows))
if (!is.null(batch_summary)) {
  out_path <- file.path(
    base_dir,
    "results/segmentation/sagui_compare",
    if (nzchar(run_suffix)) sprintf("comparison_summary_all_n%02d_%s.csv", ncomp, gsub("[^A-Za-z0-9]+", "_", tolower(run_suffix))) else sprintf("comparison_summary_all_n%02d.csv", ncomp)
  )
  utils::write.csv(batch_summary, out_path, row.names = FALSE)
  message("Batch summary written to ", out_path)
}
