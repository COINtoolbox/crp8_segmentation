suppressPackageStartupMessages({
  library(tools)
})

base_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation"
script_path <- file.path(base_dir, "src/R/pipelines/compare_sagui_voronoi_demo.R")
tag <- "sagui10"
groups <- c(10L, 15L, 20L, 30L)
suffix <- "paperexact"
summary_rows <- vector("list", length(groups))

for (i in seq_along(groups)) {
  ncomp <- groups[i]
  message("Running ", tag, " comparison for N=", ncomp)
  cmd <- sprintf(
    "env SAGUI_COMPARE_TAG=%s SAGUI_COMPARE_NCOMP=%d SAGUI_COMPARE_SUFFIX=%s Rscript %s",
    shQuote(tag),
    ncomp,
    shQuote(suffix),
    shQuote(script_path)
  )
  status <- system(cmd)
  if (!identical(status, 0L)) {
    warning("Comparison failed for ", tag, " N=", ncomp, " with exit status ", status)
    next
  }

  run_tag <- sprintf("%s_n%02d_%s", tag, ncomp, suffix)
  summary_path <- file.path(base_dir, "results/segmentation/sagui_compare", run_tag, "comparison_summary.csv")
  if (file.exists(summary_path)) {
    dat <- utils::read.csv(summary_path, stringsAsFactors = FALSE)
    dat$tag <- tag
    dat$ncomp_requested <- ncomp
    summary_rows[[i]] <- dat
  }
}

batch_summary <- do.call(rbind, Filter(Negate(is.null), summary_rows))
if (!is.null(batch_summary)) {
  out_path <- file.path(base_dir, "results/segmentation/sagui_compare", "sagui10_group_sweep_paperexact.csv")
  utils::write.csv(batch_summary, out_path, row.names = FALSE)
  message("Group sweep summary written to ", out_path)
}
