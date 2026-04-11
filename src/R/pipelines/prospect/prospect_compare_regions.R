#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

base_dir <- "/Users/rd23aag/Documents/GitHub/crp8_segmentation"

prospect_path <- file.path(
  base_dir,
  "results/sed_fitting/prospect/sagui10/raw/sagui10_prospect_results.csv"
)

current_path <- file.path(
  base_dir,
  "results/sed_fitting/sagui-10_sed_results_new.csv"
)

out_dir <- file.path(base_dir, "results/sed_fitting/prospect/sagui10/raw/comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

prospect_df <- readr::read_csv(prospect_path, show_col_types = FALSE) |>
  transmute(
    region,
    logM_prospect = .data$logMformed,
    logZ_prospect = .data$logzsol,
    Av_prospect = .data$Av,
    age_mw_prospect = .data$age_mw,
    sfr_prospect = .data$sfr
  )

current_df <- readr::read_csv(current_path, show_col_types = FALSE) |>
  transmute(
    region,
    logM_current = log10(.data$mass),
    logZ_current = .data$logzsol,
    Av_current = 1.086 * .data$dust2,
    age_mw_current = .data$age_mass_weighted,
    sfr_current = .data$sfr_100myr
  )

cmp <- left_join(current_df, prospect_df, by = "region")

property_specs <- list(
  list(cur = "logM_current", pro = "logM_prospect", label = "log M"),
  list(cur = "logZ_current", pro = "logZ_prospect", label = "log Z/Zsun"),
  list(cur = "Av_current", pro = "Av_prospect", label = "A_V"),
  list(cur = "age_mw_current", pro = "age_mw_prospect", label = "Age_MW"),
  list(cur = "sfr_current", pro = "sfr_prospect", label = "SFR")
)

summary_df <- bind_rows(lapply(property_specs, function(spec) {
  x <- cmp[[spec$cur]]
  y <- cmp[[spec$pro]]
  tibble::tibble(
    property = spec$label,
    pearson = suppressWarnings(cor(x, y, use = "pairwise.complete.obs", method = "pearson")),
    spearman = suppressWarnings(cor(x, y, use = "pairwise.complete.obs", method = "spearman")),
    median_offset = median(y - x, na.rm = TRUE),
    mad_offset = mad(y - x, na.rm = TRUE)
  )
}))

readr::write_csv(cmp, file.path(out_dir, "region_comparison.csv"))
readr::write_csv(summary_df, file.path(out_dir, "summary_metrics.csv"))

make_scatter <- function(spec) {
  ggplot(cmp, aes(x = .data[[spec$cur]], y = .data[[spec$pro]])) +
    geom_abline(slope = 1, intercept = 0, colour = "#bdbdbd", linewidth = 0.5, linetype = "dashed") +
    geom_point(size = 2.2, colour = "#245b8a", alpha = 0.85) +
    labs(x = paste("Current", spec$label), y = paste("ProSpect", spec$label)) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#ececec"),
      plot.margin = margin(6, 6, 6, 6)
    )
}

plots <- lapply(property_specs, make_scatter)
panel <- wrap_plots(plots, ncol = 2)

png_path <- file.path(out_dir, "prospect_vs_current_scatter.png")
png(png_path, width = 1800, height = 2200, res = 220)
print(panel)
dev.off()

cat("Wrote comparison outputs to:\n", out_dir, "\n")
