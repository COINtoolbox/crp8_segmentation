plot_sed_with_spec <- function(reg_long,
                               region_id,
                               spec,
                               scale_mode = c("auto", "linear", "log"),
                               band_height_frac = 0.90) {
  
  scale_mode <- match.arg(scale_mode)
  
  # ---- 1) Subset one region (column is 'region') ----
  if (!"region" %in% names(reg_long)) {
    stop("reg_long must have a 'region' column. You have: ",
         paste(names(reg_long), collapse = ", "))
  }
  
  dat <- reg_long[reg_long$region == region_id, , drop = FALSE]
  if (!nrow(dat)) stop("No rows found for region = ", region_id)
  
  # Need 'band' and 'flux'
  needed <- c("band", "flux")
  miss   <- setdiff(needed, names(dat))
  if (length(miss) > 0) {
    stop("reg_long is missing columns: ", paste(miss, collapse = ", "))
  }
  
  # ---- 2) Map band indices (0,1,2,...) -> lambda_um using spec ----
  # band is stored as character "0","1",...; convert to numeric
  band_idx <- suppressWarnings(as.numeric(dat$band))
  if (any(!is.finite(band_idx))) {
    stop("Cannot convert 'band' to numeric. Values are: ",
         paste(unique(dat$band), collapse = ", "))
  }
  
  # From spec: get band_index and lambda_um for each filter
  idx_vec <- vapply(spec$spec, function(s) {
    if (is.null(s$band_index)) NA_real_ else s$band_index
  }, numeric(1))
  
  lam_vec <- vapply(spec$spec, function(s) s$lambda_um, numeric(1))
  fil_vec <- vapply(spec$spec, function(s) s$filter,    character(1))
  
  if (any(!is.finite(idx_vec))) {
    stop("spec$spec elements must contain 'band_index'.")
  }
  
  lut_lambda <- setNames(lam_vec, idx_vec)
  lut_filter <- setNames(fil_vec, idx_vec)
  
  dat$lambda_um <- lut_lambda[as.character(band_idx)]
  dat$filter    <- lut_filter[as.character(band_idx)]
  
  if (any(!is.finite(dat$lambda_um))) {
    warning("Some bands could not be mapped to lambda_um via band_index.")
  }
  
  # ---- 3) Build band coverage rectangles from throughput curves ----
  band_df_list <- lapply(spec$spec, function(s) {
    thr <- s$throughput
    w   <- s$wave
    
    if (all(!is.finite(thr))) return(NULL)
    
    mthr <- max(thr, na.rm = TRUE)
    if (!is.finite(mthr) || mthr <= 0) return(NULL)
    
    sel <- thr > 0.05 * mthr  # 5% of peak throughput
    if (!any(sel)) sel <- rep(TRUE, length(thr))
    
    data.frame(
      filter    = s$filter,
      wave_min  = min(w[sel], na.rm = TRUE),
      wave_max  = max(w[sel], na.rm = TRUE),
      lambda_um = s$lambda_um
    )
  })
  
  band_df <- do.call(rbind, band_df_list)
  
  # keep only filters that actually appear in this region’s data
  if (!is.null(band_df) && nrow(band_df) > 0) {
    band_df <- band_df[band_df$filter %in% dat$filter, , drop = FALSE]
  } else {
    band_df <- NULL
  }
  
  # ---- 4) y-scale ----
  y <- dat$flux
  y_max <- max(y, na.rm = TRUE)
  if (!is.finite(y_max) || y_max <= 0) y_max <- 1
  
  # ---- 5) Build ggplot ----
  library(ggplot2)
  
  p <- ggplot(dat, aes(x = lambda_um, y = flux)) +
    theme_minimal(base_size = 12) +
    ggtitle(paste("Region", region_id, "- SED with NIRCam bands")) +
    xlab(expression(lambda ~ "[µm]")) +
    ylab("Flux")
  
  if (!is.null(band_df) && nrow(band_df) > 0) {
    p <- p +
      geom_rect(
        data = band_df,
        aes(xmin = wave_min, xmax = wave_max,
            ymin = 0,
            ymax = band_height_frac * y_max),
        inherit.aes = FALSE,
        alpha = 0.15
      )
  }
  
  p <- p +
    geom_line() +
    geom_point(size = 2)
  
  if (scale_mode == "log") {
    p <- p + scale_y_log10()
  }
  
  p
}

