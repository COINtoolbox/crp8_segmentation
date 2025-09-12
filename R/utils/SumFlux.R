#' @import FITSio
#' @importFrom dplyr group_by summarise ungroup mutate select tibble arrange
#' @importFrom tidyr pivot_wider
#' @export
# install.packages(c("dplyr","tidyr"))
# install.packages(c("dplyr","tidyr"))  # if needed
library(dplyr)
library(tidyr)

RegionPhotometry <- function(
    cube, labels,
    bkg = NULL,
    var_cube = NULL,
    sigma_band = NULL,
    band_values = NULL,
    digits_lambda_colnames = 6,
    return_painted_cube = FALSE
) {
  # --- helpers ---
  get_imdat <- function(x) if (is.list(x) && !is.null(x$imDat)) x$imDat else x
  to_mat <- function(A3) {  # [nx,ny,nb] -> [nx*ny, nb]
    nx <- dim(A3)[1]; ny <- dim(A3)[2]; nb <- dim(A3)[3]
    out <- matrix(NA_real_, nrow = nx*ny, ncol = nb)
    for (j in seq_len(nb)) out[, j] <- as.vector(A3[,,j])
    out
  }
  
  # --- inputs & bands ---
  M <- get_imdat(cube)
  stopifnot(is.array(M), length(dim(M)) == 3L)
  nx <- dim(M)[1]; ny <- dim(M)[2]; nb <- dim(M)[3]
  stopifnot(is.matrix(labels), all(dim(labels)[1:2] == c(nx,ny)))
  
  if (!is.null(band_values)) {
    bands <- band_values
  } else if (is.list(cube) && !is.null(cube$axDat)) {
    bands <- FITSio::axVec(3, cube$axDat)
  } else bands <- seq_len(nb)
  stopifnot(length(bands) == nb)
  band_names <- if (is.numeric(bands)) formatC(bands, digits = digits_lambda_colnames, format = "fg")
  else as.character(bands)
  band_names <- make.unique(band_names, sep = "_")
  
  # --- background (deterministic) ---
  if (!is.null(bkg)) {
    if (is.list(bkg) && !is.null(bkg$imDat)) bkg <- bkg$imDat
    if (is.array(bkg) && length(dim(bkg)) == 3L) {
      stopifnot(all(dim(bkg) == dim(M))); M <- M - bkg
    } else if (is.matrix(bkg) && all(dim(bkg) == c(nx,ny))) {
      for (j in seq_len(nb)) M[,,j] <- M[,,j] - bkg
    } else if (length(bkg) == nb) {
      for (j in seq_len(nb)) M[,,j] <- M[,,j] - as.numeric(bkg[j])
    } else if (length(bkg) == 1L) {
      M <- M - as.numeric(bkg)
    } else stop("Unsupported 'bkg' shape.")
  }
  
  # --- variance source ---
  use_var <- FALSE
  if (!is.null(var_cube)) {
    Vc <- get_imdat(var_cube); stopifnot(is.array(Vc), all(dim(Vc) == dim(M)))
    use_var <- TRUE
  } else if (!is.null(sigma_band)) {
    stopifnot(length(sigma_band) == nb)
    Vc <- array(0, dim = dim(M))
    for (j in seq_len(nb)) Vc[,,j] <- sigma_band[j]^2
    use_var <- TRUE
  } else Vc <- NULL
  
  # --- flatten to per-pixel rows ---
  X   <- to_mat(M)                # [n_pix_total, nb]
  cls <- as.vector(labels)
  valid <- is.finite(cls) & cls > 0
  Xv   <- X[valid, , drop = FALSE]
  clsv <- cls[valid]
  pix  <- seq_len(nrow(Xv))
  
  flux_df <- as_tibble(Xv, .name_repair = "minimal")
  names(flux_df) <- band_names
  flux_long <- flux_df %>%
    mutate(pix = pix) %>%
    pivot_longer(-pix, names_to = "band", values_to = "flux")
  
  if (use_var) {
    V   <- to_mat(Vc)[valid, , drop = FALSE]
    var_df <- as_tibble(V, .name_repair = "minimal")
    names(var_df) <- band_names
    var_long <- var_df %>%
      mutate(pix = pix) %>%
      pivot_longer(-pix, names_to = "band", values_to = "var")
    df <- left_join(flux_long, var_long, by = c("pix","band"))
  } else {
    df <- flux_long %>% mutate(var = NA_real_)
  }
  
  df <- df %>% mutate(region = clsv[pix])
  
  # --- integrate per region × band (single-row per group) ---
  df_sum <- df %>%
    group_by(region, band) %>%
    summarise(
      n_eff   = if (use_var)
        sum(is.finite(flux) & is.finite(var))
      else
        sum(is.finite(flux)),
      flux    = if (use_var)
        sum(ifelse(is.finite(flux) & is.finite(var), flux, 0), na.rm = TRUE)
      else
        sum(ifelse(is.finite(flux), flux, 0), na.rm = TRUE),
      var_sum = if (use_var)
        sum(ifelse(is.finite(flux) & is.finite(var), var, 0), na.rm = TRUE)
      else
        NA_real_,
      .groups = "drop"
    ) %>%
    mutate(flux_err = if (use_var) sqrt(pmax(var_sum, 0)) else NA_real_) %>%
    select(region, band, flux, flux_err, n_eff)
  
  # region-level raw pixel count (independent of band)
  n_pix_region <- tibble(region = clsv) %>%
    count(region, name = "n_pix")
  
  flux_long <- df_sum %>%
    left_join(n_pix_region, by = "region") %>%
    arrange(region, band)
  
  # --- WIDE tables (numeric only; no list-cols) ---
  flux_wide <- flux_long %>%
    select(region, n_pix, band, flux) %>%
    pivot_wider(
      names_from  = band,
      values_from = flux,
      values_fill = 0,           # numeric fill
      values_fn   = sum          # in case duplicates ever sneak in
    ) %>%
    arrange(region)
  
  neff_wide <- flux_long %>%
    select(region, band, n_eff) %>%
    pivot_wider(
      names_from  = band,
      values_from = n_eff,
      values_fill = 0,
      values_fn   = sum
    ) %>%
    rename_with(~ paste0(.x, "_n_eff"), -region)
  
  if (use_var) {
    err_wide <- flux_long %>%
      select(region, band, flux_err) %>%
      pivot_wider(
        names_from  = band,
        values_from = flux_err,
        values_fill = NA_real_,
        values_fn   = sum
      ) %>%
      rename_with(~ paste0(.x, "_err"), -region)
    
    flux_wide <- flux_wide %>%
      left_join(neff_wide, by = "region") %>%
      left_join(err_wide, by = "region")
  } else {
    flux_wide <- flux_wide %>%
      left_join(neff_wide, by = "region")
  }
  
  # --- painted cube (optional) ---
  painted_cube <- NULL
  if (return_painted_cube) {
    painted_cube <- array(NA_real_, dim = dim(M))
    # lookup region×band -> flux
    look <- flux_long %>% select(region, band, flux)
    for (j in seq_len(nb)) {
      bname <- band_names[j]
      lut <- setNames(look$flux[look$band == bname], look$region[look$band == bname])
      reg_vec <- as.vector(labels)
      layer <- ifelse(is.finite(reg_vec) & reg_vec > 0, lut[as.character(reg_vec)], NA_real_)
      painted_cube[,,j] <- matrix(layer, nrow = nx, ncol = ny)
    }
    dimnames(painted_cube) <- list(NULL, NULL, band_names)
  }
  
  list(
    flux_long    = flux_long,
    flux_wide    = flux_wide,
    painted_cube = painted_cube,
    bands        = bands,
    band_names   = band_names
  )
}
