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
    return_painted_cube = FALSE,
    error_fallback = c("none","flux_over_sqrt_n","poisson","mad_sky")
) {
  error_fallback <- match.arg(error_fallback)
  
  # --- helpers ---
  get_imdat <- function(x) if (is.list(x) && !is.null(x$imDat)) x$imDat else x
  to_mat <- function(A3) {  # [nx,ny,nb] -> [n, nb]
    nx <- dim(A3)[1]; ny <- dim(A3)[2]; nb <- dim(A3)[3]
    out <- matrix(NA_real_, nrow = nx*ny, ncol = nb)
    for (j in seq_len(nb)) out[, j] <- as.vector(A3[,,j])
    out
  }
  
  # --- inputs & bands ---
  M <- get_imdat(cube)
  stopifnot(is.array(M), length(dim(M)) == 3L)
  nx <- dim(M)[1]; ny <- dim(M)[2]; nb <- dim(M)[3]
  
  # normalize labels to [nx,ny]
  if (is.array(labels) && length(dim(labels)) == 3L && dim(labels)[3] == 1L) {
    labels <- labels[,,1]
  }
  stopifnot(is.matrix(labels), all(dim(labels)[1:2] == c(nx, ny)))
  
  if (!is.null(band_values)) {
    bands <- band_values
  } else if (is.list(cube) && !is.null(cube$axDat)) {
    bands <- FITSio::axVec(3, cube$axDat)
  } else {
    bands <- seq_len(nb)
  }
  stopifnot(length(bands) == nb)
  band_names <- if (is.numeric(bands)) formatC(bands, digits = digits_lambda_colnames, format = "fg") else as.character(bands)
  band_names <- make.unique(band_names, sep = "_")
  
  # --- background (deterministic) ---
  if (!is.null(bkg)) {
    bkgM <- get_imdat(bkg)
    if (is.array(bkgM) && length(dim(bkgM)) == 3L) {
      stopifnot(all(dim(bkgM) == dim(M))); M <- M - bkgM
    } else if (is.matrix(bkgM) && all(dim(bkgM) == c(nx, ny))) {
      for (j in seq_len(nb)) M[,,j] <- M[,,j] - bkgM
    } else if (length(bkgM) == nb) {
      for (j in seq_len(nb)) M[,,j] <- M[,,j] - as.numeric(bkgM[j])
    } else if (length(bkgM) == 1L) {
      M <- M - as.numeric(bkgM)
    } else {
      stop("Unsupported 'bkg' shape.")
    }
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
  } else {
    Vc <- NULL
  }
  
  # --- flatten ---
  X   <- to_mat(M)                 # [n_pix_total, nb]
  cls <- as.vector(labels)         # [n_pix_total]
  valid <- is.finite(cls) & (cls > 0)
  
  if (!any(valid)) {
    warning("No valid (labels > 0 & finite) pixels found. Returning empty tables.")
    empty <- tibble(region = integer(), band = character(), flux = numeric(),
                    flux_err = numeric(), n_eff = integer(), n_pix = integer())
    return(list(
      flux_long    = empty,
      flux_wide    = tibble(),
      painted_cube = if (isTRUE(return_painted_cube)) array(NA_real_, dim = dim(M)) else NULL,
      bands        = bands,
      band_names   = band_names
    ))
  }
  
  Xv   <- X[valid, , drop = FALSE]
  clsv <- cls[valid]
  pix  <- seq_len(nrow(Xv))
  
  flux_df <- as_tibble(Xv, .name_repair = "minimal"); names(flux_df) <- band_names
  flux_long0 <- flux_df %>% mutate(pix = pix) %>% pivot_longer(-pix, names_to = "band", values_to = "flux")
  
  if (use_var) {
    V   <- to_mat(Vc)[valid, , drop = FALSE]
    var_df <- as_tibble(V, .name_repair = "minimal"); names(var_df) <- band_names
    var_long <- var_df %>% mutate(pix = pix) %>% pivot_longer(-pix, names_to = "band", values_to = "var")
    df <- left_join(flux_long0, var_long, by = c("pix","band"))
  } else {
    df <- flux_long0 %>% mutate(var = NA_real_)
  }
  
  # attach region + ok flag
  df2 <- df %>%
    mutate(
      region = clsv[pix],
      ok = if (use_var) (is.finite(flux) & is.finite(var)) else is.finite(flux)
    )
  
  # --- aggregate per region × band ---
  df_sum <- df2 %>%
    group_by(region, band, .drop = FALSE) %>%
    summarise(
      n_eff   = sum(ok),
      flux    = sum(ifelse(ok, flux, 0), na.rm = TRUE),
      var_sum = if (use_var) sum(ifelse(ok, var, 0), na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
  
  # --- flux_err ---
  if (use_var) {
    df_sum <- df_sum %>%
      mutate(flux_err = sqrt(pmax(var_sum, 0)))
  } else {
    # choose fallback
    if (error_fallback == "flux_over_sqrt_n") {
      df_sum <- df_sum %>% mutate(flux_err = abs(flux) / sqrt(pmax(n_eff, 1)))
    } else if (error_fallback == "poisson") {
      df_sum <- df_sum %>% mutate(flux_err = sqrt(pmax(flux, 0)))
    } else if (error_fallback == "mad_sky") {
      # estimate per-band sigma from sky pixels (labels <= 0)
      sky <- !is.finite(cls) | (cls <= 0)
      if (!any(sky)) {
        warning("No sky pixels (labels <= 0). Falling back to flux_over_sqrt_n.")
        df_sum <- df_sum %>% mutate(flux_err = abs(flux) / sqrt(pmax(n_eff, 1)))
      } else {
        Xsky <- X[sky, , drop = FALSE]
        # robust per-band sigma
        sky_sigma <- apply(Xsky, 2L, function(v) 1.4826 * stats::mad(v, center = stats::median(v, na.rm = TRUE), na.rm = TRUE))
        names(sky_sigma) <- band_names
        df_sum <- df_sum %>%
          mutate(
            flux_err = sky_sigma[band] * sqrt(pmax(n_eff, 1))
          )
      }
    } else {
      df_sum <- df_sum %>% mutate(flux_err = NA_real_)  # "none"
    }
  }
  
  df_sum <- df_sum %>% select(region, band, flux, flux_err, n_eff)
  
  # region pixel counts
  n_pix_region <- tibble(region = clsv) %>% count(region, name = "n_pix")
  
  flux_long <- df_sum %>%
    left_join(n_pix_region, by = "region") %>%
    arrange(region, band)
  
  # --- WIDE tables ---
  flux_wide <- flux_long %>%
    select(region, n_pix, band, flux) %>%
    pivot_wider(names_from = band, values_from = flux, values_fn = dplyr::first) %>%
    arrange(region)
  
  err_wide <- flux_long %>%
    select(region, band, flux_err) %>%
    pivot_wider(names_from = band, values_from = flux_err, values_fn = dplyr::first) %>%
    rename_with(~ paste0(.x, "_err"), -region)
  
  neff_wide <- flux_long %>%
    select(region, band, n_eff) %>%
    pivot_wider(names_from = band, values_from = n_eff, values_fn = dplyr::first) %>%
    rename_with(~ paste0(.x, "_n_eff"), -region)
  
  flux_wide <- flux_wide %>%
    left_join(neff_wide, by = "region") %>%
    left_join(err_wide, by = "region")
  
  # --- painted cube (optional)
  painted_cube <- NULL
  if (isTRUE(return_painted_cube)) {
    painted_cube <- array(NA_real_, dim = dim(M))
    look <- flux_long %>% select(region, band, flux)
    reg_vec <- as.vector(labels)
    for (j in seq_len(nb)) {
      bname <- band_names[j]
      lut <- setNames(look$flux[look$band == bname], look$region[look$band == bname])
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
