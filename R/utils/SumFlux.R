#' @import FITSio
#' @importFrom dplyr group_by summarise ungroup mutate select tibble arrange
#' @importFrom tidyr pivot_wider
#' @export
RegionPhotometry <- function(cluster_result,
                             bkg = NULL,
                             var_cube = NULL,
                             sigma_band = NULL,
                             digits_lambda_colnames = 6,
                             return_painted_cube = FALSE) {
  # ---- Extrai peças ----
  cubedat     <- cluster_result$original_cube
  M           <- cubedat$imDat
  lab         <- cluster_result$cluster_map
  
  if (is.null(M) || length(dim(M)) != 3L) stop("original_cube$imDat deve ser (nx,ny,n_bandas).")
  if (is.null(lab) || !is.matrix(lab))     stop("cluster_map deve ser uma matriz (nx,ny).")
  
  dims <- dim(M); nx <- dims[1]; ny <- dims[2]; nb <- dims[3]
  if (!all(dim(lab)[1:2] == c(nx, ny))) stop("cluster_map não bate com as dimensões espaciais do cubo.")
  
  # Eixo espectral/bandas
  bands <- FITSio::axVec(3, cubedat$axDat)
  band_names <- formatC(bands, digits = digits_lambda_colnames, format = "fg")
  
  # ---- Subtração de fundo (opcional) ----
  if (!is.null(bkg)) {
    # normaliza bkg ao shape do cubo
    if (is.array(bkg) && length(dim(bkg)) == 3L) {
      if (!all(dim(bkg) == dim(M))) stop("bkg cubo deve ter as mesmas dimensões que o cubo de fluxos.")
      M <- M - bkg
    } else if (is.matrix(bkg) && all(dim(bkg) == c(nx, ny))) {
      # mesmo mapa para todas as bandas
      for (j in seq_len(nb)) M[, , j] <- M[, , j] - bkg
    } else if (length(bkg) == nb) {
      # offset escalar por banda
      for (j in seq_len(nb)) M[, , j] <- M[, , j] - bkg[j]
    } else if (length(bkg) == 1L) {
      M <- M - as.numeric(bkg)
    } else {
      stop("Formato de 'bkg' não reconhecido (use escalar, vetor por banda, matriz nx×ny, ou cubo nx×ny×nb).")
    }
  }
  
  # ---- Variância / incerteza por pixel (opcional) ----
  use_var <- FALSE
  if (!is.null(var_cube)) {
    if (!is.array(var_cube) || !all(dim(var_cube) == dim(M))) {
      stop("var_cube deve ser um array (nx,ny,nb) igual ao cubo de fluxos.")
    }
    use_var <- TRUE
  } else if (!is.null(sigma_band)) {
    if (length(sigma_band) != nb) stop("sigma_band deve ter length = n_bandas.")
    # constroi cubo de variâncias a partir de sigma_band
    var_cube <- array(0, dim = dim(M))
    for (j in seq_len(nb)) var_cube[, , j] <- sigma_band[j]^2
    use_var <- TRUE
  }
  
  # ---- Lineariza para matrizes (n_spaxel × nb) ----
  X   <- cube_to_matrix(list(imDat = M))           # fluxos
  cls <- as.vector(lab)                            # rótulos
  if (nrow(X) != length(cls)) stop("cube_to_matrix() não compatível com as dimensões espaciais.")
  
  # conta pixels por região
  # regiões NA serão ignoradas
  valid <- !is.na(cls)
  cls_valid <- cls[valid]
  X_valid   <- X[valid, , drop = FALSE]
  
  # ---- Soma por região/banda ----
  # fluxo integrado = soma dos fluxos dos pixels da região
  # erro do fluxo integrado = sqrt(soma das variâncias) se var disponível
  # (assumindo pixels independentes; se não, isso é um lower bound)
  # Preparar índice por região
  region_ids <- sort(unique(cls_valid))
  region_index <- match(cls_valid, region_ids)  # 1..n_regions
  
  # Soma eficiente via rowsum por banda
  flux_mat <- matrix(0, nrow = length(region_ids), ncol = nb)
  n_pix    <- as.numeric(rowsum(rep(1, length(cls_valid)), group = region_index))
  
  for (j in seq_len(nb)) {
    flux_mat[, j] <- rowsum(X_valid[, j], group = region_index, reorder = FALSE)
  }
  
  # Incertezas (opcional)
  flux_err_mat <- NULL
  if (use_var) {
    V <- cube_to_matrix(list(imDat = var_cube))[valid, , drop = FALSE]
    flux_var_mat <- matrix(0, nrow = length(region_ids), ncol = nb)
    for (j in seq_len(nb)) {
      flux_var_mat[, j] <- rowsum(V[, j], group = region_index, reorder = FALSE)
    }
    flux_err_mat <- sqrt(pmax(flux_var_mat, 0))
  }
  
  # ---- Tabelas de saída ----
  # wide
  flux_wide <- data.frame(region = region_ids, n_pix = n_pix, check.names = FALSE)
  colnames(flux_mat) <- band_names
  flux_wide <- cbind(flux_wide, as.data.frame(flux_mat, check.names = FALSE))
  if (!is.null(flux_err_mat)) {
    colnames(flux_err_mat) <- paste0(band_names, "_err")
    flux_wide <- cbind(flux_wide, as.data.frame(flux_err_mat, check.names = FALSE))
  }
  
  # long
  to_long <- function(M, nm) {
    df <- as.data.frame(M, check.names = FALSE)
    df$region <- region_ids
    tidyr::pivot_longer(df, -region, names_to = "band", values_to = nm)
  }
  long_flux <- to_long(flux_mat, "flux")
  if (!is.null(flux_err_mat)) {
    long_err  <- to_long(flux_err_mat, "flux_err")
    flux_long <- dplyr::left_join(long_flux, long_err, by = c("region", "band"))
  } else {
    flux_long <- long_flux
    flux_long$flux_err <- NA_real_
  }
  # adiciona n_pix
  flux_long <- dplyr::left_join(
    flux_long,
    dplyr::tibble(region = region_ids, n_pix = n_pix),
    by = "region"
  ) %>% dplyr::arrange(region, band)
  
  # ---- Cubo pintado (opcional) ----
  painted_cube <- NULL
  if (return_painted_cube) {
    # Para cada banda j, cada pixel recebe o FLUXO INTEGRADO da sua região (constante por região)
    painted_cube <- array(NA_real_, dim = dim(M))
    # dicionário: region -> linha no flux_mat
    reg_row <- match(cls, region_ids)  # inclui NAs
    for (j in seq_len(nb)) {
      layer <- rep(NA_real_, length(cls))
      ok <- !is.na(reg_row)
      layer[ok] <- flux_mat[reg_row[ok], j]
      painted_cube[, , j] <- matrix(layer, nrow = nx, ncol = ny)
    }
    dimnames(painted_cube) <- list(NULL, NULL, band_names)
  }
  
  # ---- Saída ----
  list(
    flux_long        = flux_long,
    flux_wide        = flux_wide,
    painted_cube     = painted_cube,
    bands            = bands,
    band_names       = band_names
  )
}