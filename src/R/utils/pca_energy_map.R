pca_energy_map <- function(df_mat, H, W, d = 5) {
  pca <- stats::prcomp(df_mat, center = TRUE, scale. = FALSE)
  d <- min(d, ncol(pca$x))
  E <- rowSums(pca$x[, 1:d, drop = FALSE]^2)
  matrix(E, nrow = H, ncol = W)
}