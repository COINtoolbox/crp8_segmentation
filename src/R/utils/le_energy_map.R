# Compress SEDs to a scalar anomaly score (HxW), robust to outliers.
sed_mahalanobis_map <- function(df_mat, H, W,
                                log_flux = TRUE,         # log1p to stabilize dynamic range
                                norm = c("l1","l2"),     # color-normalize per pixel
                                use_mcd = TRUE,          # robust covariance (MCD) if available
                                ridge = 1e-6) {          # tiny diagonal jitter if needed
  norm <- match.arg(norm)
  stopifnot(nrow(df_mat) == H*W)
  X <- as.matrix(df_mat)
  X[!is.finite(X)] <- 0
  if (log_flux) X <- log1p(pmax(X, 0))
  
  # per-pixel color normalization (focus on SED shape, not brightness)
  rs <- if (norm == "l1") rowSums(X) else sqrt(rowSums(X^2))
  rs[rs == 0 | !is.finite(rs)] <- 1
  Xn <- X / rs
  
  # robust center & covariance
  get_cov <- function(M) {
    if (use_mcd && requireNamespace("rrcov", quietly = TRUE)) {
      cfit <- rrcov::CovMcd(M); list(center = as.numeric(cfit@center), cov = cfit@cov)
    } else if (requireNamespace("MASS", quietly = TRUE)) {
      cfit <- MASS::cov.rob(M); list(center = as.numeric(cfit$center), cov = cfit$cov)
    } else {
      list(center = colMeans(M), cov = stats::cov(M))
    }
  }
  cfit <- get_cov(Xn)
  # stabilize covariance (in case of near-singularity)
  S <- cfit$cov
  d <- ncol(Xn)
  S <- S + diag(ridge * mean(diag(S), na.rm = TRUE), d)
  D2 <- stats::mahalanobis(Xn, center = cfit$center, cov = S)   # squared distances
  D  <- sqrt(pmax(D2, 0))
  # robust [0,1] map
  q <- stats::quantile(D, probs = c(0.01, 0.99), na.rm = TRUE)
  P <- (pmin(pmax(D, q[1]), q[2]) - q[1]) / max(q[2]-q[1], .Machine$double.eps)
  matrix(P, H, W)
}
