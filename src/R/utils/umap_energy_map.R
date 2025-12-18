umap_energy_map <- function(df_mat, H, W,
                            d = 3, n_neighbors = 30, min_dist = 0.05,
                            metric = "euclidean",
                            xy_weight = 0.25,      # spatial nudging (0–0.4 is typical)
                            center = TRUE, scale. = TRUE,
                            energy = c("z2","mahalanobis"),
                            ridge = 1e-6, seed = 42) {
  energy <- match.arg(energy)
  stopifnot(length(H) == 1, length(W) == 1, H > 0, W > 0, nrow(df_mat) == H*W)
  if (!requireNamespace("uwot", quietly = TRUE)) stop("Please install 'uwot'.")
  
  X <- as.matrix(df_mat); X[!is.finite(X)] <- 0
  
  # robust scaling (no NaNs from zero-variance cols)
  if (center || scale.) {
    mu  <- colMeans(X, na.rm = TRUE); X <- sweep(X, 2, mu, "-")
    if (scale.) {
      sdv <- apply(X, 2, stats::sd, na.rm = TRUE); sdv[!is.finite(sdv) | sdv == 0] <- 1
      X   <- sweep(X, 2, sdv, "/")
    }
  }
  
  # ---- SAFE spatial features (no scale() on constant columns!) ----
  if (xy_weight > 0) {
    xc <- rep(seq_len(W), each = H); yc <- rep(seq_len(H), times = W)
    sscale <- function(v){ s <- stats::sd(v); if (!is.finite(s) || s == 0) s <- 1; (v - mean(v)) / s }
    XY <- cbind(sscale(xc), sscale(yc)) * xy_weight
    X  <- cbind(X, XY)
  }
  X[!is.finite(X)] <- 0  # belt & suspenders
  
  set.seed(seed)
  emb <- uwot::umap(
    X,
    n_neighbors = n_neighbors,
    n_components = d,
    min_dist = min_dist,
    metric = metric,
    scale = FALSE, init = "spectral", verbose = FALSE
  )
  
  # collapse to scalar map
  if (energy == "z2") {
    mu <- colMeans(emb); Z <- sweep(emb, 2, mu, "-")
    sdv <- apply(Z, 2, stats::sd); sdv[!is.finite(sdv) | sdv == 0] <- 1
    Z <- sweep(Z, 2, sdv, "/"); E <- rowSums(Z^2)
  } else {
    mu <- colMeans(emb); S <- stats::cov(emb)
    dE <- ncol(emb); S <- S + diag(ridge * mean(diag(S), na.rm = TRUE), dE)
    E <- stats::mahalanobis(emb, center = mu, cov = S)
  }
  
  # robust [0,1]
  v <- E[is.finite(E)]; q <- stats::quantile(v, c(0.01, 0.99), na.rm = TRUE)
  P <- (pmin(pmax(E, q[1]), q[2]) - q[1]) / max(q[2]-q[1], .Machine$double.eps)
  matrix(P, nrow = H, ncol = W)
}

