copula_energy <- function(cube, eps=1e-6, positive_tail=TRUE, whiten=FALSE) {
  stopifnot(length(dim(cube))==3)
  H <- dim(cube)[1]; W <- dim(cube)[2]; K <- dim(cube)[3]
  M <- matrix(cube, H*W, K)
  # empirical CDF -> uniform -> normal scores
  to_uniform <- function(v){
    r <- rank(v, ties.method="average", na.last="keep")
    n <- sum(is.finite(v)); u <- r/(n+1)
    pmin(pmax(u, eps), 1-eps)
  }
  U <- apply(M, 2, to_uniform)
  Z <- qnorm(U)
  if (positive_tail) Z[Z < 0] <- 0
  if (!whiten) {
    Pvec <- rowSums(Z^2)
  } else {
    # optional whitening across bands (copula-Mahalanobis)
    S <- stats::cov(Z, use="pairwise.complete.obs")
    Sinv <- tryCatch(solve(S), error=function(e) MASS::ginv(S))
    Pvec <- rowSums((Z %*% Sinv) * Z)
  }
  matrix(Pvec, H, W)
}
