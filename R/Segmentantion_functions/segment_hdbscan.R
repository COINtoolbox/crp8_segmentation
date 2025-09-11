segment_hdbscan <- function(P, q_fore=0.92, scale_xy=1.0, scale_I=1.5, minPts=25) {
  H <- nrow(P); W <- ncol(P)
  bf <- build_features(P, q_fore=q_fore, scale_xy=scale_xy, scale_I=scale_I)
  idx <- bf$idx; feat <- bf$feat; if (!length(idx)) return(matrix(0L,H,W))
  hd <- dbscan::hdbscan(feat, minPts = minPts)
  lab <- hd$cluster
  L <- matrix(0L, H, W)
  if (any(lab > 0)) {
    u <- sort(unique(lab[lab>0])); re <- seq_along(u); names(re) <- u
    L[idx] <- as.integer(re[as.character(lab)])
  }
  L
}