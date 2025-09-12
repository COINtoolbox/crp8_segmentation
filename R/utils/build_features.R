build_features <- function(P, q_fore = 0.85, scale_xy = 1, scale_I = 1) {
  H <- nrow(P); W <- ncol(P)
  v <- as.numeric(P); 
  v <- v[is.finite(v)]
  a <- stats::quantile(v, 0.01, na.rm=TRUE); 
  b <- stats::quantile(v, 0.99, na.rm=TRUE)
  P01 <- (pmin(pmax(P, a), b) - a) / max(b - a, .Machine$double.eps)
  Fg  <- P01 >= stats::quantile(P01, q_fore, na.rm=TRUE)
  idx <- which(Fg); if (!length(idx)) return(list(idx=integer(0), feat=NULL, P01=P01))
  y <- ((idx-1L) %% H) + 1L; x <- ((idx-1L) %/% H) + 1L
  feat <- cbind(scale_xy * (x-1)/(W-1), scale_xy * (y-1)/(H-1), scale_I * P01[idx])
  list(idx = idx, feat = feat, P01 = P01)
}
