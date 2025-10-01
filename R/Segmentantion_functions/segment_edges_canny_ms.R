# ------------------------------------------------------------------
# Multi-scale Canny-like edge grouping with Scharr/Sobel, NMS, hysteresis
# Pure base R; uses your conv2() (zero-padded) exactly.
# Returns HxW INTEGER labels: 0=background; 1..K=edge groups (8-connected)
# ------------------------------------------------------------------
segment_edges_canny_ms <- function(
    P,
    sigmas     = c(1.0, 2.0, 3.0),  # multi-scale Gaussian sigmas
    use_scharr = TRUE,              # Scharr is more rotation-accurate than Sobel
    q_high     = 0.92,              # high quantile for hysteresis (on NMS magnitudes >0)
    q_low      = 0.60,              # low quantile (weak edges)
    scale_norm = TRUE,              # scale-normalised gradients: |∇(Gσ * P)| × σ
    drop_frame = TRUE,              # drop the artificial frame from zero-padding
    conn       = 8L,                # 8-connected edge groups
    min_size   = 0L                 # filter tiny edge groups (0 = keep all)
) {
  if (!is.matrix(P)) P <- as.matrix(P)
  P <- matrix(as.numeric(P), nrow(P), ncol(P))
  H <- nrow(P); W <- ncol(P)
  
  # --- kernels ---
  get_kernels <- function(scharr = TRUE) {
    if (scharr) {
      # 3x3 Scharr (un-normalised; scaling cancels in quantiles)
      Gx <- matrix(c(-3, 0, 3,
                     -10,0,10,
                     -3, 0, 3), 3, 3, byrow = TRUE)
      Gy <- matrix(c(-3,-10,-3,
                     0,  0,  0,
                     3, 10,  3), 3, 3, byrow = TRUE)
    } else {
      # Sobel
      Gx <- matrix(c(-1,0,1,-2,0,2,-1,0,1), 3, 3, byrow = TRUE)
      Gy <- matrix(c(-1,-2,-1,0,0,0,1,2,1), 3, 3, byrow = TRUE)
    }
    list(Gx=Gx, Gy=Gy)
  }
  GK <- get_kernels(use_scharr)
  
  gauss_kernel <- function(sigma) {
    if (sigma <= 0) return(matrix(1,1,1))
    sz <- max(3L, as.integer(2*ceiling(3*sigma)+1L))
    ax <- seq(-(sz%/%2), sz%/%2)
    g  <- exp(-0.5*(ax/sigma)^2); g <- g/sum(g)
    outer(g, g)
  }
  
  # Non-maximum suppression (4 orientation bins: 0°,45°,90°,135°)
  nms_thin <- function(M, Gx, Gy) {
    ang <- atan2(Gy, Gx) * 180/pi
    ang[ang < 0] <- ang[ang < 0] + 180
    out <- matrix(0.0, H, W)
    for (i in 2:(H-1)) {
      for (j in 2:(W-1)) {
        a <- ang[i,j]; m <- M[i,j]
        if (m <= 0) next
        if ((a < 22.5) || (a >= 157.5)) {      # 0° (horizontal)
          if (m >= M[i, j-1] && m >= M[i, j+1]) out[i,j] <- m
        } else if (a < 67.5) {                 # 45°
          if (m >= M[i-1, j+1] && m >= M[i+1, j-1]) out[i,j] <- m
        } else if (a < 112.5) {                # 90° (vertical)
          if (m >= M[i-1, j] && m >= M[i+1, j]) out[i,j] <- m
        } else {                               # 135°
          if (m >= M[i-1, j-1] && m >= M[i+1, j+1]) out[i,j] <- m
        }
      }
    }
    out
  }
  
  # Hysteresis: strong edges + weak edges connected to strong (8-conn)
  hysteresis <- function(N) {
    vals <- N[N > 0]
    if (!length(vals)) return(matrix(FALSE, H, W))
    hi <- as.numeric(quantile(vals, q_high, na.rm=TRUE))
    lo <- as.numeric(quantile(vals, q_low,  na.rm=TRUE))
    strong <- N >= hi
    weak   <- (N >= lo) & !strong
    if (drop_frame) {
      strong[c(1,H), ] <- FALSE; strong[, c(1,W)] <- FALSE
      weak[c(1,H), ]   <- FALSE; weak[, c(1,W)]   <- FALSE
    }
    if (!any(strong)) return(matrix(FALSE, H, W))
    E <- matrix(FALSE, H, W)
    q_i <- which(strong, arr.ind=TRUE)[,1]
    q_j <- which(strong, arr.ind=TRUE)[,2]
    head <- 1L
    while (head <= length(q_i)) {
      i <- q_i[head]; j <- q_j[head]; head <- head + 1L
      if (!E[i,j]) {
        E[i,j] <- TRUE
        for (di in -1:1) for (dj in -1:1) if (!(di==0 && dj==0)) {
          ii <- i + di; jj <- j + dj
          if (ii>=1L && ii<=H && jj>=1L && jj<=W) {
            if ((weak[ii,jj] || strong[ii,jj]) && !E[ii,jj]) {
              q_i <- c(q_i, ii); q_j <- c(q_j, jj)
            }
          }
        }
      }
    }
    E
  }
  
  # 8-connected component labelling
  cc_labels_conn <- function(B) {
    lab <- matrix(0L, H, W); cur <- 0L
    for (i in 1:H) for (j in 1:W) if (B[i,j] && lab[i,j]==0L) {
      cur <- cur + 1L
      q_i <- i; q_j <- j; head <- 1L; lab[i,j] <- cur
      while (head <= length(q_i)) {
        ii <- q_i[head]; jj <- q_j[head]; head <- head + 1L
        for (di in -1:1) for (dj in -1:1) if (!(di==0 && dj==0)) {
          iii <- ii + di; jjj <- jj + dj
          if (iii>=1L && iii<=H && jjj>=1L && jjj<=W) {
            if (B[iii,jjj] && lab[iii,jjj]==0L) {
              lab[iii,jjj] <- cur
              q_i <- c(q_i, iii); q_j <- c(q_j, jjj)
            }
          }
        }
      }
    }
    storage.mode(lab) <- "integer"
    lab
  }
  
  E_all <- matrix(FALSE, H, W)
  for (sigma in sigmas) {
    K  <- gauss_kernel(sigma)
    Pb <- conv2(P, K)                 # Gaussian smoothing
    Gx <- conv2(Pb, GK$Gx)
    Gy <- conv2(Pb, GK$Gy)
    M  <- sqrt(Gx^2 + Gy^2)
    if (scale_norm) M <- sigma * M    # scale-normalised derivative
    N  <- nms_thin(M, Gx, Gy)         # non-maximum suppression
    E  <- hysteresis(N)               # strong + connected weak
    E_all <- E_all | E
  }
  
  if (drop_frame) {
    E_all[c(1,H), ] <- FALSE; E_all[, c(1,W)] <- FALSE
  }
  if (!any(E_all)) return(matrix(0L, H, W))
  
  L <- cc_labels_conn(E_all)
  
  if (min_size > 0L) {
    ids <- setdiff(unique(as.integer(L)), 0L)
    keep <- ids[vapply(ids, function(id) sum(L==id), integer(1)) >= min_size]
    L[!(L %in% keep)] <- 0L
    if (length(keep)) {
      map <- setNames(seq_along(keep), keep)
      sel <- L > 0L
      L[sel] <- as.integer(map[as.character(L[sel])])
    }
  }
  L
}
