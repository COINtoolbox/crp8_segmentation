# --- tiny helpers (base R morphology via conv2 with ones) ---
box_kernel <- function(k) { if (k %% 2 == 0) k <- k + 1L; matrix(1, k, k) }
dilate_bin <- function(B, k = 3L) (conv2(B, box_kernel(k)) > 0)
erode_bin  <- function(B, k = 3L) { K <- box_kernel(k); conv2(B, K) >= sum(K) }
close_bin  <- function(B, k = 3L, iters = 1L) {
  if (iters <= 0L) return(B)
  for (t in seq_len(iters)) B <- erode_bin(dilate_bin(B, k), k)
  B
}

# 8-connected component labeller for a logical mask
cc_labels_conn <- function(B, conn = 8L) {
  H <- nrow(B); W <- ncol(B)
  lab <- matrix(0L, H, W); cur <- 0L
  nbs <- if (conn == 8L) {
    rbind(c(-1,0), c(1,0), c(0,-1), c(0,1), c(-1,-1), c(-1,1), c(1,-1), c(1,1))
  } else {
    rbind(c(-1,0), c(1,0), c(0,-1), c(0,1))
  }
  for (i in 1:H) for (j in 1:W) if (B[i,j] && lab[i,j]==0L) {
    cur <- cur + 1L
    q_i <- i; q_j <- j; head <- 1L; lab[i,j] <- cur
    while (head <= length(q_i)) {
      ii <- q_i[head]; jj <- q_j[head]; head <- head + 1L
      for (t in 1:nrow(nbs)) {
        iii <- ii + nbs[t,1]; jjj <- jj + nbs[t,2]
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

# --- main: group border/edge pixels ---
segment_sobel <- function(P,
                                 edge_q      = 0.90,  # quantile on |grad|
                                 drop_frame  = TRUE,  # remove image-frame edge
                                 close_k     = 3L,    # closing kernel (odd)
                                 close_iters = 1L,    # 0 = skip closing
                                 conn        = 8L,    # 8-connected curves
                                 min_size    = 0L) {  # filter tiny edge groups
  if (!is.matrix(P)) P <- as.matrix(P)
  P <- matrix(as.numeric(P), nrow(P), ncol(P))
  H <- nrow(P); W <- ncol(P)
  
  # Sobel (exactly as in gpa)
  Gxk <- matrix(c(-1,0,1,-2,0,2,-1,0,1), 3, 3, byrow=TRUE)
  Gyk <- matrix(c(-1,-2,-1,0,0,0,1,2,1), 3, 3, byrow=TRUE)
  Gx  <- conv2(P, Gxk)
  Gy  <- conv2(P, Gyk)
  mag <- sqrt(Gx^2 + Gy^2)
  
  # threshold edges
  thr <- as.numeric(quantile(mag, edge_q, na.rm=TRUE))
  E <- (mag >= thr) & !is.na(mag)
  
  # kill the artificial frame edge from zero-padding
  if (drop_frame) {
    E[1, ] <- FALSE; E[H, ] <- FALSE; E[, 1] <- FALSE; E[, W] <- FALSE
  }
  
  # bridge small gaps so contours stay connected
  if (close_iters > 0L) E <- close_bin(E, k = close_k, iters = close_iters)
  
  if (!any(E)) return(matrix(0L, H, W))
  
  # label edge groups
  L <- cc_labels_conn(E, conn = conn)
  
  # optional: drop tiny groups and reindex to 1..K
  if (min_size > 0L) {
    ids <- setdiff(unique(as.integer(L)), 0L)
    if (length(ids)) {
      keep <- ids[vapply(ids, function(id) sum(L==id), integer(1)) >= min_size]
      L[!(L %in% keep)] <- 0L
      if (length(keep)) {
        map <- setNames(seq_along(keep), keep)
        sel <- L > 0L
        L[sel] <- as.integer(map[as.character(L[sel])])
      }
    }
  }
  storage.mode(L) <- "integer"
  L
}

