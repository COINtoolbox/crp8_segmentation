mask_cube <- function(cube, labels, mode = c("zero","na","soft"), sigma = 2) {
  mode <- match.arg(mode)
  H <- dim(cube)[1]; W <- dim(cube)[2]; B <- dim(cube)[3]
  m <- (labels > 0)
  if (mode == "soft") {
    blur2d <- function(M, k = 5L) {
      r <- k %/% 2; out <- matrix(0, H, W)
      S <- matrix(0, H+1, W+1); S[2:(H+1),2:(W+1)] <- apply(apply(M,2,cumsum),1,cumsum)
      for (y in 1:H) for (x in 1:W) {
        y1 <- max(1,y-r); y2 <- min(H,y+r); x1 <- max(1,x-r); x2 <- min(W,x+r)
        out[y,x] <- (S[y2+1,x2+1]-S[y1,x2+1]-S[y2+1,x1]+S[y1,x1]) / ((y2-y1+1)*(x2-x1+1))
      }
      out
    }
    k <- max(3L, 2L*round(sigma) + 1L)
    ms <- blur2d(m*1, k = k); ms <- ms / max(ms)
    out <- cube; for (b in seq_len(B)) out[,,b] <- out[,,b] * ms; return(out)
  } else {
    out <- cube
    for (b in seq_len(B)) { Xb <- out[,,b]; if (mode == "zero") Xb[!m] <- 0 else Xb[!m] <- NA_real_; out[,,b] <- Xb }
    return(out)
  }
}