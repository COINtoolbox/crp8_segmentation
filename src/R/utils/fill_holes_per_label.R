fill_holes_per_label <- function(L) {
  H <- nrow(L); W <- ncol(L); out <- L
  for (k in sort(unique(L[L > 0]))) {
    idx <- which(L == k); if (!length(idx)) next
    yy <- ((idx - 1L) %% H) + 1L; xx <- ((idx - 1L) %/% H) + 1L
    y1 <- min(yy); y2 <- max(yy); x1 <- min(xx); x2 <- max(xx)
    sub <- out[y1:y2, x1:x2]; bg <- (sub == 0L)
    h <- nrow(sub); w <- ncol(sub)
    mark <- matrix(FALSE, h, w)
    q <- integer(h*w); head <- 1L; tail <- 0L
    push <- function(y,x){
      if (y>=1 && y<=h && x>=1 && x<=w && bg[y,x] && !mark[y,x]) {
        mark[y,x] <<- TRUE; tail <<- tail + 1L; q[tail] <<- (x-1L)*h + y
      }
    }
    for (x in 1:w) { if (bg[1,x]) push(1,x); if (bg[h,x]) push(h,x) }
    for (y in 1:h) { if (bg[y,1]) push(y,1); if (bg[y,w]) push(y,w) }
    while (head <= tail) {
      p <- q[head]; head <- head + 1L
      y <- ((p-1L) %% h) + 1L; x <- ((p-1L) %/% h) + 1L
      if (y>1) push(y-1,x); if (y<h) push(y+1,x)
      if (x>1) push(y,x-1); if (x<w) push(y,x+1)
    }
    holes <- bg & !mark
    sub[holes] <- k
    out[y1:y2, x1:x2] <- sub
  }
  out
}