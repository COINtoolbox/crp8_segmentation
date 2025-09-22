# Segment via the Sobel + quantile with your zero-padded conv2
# Returns HxW integer labels: 0 = background; 1..K = kept components
# mode = "interior" -> fill regions enclosed by edges (galaxy mask)
# mode = "edge"     -> label the edge pixels themselves
segment_sobel <- function(P,
                              edge_q   = 0.90,        # exactly like gpa(edge_q)
                              mode     = c("interior","edge"),
                              keep     = c("largest","all"),
                              min_size = 300L) {
  mode <- match.arg(mode)
  keep <- match.arg(keep)
  if (!is.matrix(P)) P <- as.matrix(P)
  P <- matrix(as.numeric(P), nrow(P), ncol(P))
  H <- nrow(P); W <- ncol(P)
  
  # --- Sobel identical to gpa() ---
  Gxk <- matrix(c(-1,0,1,-2,0,2,-1,0,1), 3, 3, byrow=TRUE)
  Gyk <- matrix(c(-1,-2,-1,0,0,0,1,2,1), 3, 3, byrow=TRUE)
  
  # IMPORTANT: use *your* conv2 (zero-padded), not a different one
  Gx <- conv2(P, Gxk)
  Gy <- conv2(P, Gyk)
  mag <- sqrt(Gx^2 + Gy^2)
  
  # Threshold exactly as in gpa()
  thr <- as.numeric(quantile(mag, edge_q, na.rm = TRUE))
  E   <- mag >= thr
  E[is.na(E)] <- FALSE
  
  # Choose what to label
  if (mode == "edge") {
    FG <- E  # label the edge pixels themselves
  } else {
    # "interior" (galaxy): regions enclosed by edges = non-edges NOT connected to border
    flood_from_borders <- function(passable) {
      passable[is.na(passable)] <- FALSE
      vis <- matrix(FALSE, H, W)
      q_i <- integer(0); q_j <- integer(0)
      enqueue <- function(ii, jj) {
        ok <- passable[cbind(ii, jj)] & !vis[cbind(ii, jj)]
        if (any(ok)) {
          ii <- ii[ok]; jj <- jj[ok]
          vis[cbind(ii, jj)] <<- TRUE
          q_i <<- c(q_i, ii); q_j <<- c(q_j, jj)
        }
      }
      # seed all borders
      enqueue(rep(1L, W), 1L:W)     # top
      enqueue(rep(H, W), 1L:W)     # bottom
      enqueue(1L:H, rep(1L, H))    # left
      enqueue(1L:H, rep(W, H))     # right
      
      head <- 1L
      while (head <= length(q_i)) {
        i <- q_i[head]; j <- q_j[head]; head <- head + 1L
        if (i > 1L) enqueue(i-1L, j)
        if (i < H)  enqueue(i+1L, j)
        if (j > 1L) enqueue(i, j-1L)
        if (j < W)  enqueue(i, j+1L)
      }
      vis
    }
    non_edge <- !E
    bg_vis   <- flood_from_borders(non_edge)
    FG       <- non_edge & !bg_vis  # enclosed by edges
  }
  
  if (!any(FG)) return(matrix(0L, H, W))
  
  # Connected components (4-connected), then keep largest or union >= min_size
  cc_labels <- function(B) {
    lab <- matrix(0L, H, W); cur <- 0L
    for (i in 1:H) for (j in 1:W) if (B[i,j] && lab[i,j]==0L) {
      cur <- cur + 1L
      q_i <- i; q_j <- j; head <- 1L; lab[i,j] <- cur
      while (head <= length(q_i)) {
        ii <- q_i[head]; jj <- q_j[head]; head <- head + 1L
        if (ii > 1L && B[ii-1L, jj] && lab[ii-1L, jj]==0L) {
          lab[ii-1L, jj] <- cur; q_i <- c(q_i, ii-1L); q_j <- c(q_j, jj)
        }
        if (ii < H  && B[ii+1L, jj] && lab[ii+1L, jj]==0L) {
          lab[ii+1L, jj] <- cur; q_i <- c(q_i, ii+1L); q_j <- c(q_j, jj)
        }
        if (jj > 1L && B[ii, jj-1L] && lab[ii, jj-1L]==0L) {
          lab[ii, jj-1L] <- cur; q_i <- c(q_i, ii); q_j <- c(q_j, jj-1L)
        }
        if (jj < W  && B[ii, jj+1L] && lab[ii, jj+1L]==0L) {
          lab[ii, jj+1L] <- cur; q_i <- c(q_i, ii); q_j <- c(q_j, jj+1L)
        }
      }
    }
    lab
  }
  
  labs <- cc_labels(FG)
  ids  <- setdiff(unique(as.integer(labs)), 0L)
  if (!length(ids)) return(matrix(0L, H, W))
  
  sizes <- vapply(ids, function(id) sum(labs == id), integer(1))
  keep_ids <- if (keep == "largest") ids[which.max(sizes)] else ids[sizes >= min_size]
  if (!length(keep_ids)) keep_ids <- ids[which.max(sizes)]
  
  # Reindex kept ids to 1..K; others → 0
  remap <- setNames(seq_along(keep_ids), keep_ids)
  L <- matrix(0L, H, W)
  sel <- labs %in% keep_ids
  L[sel] <- as.integer(remap[as.character(labs[sel])])
  L
}
