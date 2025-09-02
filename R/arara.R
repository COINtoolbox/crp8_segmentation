# ============================================================
# Hyperspectral clump segmentation via PCA + Density Clustering
# (DBSCAN / HDBSCAN / OPTICS), with:
#  - PCA detection image
#  - Segmentation (configurable)
#  - Size filtering
#  - Hole filling (per label)
#  - Masking the original cube
#  - Coordinate-preserving mini-cube extraction
# ============================================================

# ---- Dependencies ----
# install.packages(c("FITSio","dbscan","viridis"))
library(FITSio)
library(dbscan)

# ---- Utils ----

# Robust [0,1] scaling with quantile clipping
robust01 <- function(M, q = c(0.01, 0.99)) {
  v <- as.numeric(M); v <- v[is.finite(v)]
  a <- quantile(v, q[1], na.rm = TRUE); b <- quantile(v, q[2], na.rm = TRUE)
  (pmin(pmax(M, a), b) - a) / max(b - a, .Machine$double.eps)
}

# Linear index -> (y,x)
xy_from_idx <- function(idx, H, W) {
  cbind(y = ((idx - 1L) %% H) + 1L,
        x = ((idx - 1L) %/% H) + 1L)
}

# Relabel 1..K (keep 0 as background)
compact_labels <- function(L) {
  u <- sort(unique(L[L > 0]))
  if (!length(u)) return(L)
  re <- seq_along(u); names(re) <- u
  L[L > 0] <- as.integer(re[as.character(L[L > 0])])
  L
}

# Keep only components with size in [min_size, max_size]
filter_by_size <- function(L, min_size = 30L, max_size = Inf) {
  if (!any(L>0)) return(L)
  tb <- table(L[L>0])
  keep <- as.integer(names(tb)[tb >= min_size & tb <= max_size])
  L[!(L %in% keep)] <- 0L
  compact_labels(L)
}

# Remove components touching the image border
remove_border_components <- function(L) {
  H <- nrow(L); W <- ncol(L)
  for (k in sort(unique(L[L > 0]))) {
    if (any(L[1,] == k) || any(L[H,] == k) || any(L[,1] == k) || any(L[,W] == k)) {
      L[L == k] <- 0L
    }
  }
  compact_labels(L)
}

# ---- PCA detection image ----
# df_mat: (H*W) x B matrix of spectra (rows=pixels)
# d: number of PCs whose energy to sum into the detection image
pca_energy_map <- function(df_mat, H, W, d = 5) {
  pca <- prcomp(df_mat, center = TRUE, scale. = FALSE)
  d <- min(d, ncol(pca$x))
  E <- rowSums(pca$x[, 1:d, drop = FALSE]^2)
  matrix(E, nrow = H, ncol = W)
}

# ---- Build features for density clustering ----
build_features <- function(P, q_fore = 0.90, scale_xy = 1, scale_I = 1) {
  H <- nrow(P); W <- ncol(P)
  v <- as.numeric(P); v <- v[is.finite(v)]
  a <- quantile(v, 0.01, na.rm=TRUE); b <- quantile(v, 0.99, na.rm=TRUE)
  P01 <- (pmin(pmax(P, a), b) - a) / max(b - a, .Machine$double.eps)
  
  Fg  <- P01 >= quantile(P01, q_fore, na.rm=TRUE)
  idx <- which(Fg); if (!length(idx)) return(list(idx=integer(0), feat=NULL, P01=P01))
  
  y <- ((idx-1L) %% H) + 1L
  x <- ((idx-1L) %/% H) + 1L
  feat <- cbind(
    scale_xy * (x-1)/(W-1),
    scale_xy * (y-1)/(H-1),
    scale_I  *  P01[idx]
  )
  list(idx = idx, feat = feat, P01 = P01)
}

# ---- Segmentation: DBSCAN on detection image ----
segment_dbscan_image <- function(P,
                                 q_fore   = 0.90,
                                 scale_xy = 1.0,
                                 scale_I  = 1.0,
                                 eps      = 0.06,
                                 minPts   = 40) {
  H <- nrow(P); W <- ncol(P)
  P01 <- robust01(P)
  Fg  <- P01 >= quantile(P01, q_fore, na.rm = TRUE)
  idx <- which(Fg); if (!length(idx)) return(matrix(0L, H, W))
  
  yy <- ((idx - 1L) %% H) + 1L
  xx <- ((idx - 1L) %/% H) + 1L
  
  feat <- cbind(
    scale_xy * (xx - 1)/(W - 1),
    scale_xy * (yy - 1)/(H - 1),
    scale_I  *  P01[idx]
  )
  
  fit <- dbscan(feat, eps = eps, minPts = minPts)
  L   <- matrix(0L, H, W)
  if (any(fit$cluster > 0)) {
    u  <- sort(unique(fit$cluster[fit$cluster > 0]))
    re <- seq_along(u); names(re) <- u
    L[idx] <- as.integer(re[as.character(fit$cluster)])
  }
  L
}

# ---- Segmentation: HDBSCAN (parameter-light) ----
segment_hdbscan <- function(P, q_fore = 0.92, scale_xy = 1.0, scale_I = 1.5, minPts = 25) {
  H <- nrow(P); W <- ncol(P)
  bf <- build_features(P, q_fore=q_fore, scale_xy=scale_xy, scale_I=scale_I)
  idx <- bf$idx; feat <- bf$feat; if (!length(idx)) return(matrix(0L,H,W))
  hd  <- dbscan::hdbscan(feat, minPts = minPts)
  lab <- hd$cluster
  L <- matrix(0L, H, W)
  if (any(lab > 0)) {
    u <- sort(unique(lab[lab>0])); re <- seq_along(u); names(re) <- u
    L[idx] <- as.integer(re[as.character(lab)])
  }
  L
}

# ---- Segmentation: OPTICS + Xi (optional) ----
segment_optics_xi <- function(P, q_fore=0.90, scale_xy=1, scale_I=1.2, minPts=15, xi=0.05) {
  H <- nrow(P); W <- ncol(P)
  bf <- build_features(P, q_fore, scale_xy, scale_I)
  idx <- bf$idx; feat <- bf$feat; if (!length(idx)) return(matrix(0L,H,W))
  opt <- dbscan::optics(feat, minPts = minPts)
  lab <- dbscan::extractXi(opt, xi = xi)$cluster
  L <- matrix(0L, H, W)
  if (any(lab>0)) {
    u <- sort(unique(lab[lab>0])); re <- seq_along(u); names(re) <- u
    L[idx] <- as.integer(re[as.character(lab)])
  }
  L
}

# ---- Fill internal holes for each labeled region (base R) ----
fill_holes_per_label <- function(L) {
  H <- nrow(L); W <- ncol(L); out <- L
  for (k in sort(unique(L[L > 0]))) {
    idx <- which(L == k); if (!length(idx)) next
    yy <- ((idx - 1L) %% H) + 1L; xx <- ((idx - 1L) %/% H) + 1L
    y1 <- min(yy); y2 <- max(yy); x1 <- min(xx); x2 <- max(xx)
    sub <- out[y1:y2, x1:x2]
    bg  <- (sub == 0L)
    
    # mark border-connected background via flood fill (4-neighb.)
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

# ---- Apply binary/soft mask to cube ----
# mode = "zero" -> background set to 0
#      = "na"   -> background set to NA
#      = "soft" -> multiply by a blurred mask
mask_cube <- function(cube, labels, mode = c("zero","na","soft"), sigma = 2) {
  mode <- match.arg(mode)
  H <- dim(cube)[1]; W <- dim(cube)[2]; B <- dim(cube)[3]
  m <- (labels > 0)
  if (mode == "soft") {
    # lightweight box-blur (odd kernel)
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
    out <- cube
    for (b in seq_len(B)) out[,,b] <- out[,,b] * ms
    return(out)
  } else {
    out <- cube
    for (b in seq_len(B)) {
      Xb <- out[,,b]
      if (mode == "zero") Xb[!m] <- 0 else Xb[!m] <- NA_real_
      out[,,b] <- Xb
    }
    return(out)
  }
}

# ---- Extract tiny cubes around each compact region (preserves coords) ----
# Returns list with: id, bbox=(y1,y2,x1,x2), cube, mask (optional), npix
extract_region_cubes <- function(cube, labels, margin = 6L, min_size = 50L,
                                 return_mask = TRUE) {
  H <- dim(cube)[1]; W <- dim(cube)[2]
  ids <- sort(unique(labels[labels > 0])); if (!length(ids)) return(list())
  out <- vector("list", length(ids)); names(out) <- ids
  j <- 0L
  for (id in ids) {
    idx <- which(labels == id); if (length(idx) < min_size) next
    yx <- xy_from_idx(idx, H, W); y <- yx[,1]; x <- yx[,2]
    y1 <- max(1, min(y) - margin); y2 <- min(H, max(y) + margin)
    x1 <- max(1, min(x) - margin); x2 <- min(W, max(x) + margin)
    j <- j + 1L
    out[[j]] <- list(
      id   = id,
      bbox = c(y1, y2, x1, x2),
      cube = cube[y1:y2, x1:x2, , drop = FALSE],
      mask = if (return_mask) (labels[y1:y2, x1:x2] == id) else NULL,
      npix = length(idx)
    )
  }
  out[seq_len(j)]
}

# ---- Optional: save subcubes to FITS ----
write_region_cubes_fits <- function(regions, prefix = "cutout") {
  i <- 0L
  for (nm in names(regions)) {
    i <- i + 1L
    fn <- sprintf("%s_%03d.fits", prefix, i)
    FITSio::writeFITSim(fn, regions[[nm]]$cube)
  }
}

# ============================================================
# Example (adjust paths & params)
# ============================================================

# 1) Read cube & build PCA detection image
Xfits <- FITSio::readFITS("datacube_reg1.fits")
cube  <- Xfits$imDat                  # H x W x B
H <- dim(cube)[1]; W <- dim(cube)[2]

# Capivara helper gives (H*W) x B matrix; or reshape yourself if needed
df_mat <- if (requireNamespace("capivara", quietly = TRUE)) {
  capivara::cube_to_matrix(Xfits)
} else {
  B <- dim(cube)[3]; M <- matrix(NA_real_, H*W, B)
  for (b in seq_len(B)) M[,b] <- as.vector(cube[,,b]); M
}

# Detection image from first d PCs
d <- 2
P <- pca_energy_map(df_mat, H, W, d = d)

# 2) Segment (choose one)
# L <- segment_dbscan_image(P, q_fore=0.90, scale_xy=1.0, scale_I=1.0, eps=0.05, minPts=10)
L <- segment_hdbscan(P, q_fore=0.90, scale_xy=1.0, scale_I=2.0, minPts=50)
# L <- segment_optics_xi(P, q_fore=0.90, scale_xy=1.0, scale_I=1.2, minPts=15, xi=0.05)

# 3) Clean labels
L <- filter_by_size(L, min_size = 25)
L <- fill_holes_per_label(L)          # ← fills ringy interiors
# L <- remove_border_components(L)    # optional

# Quick look
image(L, main = "Labels")

# 4a) Mask the cube (avoid NaNs when plotting log)
cube_zero <- mask_cube(cube, L, mode = "zero")
cube_na   <- mask_cube(cube, L, mode = "na")

image(log(cube_na[,,1]),
        col = viridis::plasma(100))


# 4b) Extract mini-cubes
regions <- extract_region_cubes(cube, L, margin = 6, min_size = 50, return_mask = TRUE)
# write_region_cubes_fits(regions, prefix = "cutout")

