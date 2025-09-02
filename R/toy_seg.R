# ============================================================
# Hyperspectral clump segmentation via PCA + DBSCAN
# - Reads a FITS cube
# - Builds a PCA-based detection image
# - Segments with DBSCAN (configurable)
# - Cleans labels and (optionally) removes border clusters
# - Masks original cube OR extracts tiny subcubes around clumps
# ============================================================

# ---- Dependencies ----
# install.packages(c("FITSio","dbscan","viridis"))  # viridis is optional
library(FITSio)
library(dbscan)

# ---- Utils ----

robust01 <- function(M, q = c(0.01, 0.99)) {
  v <- as.numeric(M); v <- v[is.finite(v)]
  a <- quantile(v, q[1], na.rm = TRUE); b <- quantile(v, q[2], na.rm = TRUE)
  (pmin(pmax(M, a), b) - a) / max(b - a, .Machine$double.eps)
}

xy_from_idx <- function(idx, H, W) {
  cbind(y = ((idx - 1L) %% H) + 1L,
        x = ((idx - 1L) %/% H) + 1L)
}

compact_labels <- function(L) {
  u <- sort(unique(L[L > 0]))
  if (!length(u)) return(L)
  re <- seq_along(u); names(re) <- u
  L[L > 0] <- as.integer(re[as.character(L[L > 0])])
  L
}

drop_small <- function(L, min_size = 50L) {
  if (!any(L > 0)) return(L)
  tb   <- table(L[L > 0])
  keep <- as.integer(names(tb)[tb >= min_size])
  L[!(L %in% keep)] <- 0L
  compact_labels(L)
}

remove_border_components <- function(L) {
  H <- nrow(L); W <- ncol(L)
  for (k in sort(unique(L[L > 0]))) {
    if (any(L[1,] == k) || any(L[H,] == k) || any(L[,1] == k) || any(L[,W] == k)) {
      L[L == k] <- 0L
    }
  }
  compact_labels(L)
}

# ---- PCA energy map ----
# df_mat: (nx*ny) x B matrix, each row = spectrum at a pixel
# d: number of PCs to keep in the energy map (sum of squares of first d scores)
pca_energy_map <- function(df_mat, nx, ny, d = 5) {
  pca <- prcomp(df_mat, center = TRUE, scale. = FALSE)
  d    <- min(d, ncol(pca$x))
  E    <- rowSums(pca$x[, 1:d, drop = FALSE]^2)
  matrix(E, nrow = nx, ncol = ny)
}


# Build foreground + feature matrix for DBSCAN/OPTICS
build_features <- function(P, q_fore = 0.90, scale_xy = 1, scale_I = 1) {
  H <- nrow(P); W <- ncol(P)
  # robust [0,1]
  v <- as.numeric(P); v <- v[is.finite(v)]
  a <- quantile(v, 0.01, na.rm=TRUE); b <- quantile(v, 0.99, na.rm=TRUE)
  P01 <- (pmin(pmax(P, a), b) - a) / max(b - a, .Machine$double.eps)
  
  Fg  <- P01 >= quantile(P01, q_fore, na.rm=TRUE)
  idx <- which(Fg); if (!length(idx)) return(list(idx=integer(0), feat=NULL, P01=P01))
  
  y <- ((idx-1) %% H) + 1; x <- ((idx-1) %/% H) + 1
  feat <- cbind(scale_xy * (x-1)/(W-1),
                scale_xy * (y-1)/(H-1),
                scale_I  *  P01[idx])
  list(idx = idx, feat = feat, P01 = P01)
}


# ---- DBSCAN segmentation on a detection image (e.g. PCA energy) ----
# Returns an HxW integer label map (0=background, 1..K = clumps)
segment_dbscan_image <- function(P,
                                 q_fore   = 0.90,    # keep top-q pixels to cluster
                                 scale_xy = 1.0,     # weight for spatial coords
                                 scale_I  = 1.0,     # weight for intensity
                                 eps      = 0.06,    # DBSCAN radius (in feature space)
                                 minPts   = 40) {    # min points per clump
  H <- nrow(P); W <- ncol(P)
  
  P01 <- robust01(P)                          # stabilize scale
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

segment_hdbscan <- function(P,
                            q_fore = 0.92, scale_xy = 1.0, scale_I = 1.5,
                            minPts = 25) {
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

segment_optics_xi <- function(P, q_fore=0.90, scale_xy=1, scale_I=1.2, minPts=15, xi=0.05) {
  H <- nrow(P); W <- ncol(P)
  bf <- build_features(P, q_fore, scale_xy, scale_I)  # from previous message
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


# Keep clusters within size range; relabel compactly
filter_by_size <- function(L, min_size = 30, max_size = Inf) {
  if (!any(L>0)) return(L)
  tb <- table(L[L>0])
  keep <- as.integer(names(tb)[tb >= min_size & tb <= max_size])
  L[!(L %in% keep)] <- 0L
  u <- sort(unique(L[L>0])); re <- seq_along(u); names(re) <- u
  if (length(u)) L[L>0] <- as.integer(re[as.character(L[L>0])])
  L
}

# ---- Apply binary/soft mask to cube ----
# mode = "zero" -> set background to 0
#      = "na"   -> set background to NA
#      = "soft" -> multiply by a soft (blurred) mask
mask_cube <- function(cube, labels, mode = c("zero","na","soft"), sigma = 2) {
  mode <- match.arg(mode)
  H <- dim(cube)[1]; W <- dim(cube)[2]; B <- dim(cube)[3]
  m <- (labels > 0)
  if (mode == "soft") {
    # simple Gaussian blur implemented via repeated box smoothing:
    # for documentation clarity we keep it simple; swap for EBImage::gblur if desired
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
    ms <- blur2d(m*1, k = k)
    ms <- ms / max(ms)
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

# ---- Extract tiny cubes around each compact region ----
# margin: pixels to include around bbox (halo)
# min_size: ignore very small components
extract_region_cubes <- function(cube, labels, margin = 6L, min_size = 50L) {
  H <- dim(cube)[1]; W <- dim(cube)[2]
  regs <- split(seq_len(H*W), as.vector(labels))
  out <- list()
  for (lab in names(regs)) {
    if (lab == "0") next
    idx <- regs[[lab]]
    if (length(idx) < min_size) next
    xy <- xy_from_idx(idx, H, W)
    y1 <- max(1, min(xy[,1]) - margin)
    y2 <- min(H, max(xy[,1]) + margin)
    x1 <- max(1, min(xy[,2]) - margin)
    x2 <- min(W, max(xy[,2]) + margin)
    out[[lab]] <- list(
      cube = cube[y1:y2, x1:x2, , drop = FALSE],
      bbox = c(y1, y2, x1, x2),
      npix = length(idx)
    )
  }
  out
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

# 1) Read cube & build PCA energy image
Xfits <- FITSio::readFITS("datacube_reg1.fits")
cube  <- Xfits$imDat                   # H x W x B
nx <- dim(cube)[1]; ny <- dim(cube)[2]

# Capivara helper gives (nx*ny) x B matrix; or reshape by yourself
df_mat <- capivara::cube_to_matrix(Xfits)
# Detection image from first d PCs
d <- 2
P <- pca_energy_map(df_mat, nx, ny, d = d)

# 2) Segment with DBSCAN (pick eps/minPts to your data scale)
# Tip: use dbscan::kNNdistplot(feat, k=minPts) beforehand to pick eps.
#L <- segment_dbscan_image(
#  P,
#  q_fore = 0.90,     # fewer points → faster + less noise
#  scale_xy = 1.0,
#  scale_I  = 1.0,
#  eps    = 0.05,     # set from kNN elbow
#  minPts = 10        # ~ min clump size after foreground mask
#)

L <- segment_hdbscan(P, q_fore=0.9, scale_I=1, minPts=10)
L <- filter_by_size(L, min_size = 50)
image(L)


# 3) Clean labels (drop tiny, optional remove border touchers)
#L <- drop_small(L, min_size = 50)
# L <- remove_border_components(L)

# 4a) Mask whole cube (choose one)
#cube_zero <- mask_cube(cube, L, mode = "zero")   # background -> 0
#cube_na   <- mask_cube(cube, L, mode = "na")   # background -> NA
# cube_soft <- mask_cube(cube, L, mode = "soft", sigma = 2)

# 4b) OR extract tiny cubes around each clump
regions <- extract_region_cubes(cube, L, margin = 8, min_size = 10)
# write_region_cubes_fits(regions, prefix = "cutout")

# 5) Quick visualization (optional)
if (requireNamespace("viridis", quietly = TRUE)) {
  par(mfrow = c(1,2))
  image(t(robust01(P)[nx:1, ]), col = gray.colors(256), asp = 1,
        main = "PCA detection image"); contour(t(L[nx:1, ]), add = TRUE, drawlabels = FALSE, col = "red")
  image(log10(cube_zero[,,min(4, dim(cube)[3])]), col = viridis::magma(128),
        asp = 1, main = "Masked cube (log, band 4)")
  par(mfrow = c(1,1))
}

