# Hyperspectral clump segmentation via PCA + OPTICS‑Xi (with core→halo growth & hole fill)
# ---------------------------------------------------------------------------
# This script provides a clean, modular pipeline:
#   1) Build a PCA detection image from a FITS cube.
#   2) Segment compact cores with OPTICS‑Xi (dbscan pkg; no eps tuning).
#   3) Grow each core over a lower threshold (hysteresis) to include halos/rings.
#   4) Fill internal holes per region.
#   5) (Optional) Filter by size, mask cube, and extract padded sub‑cubes.
#
# All steps are small, documented functions so you can swap pieces easily.

# ---- Dependencies ----
# install.packages(c("FITSio","dbscan","viridis"))  # viridis optional for plotting
library(FITSio)
library(dbscan)

# ---- Utils ----

# Robust [0,1] scaling with quantile clipping
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

# Drop tiny/huge components and relabel
filter_by_size <- function(L, min_size = 30L, max_size = Inf) {
  if (!any(L > 0)) return(L)
  tb <- table(L[L > 0])
  keep <- as.integer(names(tb)[tb >= min_size & tb <= max_size])
  L[!(L %in% keep)] <- 0L
  compact_labels(L)
}

remove_border_components <- function(L) {
  H <- nrow(L); W <- ncol(L)
  for (k in sort(unique(L[L > 0]))) {
    if (any(L[1,] == k) || any(L[H,] == k) || any(L[,1] == k) || any(L[,W] == k)) L[L == k] <- 0L
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

# ---- Feature builder for density clustering ----
# Returns pixel indices used (idx), feature matrix (feat), and the scaled image (P01)
build_features <- function(P, q_fore = 0.92, scale_xy = 1, scale_I = 1) {
  H <- nrow(P); W <- ncol(P)
  P01 <- robust01(P)
  Fg  <- P01 >= quantile(P01, q_fore, na.rm = TRUE)
  idx <- which(Fg); if (!length(idx)) return(list(idx = integer(0), feat = NULL, P01 = P01))
  y <- ((idx - 1L) %% H) + 1L; x <- ((idx - 1L) %/% H) + 1L
  feat <- cbind(scale_xy * (x - 1)/(W - 1),
                scale_xy * (y - 1)/(H - 1),
                scale_I  *  P01[idx])
  list(idx = idx, feat = feat, P01 = P01)
}

# ---- OPTICS‑Xi segmentation (cores) ----
# minPts controls minimum local density; xi controls valley depth in reachability
segment_optics_xi <- function(P, q_fore = 0.92, scale_xy = 1.0, scale_I = 1.5,
                              minPts = 20, xi = 0.05) {
  H <- nrow(P); W <- ncol(P)
  bf <- build_features(P, q_fore = q_fore, scale_xy = scale_xy, scale_I = scale_I)
  idx <- bf$idx; feat <- bf$feat
  if (!length(idx)) return(matrix(0L, H, W))
  opt <- dbscan::optics(feat, minPts = minPts)
  lab <- dbscan::extractXi(opt, xi = xi)$cluster
  L <- matrix(0L, H, W)
  if (any(lab > 0)) {
    u <- sort(unique(lab[lab > 0])); re <- seq_along(u); names(re) <- u
    L[idx] <- as.integer(re[as.character(lab)])
  }
  L
}

# ---- Core → halo growth via hysteresis (per‑cluster or global threshold) ----
# P01: detection image (scaled to [0,1]); L_core: labels from OPTICS/HDBSCAN
# q_lo: per‑cluster quantile of the cluster seeds to set t_low; floor_lo: global minimum
# If per_cluster=FALSE, uses floor_lo for all clusters.
grow_labels_hysteresis <- function(P01, L_core,
                                   per_cluster = TRUE,
                                   q_lo = 0.3,      # was 0.40
                                   floor_lo = 0.4,  # was 0.48–0.50
                                   dilate_radius = 10L) {
  
  H <- nrow(P01); W <- ncol(P01)
  L_out <- matrix(0L, H, W)
  
  # 8-neighborhood flood within 'allow'
  flood <- function(allow, seeds_idx, lab_id) {
    q <- integer(length(allow)); head <- 1L; tail <- 0L
    push <- function(y,x){
      if (y>=1 && y<=H && x>=1 && x<=W && allow[y,x] && L_out[y,x]==0L) {
        L_out[y,x] <<- lab_id; tail <<- tail + 1L; q[tail] <<- (x-1L)*H + y
      }
    }
    for (s in seeds_idx) { y <- ((s-1L) %% H) + 1L; x <- ((s-1L) %/% H) + 1L; push(y,x) }
    while (head <= tail) {
      p <- q[head]; head <- head + 1L
      y <- ((p-1L) %% H) + 1L; x <- ((p-1L) %/% H) + 1L
      for (dy in -1:1) for (dx in -1:1) if (dy!=0 || dx!=0) push(y+dy, x+dx)
    }
  }
  
  ids <- sort(unique(L_core[L_core > 0]))
  for (k in ids) {
    seeds <- which(L_core == k); if (!length(seeds)) next
    thr <- if (per_cluster) max(quantile(P01[seeds], q_lo, na.rm=TRUE), floor_lo) else floor_lo
    allow <- (P01 >= thr)
    flood(allow, seeds, k)
  }
  
  # optional: light binary dilation to keep very faint outskirts
  if (dilate_radius > 0L && length(ids)) {
    r <- dilate_radius
    Ld <- L_out
    for (k in ids) {
      M <- (L_out == k)
      # simple square dilation (fast, no packages)
      Md <- M
      for (dy in -r:r) for (dx in -r:r) {
        y1 <- max(1, 1+dy); y2 <- min(H, H+dy); x1 <- max(1, 1+dx); x2 <- min(W, W+dx)
        Md[y1:y2, x1:x2] <- Md[y1:y2, x1:x2] | M[(y1-dy):(y2-dy), (x1-dx):(x2-dx)]
      }
      Ld[Md] <- k
    }
    L_out <- Ld
  }
  
  compact_labels(L_out)
}


# ---- Fill internal holes for each labeled region ----
fill_holes_per_label <- function(L) {
  H <- nrow(L); W <- ncol(L); out <- L
  for (k in sort(unique(L[L > 0]))) {
    idx <- which(L == k); if (!length(idx)) next
    yy <- ((idx - 1L) %% H) + 1L; xx <- ((idx - 1L) %/% H) + 1L
    y1 <- min(yy); y2 <- max(yy); x1 <- min(xx); x2 <- max(xx)
    sub <- out[y1:y2, x1:x2]
    bg  <- (sub == 0L)
    # mark border‑connected background via flood fill
    mark <- matrix(FALSE, nrow(sub), ncol(sub))
    q <- integer(length(bg)); head <- 1L; tail <- 0L
    push <- function(y,x){ if (bg[y,x] && !mark[y,x]) { mark[y,x] <<- TRUE; tail <<- tail + 1L; q[tail] <<- (x-1L)*nrow(sub) + y } }
    for (x in 1:ncol(sub)) { if (bg[1,x]) push(1,x); if (bg[nrow(sub),x]) push(nrow(sub),x) }
    for (y in 1:nrow(sub)) { if (bg[y,1]) push(y,1); if (bg[y,ncol(sub)]) push(y,ncol(sub)) }
    while (head <= tail) {
      p <- q[head]; head <- head + 1L
      y <- ((p-1L) %% nrow(sub)) + 1L; x <- ((p-1L) %/% nrow(sub)) + 1L
      if (y>1) push(y-1,x); if (y<nrow(sub)) push(y+1,x); if (x>1) push(y,x-1); if (x<ncol(sub)) push(y,x+1)
    }
    holes <- bg & !mark
    sub[holes] <- k
    out[y1:y2, x1:x2] <- sub
  }
  out
}

# ---- Apply binary/soft mask to cube ----
mask_cube <- function(cube, labels, mode = c("zero","na","soft"), sigma = 2) {
  mode <- match.arg(mode)
  H <- dim(cube)[1]; W <- dim(cube)[2]; B <- dim(cube)[3]
  m <- (labels > 0)
  if (mode == "soft") {
    # separable box blur approximation (odd kernel)
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
    for (b in seq_len(B)) { Xb <- out[,,b]; if (mode == "zero") Xb[!m] <- 0 else Xb[!m] <- NA_real_; out[,,b] <- Xb }
    return(out)
  }
}

# ---- Extract tiny cubes around each compact region ----
extract_region_cubes <- function(cube, labels, margin = 6L, min_size = 50L) {
  H <- dim(cube)[1]; W <- dim(cube)[2]
  regs <- split(seq_len(H*W), as.vector(labels))
  out <- list()
  for (lab in names(regs)) {
    if (lab == "0") next
    idx <- regs[[lab]]; if (length(idx) < min_size) next
    xy <- xy_from_idx(idx, H, W)
    y1 <- max(1, min(xy[,1]) - margin); y2 <- min(H, max(xy[,1]) + margin)
    x1 <- max(1, min(xy[,2]) - margin); x2 <- min(W, max(xy[,2]) + margin)
    out[[lab]] <- list(cube = cube[y1:y2, x1:x2, , drop = FALSE], bbox = c(y1,y2,x1,x2), npix = length(idx))
  }
  out
}

write_region_cubes_fits <- function(regions, prefix = "cutout") {
  i <- 0L
  for (nm in names(regions)) { i <- i + 1L; fn <- sprintf("%s_%03d.fits", prefix, i); FITSio::writeFITSim(fn, regions[[nm]]$cube) }
}












# ============================================================
# Example (adjust paths & parameters)
# ============================================================

# 1) Read cube & PCA detection image
Xfits <- FITSio::readFITS("datacube_reg1.fits")
cube  <- Xfits$imDat
H <- dim(cube)[1]; W <- dim(cube)[2]

# Capivara helper → (H*W) x B matrix
df_mat <- capivara::cube_to_matrix(Xfits)
P <- pca_energy_map(df_mat, H, W, d = 3)

# 2) Cores via OPTICS‑Xi (no eps to tune)
L_core <- segment_optics_xi(P, q_fore = 0.9, scale_xy = 1.0, scale_I = 1.8,
                            minPts = 10, xi = 0.05)

# 3) Grow cores to include halos/rings + fill internal holes
P01    <- robust01(P)
L_grow <- grow_labels_hysteresis(P01, L_core, per_cluster = TRUE, q_lo = 0.40, floor_lo = 0.48)
L_full <- fill_holes_per_label(L_grow)

# 4) Final clean‑up (size filtering, optional border removal)
L <- filter_by_size(L_full, min_size = 50)
# L <- remove_border_components(L)

# 5) Mask or extract sub‑cubes
cube_zero <- mask_cube(cube, L, mode = "na")
regions   <- extract_region_cubes(cube, L, margin = 10, min_size = 50)
# write_region_cubes_fits(regions, prefix = "cutout")

image(log(cube_zero[,,4]),col=viridis::viridis(100))












