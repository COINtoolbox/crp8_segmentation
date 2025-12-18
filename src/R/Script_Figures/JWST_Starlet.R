library(guara)

X <- FITSio::readFITS("CenA_1.fits")
cube <- X$imDat
cube_z <- (cube - mean(cube, na.rm=TRUE)) / sd(cube, na.rm=TRUE)
#cube <- cube[,,1:6]
nx <- dim(cube)[1]; ny <- dim(cube)[2]; nlam <- dim(cube)[3]
#collapsed <- collapse_white_light(cube, kclip = 2, use_weights = TRUE)



# Starlet: masked decomposition
J <- 7L
dec <- guara::starlet_mask(cube_z, J = J)

# Reconstruction from all planes + coarse should match collapsed (within tiny tol)
rec_all <- starlet_reconstruct(dec, keep_scales = c(2:3),
                                   include_coarse = FALSE,
                                   denoise_k = 3,
                                   mode="hard")



pal = c(viridis::inferno(256))
par(mfrow=c(1,2), mar=c(2,2,2,1))
rec_all[rec_all<=0] <- NA

image(asinh(cube_z) - asinh(rec_all), col=pal)
image(rec_all, col=pal)




out <- guara_segment(dec,
                     keep_scales = 2:3,      # braços
                     include_coarse = FALSE, # remove fundo (cJ fica fora)
                     per_scale_positive = TRUE,
                     k_lo = 5, k_hi = 8,
                     area_min = 20)



rec_all[rec_all>0] <- 1

out_file <- "mask_NGC300_1.fits"
FITSio::writeFITSim(rec_all, file = out_file, header = X$header, axDat = X$axDat)


library(dbscan)
library(dplyr)

# --- 1) Flatten the matrix ---
df <- which(is.finite(rec_all), arr.ind = TRUE)
df <- data.frame(y = df[,1], x = df[,2], value = rec_all[df])

# --- 2) Keep only bright pixels ---
thr <- quantile(df$value, 0.98, na.rm = TRUE)
pts <- df[df$value > thr, ]

# --- 3) Run HDBSCAN ---
pts_scaled <- scale(pts[, c("x", "y")])
hdb <- hdbscan(pts_scaled, minPts = 5)
pts$cluster <- hdb$cluster

# --- 4) Reconstruct cluster-labeled matrix ---
cluster_mat <- matrix(0L, nrow = nrow(rec_all), ncol = ncol(rec_all))
sel <- pts$cluster > 0
if (any(sel)) {
  cluster_mat[cbind(pts$y[sel], pts$x[sel])] <- pts$cluster[sel]
}

# cluster_mat: same dims as rec_all, 0 = background, >0 = cluster id
df_lab <- as.data.frame(as.table(cluster_mat))
names(df_lab) <- c("y","x","cluster")
df_lab$y <- as.integer(df_lab$y)
df_lab$x <- as.integer(df_lab$x)
df_lab$cluster <- as.integer(df_lab$cluster)
df_lab$cluster_f <- factor(ifelse(df_lab$cluster == 0, NA, df_lab$cluster))

ggplot(df_lab, aes(x = x, y = y, fill = cluster_f)) +
  geom_raster(na.rm = FALSE) +
  scale_y_reverse() + coord_equal() +
  scale_fill_viridis_d(na.value = "gray") +
  theme_minimal() +
  theme(legend.position = "none")



# Cluster a sample with hclust, then assign all pixels to nearest centroid
guara_hclust_sample_assign <- function(cube, mask = NULL, N = 50,
                                       m_limit = 8000, seed = 1) {
  stopifnot(length(dim(cube)) == 3L)
  nx <- dim(cube)[1]; ny <- dim(cube)[2]; nb <- dim(cube)[3]
  M  <- matrix(cube, nrow = nx*ny, ncol = nb)

  # build mask of usable pixels: finite & nonzero; AND optional user mask
  usable <- rowSums(is.finite(M)) == nb & rowSums(abs(M)) > 0
  if (!is.null(mask)) {
    stopifnot(all(dim(mask) == c(nx, ny)))
    usable <- usable & as.vector(mask)
  }
  idx_all <- which(usable)
  m <- length(idx_all)
  if (m == 0) stop("No usable pixels after filtering.")

  X <- M[idx_all, , drop = FALSE]
  X[!is.finite(X)] <- 0

  # sample for hierarchical
  set.seed(seed)
  if (m > m_limit) {
    idx_s <- sort(sample.int(m, m_limit))
  } else {
    idx_s <- seq_len(m)
  }
  Xs <- X[idx_s, , drop = FALSE]

  # hclust on the sample
  d  <- dist(Xs, method = "euclidean")
  hc <- hclust(d, method = "ward.D2")
  lab_s <- cutree(hc, k = N)

  # centroids in SED space (on the sample)
  centers <- vapply(1:N, function(k) {
    colMeans(Xs[lab_s == k, , drop = FALSE])
  }, numeric(ncol(Xs)))
  centers <- t(centers)  # N x nb

  # assign all usable pixels to nearest centroid
  # (fast squared Euclidean via matrix ops)
  # dist^2(a,b) = ||a||^2 + ||b||^2 - 2 a·b
  a2 <- rowSums(X * X)                 # m
  b2 <- rowSums(centers * centers)     # N
  G  <- X %*% t(centers)               # m x N
  D2 <- outer(a2, b2, "+") - 2 * G     # m x N
  lab_all <- max.col(-D2)              # nearest centroid index

  # map back to image
  labels_img <- matrix(0L, nx, ny)
  labels_img[idx_all] <- lab_all

  list(labels_img = labels_img,
       labels_vec = lab_all,
       used_mask  = matrix(usable, nx, ny),
       sample_idx = idx_all[idx_s])
}


res <- guara_hclust_sample_assign(cube, mask = rec_all, N = 10, m_limit = 8000)
pal = c(viridis::plasma(10))


# ---- Simple, robust asinh stretch for astronomy ----
asinh_stretch <- function(x, qlo = 0.001, qhi = 0.999, scale = NULL,
                          nonneg = TRUE, na.rm = TRUE) {
  # Copy to avoid modifying original
  z <- x

  # Handle NAs
  if (na.rm) z[!is.finite(z)] <- NA_real_

  # Optional non-negativity (useful when background is >= 0 in expectation)
  if (nonneg) z[z < 0] <- 0

  # Robust low/high via percentiles
  lo <- as.numeric(stats::quantile(z, probs = qlo, na.rm = TRUE, names = FALSE))
  hi <- as.numeric(stats::quantile(z, probs = qhi, na.rm = TRUE, names = FALSE))
  if (!is.finite(lo)) lo <- 0
  if (!is.finite(hi) || hi <= lo) hi <- max(z, na.rm = TRUE)

  # Shift + scale into a positive range
  z <- pmax(z - lo, 0)

  # If no scale given, set it to a fraction of dynamic range (Lupton-ish)
  if (is.null(scale)) {
    scale <- (hi - lo) / 20  # smaller => more aggressive compression
    if (!is.finite(scale) || scale == 0) scale <- 1
  }

  # Asinh compress and normalize to 0..1
  z2  <- asinh(z / scale)
  z2  <- z2 / max(z2, na.rm = TRUE)

  # Keep NAs where input was NA
  z2[!is.finite(z2)] <- NA_real_
  z2
}



par(mfrow=c(1,2), mar=c(2,2,2,1))
image(asinh_stretch(collapsed), col=pal)
#image(asinh_stretch(out$hard_rec),col=pal)
image(log((res$labels_img)),  col=pal)


SED  <- RegionPhotometry(X,res$labels_img,
                        error_fallback = "poisson")

flux_wide_pretty <- make_flux_wide(SED$flux_long, spec_nircam)






# Look at the header
names(flux_wide_pretty)

filters <- c(
  "F090W","F115W","F150W","F182M","F200W",
  "F210M","F277W","F335M","F356W","F410M",
  "F430M","F444W","F460M","F480M"
)

nircam_lambda_um <- c(
  F090W = 0.902, F115W = 1.154, F150W = 1.501, F182M = 1.842,
  F200W = 1.989, F210M = 2.099, F277W = 2.770, F335M = 3.365,
  F356W = 3.563, F410M = 4.082, F430M = 4.286, F444W = 4.421,
  F460M = 4.620,  F480M = 4.828
)


spec_nircam <- build_filter_spec(
  filters = filters,
  lambda  = nircam_lambda_um,
  throughput_input = "/Users/rd23aag/Downloads/nircam_throughputs/mean_throughputs",
  exclude = c("F150W2"),
  assume_ordered = TRUE,
  band_index = 0:(length(filters)-1)  # if your SED bands are 0..12
)

p <- plot_sed_with_spec(
  reg_long = SED$flux_long,
  region_id = 10,
  spec = spec_nircam,
  scale_mode = "auto",
  band_height_frac = 0.90
)
p



Xm <- FITSio::readFITS("/Users/rd23aag/Documents/GitHub/crp8_segmentation/data/raw/datacube_reg1.fits")
cm <- Xm$imDat


manga <- capivara::segment(Xm,N=12)
image(manga$cluster_map,col=plasma(12),asp=1)



res <- guara_sed_cluster(
  cube,
  mask    = out$hard_rec,   # logical mask from your cleaned reconstruction
  N       = 50,            # number of clusters
  method  = "kmeans",       # or "hclust" for small datasets
  max_hclust = 100000,
  seed    = 1
)


lab_img_plot <- res$labels_img
lab_img_plot[is.na(lab_img_plot) | lab_img_plot == 0] <- NA  # keep background as NA

# Plot safely
image(lab_img_plot,
      col = c(magma(100)),
      useRaster = TRUE,
      asp = 1,
      main = "SED-based Clustering")


cube_na <- guara::mask_cube(cube,out$hard_rec, mode = "na")
cub_cut <- cube_na
cube_cap <- list(imDat = cub_cut)   # cube_2 is your [nx,ny,nb] array
seg <- capivara::segment(cube_cap,N=20)

image((seg$cluster_map),col=seg$cluster_snr)

