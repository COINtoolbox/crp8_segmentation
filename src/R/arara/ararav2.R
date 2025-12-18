# ============================================================
# Hyperspectral clump segmentation via PCA + Density Clustering
# - PCA detection image
# - Child masks via DBSCAN/HDBSCAN/OPTICS
# - Hole fill (per label)
# - Merge children -> PARENT super-regions (no nested groups)
# - Cutouts
# ============================================================

# ---- Dependencies ----
# install.packages(c("FITSio","dbscan","viridis"))    # core
# install.packages("mclust")                          # optional for SED-GMM
library(FITSio)
library(dbscan)

# -------------------------- Utils ---------------------------
robust01 <- function(M, q = c(0.01, 0.99)) {
  v <- as.numeric(M); 
  v <- v[is.finite(v)]
  a <- stats::quantile(v, q[1], na.rm = TRUE);
  b <- stats::quantile(v, q[2], na.rm = TRUE)
  (pmin(pmax(M, a), b) - a) / max(b - a, .Machine$double.eps)
}

xy_from_idx <- function(idx, H, W) {
  cbind(y = ((idx - 1L) %% H) + 1L,
        x = ((idx - 1L) %/% H) + 1L)
}

compact_labels <- function(L) {
  u <- sort(unique(L[L > 0])); if (!length(u)) return(L)
  re <- seq_along(u); names(re) <- u
  L[L > 0] <- as.integer(re[as.character(L[L > 0])]); L
}

filter_by_size <- function(L, min_size = 30L, max_size = Inf) {
  if (!any(L>0)) return(L)
  tb <- table(L[L>0])
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

# ------------------ PCA detection image --------------------
pca_energy_map <- function(df_mat, H, W, d = 5) {
  pca <- stats::prcomp(df_mat, center = TRUE, scale. = FALSE)
  d <- min(d, ncol(pca$x))
  E <- rowSums(pca$x[, 1:d, drop = FALSE]^2)
  matrix(E, nrow = H, ncol = W)
}

# ---- Build features for density clustering on detection image
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

# ----------------- Segmentation methods --------------------
segment_dbscan <- function(P, q_fore=0.90, scale_xy=1.0, scale_I=1.0, eps=0.06, minPts=30) {
  H <- nrow(P); W <- ncol(P)
  P01 <- robust01(P); Fg <- P01 >= stats::quantile(P01, q_fore, na.rm=TRUE)
  idx <- which(Fg); if (!length(idx)) return(matrix(0L, H, W))
  yy <- ((idx - 1L) %% H) + 1L; xx <- ((idx - 1L) %/% H) + 1L
  feat <- cbind(scale_xy * (xx - 1)/(W - 1), scale_xy * (yy - 1)/(H - 1), scale_I * P01[idx])
  fit <- dbscan(feat, eps = eps, minPts = minPts)
  L <- matrix(0L, H, W)
  if (any(fit$cluster > 0)) {
    u <- sort(unique(fit$cluster[fit$cluster > 0])); re <- seq_along(u); names(re) <- u
    L[idx] <- as.integer(re[as.character(fit$cluster)])
  }
  L
}

segment_hdbscan <- function(P, q_fore=0.92, scale_xy=1.0, scale_I=1.5, minPts=25) {
  H <- nrow(P); W <- ncol(P)
  bf <- build_features(P, q_fore=q_fore, scale_xy=scale_xy, scale_I=scale_I)
  idx <- bf$idx; feat <- bf$feat; if (!length(idx)) return(matrix(0L,H,W))
  hd <- dbscan::hdbscan(feat, minPts = minPts)
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

# --------------- Fill internal holes (per label) --------------
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

# ---- Merge children → Parent super-regions (no nested groups)
merge_to_superregions <- function(L, bridge = 2L, min_size = 200L) {
  H <- nrow(L); W <- ncol(L)
  M <- (L > 0)
  
  if (bridge > 0L) {                       # small dilation bridges gaps
    Md <- M
    for (dy in -bridge:bridge) for (dx in -bridge:bridge) {
      y1 <- max(1, 1+dy); y2 <- min(H, H+dy)
      x1 <- max(1, 1+dx); x2 <- min(W, W+dx)
      Md[y1:y2, x1:x2] <- Md[y1:y2, x1:x2] | M[(y1-dy):(y2-dy), (x1-dx):(x2-dx)]
    }
    M <- Md
  }
  
  # connected components on M (8-neighb.)
  Lp <- matrix(0L, H, W)
  seen <- matrix(FALSE, H, W); id <- 0L
  for (y in 1:H) for (x in 1:W) {
    if (!M[y,x] || seen[y,x]) next
    id <- id + 1L
    qy <- y; qx <- x; head <- 1L; tail <- 1L
    seen[y,x] <- TRUE; Lp[y,x] <- id
    while (head <= tail) {
      cy <- qy[head]; cx <- qx[head]; head <- 1L + head
      for (dy in -1:1) for (dx in -1:1) if (dy!=0 || dx!=0) {
        ny <- cy+dy; nx <- cx+dx
        if (ny>=1 && ny<=H && nx>=1 && nx<=W && M[ny,nx] && !seen[ny,nx]) {
          seen[ny,nx] <- TRUE; Lp[ny,nx] <- id
          tail <- tail + 1L; qy[tail] <- ny; qx[tail] <- nx
        }
      }
    }
  }
  
  Lp <- filter_by_size(Lp, min_size = min_size)
  
  # map each child to a parent (majority overlap)
  child_ids <- sort(unique(L[L>0])); child2parent <- integer(max(c(child_ids, 0L)))
  for (k in child_ids) {
    ov <- table(Lp[L == k]); pid <- as.integer(names(ov)[which.max(ov)])
    child2parent[k] <- ifelse(length(pid), pid, 0L)
  }
  list(L_parent = Lp, child2parent = child2parent)
}

# -------- Cutouts for parent regions (preserve coordinates)
# fast square dilation (no extra deps)
dilate_mask_square <- function(M, r = 0L) {
  if (is.null(M) || r <= 0L) return(M)
  H <- nrow(M); W <- ncol(M)
  out <- M
  for (dy in -r:r) for (dx in -r:r) {
    y1 <- max(1, 1+dy); y2 <- min(H, H+dy)
    x1 <- max(1, 1+dx); x2 <- min(W, W+dx)
    out[y1:y2, x1:x2] <- out[y1:y2, x1:x2] | M[(y1-dy):(y2-dy), (x1-dx):(x2-dx)]
  }
  out
}

# Improved extractor
extract_region_cubes <- function(cube, labels,
                                 min_size = 50L,
                                 pad      = c("margin","none"),
                                 margin   = 6L,
                                 dilate_r = 0L,
                                 return_mask = TRUE,
                                 return_masked_cube = FALSE) {
  pad <- match.arg(pad)
  H <- dim(cube)[1]; W <- dim(cube)[2]; B <- dim(cube)[3]
  regs <- split(seq_len(H*W), as.vector(labels))
  out <- list()
  for (lab in names(regs)) {
    if (lab == "0") next
    idx <- regs[[lab]]
    if (length(idx) < min_size) next
    
    yy <- ((idx - 1L) %% H) + 1L
    xx <- ((idx - 1L) %/% H) + 1L
    # build a tight mask for this label
    Mfull <- matrix(FALSE, H, W); Mfull[cbind(yy, xx)] <- TRUE
    if (dilate_r > 0L) Mfull <- dilate_mask_square(Mfull, r = dilate_r)
    
    # bbox choice
    y1 <- min(yy); y2 <- max(yy); x1 <- min(xx); x2 <- max(xx)
    if (pad == "margin") {
      y1 <- max(1, y1 - margin); y2 <- min(H, y2 + margin)
      x1 <- max(1, x1 - margin); x2 <- min(W, x2 + margin)
    }
    cut_cube <- cube[y1:y2, x1:x2, , drop = FALSE]
    cut_mask <- Mfull[y1:y2, x1:x2, drop = FALSE]
    
    if (return_masked_cube) {
      for (b in seq_len(B)) {
        Xb <- cut_cube[,,b]; Xb[!cut_mask] <- NA_real_
        cut_cube[,,b] <- Xb
      }
    }
    
    out[[lab]] <- list(
      parent_id = as.integer(lab),
      cube = cut_cube,
      mask = if (return_mask) cut_mask else NULL,
      bbox = c(y1, y2, x1, x2),
      npix = length(idx)
    )
  }
  out
}





# --------------- Parent-level catalogue --------------------
catalog_parent <- function(cube, L_parent) {
  B <- dim(cube)[3]; H <- dim(cube)[1]; W <- dim(cube)[2]
  ids <- sort(unique(L_parent[L_parent > 0])); if (!length(ids)) return(list(meta=NULL, flux=NULL))
  meta <- data.frame(parent_id = ids, npix = 0L, ycen = 0.0, xcen = 0.0)
  flux <- matrix(0.0, nrow = length(ids), ncol = B)
  for (i in seq_along(ids)) {
    id <- ids[i]; mask <- (L_parent == id); idx <- which(mask)
    meta$npix[i] <- length(idx)
    xy <- xy_from_idx(idx, H, W)
    meta$ycen[i] <- mean(xy[,1]); meta$xcen[i] <- mean(xy[,2])
    for (b in seq_len(B)) flux[i,b] <- sum(cube[,,b][mask], na.rm = TRUE)
  }
  colnames(flux) <- paste0("flux_b", seq_len(B))
  list(meta = meta, flux = flux)
}

# ---------------- Optional SED segmentation (GMM) ------------
segment_sed_gmm <- function(cube_reg, mask_reg = NULL, d_spec = 5, kmax = 6, min_size = 20) {
  if (!requireNamespace("mclust", quietly = TRUE)) stop("Please install 'mclust'")
  H <- dim(cube_reg)[1]; W <- dim(cube_reg)[2]; B <- dim(cube_reg)[3]
  if (is.null(mask_reg)) mask_reg <- matrix(TRUE, H, W)
  idx <- which(mask_reg); if (!length(idx)) return(list(L_sed = matrix(0L,H,W), model = NULL, pca = NULL))
  X <- matrix(NA_real_, length(idx), B)
  for (b in seq_len(B)) X[,b] <- cube_reg[,,b][idx]
  X <- log1p(pmax(X, 0))
  X <- scale(X, center = TRUE, scale = TRUE)
  pc <- stats::prcomp(X, center = TRUE, scale. = FALSE)
  d <- min(d_spec, ncol(pc$x)); Z <- pc$x[, 1:d, drop = FALSE]
  fit <- mclust::Mclust(Z, G = 1:kmax, verbose = FALSE)
  lab <- integer(H * W); lab[idx] <- as.integer(fit$classification)
  Ls <- matrix(lab, H, W)
  Ls <- filter_by_size(Ls, min_size = min_size)
  list(L_sed = Ls, model = fit, pca = pc)
}

# ----------------- Mask whole cube (helper) -----------------
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

# ========================= Example ==========================
# 1) Read cube & build PCA detection image
Xfits <- FITSio::readFITS("datacube_reg1.fits")
cube  <- Xfits$imDat
H <- dim(cube)[1]; W <- dim(cube)[2]

df_mat <- if (requireNamespace("capivara", quietly = TRUE)) {
  capivara::cube_to_matrix(Xfits)
} else {
  B <- dim(cube)[3]; M <- matrix(NA_real_, H*W, B)
  for (b in seq_len(B)) M[,b] <- as.vector(cube[,,b]); M
}

P  <- pca_energy_map(df_mat, H, W, d = 2)

# 2) CHILD segmentation (choose one)
# L_child <- segment_dbscan(P, q_fore=0.90, scale_xy=1.0, scale_I=1.0, eps=0.05, minPts=10)
L_child <- segment_hdbscan(P, q_fore=0.85, scale_xy=1.0, scale_I=2.0, minPts=30)
# L_child <- segment_optics_xi(P, q_fore=0.90, scale_xy=1.0, scale_I=1.2, minPts=15, xi=0.05)

# 3) Clean children + fill interiors
L_child <- filter_by_size(L_child, min_size = 25)
L_child <- fill_holes_per_label(L_child)

# 4) Merge to PARENT super-regions (no nested groups)
super    <- merge_to_superregions(L_child, bridge = 2, min_size = 300)
L_parent <- super$L_parent

# 5) Catalogue + cutouts (per parent)
catalog  <- catalog_parent(cube, L_parent)


cube_2 <- mask_cube(cube,L_parent)

regions  <- extract_region_cubes(cube, L_parent,
                                 min_size = 100,
                                 pad = "margin", margin = 8,
                                 return_masked_cube = TRUE)





save_regions_as_fits <- function(regions,
                                 out_dir  = "regions_fits",
                                 prefix   = "region",
                                 write_masks = TRUE,
                                 catalog_csv = "regions_catalog.csv",
                                 overwrite   = TRUE) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  .safe_write <- function(expr, path) {
    if (file.exists(path) && !overwrite) {
      message("Skipping existing file: ", path)
      return(invisible(FALSE))
    }
    eval(expr)
    invisible(TRUE)
  }
  
  meta <- list()
  labs <- names(regions)
  ord  <- suppressWarnings(order(as.integer(labs)))
  labs <- labs[if (all(is.finite(ord))) ord else seq_along(labs)]
  
  for (lab in labs) {
    reg <- regions[[lab]]
    if (is.null(reg) || is.null(reg$cube)) next
    
    id   <- if (nzchar(lab)) as.integer(lab) else NA_integer_
    cube <- reg$cube               # 3D array
    mask <- reg$mask               # 2D logical
    bbox <- reg$bbox
    npix <- reg$npix
    
    # --- Replace NAs or out-of-mask with zeros ---
    if (!is.null(mask)) {
      for (b in seq_len(dim(cube)[3])) {
        Xb <- cube[,,b]
        Xb[!mask | !is.finite(Xb)] <- 0
        cube[,,b] <- Xb
      }
    } else {
      cube[!is.finite(cube)] <- 0
    }
    
    storage.mode(cube) <- "double"
    
    f_cube <- file.path(out_dir, sprintf("%s_%03d.fits", prefix, id))
    f_mask <- file.path(out_dir, sprintf("%s_%03d_mask.fits", prefix, id))
    
    .safe_write(quote(FITSio::writeFITSim(cube, f_cube)), f_cube)
    
    if (isTRUE(write_masks) && !is.null(mask)) {
      mask_im <- matrix(as.integer(mask), nrow(mask), ncol(mask))
      .safe_write(quote(FITSio::writeFITSim(mask_im, f_mask)), f_mask)
    }
    
    meta[[length(meta) + 1L]] <- data.frame(
      parent_id = id,
      file_cube = f_cube,
      file_mask = if (isTRUE(write_masks) && !is.null(mask)) f_mask else NA_character_,
      y1 = bbox[1], y2 = bbox[2], x1 = bbox[3], x2 = bbox[4],
      npix = npix,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(meta)) {
    catdf <- do.call(rbind, meta)
    utils::write.csv(catdf, file.path(out_dir, catalog_csv), row.names = FALSE)
    message("Wrote catalog: ", file.path(out_dir, catalog_csv))
    return(invisible(catdf))
  } else {
    warning("No regions found to write.")
    invisible(NULL)
  }
}


catdf <- save_regions_as_fits(
  regions,
  out_dir = "regions_fits",   # <- inside Region2
  prefix  = "parent",
  write_masks = TRUE,
  catalog_csv = "parent_regions_catalog.csv",
  overwrite = TRUE
)



# 7) Quick looks
image(L_parent, main = "Parent labels")
cube_na <- mask_cube(cube, L_parent, mode = "na")

image(log(cube_na[,,5]), col = viridis::plasma(100))




# 5) Catalogue + cutouts (per parent)
catalog  <- catalog_parent(cube, L_parent)
# 7) Quick looks
image(L_parent, main = "Parent labels")









# install.packages(c("ggplot2","viridis","patchwork"))  # if needed
library(ggplot2)
library(viridis)
library(patchwork)

# ---------------------- helpers ----------------------

# Convert a 2D matrix to a data.frame for ggplot (origin at top-left)
# helper from earlier
as_image_df <- function(M) {
  reshape2::melt(M)
  }

# pick region "1", band 1
Xp <-  regions$`2`$cube[,,7]                # or regions[[1]] if you want by index
# optional safe log:
Xp <- log(Xp)

df <- as_image_df(Xp)


mask_contours_df <- function(Mask) as_image_df(Mask * 1)

df_mask <- mask_contours_df(regions[["2"]][["mask"]])


ggplot(df, aes(Var1, Var2)) +
  geom_raster(interpolate = FALSE,aes(fill = value)) +
  geom_contour(data = df_mask, aes(Var1, Var2,z = value), color = "red3", size = 0.1) +
  coord_fixed() +
  scale_fill_continuous_sequential("grays",na.value = "white") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  xlab("px") + ylab("px")




# df: columns Var1 (x), Var2 (y), value
nx <- max(df$Var1); ny <- max(df$Var2)

# put values into a matrix indexed by (row=y, col=x)
M <- matrix(NA_real_, nrow = ny, ncol = nx)
M[cbind(df$Var2, df$Var1)] <- df$value

# terra raster (flip rows to match raster’s lower-left origin)
r <- rast(M[ny:1, , drop = FALSE])
ext(r) <- ext(1, nx, 1, ny)       # cell centers at integer coords
crs(r) <- ""                      # image/grid space

# -> sf pixel polygons with the value attribute
sf_pixels <- st_as_sf(as.polygons(r, values = TRUE, dissolve = FALSE))
names(sf_pixels)[names(sf_pixels) == "lyr.1"] <- "value"

# save if you want
# st_write(sf_pixels, "flux_pixels.gpkg", delete_dsn = TRUE)

# df_mask: Var1, Var2, value (1 inside; 0/NA outside)
Mmask <- matrix(0L, nrow = ny, ncol = nx)
Mmask[cbind(df_mask$Var2, df_mask$Var1)] <- as.integer(df_mask$value > 0)

rmask <- rast(Mmask[ny:1, , drop = FALSE])
ext(rmask) <- ext(1, nx, 1, ny)

# polygonize: one (or more) polygons for the mask == 1
polys <- as.polygons(rmask, values = TRUE, dissolve = TRUE, na.rm = TRUE)
mask_sf <- st_as_sf(polys[polys$lyr.1 == 1, ])
st_crs(mask_sf) <- NA_character_




# Example kernel: 3x3 mean
K <- matrix(1/9, 3, 3)

# Option A: mask first, then convolve (outside is NA; omit NAs at edges)
r_in  <- mask(r, vect(mask_sf))
r_flt <- focal(r_in, w = K, na.policy = "omit")

# (Alternative: convolve then mask — similar result; pick what you prefer)
# r_flt <- mask(focal(r, w = K, na.policy = "omit"), vect(mask_sf))



sf_pixels_flt <- st_as_sf(as.polygons(r_flt, values = TRUE, dissolve = FALSE))
names(sf_pixels_flt)[names(sf_pixels_flt) == "lyr.1"] <- "value"

# save/plot
# st_write(sf_pixels_flt, "flux_filtered_pixels.gpkg", delete_dsn = TRUE)
# plot(sf_pixels_flt["value"])


ggplot() +
  geom_sf(data = sf_pixels, aes(fill = value), color = NA) +
  geom_sf(data = mask_sf, fill = NA, color = "white", linewidth = 0.6) +
  coord_sf(expand = FALSE) +
  scale_fill_viridis_c(na.value = "black") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), legend.position = "none")





# --- Gaussian kernel (σ in pixels) ------------------------------------------
gauss_kernel <- function(sigma, radius = ceiling(3*sigma)) {
  x <- -radius:radius
  K <- exp(-(outer(x,x,function(a,b) (a^2+b^2)))/(2*sigma^2))
  K / sum(K)
}

# --- Find flux peaks (centroids) --------------------------------------------
# r: SpatRaster with your flux (has values!)
# mask_sf: ROI polygon(s)
# fwhm_px: PSF FWHM in pixels  (sigma ≈ fwhm/2.355)
# thr_q: keep only peaks above this quantile of the smoothed flux inside mask
# min_sep: minimum separation between peaks (pixels) ~ PSF
# refine_radius: subpixel COM refinement window radius (pixels)
find_flux_centroids <- function(r, mask_sf, fwhm_px = 3,
                                thr_q = 0.90, min_sep = NULL, refine_radius = 0) {
  stopifnot(terra::hasValues(r))
  if (is.null(min_sep)) min_sep <- fwhm_px
  sigma <- fwhm_px / 2.355
  
  # 1) restrict to mask
  rin <- terra::mask(r, terra::vect(mask_sf))
  
  # 2) smooth ~ PSF
  K  <- gauss_kernel(sigma)
  rs <- terra::focal(rin, w = K, fun = sum, na.policy = "omit")  # weighted sum
  rs <- rs / sum(K)
  
  # 3) local maxima by max-filter
  win <- matrix(1, nrow = 2*ceiling(min_sep)+1, ncol = 2*ceiling(min_sep)+1)
  rmax <- terra::focal(rs, w = win, fun = max, na.policy = "omit")
  vals <- terra::values(rs)
  keep <- which(!is.na(vals) & abs(vals - terra::values(rmax)) < 1e-12)
  
  # 4) threshold
  thr  <- stats::quantile(vals, thr_q, na.rm = TRUE)
  keep <- keep[ vals[keep] >= thr ]
  if (!length(keep)) return(st_sf(flux = numeric(0), geometry = st_sfc(), crs = NA))
  
  # 5) candidate coordinates
  xy   <- terra::xyFromCell(rs, keep)
  v    <- as.numeric(vals[keep])
  
  # 6) suppress nearby duplicates (keep brightest within min_sep)
  if (length(keep) > 1L && min_sep > 0) {
    d  <- as.matrix(dist(xy))
    hc <- hclust(as.dist(d))
    grp <- cutree(hc, h = min_sep)
    keep_idx <- unlist(tapply(seq_along(grp), grp, function(ii) ii[which.max(v[ii])]))
    xy <- xy[keep_idx, , drop = FALSE]
    v  <- v[keep_idx]
  }
  
  # 7) optional sub-pixel refinement by center-of-mass in a small window
  if (refine_radius > 0) {
    refine_one <- function(px, py) {
      ex <- terra::ext(px - refine_radius, px + refine_radius,
                       py - refine_radius, py + refine_radius)
      cwin <- terra::crop(rin, ex)
      idx <- which(!is.na(terra::values(cwin)))
      if (!length(idx)) return(c(px, py))
      xyw <- terra::xyFromCell(cwin, idx)
      fw  <- as.numeric(terra::values(cwin)[idx])
      w   <- fw / sum(fw)
      c(sum(w*xyw[,1]), sum(w*xyw[,2]))
    }
    xy <- t(apply(xy, 1, function(p) refine_one(p[1], p[2])))
  }
  
  st_as_sf(data.frame(flux = v, x = xy[,1], y = xy[,2]),
           coords = c("x","y"), crs = NA)
}


# r: SpatRaster built from your df (has values)
# mask_sf: your ROI polygon(s)
cent <- find_flux_centroids(r, mask_sf,
                            fwhm_px = 3,   # set to your PSF in pixels
                            thr_q   = 0.90,
                            min_sep = 3,   # ~ PSF
                            refine_radius = 2)

# inspect numeric coordinates
st_coordinates(cent)
# plot on your ggplot
#library(ggplot2)
ggplot() +
  geom_sf(data = sf_pixels, aes(fill = value), color = NA) +
  geom_sf(data = mask_sf, fill = NA, color = "white", linewidth = 0.6) +
  geom_sf(data = cent, color = "red", shape = 1, size = 5, stroke = 1) +
  coord_sf(expand = FALSE) +
  scale_fill_continuous_sequential("grays",na.value = "white") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), legend.position = "none")


