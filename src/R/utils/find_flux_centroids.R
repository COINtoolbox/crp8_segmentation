library(terra)
library(sf)

# ---- SpatRaster extent -> polygon (sf) ----
raster_extent_polygon <- function(r) {
  e <- terra::ext(r)
  cr <- terra::crs(r, proj = TRUE)
  coords <- matrix(c(e$xmin, e$ymin,
                     e$xmax, e$ymin,
                     e$xmax, e$ymax,
                     e$xmin, e$ymax,
                     e$xmin, e$ymin), ncol = 2, byrow = TRUE)
  st_sf(geometry = st_sfc(st_polygon(list(coords)), crs = cr))
}

# ---- Matrix -> SpatRaster (your version, kept) ----
# ---- Matrix -> SpatRaster (com controle de transposição) ----
matrix_to_rast <- function(M,
                           template = NULL,
                           extent   = NULL,
                           crs      = NA,
                           xres = 1, yres = 1, x0 = 0, y0 = 0,
                           origin = c("upper","lower"),
                           transpose = FALSE) {       # <- NOVO
  stopifnot(is.matrix(M))
  origin <- match.arg(origin)
  if (transpose) M <- t(M)                           # <- NOVO
  
  r <- terra::rast(M)                                # row 1 = TOP (em terra)
  if (origin == "lower") r <- terra::flip(r, direction = "vertical")
  
  if (!is.null(template)) {
    terra::ext(r) <- terra::ext(template)
    terra::crs(r) <- terra::crs(template)
  } else if (!is.null(extent)) {
    terra::ext(r) <- extent
    terra::crs(r) <- crs
  } else {
    terra::ext(r) <- terra::ext(x0, x0 + ncol(r)*xres, y0, y0 + nrow(r)*yres)
    terra::crs(r) <- crs
  }
  r
}


# ---- SpatRaster -> Matrix (your version, kept) ----
rast_to_matrix <- function(r, origin = c("upper","lower")) {
  origin <- match.arg(origin)
  M <- terra::as.matrix(r, wide = TRUE)  # rows from TOP to BOTTOM
  if (origin == "lower") M <- M[nrow(M):1, , drop = FALSE]
  M
}

gauss_kernel <- function(sigma, radius = ceiling(3*sigma)) {
  stopifnot(sigma > 0)
  x <- -radius:radius
  K <- exp(-(outer(x, x, function(a,b) (a^2 + b^2))) / (2*sigma^2))
  K / sum(K)
}

# ---- Peak finder with pixel-aware 'min_sep' and CRS propagation ----
find_flux_centroids <- function(r, mask_sf,
                                fwhm_px = 3, thr_q = 0.90,
                                min_sep = NULL, refine_radius = 0) {
  stopifnot(terra::hasValues(r))
  if (terra::nlyr(r) != 1L) r <- r[[1]]            # ensure single layer
  if (is.null(min_sep)) min_sep <- fwhm_px
  sigma <- fwhm_px / 2.355
  
  # 1) restrict to mask
  rin <- terra::mask(r, terra::vect(mask_sf))
  
  # 2) smooth ~ PSF
  K  <- gauss_kernel(sigma)
  rs <- terra::focal(rin, w = K, fun = "sum", na.rm = TRUE, pad = TRUE)
  rs <- rs / sum(K)
  
  # 3) local maxima (moving max with window ~ min_sep px)
  win_sz <- 2L * ceiling(min_sep) + 1L
  win  <- matrix(1, nrow = win_sz, ncol = win_sz)
  rmax <- terra::focal(rs, w = win, fun = "max", na.rm = TRUE, pad = TRUE)
  
  vals <- as.numeric(terra::values(rs))
  vmax <- as.numeric(terra::values(rmax))
  
  # avoid float issues
  keep <- which(!is.na(vals) & !is.na(vmax) & (vals >= (vmax - .Machine$double.eps^0.5)))
  
  # 4) threshold inside mask (quantile over in-mask values)
  thr  <- stats::quantile(vals, thr_q, na.rm = TRUE, names = FALSE)
  keep <- keep[ vals[keep] >= thr ]
  if (!length(keep)) {
    return(sf::st_sf(flux = numeric(0), geometry = sf::st_sfc(), crs = terra::crs(r, proj = TRUE)))
  }
  
  # 5) candidate coords
  xy <- terra::xyFromCell(rs, keep)
  v  <- vals[keep]
  
  # 6) suppress neighbors within min_sep *pixels* (convert to map units)
  res_xy <- terra::res(r)  # (xres, yres) in map units
  # use average pixel size to build an isotropic distance in map units
  px_unit <- mean(res_xy)
  h <- min_sep * px_unit
  if (length(keep) > 1L && h > 0) {
    d  <- as.matrix(stats::dist(xy))
    hc <- stats::hclust(stats::as.dist(d), method = "complete")
    grp <- stats::cutree(hc, h = h)
    keep_idx <- unlist(tapply(seq_along(grp), grp, function(ii) ii[which.max(v[ii])]))
    xy <- xy[keep_idx, , drop = FALSE]; v <- v[keep_idx]
  }
  
  # 7) optional sub-pixel COM refinement (radius in *pixels*)
  if (refine_radius > 0) {
    refine_one <- function(px, py) {
      rx <- refine_radius * res_xy[1]
      ry <- refine_radius * res_xy[2]
      ex <- terra::ext(px - rx, px + rx, py - ry, py + ry)
      cwin <- terra::crop(rin, ex)
      vv   <- terra::values(cwin)
      idx  <- which(!is.na(vv))
      if (!length(idx)) return(c(px, py))
      xyw  <- terra::xyFromCell(cwin, idx)
      w    <- as.numeric(vv[idx])
      w    <- w / sum(w)
      c(sum(w * xyw[,1]), sum(w * xyw[,2]))
    }
    xy <- t(apply(xy, 1, function(p) refine_one(p[1], p[2])))
  }
  
  # propagate CRS from raster
  cr <- terra::crs(r, proj = TRUE)
  sf::st_as_sf(data.frame(flux = v, x = xy[,1], y = xy[,2]),
               coords = c("x","y"), crs = cr)
}
