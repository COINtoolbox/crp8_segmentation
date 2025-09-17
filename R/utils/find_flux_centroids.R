# ---- Matrix -> SpatRaster (robust) ----
matrix_to_rast <- function(M,
                           template = NULL,          # optional SpatRaster to copy ext+CRS
                           extent   = NULL,          # or terra::ext(xmin,xmax,ymin,ymax)
                           crs      = NA,            # set CRS if not using template
                           xres = 1, yres = 1, x0 = 0, y0 = 0,  # used if extent=NULL & no template
                           origin = c("upper","lower")) {        # does M[1,] represent TOP or BOTTOM row?
  stopifnot(is.matrix(M))
  origin <- match.arg(origin)
  
  # start from terra's native constructor (handles row-major correctly)
  r <- terra::rast(M)  # by default row 1 -> TOP row
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

# ---- SpatRaster -> Matrix (inverse; matches 'origin') ----
rast_to_matrix <- function(r, origin = c("upper","lower")) {
  origin <- match.arg(origin)
  M <- terra::as.matrix(r, wide = TRUE)  # rows from TOP to BOTTOM
  if (origin == "lower") M <- M[nrow(M):1, , drop = FALSE]
  M
}

gauss_kernel <- function(sigma, radius = ceiling(3*sigma)) {
  x <- -radius:radius
  K <- exp(-(outer(x, x, function(a,b) (a^2 + b^2))) / (2*sigma^2))
  K / sum(K)
}

find_flux_centroids <- function(r, mask_sf,
                                fwhm_px = 3, thr_q = 0.90,
                                min_sep = NULL, refine_radius = 0) {
  stopifnot(terra::hasValues(r))
  if (is.null(min_sep)) min_sep <- fwhm_px
  sigma <- fwhm_px / 2.355
  
  # 1) restrict to mask
  rin <- terra::mask(r, terra::vect(mask_sf))
  
  # 2) smooth ~ PSF
  K  <- gauss_kernel(sigma)
  rs <- terra::focal(rin, w = K, fun = "sum", na.rm = TRUE, pad = TRUE)
  rs <- rs / sum(K)
  
  # 3) local maxima (moving max with window ~ min_sep)
  win  <- matrix(1, nrow = 2*ceiling(min_sep)+1, ncol = 2*ceiling(min_sep)+1)
  rmax <- terra::focal(rs, w = win, fun = "max", na.rm = TRUE, pad = TRUE)
  
  vals <- terra::values(rs)
  vmax <- terra::values(rmax)
  keep <- which(!is.na(vals) & (abs(vals - vmax) < 1e-12))
  
  # 4) threshold inside mask
  thr  <- stats::quantile(vals, thr_q, na.rm = TRUE)
  keep <- keep[ vals[keep] >= thr ]
  if (!length(keep)) return(sf::st_sf(flux = numeric(0), geometry = sf::st_sfc(), crs = NA))
  
  # 5) candidate coords
  xy <- terra::xyFromCell(rs, keep)
  v  <- as.numeric(vals[keep])
  
  # 6) suppress neighbors within min_sep (keep brightest)
  if (length(keep) > 1L && min_sep > 0) {
    d  <- as.matrix(dist(xy))
    hc <- hclust(as.dist(d))
    grp <- cutree(hc, h = min_sep)
    keep_idx <- unlist(tapply(seq_along(grp), grp, function(ii) ii[which.max(v[ii])]))
    xy <- xy[keep_idx, , drop = FALSE]; v <- v[keep_idx]
  }
  
  # 7) optional sub-pixel COM refinement
  if (refine_radius > 0) {
    refine_one <- function(px, py) {
      ex <- terra::ext(px - refine_radius, px + refine_radius,
                       py - refine_radius, py + refine_radius)
      cwin <- terra::crop(rin, ex)
      vv   <- terra::values(cwin)
      idx  <- which(!is.na(vv))
      if (!length(idx)) return(c(px, py))
      xyw  <- terra::xyFromCell(cwin, idx)
      fw   <- as.numeric(vv[idx]); w <- fw / sum(fw)
      c(sum(w * xyw[,1]), sum(w * xyw[,2]))
    }
    xy <- t(apply(xy, 1, function(p) refine_one(p[1], p[2])))
  }
  
  sf::st_as_sf(data.frame(flux = v, x = xy[,1], y = xy[,2]),
               coords = c("x","y"), crs = sf::st_crs(NA))
}
