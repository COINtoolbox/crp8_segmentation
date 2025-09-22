# P is an HxW matrix, e.g. from PCA/UMAP/LE/NMF/RPCA
res <- find_clumps(
  input = P,
  origin = "upper",          # or "lower" if row 1 is bottom
  fwhm_px = 2,
  thr_q   = 0.875,
  min_sep_px = 4,
  refine_radius_px = 2,
  make_plot = TRUE,
  sample_points = 0          # set >0 to subsample plot points
)


# quick plot
res$plot


wcs <- make_wcs_from_axDat(Xfits$axDat)

# Sanity check: CRPIX should map to CRVAL exactly
wcs_pix2radec_tan(wcs$crpix[1], wcs$crpix[2], wcs)
cent <- res$centroids


M <- cube[,,4]
H <- nrow(M); 
W <- ncol(M)
pts <- as.data.frame(sf::st_coordinates(cent))  # columns X, Y

# Image with pixel-index axes (and the usual orientation fix for image()):
image(x = 1:W, y = 1:H, z = t(M)[, H:1],
      col = hcl.colors(100), useRaster = TRUE, asp = 1,
      xlab = "x (pixels)", ylab = "y (pixels)",
      zlim = range(M, na.rm = TRUE))  # ignore NAs for color scale

# Overlay the points (no swaps, no flips)
points(pts$X, pts$Y, pch = 1, col = "red", cex = 2, lwd = 1.2)


