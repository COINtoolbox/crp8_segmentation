# ------------------------------------------------------------------------------
# Region Photometry on Labelled Data Cubes
# ------------------------------------------------------------------------------
# Overview
#   Reads (i) an original photometric/spectral cube and (ii) a labelled mask,
#   integrates the flux per region and per band, and writes tidy tables plus an
#   optional “painted” cube to an output folder.
#
# Inputs
#   - cube   : 3-D array [nx, ny, nb], or FITSio-style list with $imDat[,,], and
#              optionally $axDat (used to derive band/wavelength names).
#   - labels : integer matrix [nx, ny]; pixels with value > 0 belong to a region;
#              values <= 0 are treated as sky/ignored.
#
# Optional parameters (if your script forwards them to RegionPhotometry)
#   - bkg          : background to subtract. One of:
#                      * scalar (single value)
#                      * numeric vector length nb (per band)
#                      * matrix [nx, ny]
#                      * cube [nx, ny, nb]
#   - var_cube     : variance cube [nx, ny, nb] (same shape as cube) for errors.
#   - sigma_band   : numeric vector length nb with per-band σ; used if var_cube
#                    is missing (errors = sqrt(sum σ^2 over pixels)).
#   - error_fallback : when no variance is available, estimate uncertainties:
#                      "none" | "flux_over_sqrt_n" | "poisson" | "mad_sky".
#                      (See “Uncertainties” below.)
#   - band_values  : custom band/wavelength labels (length nb).
#   - return_painted_cube : if TRUE, returns a cube where each layer contains
#                           the integrated flux of the corresponding region.
#
# Outputs (written to <outdir>/)
#   - flux_long.csv : tidy table with columns:
#         region | band | flux | flux_err | n_eff | n_pix
#     where:
#       * flux    = sum of pixel values in the region for that band
#       * n_eff   = number of finite pixels used in the sum
#       * n_pix   = total pixels in the region (band-independent)
#   - flux_wide.csv : one row per region; columns per band (plus _n_eff and _err).
#   - labels_<name>.fits (optional) : painted cube with region-integrated fluxes
#                                     replicated across each region.
#
# Method (brief)
#   1) Optional background subtraction.
#   2) Flatten cube to per-pixel spectra; ignore non-finite/label<=0 pixels.
#   3) Per (region × band) sum of flux; track n_eff and, if provided, variance.
#   4) Compute flux_err from variance or the chosen fallback (see below).
#   5) Emit tidy and wide tables; optionally emit painted cube.
#
# Uncertainties
#   - With var_cube or sigma_band: flux_err = sqrt(sum of variances).
#   - error_fallback choices when no variance is supplied:
#       * "flux_over_sqrt_n": flux_err = |flux| / sqrt(n_eff)
#       * "poisson"         : flux_err = sqrt(max(flux, 0))   # counts-like
#       * "mad_sky"         : per-band σ from sky (labels<=0) via MAD, then
#                             flux_err = σ_band * sqrt(n_eff)
#
# Assumptions & caveats
#   - Orientation: labels must align with cube axes [nx, ny]. If your mask looks
#     rotated/flipped w.r.t. the image, rotate/flip the label matrix beforehand.
#   - Units: this script preserves the cube’s native units; ensure background and
#     variances are provided in consistent units.
#   - Missing data: non-finite pixels are ignored; regions with n_eff = 0 yield
#     flux = 0 and flux_err = NA (or fallback estimate if selected).
#



Xfits <- FITSio::readFITS("..//data/datacube_reg1_Jy.fits")
cube  <- Xfits$imDat
Seg_wavelets <- FITSio::readFITS("..//Segmented_regions/Wavelets/datacube_reg1_segmentationmap.fits") 



# ---------- 2-D helpers ----------
rot90cw  <- function(M) t(M[nrow(M):1, , drop = FALSE])          # 90° clockwise
rot90ccw <- function(M) t(M)[, nrow(M):1, drop = FALSE]          # 90° counter-clockwise
rot180   <- function(M) M[nrow(M):1, ncol(M):1, drop = FALSE]
flipV    <- function(M) M[nrow(M):1, , drop = FALSE]             # mirror bottom↔top
flipH    <- function(M) M[, ncol(M):1, drop = FALSE]             # mirror left↔right

# Compose: rotate then optional flips (order: rotate → flips)
# --- 2-D transforms (labels only) ---
rot90cw  <- function(M) t(M[nrow(M):1, , drop = FALSE])          # 90° clockwise
rot90ccw <- function(M) t(M)[, nrow(M):1, drop = FALSE]          # 90° counter-clockwise
rot180   <- function(M) M[nrow(M):1, ncol(M):1, drop = FALSE]
flipV    <- function(M) M[nrow(M):1, , drop = FALSE]             # mirror bottom↔top
flipH    <- function(M) M[, ncol(M):1, drop = FALSE]             # mirror left↔right

# Compose: rotate → optional flips
orient2D <- function(M, rotate = c(0,90,180,270), flip_vert = FALSE, flip_horiz = FALSE) {
  rotate <- as.character(match.arg(as.character(rotate), c("0","90","180","270")))
  R <- switch(rotate, "0"=M, "90"=rot90cw(M), "180"=rot180(M), "270"=rot90ccw(M))
  if (flip_vert)  R <- flipV(R)
  if (flip_horiz) R <- flipH(R)
  R
}

labels_fix <- Seg_wavelets$imDat |> flipV() |> rot90ccw()

cub_cut_wav <- mask_cube(cube,labels_fix,mode="na")



reg_Wavelet <- RegionPhotometry(cube, labels_fix, 
                        error_fallback = "flux_over_sqrt_n")

flux_wide_out_wavelets <- reorder_by_band(reg_Wavelet$flux_wide, filters)

readr::write_csv(flux_wide_out_wavelets,
                 "/Users/rd23aag/Documents/GitHub/crp8_segmentation/outdir/flux_wide_wavelets_reg1.csv")




