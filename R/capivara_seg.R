# ==============================================================================
# Hyperspectral Clump Segmentation & Regional Photometry
# ==============================================================================
# Overview
#   End-to-end pipeline that (1) builds a PCA-based detection image from a
#   hyperspectral/photometric cube, (2) extracts “child” clumps via HDBSCAN,
#   (3) cleans masks (hole fill per label), (4) merges children into PARENT
#   super-regions (no nesting), (5) creates cutouts, (6) refines labels with
#   Capivara (spectro-spatial segmentation), and (7) summarizes fluxes per
#   region and band.
#
# Inputs
#   - cube   : 3-D array [nx, ny, nb] or FITSio-style list ($imDat, $axDat)
#   - (optional) initial_mask : logical/int matrix [nx, ny] to restrict search
#
# Outputs (written to <outdir>/)
#   - labels_children.fits : child clumps (dense substructures)
#   - labels_parent.fits   : merged super-regions (non-overlapping)
#   - labels_capivara.fits : final label map after Capivara refinement
#   - cutouts/             : per-region FITS/PNG cutouts for QA
#   - flux_long.csv        : tidy photometry (region, band, flux, flux_err, n_eff)
#   - flux_wide.csv        : wide photometry (one row per region)
#
# Pipeline (high level)
#   1) PCA detection image
#      - Flatten cube -> [n_pix × nb]; robust-scale bands
#      - PCA (e.g., prcomp); detection = positive combo of top PCs (or PC1)
#   2) Density clustering (HDBSCAN)
#      - dbscan::hdbscan on (x, y, detection) features
#      - minPts controls minimum clump size & density sensitivity
#   3) Morphological cleanup
#      - Fill holes per label; remove single-pixel islands; optional opening/closing
#   4) Merge -> PARENT super-regions
#      - Union adjacent/overlapping children; prevent nested labels
#   5) Cutouts
#      - Write per-region cutouts for rapid visual QC (PNG/FITS)
#   6) Capivara refinement
#      - Run Capivara on the cube to align spectro-spatial structure; reconcile
#        with parent map to produce the final label field
#   7) Regional photometry
#      - RegionPhotometry(cube, labels_final, …) → flux tables + painted cube
#      - Uncertainties: variance cube/sigma per band, or fallback estimators
#
# Key parameters (tune per dataset)
#   - n_pca_comp     : 1–3 usually suffice for detection
#   - minPts         : HDBSCAN minimum cluster size (controls fragmentation)
#   - min_region_pix : post-processing floor on region area
#   - smooth_sigma   : optional Gaussian smoothing on detection image
#
# Conventions & caveats
#   - Orientation: ensure label matrix aligns with cube; rotate/flip if needed.
#   - Pixels with label <= 0 are treated as sky/ignored by photometry.
#   - Units: background/variance must match cube units.
#
# Quick start
#   labels_final <- run_pipeline(cube, outdir = "out", minPts = 20, n_pca_comp = 2)
#   reg <- RegionPhotometry(cube, labels_final, error_fallback = "flux_over_sqrt_n")
#   readr::write_csv(reg$flux_long, file.path("out", "flux_long.csv"))
#   readr::write_csv(reg$flux_wide, file.path("out", "flux_wide.csv"))
# ==============================================================================

library(FITSio)
library(dbscan)
source("utils/fill_holes_per_label.R")
source("utils/filter_by_size.R")
source("utils/mask_cube.R")
source("utils/RegionPhotometry.R")
source("utils/pca_energy_map.R")
source("utils/pca_energy_map.R")
source("utils/build_features.R")
source("utils/compact_labels.R")
source("utils/cube_to_matrix.R")
source("Segmentantion_functions/segment_hdbscan.R")

Xfits <- FITSio::readFITS("..//data/datacube_reg1_Jy.fits")
cube  <- Xfits$imDat
H <- dim(cube)[1]; W <- dim(cube)[2]

df_mat <- cube_to_matrix(Xfits)
P  <- pca_energy_map(df_mat, H, W, d = 2)
# Child segmentation → clean → fill
L_child <- segment_hdbscan(P, q_fore = 0.9, scale_xy = 1.0, 
                           scale_I = 2.0, minPts = 30)
L_child2  <- L_child |> filter_by_size(min_size = 25) |> fill_holes_per_label()
image(L_child2,col=viridis(100))
cub_cut <- mask_cube(cube,L_child2,mode="na")
image(cub_cut[,,1],col=viridis(100))
cube_cap <- list(imDat = cub_cute)   # cube_2 is your [nx,ny,nb] array
seg <- capivara::segment(cube_cap,N=30)
image(seg$cluster_map,col=viridis(30))
reg <- RegionPhotometry(cube, seg$cluster_map, 
                        error_fallback = "flux_over_sqrt_n")



# ---- mapping 1..13 -> JWST filters ----
filters <- c("F090W","F115W","F150W","F182M","F200W",
             "F210M","F277W","F335M","F356W","F410M",
             "F430M","F444W","F480M")
map_num2filt <- setNames(filters, as.character(seq_along(filters)))

# LONG: rename numeric bands -> JWST filters
reg$flux_long <- reg$flux_long %>%
  dplyr::mutate(band = {
    b   <- trimws(as.character(band))
    idx <- suppressWarnings(as.integer(b))
    keep <- !is.na(idx) & idx >= 1 & idx <= length(filters)
    out  <- b
    out[keep] <- map_num2filt[as.character(idx[keep])]
    out
  })

# WIDE: rename columns "1", "1_err", "1_n_eff" -> "F090W", "F090W_err", ...
rename_wide_by_filters <- function(df, filters) {
  old <- colnames(df)
  m   <- regexpr("(_n_eff|_err)$", old)
  suf  <- ifelse(m > 0, regmatches(old, m), "")
  base <- trimws(ifelse(m > 0, substring(old, 1, m - 1), old))
  
  idx     <- suppressWarnings(as.integer(base))
  is_band <- !is.na(idx) & idx >= 1 & idx <= length(filters)
  newbase <- base
  newbase[is_band] <- filters[idx[is_band]]
  
  new <- paste0(newbase, suf)
  if (any(duplicated(new))) new <- make.unique(new, sep = "_")
  colnames(df) <- new
  df
}
reg$flux_wide <- rename_wide_by_filters(reg$flux_wide, filters)

# Recompute errors as |flux|/sqrt(n_pix) and DROP *_n_eff columns
recompute_err_drop_neff <- function(df, filters) {
  stopifnot("n_pix" %in% names(df))
  for (f in filters) {
    if (f %in% names(df)) {
      df[[paste0(f, "_err")]] <- abs(df[[f]]) / sqrt(pmax(df$n_pix, 1))
    }
  }
  df[, !grepl("_n_eff$", names(df)), drop = FALSE]
}
reg$flux_wide <- recompute_err_drop_neff(reg$flux_wide, filters)

# Order nicely: per band (flux, err) only
reorder_by_band <- function(df, filters) {
  parts <- unlist(lapply(filters, function(f) c(f, paste0(f, "_err"))))
  keep  <- c("region", "n_pix", parts)
  dplyr::select(df, dplyr::any_of(keep))
}
flux_wide_out <- reorder_by_band(reg$flux_wide, filters)

# ---- salvar ----
outdir <- normalizePath("./../outdir", mustWork = FALSE)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
readr::write_csv(flux_wide_out, file.path(outdir, "flux_wide_capivara_reg1.csv"))






