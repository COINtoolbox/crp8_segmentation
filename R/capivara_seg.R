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
source("utils/le_energy_map.R")
source("utils/umap_energy_map.R")
source("utils/rgb_mask_to_df.R")
source("utils/smart_sum.R")
source("utils/make_rgb.R")

conv2 <- function(x, k) {
  H <- nrow(x); W <- ncol(x)
  kh <- nrow(k); kw <- ncol(k)
  ph <- floor(kh/2); pw <- floor(kw/2)
  y <- matrix(0, H, W)
  k <- k[seq(kh,1), seq(kw,1)]  # flip kernel
  xpad <- matrix(0, H+2*ph, W+2*pw)
  xpad[(ph+1):(ph+H), (pw+1):(pw+W)] <- x
  for (i in 1:H)
    for (j in 1:W)
      y[i,j] <- sum(xpad[i:(i+kh-1), j:(j+kw-1)] * k)
  y
}

source("Segmentantion_functions/segment_sobel.R")
source("Segmentantion_functions/segment_hdbscan.R")

Xfits <- FITSio::readFITS("..//data/raw/datacube_reg4.fits")
cube  <- Xfits$imDat
H <- dim(cube)[1]; W <- dim(cube)[2]

df_mat <- cube_to_matrix(Xfits)
#P  <- pca_energy_map(df_mat, H, W, d = 10)
P <- smart_sum(cube, method="localivar", bg_q=0.5)

#P <- copula_energy(cube)
#P <- umap_energy_map(df_mat, H, W, d = 5, n_neighbors = 20, min_dist = 0.1,
#                     xy_weight = 0.25, energy = "mahalanobis")

# Child segmentation → clean → fill
L_edge_groups <- segment_sobel(
  P, edge_q = 0.45, drop_frame = TRUE,
  close_k = 3L, close_iters = 1L,
  conn = 8L, min_size = 10L
)



#image(L_edge_groups > 0L, col = c("#fff7c7", "#6a0027"))  # quick visual

cube_rgb <- make_rgb(cube, r=7, g=4, b=2,
                     pansharpen=0.5, guide_band=2,
                     upscale=1, unsharp_sigma=1.1, unsharp_amount=0.7,
                     sat=0.9, gamma=1.0)


df <- rgb_mask_to_df(cube_rgb, L_edge_groups, rotate_ccw = TRUE)
ggplot(df, aes(x = -x_plot, y = y_plot, fill = hex)) +
  geom_raster() +
  scale_fill_identity() +
  coord_equal() +
  theme_void()



# C) Apply before your HDBSCAN
L_child <- segment_hdbscan(P, q_fore = 0.875, scale_xy = 1, scale_I = 1, minPts = 25)



cube_na <- mask_cube(cube,L_edge_groups , mode = "na")

cub_cut <- cube_na

cube_cap <- list(imDat = cub_cut)   # cube_2 is your [nx,ny,nb] array
seg <- capivara::segment(cube_cap,N=50)
pal <- viridis(50)

# make group 1 white
pal[1:3] <- "#FFFFFF"

image(seg$cluster_map,col=pal)


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



# Define the directory where your throughput files are located
# Example: local path
library(dplyr)
library(purrr)
library(readr)
library(ggplot2)
library(stringr)

path <- "/Users/rd23aag/Downloads/nircam_throughputs/mean_throughputs/"

# List all files in the directory
files <- list.files(path, pattern = "_mean_system_throughput\\.txt$", full.names = TRUE)

# Filter only the ones that match your chosen filters
files_filtered <- files[str_detect(files, paste(filters, collapse="|"))]

# Function to read and clean each file
read_filter <- function(f) {
  df <- read.table(f, header = TRUE)
  filter_name <- str_extract(basename(f), "^F[0-9]+[A-Z]*[0-9]*[A-Z]*")
  df <- df %>%
    rename(wavelength = 1, throughput = 2) %>%
    mutate(filter = filter_name,
           lambda_c = nircam_lambda_um[filter_name])
  return(df)
}

# Combine into one dataframe
throughput_all <- map_dfr(files_filtered, read_filter)
throughput_all  <- filter(throughput_all,filter!="F150W2")

# Inspect
head(throughput_all)

ggplot(throughput_all, aes(x = wavelength, y = throughput, color = filter)) +
  geom_line(size = 0.8) +
  theme_minimal(base_size = 14) +
  labs(x = "Wavelength (μm)",
       y = "System Throughput",
       color = "Filter") +
  scale_color_viridis_d(option = "turbo")

region_id <- 6  # <- choose the group/region you want

nircam_lambda_um <- c(
  F090W = 0.902, F115W = 1.154, F150W = 1.501, F182M = 1.842,
  F200W = 1.989, F210M = 2.099, F277W = 2.770, F335M = 3.365,
  F356W = 3.563, F410M = 4.082, F430M = 4.286, F444W = 4.421,
  F480M = 4.828
)

sed <- reg$flux_long %>%
  filter(region == region_id) %>%
  mutate(
    lambda_um = nircam_lambda_um[band],
    band = factor(band, levels = names(sort(nircam_lambda_um[names(nircam_lambda_um) %in% band])))
  ) %>%
  arrange(lambda_um)

ggplot(sed, aes(x = lambda_um, y = flux)) +
#  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(flux - flux_err, 0), ymax = flux + flux_err), width = 0) +
#  scale_x_continuous("Wavelength (µm)", breaks = c(0.9, 1.15, 1.5, 1.84, 2.0, 2.1, 2.77, 3.37, 3.56, 4.08, 4.29, 4.42, 4.83)) +
  scale_y_continuous("Flux", labels = comma) +
  theme_bw() 




throughput_keep <- throughput_all %>%
  filter(filter %in% filters)

# ---- normalize the SED flux
sed_norm <- sed %>%
  mutate(flux_norm = flux / max(flux, na.rm = TRUE),
         flux_err_norm = flux_err / max(flux, na.rm = TRUE))

# ---- normalize throughput curves (they usually already go 0–1)
throughput_norm <- throughput_all %>%
  filter(filter %in% names(nircam_lambda_um)) %>%
  mutate(throughput_scaled = throughput )   # scale them to sit at bottom (15% height)

# ---- plot normalized SED with filter bands
ggplot() +
  # filters as filled curves slightly below the SED
  geom_area(
    data = throughput_norm,
    aes(x = wavelength, y = 100*throughput_scaled, fill = filter),
    alpha = 0.3,                      # transparency of filled region
    position = "identity",  
    color = NA                        # no border for the filled area
  ) +
#  geom_line(data = throughput_norm,
#            aes(x = wavelength, y = throughput_scaled, color = filter),
#            linewidth = 0.8, inherit.aes = FALSE) +
  
  # normalized SED points + errors
 geom_errorbar(data = sed_norm,
                aes(x = lambda_um,
                    ymin = pmax(flux - flux_err, 0),
                    ymax = flux + flux_err),
                width = 0.05) +
  geom_point(data = sed_norm,
             aes(x = lambda_um, y = flux),
             size = 2,shape=21,fill="orange2") +
  
  scale_y_continuous("Flux") +
  scale_color_viridis_d(option = "turbo", name = "Filter") +
  scale_fill_viridis_d(option = "turbo", name = "Filter") +
  labs(x = "Wavelength (μm)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none",
        legend.title = element_text(face = "bold"),
        panel.grid = element_blank(),  )
