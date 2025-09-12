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




image(cube_2[,,10])


cube_2

cube_cap <- list(imDat = cube_2)   # cube_2 is your [nx,ny,nb] array
seg <- capivara::segment(cube_cap,N=15)




regflux <- RegionPhotometry(seg$original_cube,seg$cluster_map,error_fallback = "flux_over_sqrt_n")