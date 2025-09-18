library(FITSio)
library(imager)
library(ggplot2)
library(png)

# ---------- Helpers ----------
safe_q <- function(x, p){
  x <- x[is.finite(x)]
  if (!length(x)) return(0)
  as.numeric(quantile(x, p, na.rm = TRUE))
}
clp <- function(z) pmax(pmin(z, 1), 0)

stretch_builder <- function(qlo, qhi, Q = 6, eps = 1e-12){
  rng <- max(qhi - qlo, eps)
  function(x){
    v <- (pmin(pmax(x, qlo), qhi) - qlo) / rng
    asinh(Q * v) / asinh(Q)
  }
}

gray_world <- function(R,G,B, mask, eps = 1e-12){
  m <- function(z) mean(z[mask], na.rm = TRUE)
  sR <- m(R) + eps; sG <- m(G) + eps; sB <- m(B) + eps
  s  <- (sR + sG + sB) / 3
  list(R = R * (s/sR), G = G * (s/sG), B = B * (s/sB))
}

psf_match_optional <- function(img, sigma = 0){
  if (sigma <= 0) return(img)
  if (!requireNamespace("imager", quietly = TRUE)) stop("imager not installed.")
  as.matrix(imager::isoblur(imager::as.cimg(img), sigma = sigma))
}

make_mask <- function(ref, q = 0.30){
  th <- safe_q(ref, q)
  is.finite(ref) & ref > th
}

# ---------- Main ----------
make_rgb <- function(cube,
                     r = 10, g = 7, b = 2,                  # bands if method="bands"
                     method = c("bands","nmf"),
                     nmf_bands = NULL, nmf_rank = 3,        # NMF params
                     ref_band = NULL,                       # mask reference band
                     mask_q = 0.30, sky_q = 0.05,           # mask & bg percentiles
                     Q = 6, qlo_p = 0.02, qhi_p = 0.99,     # stretch params
                     psf_sigma = 0,                         # optional PSF equalization
                     sat = 0.95, gamma = 1.0,               # look controls
                     wavelengths = NULL, seed = 1L){
  
  stopifnot(length(dim(cube)) == 3)
  H <- dim(cube)[1]; W <- dim(cube)[2]; Lbands <- dim(cube)[3]
  
  method <- match.arg(method)
  
  # --- choose reference for mask ---
  if (is.null(ref_band)) ref_band <- g
  ref  <- cube[,,ref_band]
  
  # optional PSF match (applied to all 3 for "bands" path)
  get_band <- function(idx){
    x <- cube[,,idx]
    if (psf_sigma > 0) x <- psf_match_optional(x, psf_sigma)
    x
  }
  
  # --- Build RGB sources (pre-stretch) ---
  if (method == "bands"){
    R0 <- get_band(r); G0 <- get_band(g); B0 <- get_band(b)
  } else {
    # ------- NMF path -------
    if (is.null(nmf_bands)) nmf_bands <- seq_len(Lbands)
    Vraw <- apply(cube[,,nmf_bands], 3, c)        # (H*W) x n_bands
    
    # Mask from reference (not PSF-matched for speed)
    mask <- make_mask(ref, q = mask_q)
    sky  <- is.finite(ref) & !mask
    
    # sky subtraction per band
    for (j in seq_along(nmf_bands)){
      bj <- safe_q(cube[,,nmf_bands[j]][sky], sky_q)
      Vraw[, j] <- pmax(Vraw[, j] - bj, 0)
    }
    
    # optional per-pixel normalization to focus on chroma
    Lum_pix <- rowSums(Vraw) + 1e-12
    V <- Vraw / Lum_pix
    
    # NMF
    if (!requireNamespace("NMF", quietly = TRUE)) stop("NMF package not installed.")
    set.seed(seed)
    fit <- NMF::nmf(V, rank = nmf_rank, method = "brunet", nrun = 5)
    Wmat <- NMF::basis(fit)    # (H*W) x k
    Hmat <- NMF::coef(fit)     # k x n_bands
    
    # Order components by wavelength centroid (blue->red)
    if (is.null(wavelengths)) {
      λ <- seq_along(nmf_bands)
    } else {
      λ <- wavelengths[ nmf_bands ]
    }
    # safeguard against zeros
    den <- rowSums(Hmat) + 1e-12
    centroid <- as.numeric((Hmat %*% λ) / den)  # length k
    ord <- order(centroid)
    Wmat <- Wmat[, ord, drop = FALSE]
    
    # Choose 3 components (first/median/last by centroid if k>=3)
    k <- ncol(Wmat)
    pick <- if (k >= 3) c(1, floor((k+1)/2), k) else seq_len(k)
    W3 <- Wmat[, pick, drop = FALSE]
    
    # Back to images (B,G,R from low->mid->high)
    B0 <- matrix(W3[,1], H, W)
    G0 <- if (ncol(W3) >= 2) matrix(W3[,2], H, W) else matrix(0, H, W)
    R0 <- if (ncol(W3) >= 3) matrix(W3[,3], H, W) else matrix(0, H, W)
    
    # Equalize PSF if requested
    if (psf_sigma > 0){
      B0 <- psf_match_optional(B0, psf_sigma)
      G0 <- psf_match_optional(G0, psf_sigma)
      R0 <- psf_match_optional(R0, psf_sigma)
    }
  }
  
  # --- Mask & Sky using (possibly PSF-matched) reference ---
  mask <- make_mask(ref, q = mask_q)
  sky  <- is.finite(ref) & !mask
  
  # --- Background neutralization (use sky!) ---
  bgR <- safe_q(R0[sky], sky_q)
  bgG <- safe_q(G0[sky], sky_q)
  bgB <- safe_q(B0[sky], sky_q)
  R1 <- pmax(R0 - bgR, 0); G1 <- pmax(G0 - bgG, 0); B1 <- pmax(B0 - bgB, 0)
  
  # --- Common stretch from luminance (same qlo/qhi) ---
  L0   <- 0.30*R1 + 0.59*G1 + 0.11*B1
  qlo  <- safe_q(L0[mask], qlo_p)
  qhi  <- safe_q(L0[mask], qhi_p)
  S    <- stretch_builder(qlo, qhi, Q = Q)
  R <- S(R1); G <- S(G1); B <- S(B1)
  
  # --- Gray-world WB ---
  WB <- gray_world(R, G, B, mask = mask); R <- WB$R; G <- WB$G; B <- WB$B
  
  # --- Optional global desaturation + gamma ---
  Lum <- 0.299*R + 0.587*G + 0.114*B
  R <- clp((sat*R + (1 - sat)*Lum)^gamma)
  G <- clp((sat*G + (1 - sat)*Lum)^gamma)
  B <- clp((sat*B + (1 - sat)*Lum)^gamma)
  
  list(R = R, G = G, B = B, mask = mask)
}


make_rgb <- function(cube,
                     r = 12, g = 7, b = 2,                  # bands (can be vectors)
                     ref_band = NULL,                       # mask reference band
                     mask_q = 0.30, sky_q = 0.05,           # mask & bg percentiles
                     Q = 6, qlo_p = 0.02, qhi_p = 0.99,     # stretch params
                     psf_sigma = 0,                         # optional PSF blur (all 3)
                     sat = 0.95, gamma = 1.0,               # look controls
                     pansharpen = 0.0, guide_band = 2,      # SR-like luminance inject
                     upscale = 1L,                          # 1=no upsample; 2,3,...
                     unsharp_sigma = 1.2, unsharp_amount = 0.0) {  # 0 disables sharpen
  
  stopifnot(length(dim(cube)) == 3)
  H <- dim(cube)[1]; W <- dim(cube)[2]
  
  # ---- helpers ----
  safe_q <- function(x, p){
    x <- x[is.finite(x)]; if (!length(x)) return(0)
    as.numeric(stats::quantile(x, p, na.rm=TRUE))
  }
  # PSF blur that preserves HxW
  psf_match_optional <- function(img, sigma=0) {
    if (sigma <= 0) return(img)
    cim <- imager::as.cimg(aperm(array(img, dim=c(H,W,1)), c(2,1,3)))  # x,y,cc
    bl  <- imager::isoblur(cim, sigma = sigma)
    arr <- as.array(bl)                             # x,y,z=1,c=1
    aperm(arr, c(2,1,3,4))[,,1,1, drop=TRUE]        # back to HxW matrix
  }
  make_mask <- function(ref, q=0.30){
    thr <- safe_q(ref, q); is.finite(ref) & ref > thr
  }
  stretch_builder <- function(qlo, qhi, Q=6){
    force(qlo); force(qhi); force(Q)
    function(x){
      v <- pmin(pmax(x, qlo), qhi); v <- (v - qlo)/(qhi - qlo + 1e-12)
      asinh(Q*v)/asinh(Q)
    }
  }
  gray_world <- function(R,G,B, mask){
    mw <- function(z) mean(z[mask], na.rm=TRUE)
    sR <- mw(R); sG <- mw(G); sB <- mw(B); s <- (sR+sG+sB)/3
    list(R = R*(s/(sR+1e-12)),
         G = G*(s/(sG+1e-12)),
         B = B*(s/(sB+1e-12)))
  }
  clp <- function(z) pmax(pmin(z,1),0)
  
  # average multiple bands if vectors are provided
  mix_bands <- function(idx_vec){
    idx_vec <- as.integer(idx_vec)
    if (length(idx_vec)==1L) return(cube[,,idx_vec])
    rowMeans(aperm(cube[,,idx_vec, drop=FALSE], c(3,1,2)), dims=1)
  }
  
  # pick bands (+ optional PSF equalization)
  R0 <- mix_bands(r); G0 <- mix_bands(g); B0 <- mix_bands(b)
  if (psf_sigma > 0) {
    R0 <- psf_match_optional(R0, psf_sigma)
    G0 <- psf_match_optional(G0, psf_sigma)
    B0 <- psf_match_optional(B0, psf_sigma)
  }
  
  # reference for mask
  if (is.null(ref_band)) ref_band <- if (length(g)) g[1] else 1
  ref <- mix_bands(ref_band)
  
  # mask & sky
  mask <- make_mask(ref, q = mask_q)
  sky  <- is.finite(ref) & !mask
  
  # background neutralization (use sky)
  bgR <- safe_q(R0[sky], sky_q)
  bgG <- safe_q(G0[sky], sky_q)
  bgB <- safe_q(B0[sky], sky_q)
  R1 <- pmax(R0 - bgR, 0); G1 <- pmax(G0 - bgG, 0); B1 <- pmax(B0 - bgB, 0)
  
  # common stretch (same qlo/qhi for all)
  L0   <- 0.30*R1 + 0.59*G1 + 0.11*B1
  qlo  <- safe_q(L0[mask], qlo_p)
  qhi  <- safe_q(L0[mask], qhi_p); if (qhi <= qlo) qhi <- qlo + 1e-6
  S    <- stretch_builder(qlo, qhi, Q = Q)
  R <- S(R1); G <- S(G1); B <- S(B1)
  
  # gray-world white balance
  WB <- gray_world(R,G,B, mask); R <- WB$R; G <- WB$G; B <- WB$B
  
  # SR-like pan-sharpen (inject short-λ detail into luminance), color-safe
  if (pansharpen > 0) {
    guide0 <- mix_bands(guide_band)
    guide0 <- pmax(guide0 - safe_q(guide0[sky], sky_q), 0)
    Gd <- S(guide0)
    Lum <- 0.299*R + 0.587*G + 0.114*B
    Lum_new <- (1 - pansharpen)*Lum + pansharpen*Gd
    s <- Lum_new / (Lum + 1e-9)
    R <- R*s; G <- G*s; B <- B*s
  }
  
  # desaturation + gamma
  Lum <- 0.299*R + 0.587*G + 0.114*B
  R <- clp((sat*R + (1 - sat)*Lum)^gamma)
  G <- clp((sat*G + (1 - sat)*Lum)^gamma)
  B <- clp((sat*B + (1 - sat)*Lum)^gamma)
  
  # compose RGB
  out <- array(0, dim=c(H,W,3)); out[,,1] <- R; out[,,2] <- G; out[,,3] <- B
  
  # ---- Fast path: no upsample & no unsharp ----
  if (upscale <= 1L && unsharp_amount <= 0) {
    return(list(RGB = out, R = R, G = G, B = B,
                mask = mask, qlo = qlo, qhi = qhi))
  }
  
  # ---- Upscale / Unsharp path ----
  img  <- imager::as.cimg(aperm(out, c(2,1,3)))   # x,y,cc
  if (upscale > 1L) img <- imager::imresize(img, scale = upscale)
  if (unsharp_amount > 0) {
    blur <- imager::isoblur(img, sigma = unsharp_sigma)
    img  <- img + unsharp_amount*(img - blur)      # sharpen
  }
  # clamp to [0,1] without imclip()
  arr0 <- as.array(img)
  arr0 <- pmax(pmin(arr0, 1), 0)
  
  # handle 3D or 4D array return from imager
  if (length(dim(arr0)) == 4) {
    # x,y,z,c -> H,W,c
    arr <- aperm(arr0, c(2,1,4,3))[,,1:3,1, drop=TRUE]
  } else if (length(dim(arr0)) == 3) {
    # x,y,c -> H,W,c
    arr <- aperm(arr0, c(2,1,3))
  } else stop("Unexpected number of dimensions from imager.")
  
  list(RGB = arr, R = arr[,,1], G = arr[,,2], B = arr[,,3],
       mask = mask, qlo = qlo, qhi = qhi)
}





# Load cube
# Paths
files <- sprintf("..//..//data/raw/datacube_reg%d.fits", 1:4)

# Read all FITS and extract the image data
cubes <- lapply(files, function(f) FITSio::readFITS(f)$imDat)


mk <- function(cube) make_rgb(cube, r=9, g=6, b=1,
                              pansharpen=0.35, guide_band=2,
                              upscale=10, unsharp_sigma=1.1, unsharp_amount=0.6,
                              sat=0.8, gamma=1.0)
rgbs <- lapply(cubes, mk)

p1 <- plot_rgb_gg(rgbs[[1]])
p2 <- plot_rgb_gg(rgbs[[2]])
p3 <- plot_rgb_gg(rgbs[[3]])
p4 <- plot_rgb_gg(rgbs[[4]])
p5 <- plot_rgb_gg(rgbs[[5]])
p6 <- plot_rgb_gg(rgbs[[6]])
p7 <- plot_rgb_gg(rgbs[[7]])


mosaic <- (p1 | p2 | p3 ) /
  ( p4| p5 | p6 | p7 ) & theme(
    panel.spacing = grid::unit(0, "pt"),           # no gaps
    plot.margin   = grid::unit(c(0,0,0,0), "pt")   # no outer margin
  )

ragg::agg_png("mosaic_8000x5000.png", width = 8000, height = 5000, units = "px", background = "white")
print(mosaic)
dev.off()
