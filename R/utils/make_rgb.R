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
