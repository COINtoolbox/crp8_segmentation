# ----------------------------
# Gradient Pattern Analysis (GPA) in base R
# ----------------------------

# --- small helpers ---
conv2 <- function(img, ker) {
  # zero-padded 2D convolution (valid for modest image sizes)
  kx <- nrow(ker); ky <- ncol(ker)
  px <- floor(kx/2); py <- floor(ky/2)
  nr <- nrow(img); nc <- ncol(img)
  pad <- matrix(0, nr + 2*px, nc + 2*py)
  pad[(px+1):(px+nr), (py+1):(py+nc)] <- img
  out <- matrix(0, nr, nc)
  for (i in 1:nr) {
    for (j in 1:nc) {
      region <- pad[i:(i+kx-1), j:(j+ky-1)]
      out[i, j] <- sum(region * ker)
    }
  }
  out
}

box_blur <- function(x, k = 3L) {
  k <- max(1L, as.integer(k))
  ker <- matrix(1, k, k) / (k*k)
  conv2(x, ker)
}

entropy <- function(p) {
  p <- p[p > 0]
  -sum(p * log(p))
}

# --- main GPA function ---
gpa <- function(img,
                bins_theta = 36L,        # orientation histogram bins
                blur_k = 5L,             # smoothing for structure tensor
                edge_q = 0.90,           # quantile threshold for "strong" edges
                fract_scales = c(2,3,4,6,8,12)) {
  if (!is.matrix(img)) img <- as.matrix(img)
  img <- matrix(as.numeric(img), nrow(img), ncol(img))
  
  # Sobel kernels
  Gxk <- matrix(c(-1,0,1,-2,0,2,-1,0,1), 3, 3, byrow = TRUE)
  Gyk <- matrix(c(-1,-2,-1,0,0,0,1,2,1), 3, 3, byrow = TRUE)
  
  Gx <- conv2(img, Gxk)
  Gy <- conv2(img, Gyk)
  
  mag   <- sqrt(Gx^2 + Gy^2)
  theta <- atan2(Gy, Gx)                 # range [-pi, pi]
  
  # Orientation histogram + entropy (Shannon)
  brks <- seq(-pi, pi, length.out = bins_theta + 1L)
  h    <- hist(theta, breaks = brks, plot = FALSE)
  p    <- h$counts / sum(h$counts)
  Htheta <- entropy(p)                   # in nats
  Htheta_bits <- Htheta / log(2)
  
  # Structure tensor J = [[<Gx^2>, <GxGy>], [<GxGy>, <Gy^2>]]
  J11 <- box_blur(Gx^2, blur_k)
  J22 <- box_blur(Gy^2, blur_k)
  J12 <- box_blur(Gx*Gy, blur_k)
  
  # eigenvalues per-pixel (closed form for 2x2)
  tr  <- J11 + J22
  det <- J11*J22 - J12^2
  disc <- pmax(0, tr^2/4 - det)
  lam1 <- tr/2 + sqrt(disc)
  lam2 <- tr/2 - sqrt(disc)
  
  # Coherence (anisotropy) per-pixel and global
  coh <- (lam1 - lam2) / pmax(1e-12, lam1 + lam2)
  coh[!is.finite(coh)] <- 0
  mean_coh <- mean(coh, na.rm = TRUE)
  median_coh <- median(coh, na.rm = TRUE)
  
  # Basic gradient stats
  mag_stats <- c(mean = mean(mag), sd = sd(mag), median = median(mag),
                 q90 = as.numeric(quantile(mag, 0.90)),
                 q95 = as.numeric(quantile(mag, 0.95)),
                 max = max(mag))
  
  # Edge density above quantile threshold
  thr <- as.numeric(quantile(mag, edge_q))
  strong_edges <- mag >= thr
  edge_density <- mean(strong_edges)
  
  # Simple directional balance (how uniform orientations are)
  # resultant vector length R = |E[exp(i*theta)]|
  R <- sqrt(mean(cos(theta))^2 + mean(sin(theta))^2)
  
  # Very simple box-count fractal proxy on strong edge mask:
  # count occupied boxes vs scale and fit slope ~ 1/D
  box_count <- function(mask, s) {
    nr <- nrow(mask); nc <- ncol(mask)
    r <- floor(nr / s) * s
    c <- floor(nc / s) * s
    M <- mask[1:r, 1:c]
    # reshape into s x s blocks and check any TRUE in each block
    B <- matrix(M, nrow = s)
    # Rearrangement trick: aggregate by block via max
    # We'll loop blocks along rows and cols:
    cnt <- 0L
    for (i in seq(1, r, by = s)) {
      for (j in seq(1, c, by = s)) {
        if (any(M[i:(i+s-1), j:(j+s-1)])) cnt <- cnt + 1L
      }
    }
    cnt
  }
  
  nboxes <- sapply(fract_scales, function(s) box_count(strong_edges, s))
  # Fit log N(s) vs log(1/s) -> slope ~ D_box (proxy on edge set)
  xs <- log(1 / fract_scales)
  ys <- log(pmax(1, nboxes))
  fit <- lm(ys ~ xs)
  D_box_proxy <- unname(coef(fit)[2])
  
  list(
    gradients = list(Gx = Gx, Gy = Gy, magnitude = mag, theta = theta),
    orientation_hist = list(breaks = brks, counts = h$counts, prob = p,
                            H_nats = Htheta, H_bits = Htheta_bits),
    coherence = list(mean = mean_coh, median = median_coh, field = coh),
    magnitude_stats = mag_stats,
    edge = list(threshold = thr, density = edge_density, mask = strong_edges),
    directional_resultant_R = R,
    fractal_proxy = list(scales = fract_scales,
                         nboxes = nboxes,
                         D_box_proxy = D_box_proxy,
                         fit = fit)
  )
}

# ----------------------------
# Minimal demo
# ----------------------------
if (sys.nframe() == 0) {
  set.seed(42)
  # synthetic field: smooth blob + oriented ridge + noise
  n <- 128
  x <- matrix(0, n, n)
  # radial Gaussian
  cx <- cy <- n/2
  for (i in 1:n) for (j in 1:n) {
    r2 <- (i - cx)^2 + (j - cy)^2
    x[i, j] <- exp(-r2/(2*(n/6)^2))
  }
  # add an oriented line
  for (k in 1:n) {
    i <- k
    j <- floor(0.6*k + 20)
    if (j >= 1 && j <= n) x[i, j] <- x[i, j] + 0.8
  }
  # smooth and add noise
  x <- box_blur(x, 7)
  x <- x + 0.10 * matrix(rnorm(n*n), n, n)
  
  res <- gpa(x)
  
  cat("\n--- GPA summary ---\n")
  print(res$magnitude_stats)
  cat(sprintf("Orientation entropy: %.3f bits\n", res$orientation_hist$H_bits))
  cat(sprintf("Global coherence (mean): %.3f\n", res$coherence$mean))
  cat(sprintf("Edge density (q=%.2f): %.3f\n", 0.90, res$edge$density))
  cat(sprintf("Directional resultant R: %.3f (1=aligned, 0=uniform)\n", res$directional_resultant_R))
  cat(sprintf("Fractal proxy (box-count slope): %.3f\n", res$fractal_proxy$D_box_proxy))
}

# ----------------------------
# Optional quick plots (base)
# ----------------------------
plot_gpa_quick <- function(res) {
  oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar))
  par(mfrow = c(2,3), mar = c(3,3,2,1))
  image(t(apply(res$gradients$magnitude, 2, rev)), axes = FALSE, main = "Grad magnitude")
  image(t(apply(res$coherence$field, 2, rev)), axes = FALSE, main = "Coherence")
  # Orientation histogram
  mids <- head(res$orientation_hist$breaks, -1) + diff(res$orientation_hist$breaks)/2
  plot(mids, res$orientation_hist$prob, type = "h", xlab = "theta (rad)", ylab = "p", main = "Orientation histogram")
  # Strong edge mask
  image(t(apply(res$edge$mask, 2, rev)), axes = FALSE, main = "Strong edge mask")
  # Box-count fit
  xs <- log(1 / res$fractal_proxy$scales); ys <- log(pmax(1, res$fractal_proxy$nboxes))
  plot(xs, ys, pch = 19, main = "Box-count fit", xlab = "log(1/scale)", ylab = "log(N)")
  abline(res$fractal_proxy$fit, lty = 2)
  # Text panel
  plot.new(); title("Summary")
  text(0.05, 0.9, sprintf("H(theta): %.2f bits", res$orientation_hist$H_bits), adj = 0)
  text(0.05, 0.8, sprintf("Coherence: %.2f", res$coherence$mean), adj = 0)
  text(0.05, 0.7, sprintf("Edge density: %.2f", res$edge$density), adj = 0)
  text(0.05, 0.6, sprintf("R (dir. resultant): %.2f", res$directional_resultant_R), adj = 0)
  text(0.05, 0.5, sprintf("Fractal proxy slope: %.2f", res$fractal_proxy$D_box_proxy), adj = 0)
}


Xfits <- FITSio::readFITS("..//..//data/raw/datacube_reg1.fits")
cube  <- Xfits$imDat

# img is a numeric matrix (e.g., grayscale image or 2-D field)
res <- gpa(cube[,,1])

# quick summary
res$magnitude_stats
res$orientation_hist$H_bits
res$coherence$mean
res$edge$density
res$directional_resultant_R
res$fractal_proxy$D_box_proxy

# optional plots
plot_gpa_quick(res)

plot_gradient_field <- function(res, step = 5, scale = 0.5) {
  Gx <- res$gradients$Gx
  Gy <- res$gradients$Gy
  mag <- res$gradients$magnitude
  
  nr <- nrow(Gx); nc <- ncol(Gx)
  xs <- seq(1, nc, by = step)
  ys <- seq(1, nr, by = step)
  
  # downsample field
  gx <- Gx[ys, xs]
  gy <- Gy[ys, xs]
  m  <- mag[ys, xs]
  
  # normalize length for plotting
  gx <- gx / (max(m) + 1e-9) * scale * step
  gy <- gy / (max(m) + 1e-9) * scale * step
  
  # base image of gradient magnitude
  image(t(apply(mag, 2, rev)), col = gray.colors(100),
        axes = FALSE, main = "Gradient field with arrows")
  
  # overlay quiver arrows
  for (i in seq_along(xs)) {
    for (j in seq_along(ys)) {
      x <- xs[i]; y <- ys[j]
      dx <- gx[j, i]; dy <- gy[j, i]
      # convert to image coordinates (flip y for plotting)
      x0 <- (x-1)/(nc-1)
      y0 <- 1 - (y-1)/(nr-1)
      arrows(x0, y0, x0 + dx/nc, y0 - dy/nr,
             length = 0.05, col = "red", lwd = 0.7)
    }
  }
}

# after running res <- gpa(img)
plot_gradient_field(res, step = 6, scale = 2)

