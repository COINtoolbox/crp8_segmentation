# ============================
# Gradient Pattern Analysis (fast) + Quiver plots
# ============================

suppressPackageStartupMessages({
  library(ggplot2)
  # viridis is included via ggplot2's scale_*_viridis_c;
  # if missing, install.packages("viridis")
})

# -----------------------------
# FAST GPA (asymmetric) — vectorized + binned pruning
# -----------------------------
gpa_asym_fast <- function(mat,
                          cx = NULL, cy = NULL, r = NULL,
                          mtol = 0.05,        # modulus tol (fraction of max |∇|)
                          ftol = 0.05,        # orientation tolerance (radians) around pi
                          ptol = 0L,          # shell thickness tolerance (integer dist)
                          mask = NULL,        # 0 valid; nonzero ignored
                          geom = TRUE,        # TRUE: geometric Ga; FALSE: triangulation Ga
                          tri_tol = 1e-6) {
  mat <- as.matrix(mat)
  nr <- nrow(mat); nc <- ncol(mat)
  if (is.null(cx)) cx <- (nc+1)/2
  if (is.null(cy)) cy <- (nr+1)/2
  if (is.null(r))  r  <- max(nc, nr)/2
  if (is.null(mask)) mask <- matrix(0, nr, nc)
  
  # ---- Vectorized central differences (fast) ----
  dx <- matrix(0, nr, nc); dy <- matrix(0, nr, nc)
  # x-gradient
  if (nc >= 3) dx[, 2:(nc-1)] <- (mat[, 3:nc] - mat[, 1:(nc-2)]) / 2
  dx[, 1]  <- (mat[, min(2, nc)] - mat[, 1])
  dx[, nc] <- (mat[, nc] - mat[, max(1, nc-1)])
  # y-gradient
  if (nr >= 3) dy[2:(nr-1), ] <- (mat[3:nr, ] - mat[1:(nr-2), ]) / 2
  dy[1, ]  <- (mat[min(2, nr), ] - mat[1, ])
  dy[nr, ] <- (mat[nr, ] - mat[max(1, nr-1), ])
  
  mods   <- sqrt(dx^2 + dy^2)
  phases <- atan2(dy, dx); phases[phases < 0] <- phases[phases < 0] + 2*pi
  maxGrad <- max(mods)
  
  # ---- Shells (integer distance) ----
  xs <- matrix(rep(1:nc, each = nr), nr, nc)
  ys <- matrix(rep(1:nr, times = nc), nr, nc)
  dists <- floor(sqrt((xs - cx)^2 + (ys - cy)^2))
  
  # ---- Initialize asymmetric copies ----
  adx <- dx; ady <- dy
  
  # ---- Zero weak vectors first ----
  weak <- (maxGrad <= 0) | ((mods / maxGrad) <= mtol)
  adx[weak] <- 0; ady[weak] <- 0
  
  # ---- Binned asymmetric pruning (vectorized) ----
  # orientation bins (even number so antipodal = +B/2)
  Btheta <- max(8L, 2L * as.integer(round(pi / max(ftol, 1e-6))))
  theta_bin <- as.integer(floor((phases / (2*pi)) * Btheta)) %% Btheta
  anti_bin  <- (theta_bin + Btheta %/% 2L) %% Btheta
  
  # modulus bins of width band = mtol * maxGrad
  band <- pmax(1e-12, mtol * maxGrad)
  mod_bin <- as.integer(floor(mods / band))
  
  # valid pixels (unmasked & not weak)
  valid <- which(mask == 0, arr.ind = TRUE)
  if (length(valid)) {
    keep0 <- !weak[valid]
    valid <- valid[keep0, , drop = FALSE]
  }
  if (nrow(valid) > 0) {
    i <- valid[,1]; j <- valid[,2]
    shell_b <- as.integer(floor(dists[valid] / max(1L, as.integer(ptol) + 1L)))
    mbin    <- mod_bin[valid]
    tbin    <- theta_bin[valid]
    tbin_a  <- anti_bin[valid]
    
    id_self <- paste(shell_b, mbin, tbin,   sep = "/")
    id_anti <- paste(shell_b, mbin, tbin_a, sep = "/")
    
    present <- table(id_self)
    has_self <- as.logical(present[id_self]); has_self[is.na(has_self)] <- FALSE
    has_anti <- (id_anti %in% names(present))
    
    mark_rm <- has_self & has_anti
    if (any(mark_rm)) {
      adx[cbind(i[mark_rm], j[mark_rm])] <- 0
      ady[cbind(i[mark_rm], j[mark_rm])] <- 0
    }
  }
  
  # ---- Counts over unmasked ----
  valid_all <- which(mask == 0, arr.ind = TRUE)
  totalVet <- nrow(valid_all)
  kept <- which((adx != 0 | ady != 0) & mask == 0, arr.ind = TRUE)
  totalAssimetric <- nrow(kept)
  
  # ---- phaseDiversity over kept set ----
  if (totalAssimetric >= 1) {
    ii <- kept[,1]; jj <- kept[,2]
    sumx <- sum(dx[ii, jj]); sumy <- sum(dy[ii, jj])
    smod <- sum(mods[ii, jj])
    phaseDiversity <- if (smod > 0) sqrt(sumx^2 + sumy^2) / smod else 0
  } else phaseDiversity <- 0
  
  Ga_geom <- (totalAssimetric / max(1, totalVet)) * (2 - phaseDiversity)
  
  # ---- Optional triangulation Ga (old) ----
  Ga_old <- NA_real_; n_points <- 0L; n_edges <- 0L; tri_pts <- NULL
  if (!geom) {
    amods <- sqrt(adx^2 + ady^2)
    sel <- which((amods > tri_tol) & (mask == 0), arr.ind = TRUE)
    if (nrow(sel) >= 3) {
      ii <- sel[,1]; jj <- sel[,2]
      px <- jj + 0.5 * adx[ii, jj]
      py <- ii + 0.5 * ady[ii, jj]
      tri_pts <- cbind(px, py)
      n_points <- nrow(tri_pts)
      if (!requireNamespace("geometry", quietly = TRUE))
        stop("Please install.packages('geometry') for triangulation Ga.")
      tri <- geometry::delaunayn(tri_pts)
      edges <- rbind(tri[,c(1,2), drop=FALSE],
                     tri[,c(2,3), drop=FALSE],
                     tri[,c(1,3), drop=FALSE])
      edges <- t(apply(edges, 1, sort))
      edges <- unique(edges)
      n_edges <- nrow(edges)
      Ga_old <- (n_edges - n_points) / n_points
    } else {
      n_points <- nrow(sel); n_edges <- 0L; Ga_old <- 0
    }
  }
  
  list(
    params = list(cx=cx, cy=cy, r=r, mtol=mtol, ftol=ftol, ptol=ptol, geom=geom),
    fields = list(dx=dx, dy=dy, mods=mods, phases=phases, maxGrad=maxGrad,
                  adx=adx, ady=ady, mask=mask),
    counts = list(totalVet=totalVet, totalAssimetric=totalAssimetric),
    phaseDiversity = phaseDiversity,
    Ga = if (geom) Ga_geom else Ga_old,
    Ga_geom = Ga_geom,
    Ga_old = Ga_old,
    triangulation = list(points=tri_pts, n_points=n_points, n_edges=n_edges)
  )
}

# -----------------------------
# BASE R QUIVER (quick look)
# -----------------------------
plot_quiver <- function(res, step = 6, scale = 1.6, asymmetric = TRUE) {
  dx <- if (asymmetric) res$fields$adx else res$fields$dx
  dy <- if (asymmetric) res$fields$ady else res$fields$dy
  mag <- res$fields$mods
  nr <- nrow(dx); nc <- ncol(dx)
  
  xs <- seq(1, nc, by = step)
  ys <- seq(1, nr, by = step)
  gx <- dx[ys, xs, drop = FALSE]
  gy <- dy[ys, xs, drop = FALSE]
  
  maxm <- max(res$fields$maxGrad, 1e-9)
  gx <- gx / maxm * scale * step
  gy <- gy / maxm * scale * step
  
  image(t(apply(mag, 2, rev)), col = gray.colors(128), axes = FALSE,
        main = if (asymmetric) "Asymmetric gradient field" else "Full gradient field")
  
  for (ii in seq_along(xs)) {
    for (jj in seq_along(ys)) {
      x <- xs[ii]; y <- ys[jj]
      x0 <- (x-1)/(nc-1); y0 <- 1 - (y-1)/(nr-1)
      arrows(x0, y0, x0 + gx[jj, ii]/nc, y0 - gy[jj, ii]/nr,
             length = 0.05, col = "red", lwd = 0.7)
    }
  }
}

# -----------------------------
# ggplot2 QUIVER (publication)
# -----------------------------
.mat_df <- function(M, name = "val") {
  nr <- nrow(M); nc <- ncol(M)
  df <- expand.grid(y = seq_len(nr), x = seq_len(nc))
  df[[name]] <- as.vector(M)
  df
}

build_quiver_df <- function(res, step = 6, asymmetric = TRUE, scale = 1.6) {
  dx <- if (asymmetric) res$fields$adx else res$fields$dx
  dy <- if (asymmetric) res$fields$ady else res$fields$dy
  mag <- res$fields$mods
  
  nr <- nrow(dx); nc <- ncol(dx)
  xs <- seq(1, nc, by = step); ys <- seq(1, nr, by = step)
  
  dxs <- dx[ys, xs, drop = FALSE]; dys <- dy[ys, xs, drop = FALSE]
  mds <- mag[ys, xs, drop = FALSE]
  
  maxm <- max(res$fields$maxGrad, 1e-9)
  ux <- dxs / maxm * scale * step
  uy <- dys / maxm * scale * step
  
  df <- expand.grid(y = ys, x = xs)
  df$u <- as.vector(ux); df$v <- as.vector(uy); df$mag <- as.vector(mds)
  df$alen <- sqrt(df$u^2 + df$v^2)
  df <- df[df$alen > 1e-6 & is.finite(df$mag), , drop = FALSE]  # filter here
  
  df$theta <- atan2(df$v, df$u)
  df$theta_deg <- (df$theta %% (2*pi)) * 180/pi
  df$xend <- df$x + df$u; df$yend <- df$y + df$v
  df
}


build_bg_df <- function(res, normalize = TRUE) {
  mag <- res$fields$mods
  if (normalize) {
    rng <- range(mag, finite = TRUE)
    if (diff(rng) > 0) mag <- (mag - rng[1]) / diff(rng)
  }
  .mat_df(mag, "bg")
}

plot_quiver_gg <- function(res,
                           step = 6,
                           asymmetric = TRUE,
                           scale = 1.6,
                           color_by = c("orientation", "magnitude"),
                           bg_alpha = 0.85,
                           show_legend = TRUE,
                           title = NULL,
                           subtitle = NULL) {
  color_by <- match.arg(color_by)
  bg <- build_bg_df(res)
  qv <- build_quiver_df(res, step = step, asymmetric = asymmetric, scale = scale)
  
  if (is.null(title))
    title <- if (asymmetric) "Asymmetric Gradient Field" else "Full Gradient Field"
  if (is.null(subtitle))
    subtitle <- if (color_by == "orientation")
      "Arrows colored by local orientation" else "Arrows colored by |∇|"
  
  p <- ggplot() +
    geom_raster(data = bg, aes(x = x, y = y, fill = bg),
                alpha = bg_alpha, interpolate = TRUE) +
    scale_fill_gradient(name = "|∇| (norm.)", low = "white", high = "black") +
    {
      if (color_by == "orientation") {
        geom_segment(
          data = qv,
          aes(x = x, y = y, xend = xend, yend = yend, color = theta_deg),
          arrow = arrow(length = unit(2.2, "pt"), type = "closed"),
          linewidth = 0.35
        )
      } else {
        geom_segment(
          data = qv,
          aes(x = x, y = y, xend = xend, yend = yend, color = mag),
          arrow = arrow(length = unit(2.2, "pt"), type = "closed"),
          linewidth = 0.35
        )
      }
    } +
    {
      if (color_by == "orientation")
        scale_color_viridis_c(name = "Orientation (deg)", option = "A", direction = 1)
      else
        scale_color_viridis_c(name = "|∇|", option = "C", trans = "sqrt")
    } +
    coord_fixed(expand = FALSE) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_void(base_size = 11) +
    theme(
      legend.position = if (show_legend) "right" else "none",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(margin = margin(b = 6))
    )
  
  p
}

# ============================
# DEMO
# ============================
# If you already have cube, use m <- cube[,,1]. Otherwise synthesize a test image.
if (exists("cube")) {
  m <- cube[,,1]
} else {
  set.seed(7)
  n <- 192
  m <- matrix(0, n, n)
  cx <- cy <- n/2
  for (i in 1:n) for (j in 1:n) {
    r2 <- (i-cy)^2 + (j-cx)^2
    m[i,j] <- exp(-r2/(2*(n/5)^2))
  }
  # add an oriented ridge
  for (k in 1:n) {
    i <- k; j <- floor(0.6*k + 18)
    if (j >= 1 && j <= n) m[i, j] <- m[i, j] + 0.7
  }
  # smooth + noise
  k <- 7
  ker <- matrix(1, k, k) / (k*k)
  conv2 <- function(img, ker) {
    px <- floor(nrow(ker)/2); py <- floor(ncol(ker)/2)
    pad <- matrix(0, nrow(img)+2*px, ncol(img)+2*py)
    pad[(px+1):(px+nrow(img)), (py+1):(py+ncol(img))] <- img
    out <- matrix(0, nrow(img), ncol(img))
    for (i in 1:nrow(img)) for (j in 1:ncol(img)) {
      out[i,j] <- sum(pad[i:(i+nrow(ker)-1), j:(j+ncol(ker)-1)] * ker)
    }
    out
  }
  m <- conv2(m, ker)
  m <- m + 0.08 * matrix(rnorm(n*n), n, n)
}

# Run GPA
m <- cube
res <- gpa_asym_fast(m, mtol = 0.02, ftol = 0.1, ptol = 0)

# Quick summary
print(list(Ga_geom = res$Ga_geom, Ga_old = res$Ga_old, phaseDiversity = res$phaseDiversity))

# --- BASE plot (fast glance)
plot_quiver(res, step = 4, scale = 5, asymmetric = TRUE)

# --- ggplot2 plots (publication)
p_orient <- plot_quiver_gg(res, step = 3, asymmetric = TRUE,  scale = 2.2,
                           color_by = "orientation")
p_mag    <- plot_quiver_gg(res, step = 2, asymmetric = TRUE, scale = 2.2, color_by = "magnitude")

# Show (if interactive) and save
print(p_orient)
print(p_mag)


ggsave("gpa_quiver_orientation.png", p_orient, width = 6, height = 6, dpi = 450)
ggsave("gpa_quiver_magnitude.png",  p_mag,    width = 6, height = 6, dpi = 450)



# Compute G2 from gpa_asym_fast() output (asymmetric field already pruned)
g2_from_res <- function(res) {
  adx <- res$fields$adx; ady <- res$fields$ady
  mask <- res$fields$mask
  V  <- sum(mask == 0)                        # total valid vectors
  keep <- which((adx != 0 | ady != 0) & mask == 0, arr.ind = TRUE)
  VA <- nrow(keep)
  if (VA < 1 || V < 1) return(0)
  
  ii <- keep[,1]; jj <- keep[,2]
  # vector sum and sum of moduli
  vx <- adx[ii, jj]; vy <- ady[ii, jj]
  num <- sqrt(sum(vx)^2 + sum(vy)^2)
  den <- sum(sqrt(vx^2 + vy^2))
  align <- if (den > 0) num / den else 0
  G2 <- (VA / V) * (2 - align)
  G2
}

# Example:
 res <- gpa_asym_fast(m, mtol=0.02, ftol=0.10, ptol=0)
 g2 <- g2_from_res(res)
# cat(sprintf("G2 = %.3f\n", g2))


 plot_orientation_rose <- function(res, asymmetric = TRUE, bins = 36, normalize = FALSE) {
   stopifnot(bins >= 4)
   dx <- if (asymmetric) res$fields$adx else res$fields$dx
   dy <- if (asymmetric) res$fields$ady else res$fields$dy
   mask <- res$fields$mask
   
   # angles in [0, 2pi), weights = |∇|
   theta <- (atan2(dy, dx) + 2*pi) %% (2*pi)
   w     <- sqrt(dx^2 + dy^2)
   
   # valid pixels only
   ok <- (mask == 0) & is.finite(theta) & is.finite(w)
   theta <- as.vector(theta[ok]); w <- as.vector(w[ok])
   
   # define fixed bin edges and centers
   brks <- seq(0, 2*pi, length.out = bins + 1L)
   mids <- head(brks, -1) + diff(brks)/2
   
   # bin as factor with all levels present
   bin <- cut(theta, breaks = brks, include.lowest = TRUE, labels = FALSE)
   bin <- factor(bin, levels = seq_len(bins))
   
   # sum weights per bin, preserving empty bins
   wsum <- tapply(w, bin, sum, default = 0)
   wsum[is.na(wsum)] <- 0
   if (normalize && sum(wsum) > 0) wsum <- wsum / sum(wsum)
   
   df <- data.frame(theta_mid = mids, weight = as.numeric(wsum), bin = seq_len(bins))
   
   # plot
   ggplot(df, aes(x = theta_mid, y = weight)) +
     geom_col(width = 2*pi / bins) +
     coord_polar(start = 0) +
     labs(
       title = if (asymmetric) "Orientation rose (asymmetric field)"
       else "Orientation rose (full field)",
       x = NULL, y = if (normalize) "Fractional weight" else "Weighted count (|∇|)"
     ) +
     theme_minimal(base_size = 11) +
     theme(
       axis.text.x = element_blank(),
       panel.grid.minor = element_blank()
     )
 }
 
 
plot_orientation_rose(res, asymmetric = TRUE, bins = 36)



plot_radial_G2 <- function(res, cx = NULL, cy = NULL, bins = 20) {
  adx <- res$fields$adx; ady <- res$fields$ady; mask <- res$fields$mask
  nr <- nrow(adx); nc <- ncol(adx)
  if (is.null(cx)) cx <- (nc+1)/2
  if (is.null(cy)) cy <- (nr+1)/2
  
  xs <- matrix(rep(1:nc, each = nr), nr, nc)
  ys <- matrix(rep(1:nr, times = nc), nr, nc)
  r  <- sqrt((xs - cx)^2 + (ys - cy)^2)
  rvec <- as.vector(r); mvec <- as.vector(mask)
  vx <- as.vector(adx); vy <- as.vector(ady)
  
  df <- data.frame(r = rvec, vx = vx, vy = vy, mask = mvec)
  df <- df[df$mask == 0, ]
  # define radial bins
  br <- quantile(df$r, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE)
  df$bin <- cut(df$r, br, include.lowest = TRUE, labels = FALSE)
  
  # compute G2 per bin
  per <- lapply(split(df, df$bin), function(d) {
    V <- nrow(d)
    keep <- (d$vx != 0 | d$vy != 0)
    VA <- sum(keep)
    if (VA < 1 || V < 1) return(data.frame(r_mid = median(d$r), G2 = 0))
    sx <- sum(d$vx[keep]); sy <- sum(d$vy[keep])
    num <- sqrt(sx^2 + sy^2)
    den <- sum(sqrt(d$vx[keep]^2 + d$vy[keep]^2))
    align <- if (den > 0) num/den else 0
    data.frame(r_mid = median(d$r), G2 = (VA/V) * (2 - align))
  })
  out <- do.call(rbind, per)
  
  ggplot(out, aes(r_mid, G2)) +
    geom_line() + geom_point(size = 1.8) +
    labs(title = "Radial G2 profile (asymmetric field)",
         x = "Radius (pixels)", y = "G2") +
    theme_minimal(base_size = 11)
}
# Usage:
 print(plot_radial_G2(res))
 
 gpa_field_df <- function(res, asymmetric = TRUE) {
   dx <- if (asymmetric) res$fields$adx else res$fields$dx
   dy <- if (asymmetric) res$fields$ady else res$fields$dy
   mag <- res$fields$mods
   mask <- res$fields$mask
   
   nr <- nrow(dx); nc <- ncol(dx)
   df <- expand.grid(y = seq_len(nr), x = seq_len(nc))
   df$u <- as.vector(dx)
   df$v <- as.vector(dy)
   df$mag <- as.vector(mag)
   df$mask <- as.vector(mask)
   df <- df[df$mask == 0 & is.finite(df$u) & is.finite(df$v) & is.finite(df$mag), , drop = FALSE]
   df
 }
 
 # 1) STREAMLINES (requires metR)
 # install.packages("metR")  # if needed
 plot_streamlines_gg <- function(res,
                                 asymmetric = TRUE,
                                 bg_alpha = 0.85,
                                 L = 20,            # streamline integration length (steps)
                                 n = 2000,          # number of seeds (density)
                                 resample = 1,      # resampling step (pixels)
                                 arrow_length = 2.2 # arrowhead size (pt)
 ) {
   if (!requireNamespace("metR", quietly = TRUE)) {
     stop("Please install.packages('metR') for streamlines.")
   }
   df <- gpa_field_df(res, asymmetric = asymmetric)
   
   # Background = normalized |∇|
   rng <- range(df$mag, finite = TRUE)
   df$bg <- if (diff(rng) > 0) (df$mag - rng[1]) / diff(rng) else 0
   
   # Title
   title <- if (asymmetric) "Asymmetric gradient field — streamlines" else "Full gradient field — streamlines"
   
   ggplot(df, aes(x = x, y = y)) +
     geom_raster(aes(fill = bg), alpha = bg_alpha, interpolate = TRUE) +
     scale_fill_gradient(name = "|∇| (norm.)", low = "white", high = "black") +
     metR::geom_streamline(aes(dx = u, dy = v, color = mag),
                           L = L, res = resample, n = n,
                           arrow = arrow(length = unit(arrow_length, "pt"), type = "closed"),
                           lineend = "round", linejoin = "round", size = 0.35) +
     scale_color_viridis_c(name = "|∇|", trans = "sqrt") +
     coord_fixed(expand = FALSE) +
     labs(title = title, x = NULL, y = NULL) +
     theme_void(base_size = 11) +
     theme(legend.position = "right",
           plot.title = element_text(face = "bold"))
 }
 
 # 2) DENSE QUIVER (fast) with geom_spoke
 plot_quiver_dense_gg <- function(res,
                                  asymmetric = TRUE,
                                  bg_alpha = 0.85,
                                  step = 3,        # arrow grid spacing (smaller = denser)
                                  scale = 2.0,     # arrow length multiplier
                                  color_by = c("orientation","magnitude")) {
   color_by <- match.arg(color_by)
   df <- gpa_field_df(res, asymmetric = asymmetric)
   
   # Subsample on a grid (keep every 'step' pixel)
   df <- df[df$x %% step == 0 & df$y %% step == 0, , drop = FALSE]
   
   # Normalize vector lengths for display
   maxm <- max(df$mag, na.rm = TRUE)
   df$len <- if (maxm > 0) (sqrt(df$u^2 + df$v^2) / maxm) * scale else 0
   
   # Orientation for coloring (degrees in [0,360))
   df$theta <- (atan2(df$v, df$u) %% (2*pi)) * 180/pi
   
   # Background = normalized |∇|
   rng <- range(df$mag, finite = TRUE)
   df$bg <- if (diff(rng) > 0) (df$mag - rng[1]) / diff(rng) else 0
   
   title <- if (asymmetric) "Asymmetric gradient field — dense quiver" else "Full gradient field — dense quiver"
   
   p <- ggplot(df, aes(x = x, y = y)) +
     geom_raster(aes(fill = bg), alpha = bg_alpha, interpolate = TRUE) +
     scale_fill_gradient(name = "|∇| (norm.)", low = "white", high = "black") +
     {
       if (color_by == "orientation") {
         geom_spoke(aes(angle = theta * pi/180, radius = len, color = theta),
                    linewidth = 0.35, lineend = "round")
       } else {
         geom_spoke(aes(angle = atan2(v, u), radius = len, color = mag),
                    linewidth = 0.35, lineend = "round")
       }
     } +
     {
       if (color_by == "orientation")
         scale_color_viridis_c(name = "Orientation (deg)", option = "A")
       else
         scale_color_viridis_c(name = "|∇|", option = "C", trans = "sqrt")
     } +
     coord_fixed(expand = FALSE) +
     labs(title = title, x = NULL, y = NULL) +
     theme_void(base_size = 11) +
     theme(legend.position = "right",
           plot.title = element_text(face = "bold"))
   p
 }

 
 p_stream <- plot_streamlines_gg(res, asymmetric = TRUE, L = 24, n = 3000, resample = 1)
 print(p_stream)
