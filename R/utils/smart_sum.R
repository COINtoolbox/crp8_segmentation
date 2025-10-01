smart_sum <- function(cube, method=c("ivar","equal","trimmed","quantile","template","localivar"),
                      bg_q=0.2, trim=0.1, q=0.75, template=NULL,
                      bg_k=31L, noise_k=9L) {
  method <- match.arg(method)
  stopifnot(length(dim(cube))==3)
  H <- dim(cube)[1]; W <- dim(cube)[2]; K <- dim(cube)[3]
  
  # per-band robust sigma from faint tail
  band_sigma <- function(cube, bg_q){
    s <- numeric(K)
    for (k in seq_len(K)) {
      Ik <- cube[,,k]
      qv <- as.numeric(quantile(Ik, bg_q, na.rm=TRUE))
      bg <- Ik[Ik <= qv]
      s[k] <- stats::mad(bg, constant=1.4826, na.rm=TRUE) + 1e-12
    }
    s
  }
  
  if (method=="equal") {
    return(apply(cube, c(1,2), mean, na.rm=TRUE))
  }
  
  if (method=="ivar") {
    sig <- band_sigma(cube, bg_q)
    w <- 1/(sig^2); den <- sum(w)
    num <- matrix(0, H, W)
    for (k in seq_len(K)) num <- num + w[k]*cube[,,k]
    return(num/den)
  }
  
  if (method=="trimmed") {  # winsorised/trimmed mean across bands
    return(apply(cube, c(1,2), function(v){
      v <- sort(v); n <- length(v); k <- floor(trim*n)
      if (k>0) mean(v[(k+1):(n-k)], na.rm=TRUE) else mean(v, na.rm=TRUE)
    }))
  }
  
  if (method=="quantile") {  # robust stack (e.g., 75th percentile)
    return(apply(cube, c(1,2), stats::quantile, probs=q, na.rm=TRUE, names=FALSE))
  }
  
  if (method=="template") {  # matched filter to a known SED template
    stopifnot(length(template)==K)
    sig <- band_sigma(cube, bg_q); w <- 1/(sig^2)
    alpha <- template*w
    alpha <- alpha / sum(template*template*w)  # normalised matched filter
    M <- matrix(cube, H*W, K)
    P <- matrix(M %*% alpha, H, W)
    return(P)
  }
  
  if (method=="localivar") {  # per-pixel local IVAR using background-flattened noise
    # simple box blur using your conv2():
    box_blur <- function(x, k){ ker <- matrix(1,k,k)/(k*k); conv2(x,ker) }
    num <- matrix(0, H, W); den <- matrix(0, H, W)
    for (k in seq_len(K)) {
      Ik <- cube[,,k]
      bg  <- box_blur(Ik, bg_k)
      res <- Ik - bg
      # robust local sigma ~ MAD via local absolute deviations
      sigma_loc <- 1.4826 * box_blur(abs(res - box_blur(res, noise_k)), noise_k) + 1e-12
      w <- 1/(sigma_loc^2)
      num <- num + w*Ik
      den <- den + w
    }
    return(num/den)
  }
}
