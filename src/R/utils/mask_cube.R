mask_cube <- function(cube, labels, mode = c("zero","na")) {
  mode <- match.arg(mode)
  H <- dim(cube)[1]; W <- dim(cube)[2]; B <- dim(cube)[3]
  m <- (labels > 0)
  out <- cube
  for (b in seq_len(B)) 
  { Xb <- out[,,b];
   if (mode == "zero") Xb[!m] <- 0 else 
   Xb[!m] <- NA_real_; out[,,b] <- Xb }
    return(out)
  }



