rgb_mask_to_df <- function(cube_rgb, L, rotate_ccw = TRUE) {
  R <- cube_rgb$R; G <- cube_rgb$G; B <- cube_rgb$B
  stopifnot(is.matrix(R), is.matrix(G), is.matrix(B), is.matrix(L))
  stopifnot(identical(dim(R), dim(G)), identical(dim(R), dim(B)),
            identical(dim(R), dim(L)))
  
  # Melt RGB channels -> long -> wide
  df_long <- reshape2::melt(list(R = R, G = G, B = B),
                            varnames = c("y","x"), value.name = "val")
  names(df_long)[names(df_long) == "L1"] <- "ch"
  df <- reshape2::dcast(df_long, y + x ~ ch, value.var = "val")
  df$y <- as.integer(df$y); df$x <- as.integer(df$x)
  
  # Melt mask to get alpha per (y,x), then merge
  alpha_df <- reshape2::melt(L > 0L, varnames = c("y","x"), value.name = "alpha")
  alpha_df$y <- as.integer(alpha_df$y); alpha_df$x <- as.integer(alpha_df$x)
  df <- merge(df, alpha_df, by = c("y","x"), all.x = TRUE, sort = FALSE)
  df$alpha[is.na(df$alpha)] <- FALSE
  
  # Replace NA channels (background) with 0 so rgb() is happy
  df$R[is.na(df$R)] <- 0
  df$G[is.na(df$G)] <- 0
  df$B[is.na(df$B)] <- 0
  
  # Pick maxColorValue automatically (supports [0,1], 0..255, 0..65535)
  detect_mcv <- function(R,G,B){
    m <- suppressWarnings(max(R, G, B, na.rm = TRUE))
    if (!is.finite(m) || m <= 1 + 1e-8) 1 else
      if (m <= 255 + 1e-8) 255 else
        if (m <= 65535 + 1e-8) 65535 else m
  }
  mcv <- detect_mcv(df$R, df$G, df$B)
  
  # Encode mask in alpha channel (0/1 * mcv)
  df$hex <- grDevices::rgb(df$R, df$G, df$B,
                           alpha = as.numeric(df$alpha) * mcv,
                           maxColorValue = mcv)
  
  # Plot coords (90° CCW like your example)
  if (rotate_ccw) {
    df$x_plot <- -df$y; df$y_plot <- df$x
  } else {
    df$x_plot <-  df$x; df$y_plot <- df$y
  }
  df
}
