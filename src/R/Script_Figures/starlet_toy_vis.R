## ============================================================
## Undecimated (à trous) demo em R — plots separados por linha
## ============================================================

## --- Filtro base (binomial 3-tap ~ Gauss suave) ---
h0 <- c(1, 2, 1) / 4  # soma = 1

## --- À trous: insere zeros entre taps para o nível j (espaço = 2^j - 1) ---
insert_holes <- function(h, j) {
  if (j == 0) return(h)
  spacing <- 2^j - 1
  out <- numeric(0)
  for (i in seq_along(h)) {
    out <- c(out, h[i])
    if (i < length(h)) out <- c(out, rep(0, spacing))
  }
  out
}

## --- Convolução 'same' (sem NAs): faz conv completa e recorta centralizado ---
conv_same <- function(x, h) {
  y <- convolve(x, rev(h), type = "open")    # len = length(x)+length(h)-1
  L <- length(h); N <- length(x)
  start <- floor(L/2) + 1
  end   <- start + N - 1
  y[start:end]
}

## --- Sinal 1-D: mistura de senoides + ruído suave ---
set.seed(42)
n <- 256
x <- seq(-pi, 4*pi, length.out = n)

# Define y as a function of x
y <-  sin(x)*exp(-1.5*x^2) + 0.25*cos(2*x) + rnorm(length(x), 0, 0.05)


dfg <- data.frame(x=x,y=y)

pdf("curve_toy.pdf",width=8,height = 4)
ggplot(dfg,aes(x=x,y=y))+
  geom_line(color="#a88565",linewidth = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=16),
        axis.text =  element_text(size=12))
dev.off()

J <- 4  # níveis j = 0,1,2,3

## --- Decomposição multiescala: a_j e w_j ---
a <- vector("list", J + 1)   # a_0...a_J
w <- vector("list", J)       # w_1...w_J
a[[1]] <- y                  # a_0

for (j in 0:(J-1)) {
  hj <- insert_holes(h0, j)
  a[[j+2]] <- conv_same(a[[j+1]], hj)  # a_{j+1} = a_j * h_j
  w[[j+1]] <- a[[j+1]] - a[[j+2]]      # w_{j+1} = a_j - a_{j+1}
}

## --- Reconstrução: x = a_J + sum_j w_j ---
a_coarse <- a[[J+1]]
y_recon  <- a_coarse
for (jj in seq_along(w)) y_recon <- y_recon + w[[jj]]
recon_err <- max(abs(y_recon - y))
cat(sprintf("Erro máx. de reconstrução: %.2e\n", recon_err))

## =======================
## Paleta/cores auxiliares
## =======================
cols <- c("#d9bea1", "#a88565", "#3591a0", "#255e65")

## =============
## FIGURA 1: h_j
## =============
op <- par(no.readonly = TRUE)
par(mar = c(4,4,3,1),
    oma = c(1,0,1,0),
    cex.main = 1.6,     # title size
    cex.lab  = 1.4,     # axis labels
    cex.axis = 1.2
)
plot(0, 0, type = "n", xlab = "índice", ylab = "amplitude",
     xlim = c(0, 42), ylim = c(0, 0.5),
     main = "Filters à trous h_j")
xpos <- 1
for (j in 0:(J-1)) {
  hj <- insert_holes(h0, j)
  idx <- seq(xpos, by = 1, length.out = length(hj))
  segments(idx, 0, idx, hj, lwd = 3, col = cols[j+1])
  points(idx, hj, pch = 16, col = cols[j+1])
  text(mean(idx), max(hj)*1.08, labels = paste0("j=", j), col = cols[j+1])
  xpos <- max(idx) + 2
}
legend("topright", legend = paste0("j=", 0:(J-1)), lwd = 4, col = cols, bty = "n")



# Create a tibble of kernels for j = 0...(J-1)
df_list <- vector("list", J)
xpos <- 1

for (j in 0:(J - 1)) {
  hj  <- insert_holes(h0, j)
  idx <- seq(xpos, by = 1, length.out = length(hj))
  
  df_list[[j + 1]] <- data.frame(
    j     = j,
    idx   = idx,
    value = hj
  )
  
  xpos <- max(idx) + 2
}

kernels_df <- do.call(rbind, df_list)

# compute label positions
kernel_labels <- kernels_df |>
  dplyr::group_by(j) |>
  dplyr::summarise(
    xpos = mean(idx),
    ypos = -0.025   # slightly below zero for readability
  )

pdf("starlet_kernel.pdf",width=8,height = 4)
ggplot(kernels_df, aes(x = idx, y = value, colour = factor(j))) +
  geom_segment(aes(xend = idx, yend = 0), size = 1.5) +
  geom_point(size = 3) +
  geom_text(
    data = kernel_labels,
    aes(
      x = xpos,
      y = ypos,
      label = paste0("h[", j, "]")
    ),
    parse = TRUE,      # <-- This is essential
    size = 5,
    vjust = 1
  ) +
  scale_colour_manual(values = cols, name = "j") +
  labs(
    x = NULL,
    y = "amplitude",
    title = expression("À trous filters"~ h[j])
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.05)))+
  theme_void(base_size = 16) +
  theme(
    panel.grid = element_blank(),       # remove grid
    axis.text.x = element_blank(),      # remove x numbers
    axis.ticks.x = element_blank(),     # remove x ticks
    legend.position = "none",           # labels now printed under kernels
    plot.title = element_text(size = 16)
  )
dev.off()



## ===========================================
## FIGURA 2: coluna com a_0, a_1, ..., a_J (J+1 linhas)
## ===========================================
par(mfrow = c(J+1, 1),
    mar = c(3,4,2,1), 
    oma = c(1,0,1,0),
    cex.main = 1.6,     # title size
    cex.lab  = 1.4,     # axis labels
    cex.axis = 1.2
    )
yr <- range(unlist(a))
for (i in 1:(J+1)) {
  plot(a[[i]], type = "l", lwd = 2, col = if(i==1) "black" else cols[i-1],
       ylim = yr, xlab = "", ylab = "value",
#       main = sprintf("a_%d (mesmo tamanho, mais suave a cada nível)", i-1))
       main = bquote(a[.(i-1)])
)
#  abline(h = 0, col = "#99999955")
}
#mtext("a_j: versões progressivamente mais suaves (sem decimação)", outer = TRUE, line = 0.2)

a_df <- do.call(rbind, lapply(seq_along(a), function(i) {
  j <- i - 1L  # so a[[1]] -> j = 0, ..., a[[J+1]] -> j = J
  data.frame(
    j     = j,
    label = paste0("a[", j, "]"),  # for facet labels with math
    idx   = seq_along(a[[i]]),
    value = a[[i]]
  )
}))

# common y-range, like yr <- range(unlist(a))
yr <- range(a_df$value)

# j goes from 0 to J
j_vals <- 0:J

col_vec <- c(
  "0" = "black",                          # a_0 in black
  setNames(cols[1:J], as.character(1:J))  # a_1..a_J from your palette
)

#pdf("a_starlet.pdf",width  = 6,height = 10)
g1 <- ggplot(a_df, aes(x = idx, y = value, colour = factor(j))) +
  geom_line(linewidth = 1) +
  facet_wrap(label ~ ., labeller = label_parsed,ncol=1) +  # shows a[0], a[1], ... in strips
  scale_colour_manual(values = col_vec, guide = "none") +
  scale_y_continuous(limits = yr) +
  labs(
    x = "x",
    y = "y"
  ) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill="#d9bea1"),
        strip.text = element_text(size = 16),
    panel.grid = element_blank(),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    axis.title = element_text(size = 16)
  )
#dev.off()

## ======================================
## FIGURA 3: coluna com w_1, ..., w_J (J linhas)
## ======================================
par(mfrow = c(J, 1), 
    mar = c(3,4,2,1), 
    oma = c(1,0,1,0),
    cex.main = 1.6,     # title size
    cex.lab  = 1.4,     # axis labels
    cex.axis = 1.2      # tick labels
)

w_df <- do.call(rbind, lapply(seq_along(w), function(i) {
  j <- i  # w[[1]] → j = 1 corresponds to w_1
  
  data.frame(
    j     = j,
    label = paste0("w[", j, "] == a[", j, "] - a[", j-1, "]"),
    idx   = seq_along(w[[i]]),
    value = w[[i]]
  )
}))


yrw <- range(w_df$value)

j_vals <- 1:J

col_vec <- c(
  "0" = "black",                          # a_0 in black
  setNames(cols[1:J], as.character(1:J))  # a_1..a_J from your palette
)

g2 <- ggplot(w_df, aes(x = idx, y = value, colour = factor(j))) +
  geom_line(linewidth = 1) +
  facet_wrap(label ~ ., labeller = label_parsed,ncol=1) +  # shows a[0], a[1], ... in strips
  scale_colour_manual(values = col_vec, guide = "none") +
  scale_y_continuous(limits = yrw ) +
  labs(
    x = "x",
    y = "y"
  ) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill="#d9bea1"),
        strip.text = element_text(size = 16),
        panel.grid = element_blank(),
        axis.text.x  = element_text(size = 12),
        axis.text.y  = element_text(size = 12),
        axis.title = element_text(size = 16)
  )


pdf("aw_starlet.pdf",width  = 12,height = 10)
g1+g2
dev.off()
yrw <- range(do.call(c, w))
for (i in 1:J) {
  plot(w[[i]], type = "l", lwd = 2, col = cols[i],
       ylim = yrw, xlab = "", ylab = "scale-j coefficients",
       main = bquote(w[.(i)] == a[.(i-1)] - a[.(i)])
       
       )
#  abline(h = 0, col = "#99999955")
}
#mtext("w_j: detail coefficients per scale (same dimension as the signal)", 
#      outer = TRUE, line = 0.2)

## ==========================================
## FIGURA 4: checagem de reconstrução (overlay)
## ==========================================
par(mfrow = c(1,1), mar = c(4,4,3,1))
plot(y, type = "l", lwd = 2, col = "black",
#     main = sprintf("Reconstrução (erro máx = %.2e)", recon_err),
     xlab = "x", ylab = "y")
lines(y_recon, col = "#1b9e77", lwd = 2, lty = 2)
#legend("topright", c("original", "reconstruído"),
#       col = c("black", "#1b9e77"), lwd = 2, lty = c(1,2), bty = "n")

par(op)

## ===========================================================
## Intuitive reading:
##  - Top-left: filters widen by inserting zeros (holes).
##  - Top-right: a_j are same-size blurred versions (no decimation).
##  - Bottom-left: w_j are detail bands capturing oscillations
##                 at progressively larger scales.
##  - Bottom-right: reconstruction from all w_j + coarse a_J
##                  matches the original sinusoid exactly.
## ===========================================================

## ==============================
## TL;DR intuition for the plots:
## ==============================
## - Top-left: the filters h_j keep the same weights but develop "holes" (zeros)
##   between taps as j grows → effectively wider smoothing without downsampling.
## - Top-right: a_j are same-sized signals, progressively blurrier each level.
## - Bottom-left: w_j capture what is LOST from a_j to a_{j+1} → scale-specific details.
## - Bottom-right: adding all w_j back to the coarsest a_J reconstructs x exactly.
