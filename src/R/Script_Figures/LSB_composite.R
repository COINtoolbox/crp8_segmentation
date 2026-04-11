library(FITSio)
library(imager)
library(ggplot2)
library(png)
source("..//utils/make_rgb.R")

file <- "..//..//..//data/raw/datacube_sagui12.fits"

# Read all FITS and extract the image data
cube <- FITSio::readFITS(file)$imDat


mk <-  make_rgb(cube, r=7, g=4, b=2,
                              pansharpen=0.5, guide_band=2,
                              upscale=2, unsharp_sigma=1.1, unsharp_amount=0.7,
                              sat=0.9, gamma=1.0)



p <- plot_rgb_gg(mk)


ragg::agg_png("lsb_composite.png", width = 5000, height = 5000, units = "px", background = "black")
print(p)
dev.off()
