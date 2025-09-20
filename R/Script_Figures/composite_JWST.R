library(FITSio)
library(imager)
library(ggplot2)
library(png)
source("..//utils/make_rgb.R")

files <- sprintf("..//..//data/raw/datacube_reg%d.fits", c(1,2,4,5,6,7,8))

# Read all FITS and extract the image data
cubes <- lapply(files, function(f) FITSio::readFITS(f)$imDat)


mk <- function(cube) make_rgb(cube, r=7, g=4, b=2,
                              pansharpen=0.5, guide_band=2,
                              upscale=2, unsharp_sigma=1.1, unsharp_amount=0.7,
                              sat=0.9, gamma=1.0)
rgbs <- lapply(cubes, mk)


p1 <- plot_rgb_gg(rgbs[[1]])
p2 <- plot_rgb_gg(rgbs[[2]])
p3 <- plot_rgb_gg(rgbs[[3]])
p4 <- plot_rgb_gg(rgbs[[4]])
p5 <- plot_rgb_gg(rgbs[[5]])
p6 <- plot_rgb_gg(rgbs[[6]])
p7 <- plot_rgb_gg(rgbs[[7]])




mosaic <- (p1 | p4|  p7) /
  ( p2| p3 |  p5 | p6) & theme(
    panel.spacing = grid::unit(0, "pt"),           # no gaps
    plot.margin   = grid::unit(c(0,0,0,0), "pt")   # no outer margin
  )

mosaic <- mosaic +
  plot_annotation(
    title = "**<span style='color:black;'>JWST/NIRCam</span>**<br>
             <span style='font-size:180pt'>
             (<span style='color:red;'>F277W</span> / 
              <span style='color:green;'>F182W</span> / 
              <span style='color:blue;'>F115W</span>)
             </span>"
  )  &
  theme(
    # Rich text rendering
    plot.title = element_markdown(
      size = 180, face = "bold", hjust = 0.5,
      lineheight = 1.1,
      margin = margin(b = 6)   # pull title closer to plots
    ),
    plot.title.position = "plot"
  )

ragg::agg_png("mosaic.png", width = 8000, height = 4500, units = "px", background = "black")
print(mosaic)
dev.off()
