# Workflow 

## Spectro-Photometric Segmentation

I will list below some possible scenarios we tan take to address the segmentation problem on spectro-photometric data. 


Data link: [JADES](https://archive.stsci.edu/hlsp/jades)

### Context:

In general, segmentation models—including those commonly used in astrophysics and computer vision—tend to focus on grayscale images. Even when applied to RGB data, they rarely take into account the idea of a spectral energy distribution (SED), since it is not very meaningful to define a “typical” SED for, say, a tree or a dog.

By contrast, my previous code Capivara was designed to segment IFU datacubes in MaNGA, where the field typically contains only one galaxy. Even in the presence of foreground stars or mergers, I did not need to worry much about multiple independent objects in the same field. My approach was therefore based primarily on grouping pixels by spectral similarity while ignoring spatial information. Yet it performed quite well: after all, each spectrum in a datacube implicitly encodes some degree of spatial information, even if not 
explicitly.

The challenge now is to combine both strategies in an astrophysically coherent way. While Capivara can be extended to multiband images, its performance is suboptimal when used in isolation.
