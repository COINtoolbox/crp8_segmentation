# Possible Workflow 

## Spectro-Photometric Segmentation

I will list below some possible scenarios we tan take to address the segmentation problem on spectro-photometric data. 


Data link: [JADES](https://archive.stsci.edu/hlsp/jades)

### Context:

In general, segmentation models—including those commonly used in astrophysics and computer vision—tend to focus on grayscale images. Even when applied to RGB data, they rarely take into account the idea of a spectral energy distribution (SED), since it is not very meaningful to define a “typical” SED for, say, a tree or a dog.

By contrast, my previous code Capivara was designed to segment IFU datacubes in MaNGA, where the field typically contains only one galaxy. Even in the presence of foreground stars or mergers, I did not need to worry much about multiple independent objects in the same field. My approach was therefore based primarily on grouping pixels by spectral similarity while ignoring spatial information. Yet it performed quite well: after all, each spectrum in a datacube implicitly encodes some degree of spatial information, even if not 
explicitly.

The challenge now is to combine both strategies in an astrophysically coherent way. While Capivara can be extended to multiband images, its performance is suboptimal when used in isolation. 

Below, I will include some possible, though not exhaustive, steps we need to develop. 

### Possible steps:

1. Create our own customized version of something similar to [astrodendro](https://dendrograms.readthedocs.io/en/stable/). 
The method sounds simple enough to be reproducible, there is nothing wrong with this package, it's just for the sake of a self-consistent pipeline. Here is the first issue, the method works on each filter alone, so there is no consistency across bands. On the other hands, I think it's unfeasible to work directly on the SED space for high-resolution images. So we may need a two-step process, unless someone has a better idea. Here is a opportunity for a quick comparison against  [Segment Anything Models](https://ai.meta.com/sam2/). This first step won't be the final result, but rather a series of cutouts of candidate regions for further scrutinity on the SED/Spectral space. So ideally we want widows around target regions that will not split regions that are spatially connected, but also not target noisy/sky.

2. Once we have the cutouts, we run some variation of [capivara](https://github.com/RafaelSdeSouza/capivara), so we can create some sort of spectral tesselation within each target region. So we can assign each pixel to a given group and create some sort of "representative" SED/SPECTRA. Note that we could in principle run a pixel-wise SED fitting, but this would be likely sensitive to noisy. How we do this part is open to debate. Besides defining average is more complex than one may think in a first glimpse, in particular for order-sensitive structures. Naively we can just take the mean/median, but potentially more shape preserving summaries also exist such as Wasserstein barycenter.

3. At this point, our data product consists of a series of “tiny” data cubes centered on each region of interest (e.g., star-forming clumps). In each case, the set of pixels is replaced by a “representative” SED. We then run an SED-fitting code to map these pixels into astrophysical quantities of interest (e.g., metallicity, age, star-formation history). This step involves many possible sources of uncertainty, but we will not address those here. Instead, we will rely on standard approaches for now. While not strictly necessary at this stage, if time permits, one could explore accelerating the SED-fitting process via surrogate models—provided that accuracy is not compromised.
We therefore have three independent steps, which will ultimately need to be bridged coherently by someone with more software engineering–oriented skills.

4. How to save the information? One possibility is to borrow the format used in spatial models in geostatistics, see e.g. [Spatial data in R](https://r-spatial.github.io/sf/) for the use of the Geospatial Data Abstraction Library (GDAL). This approach allows data to be manipulated in vector geospatial formats (polygons, lines, etc.) and operated within a self-consistent algebraic framework.

5. A long shot for the LLM enthusiasts. An alternative avenue would be to integrate a large language model (LLM) capable of generating textual descriptions of each region based on its representative SED/spectrum. For instance, if our final product is a dashboard where users can upload a field and k-bands, the system could output a report summarizing properties such as the fraction of young stellar populations, metallicity profiles, and so on. This idea, however, is more of a future spin-off to be considered if steps 1–4 are successfully implemented. 




