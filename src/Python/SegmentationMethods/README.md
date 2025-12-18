## SEREIAS: Spectral Energy-aware Region Extraction via Interscale Aggregation with Starlets.


COIN CRP2025. Alberto, Celine, Rafael, Aarya, Shravya, Reinaldo, Ana, Andressa, Thallis, Lilianne, Rupesh, Emille, Kristen, Darc [for this specific package only!].


TLDR: SEREIAS is a spatial segmentation tool that is based on a SED-aware analysis that segments hyperspectral cubes at multiple starlet wavelet scales, then fuses and refines the thresholded regions using interscale connection evidence and adjacency-based spectral similarity.


## Installation

To install the package, first clone the repository and navigate to the `SegmentationMethods` directory:

```bash
git clone https://github.com/COINtoolbox/crp8_segmentation.git
cd crp8_segmentation/SegmentationMethods
pip install .
```

## A bit more details

SEREIAS is an acronym for "Spectral Energy-aware Region Extraction via Interscale Aggregation with Starlets". The method is an SED-aware spatial segmentation strategy for hyperspectral data cubes (x, y, λ). 

It decomposes each cube with a 2D starlet transform acrss wavelengths, segments each scale via saliency-driven, neighborhood-constrained region growing that enforces spectral coherence using a simple cosine similarity between spectra, then fuses all scales into a single, unique global label map using interscale evidence (saliency or scale detection priority rules). A final (optional) adjacency-based refinement compares region SEDs in the original cube to merge spectrally indistinguishable neighbors while preserving genuinely distinct structures. 

SEREIAS can output per-scale maps, a fused segmentation, and optional diagnostics (mosaics, wavelet power metrics), and is highly configurable (e.g., seed quantiles, cosine-similarity thresholds, 4 or 8 pixel spatial connectivity, minimum detectable region size). Designed for astronomy and remote sensing alike, SEREIA emphasizes spectral fidelity and spatial contiguity to produce clean, interpretable object maps from complex multi-scale data.

To use SEREIAS, you just need to import the SEREIAS module (add import sereias to the top of your Python code) and then call its functions. For a primer about what functions to call and what order they need to be called, check the two notebooks that we provide in this folder (demo_sereias.ipynb and development_sereias.ipynb).

Enjoy your segmentation thanks to SEREIAS! But be careful and remember Undine... Resist their temptation... SEREIAS's chants are powerful and irresistibly alluring sounds, and if you are not careful and respectful with them, they may lure you to your destruction and death.

