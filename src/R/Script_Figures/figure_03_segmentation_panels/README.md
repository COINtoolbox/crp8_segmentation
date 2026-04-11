# Figure 03 Segmentation Panels

Standalone paper-reproduction scripts for the Sagui segmentation mosaic used in
the CRP8 paper.

These scripts reproduce the per-galaxy PNG panels that are later assembled in
Overleaf. Each entry script targets one variant of the same figure:

- `figure_03_segmentation_panels_original.R`
- `figure_03_segmentation_panels_psfmatched.R`
- `figure_03_segmentation_panels_copula.R`
- `figure_03_segmentation_panels_contourlet_psfmatched.R`

The shared implementation lives in:

- `figure_03_segmentation_panels_common.R`

Design notes:

- The `original` variant mirrors the key visual logic of
  `crp8_segmentation/src/R/pipelines/JWST_main.R`.
- The `psfmatched` variant applies the same pipeline to the PSF-homogenized
  cubes in `data/PSF_matched`.
- The `copula` variant builds the starlet mask on a Gaussian-copula transformed
  cube, but still segments the original cube masked by that support.
- The `contourlet_psfmatched` variant replaces the starlet support by a
  contourlet-like directional support grown from the starlet seed.

Each script writes:

- individual `sagui*.png` panels for Overleaf,
- segmentation FITS maps,
- per-region photometry CSV tables measured on the original cube, and
- a manifest CSV describing the generated assets.

For the PSF-matched variants, the segmentation itself is computed on the
PSF-homogenized cubes, while the region photometry is extracted from the
corresponding original raw cubes.
