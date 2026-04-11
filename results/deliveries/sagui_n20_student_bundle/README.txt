Sagui N=20 student delivery bundle

Contents
- segmentation_fits/: region-ID FITS maps used for the N=20 paper rerun
- flux_per_region_mad_sky/: per-region photometry CSVs using mad_sky uncertainties
- flux_per_region_flux_over_sqrt_n/: per-region photometry CSVs using flux/sqrt(N) uncertainties

Important conventions
- Segmentation was computed on the PSF-matched cubes.
- Photometry was summed on the original/raw cubes.
- Fluxes and flux uncertainties in the CSVs are in Jy.
- The raw and PSF-matched cube pixel units were treated as 10 * nanoJy before conversion.
- The recommended tables for external SED fitting are the mad_sky CSVs unless you explicitly want the alternative error model.

Systems included
- sagui1_2_3
- sagui4
- sagui5_6
- sagui7
- sagui8
- sagui9
- sagui10
- sagui11
