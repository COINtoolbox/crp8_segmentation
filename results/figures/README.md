# Figures Layout

This folder now has a clearer split between intermediate products and paper-style composites.

## Existing folders kept as-is

- `sagui_seg/`
  Individual segmentation PNGs and benchmark outputs.
- `sed_maps/<tag>/`
  Per-galaxy property maps and FITS rasters.
- `starlet/` and `starlet_png/`
  Starlet decomposition products.

## New folders for assembled figures

- `paper_panels/segmentation/`
  Rebuilt multi-object segmentation mosaics.
- `paper_panels/sed_properties/`
  Rebuilt 2x2 SED-property panels from the single-map PNGs.

## New surface-SFR products

The new surface-SFR maps are written to a dedicated folder so the older map set stays untouched:

- `sed_maps_surface_sfr/<tag>/logZ_smooth.png`
- `sed_maps_surface_sfr/<tag>/Av_smooth.png`
- `sed_maps_surface_sfr/<tag>/age_mw_smooth.png`
- `sed_maps_surface_sfr/<tag>/sfr_smooth.png`
- `sed_maps_surface_sfr/<tag>/surface_sfr_regions.csv`

These values come directly from the new SED result tables:

- `sigma_sfr_10myr`
- `sigma_sfr_100myr`

The current `sfr_smooth.*` files in `sed_maps_surface_sfr/` are mapped from `sigma_sfr_100myr` so the file naming stays parallel to the older `sed_maps/` tree.
