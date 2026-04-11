# ProSpect Workspace

This folder is the dedicated workspace for region-based SED fitting with
`ProSpect` inside `crp8_segmentation`.

Current science choice:

- segmentation is defined on the **PSF-matched** cubes
- region photometry is measured by default on the **original/raw** cubes
- fitted properties are then painted back onto the **segmentation map**

Why this default:

- the PSF-matched cubes make the segmentation more coherent across bands
- the original cubes preserve the native regional photometry used for fitting

Important caveat:

- using raw photometry reintroduces the native band-to-band PSF differences at
  the photometry stage
- for that reason, every script in this folder keeps a
  `photometry_mode = c("raw", "psfmatched")` switch so we can compare both
  choices later

Suggested workflow:

1. `prospect_fit_regions.R`
   - reads a region photometry table
   - builds a ProSpect-ready input bundle
   - writes per-region fitting outputs to `results/sed_fitting/prospect`
2. `prospect_paint_maps.R`
   - reads fitted region properties
   - reads the segmentation FITS
   - paints region-level quantities back onto the pixel grid
   - writes FITS and PNG products to `results/sed_maps/prospect`

Default paths:

- code:
  - `src/R/pipelines/prospect`
- fitting outputs:
  - `results/sed_fitting/prospect`
- map outputs:
  - `results/sed_maps/prospect`

Initial implementation status:

- the folder structure and reusable helpers are in place
- the scripts currently prepare and validate inputs
- the actual fitting loop is the next step

Learning workflow:

1. `prospect_learn_one_region.R`
   - fits one segmented region with a conservative ProSpect model
   - fixes redshift to the spectroscopic value
   - fixes metallicity by default to reduce degeneracy while learning
   - inflates the reported errors with a small systematic floor
   - writes:
     - a compact summary CSV
     - a model-vs-data photometry CSV
     - a diagnostic PNG/PDF

2. `prospect_batch_simple.R`
   - runs the same simple setup over every region in a target galaxy
   - writes a single summary table for quick inspection

3. `prospect_learning_helpers.R`
   - shared code used by both learning scripts
   - contains the simple model definition and plotting helper

Recommended first steps:

1. Run `prospect_learn_one_region.R` on one region that you already understand
   from the image and the current Prospector output.
2. Inspect the diagnostic plot and the model-vs-data CSV.
3. Only after that looks reasonable, run `prospect_batch_simple.R`.

Important practical notes:

- `ProSpect` expects photometry in Jy.
- The current default is still:
  - segmentation on PSF-matched cubes
  - photometry on the original/raw cubes
- The learning scripts default to a deliberately simple model because the
  broad-band JWST region SEDs do not strongly constrain dust, metallicity, and
  flexible SFHs simultaneously.
