# Figure 04 Strategy Compare

Standalone paper-reproduction script for the three-panel segmentation strategy
comparison.

Entry point:

- `figure_04_strategy_compare_psfmatched.R`

This script reruns the underlying comparison pipelines on the PSF-matched cube
for `sagui10` with `N = 10`, then copies the three paper-facing PNGs into a
flat, Overleaf-friendly folder:

- `regions_plot_voronoi_nomask.png`
- `regions_plot_pixedfit_r.png`
- `regions_plot_sagui.png`

The output folder also contains a `manifest.csv` with the originating run tags
and source PNG paths.
