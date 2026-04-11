# Figure 05 Benchmark Grid

Standalone paper-reproduction script for the 3 x 4 benchmark grid.

Entry point:

- `figure_05_benchmark_grid_psfmatched.R`

By default the script uses `sagui8` as the benchmark example and reruns the
three comparison pipelines on the PSF-matched cube for `N = 10, 15, 20, 30`.

The final output folder contains flat PNG filenames that are easy to assemble in
Overleaf:

- `sagui_n10.png`, `sagui_n15.png`, `sagui_n20.png`, `sagui_n30.png`
- `voronoi_n10.png`, ..., `voronoi_n30.png`
- `pixedfit_r_n10.png`, ..., `pixedfit_r_n30.png`

Set `SAGUI_FIG5_TAG` to another galaxy tag (for example `sagui7`) to rerun the
same benchmark grid for a different system.
