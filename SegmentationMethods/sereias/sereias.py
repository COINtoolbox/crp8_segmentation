import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from collections import deque

"""
SEREIAS: Spectral Energy-aware Region Extraction via Interscale Aggregation with Starlets.

TLDR: SEREIAS is spatial segmentation tool that is based on a SED-aware analysis that segments hyperspectral cubes at multiple starlet wavelet scales, then fuses and refines the thresholded regions using interscale connection evidence and adjacency-based spectral similarity.

A bit more details: SEREIAS is an acronym for "Spectral Energy-aware Region Extraction via Interscale Aggregation with Starlets". The method is a SED-aware spatial segmentation strategy for hyperspectral data cubes (x, y, λ). It decomposes each cube with a 2D starlet transform acrss wavelengths, segments each scale via saliency-driven, neighborhood-constrained region growing that enforces spectral coherence, then fuses all scales into a single, unique global label map using interscale evidence (saliency or scale detection priority rules). A final (optional) adjacency-based refinement compares region SEDs in the original cube to merge spectrally indistinguishable neighbors while preserving genuinely distinct structures. SEREIA can output per-scale maps, a fused segmentation, and optional diagnostics (mosaics, wavelet power metrics), and is highly configurable (e.g., seed quantiles, cosine-similarity thresholds, 4 or 8 pixel spatial connectivity, minimum detectable region size). Designed for astronomy and remote sensing alike, SEREIA emphasizes spectral fidelity and spatial contiguity to produce clean, interpretable object maps from complex multi-scale data.

COIN CRP2025. Alberto, Celine, Rafael, Aarya, Shravya, Reinaldo, Lilianne, Ana, Andressa, Emille, Kristen, Darc, Rupesh, Thallis [for this specific package only!].
"""


def starlet_filter_1d(scale: int) -> np.ndarray:
    """
    Build the 1D a trous (B3-spline) filter for a given scale.

    Base (undilated) filter: [1, 4, 6, 4, 1] / 16

    For scale j (1-indexed), we insert (2**(j-1) - 1) zeros between taps.
    The resulting length is 1 + 4 * 2**(j-1).
    """
    if scale < 1:
        raise ValueError("scale must be >= 1 (1-indexed)")
    base = np.array([1., 4., 6., 4., 1.], dtype=np.float64) / 16.0
    if scale == 1:
        return base
    step = 2 ** (scale - 1)
    L = 1 + 4 * step
    h = np.zeros(L, dtype=np.float64)
    idxs = np.arange(0, 5 * step, step, dtype=int)
    h[idxs] = base
    return h
    

def _conv1d_reflect(a: np.ndarray, k: np.ndarray, axis: int) -> np.ndarray:
    """1D convolution along `axis` with 'reflect' boundary conditions.

    Robust implementation: move the chosen axis to the end, pad/slide there,
    do a tensordot over the window axis, then move the axis back.
    """
    from numpy.lib.stride_tricks import sliding_window_view

    a = np.asarray(a, dtype=np.float64)
    k = np.asarray(k, dtype=np.float64)
    r = (k.size - 1) // 2

    # Move the target axis to the last position
    a_last = np.moveaxis(a, axis, -1)

    # Reflect-pad along the last axis
    pad = [(0, 0)] * a_last.ndim
    pad[-1] = (r, r)
    ap = np.pad(a_last, pad, mode="reflect")

    # Build sliding windows over the last axis -> shape (..., N, k)
    win = sliding_window_view(ap, window_shape=k.size, axis=-1)

    # Convolution (reverse kernel to match true convolution, not correlation)
    out_last = np.tensordot(win, k[::-1], axes=([-1], [0]))

    # Move the last axis back to its original position
    out = np.moveaxis(out_last, -1, axis)
    return out


def _conv2d_separable_xy(cube_xyl: np.ndarray, h1d: np.ndarray) -> np.ndarray:
    """Apply separable 2D convolution with the same 1D kernel along X then Y."""
    tmp = _conv1d_reflect(cube_xyl, h1d, axis=0)  # along x
    smooth = _conv1d_reflect(tmp, h1d, axis=1)    # along y
    return smooth
    

def starlet_decompose_cube(
    cube_xyl: np.ndarray,
    n_scales: int,
    *,
    positive: bool = True,
    k_sigma: float | None = None,
    return_residual: bool = False,
):
    """Sparse positive starlet decomposition for (x, y, λ) data cubes (no spectral smoothing)."""
    if not isinstance(n_scales, int) or n_scales < 1:
        raise ValueError("n_scales must be a positive integer (>= 1).")
    cube = np.asarray(cube_xyl, dtype=np.float64)
    if cube.ndim != 3:
        raise ValueError("Input cube must be 3D with shape (x, y, lambda).")
    nx, ny, nl = cube.shape
    if nx < 5 or ny < 5:
        raise ValueError("Starlet requires spatial dims >= 5 in (x, y).")

    c = cube.copy()
    scales: list[np.ndarray] = []

    for j in range(1, n_scales + 1):
        h = starlet_filter_1d(j)
        smooth = _conv2d_separable_xy(c, h)
        w = c - smooth

        if positive:
            np.maximum(w, 0, out=w)

        if k_sigma is not None and k_sigma > 0:
            med = np.median(w, axis=(0, 1), keepdims=True)
            mad = np.median(np.abs(w - med), axis=(0, 1), keepdims=True)
            sigma = mad / 0.67448975 + 1e-12
            w = np.maximum(w - k_sigma * sigma, 0.0)

        scales.append(w)
        c = smooth

    if return_residual:
        return scales, c
    return scales


def _norm01(img2d, p_lo=1.0, p_hi=99.5, *, lo=None, hi=None):
    """Percentile-based normalization to [0,1]. If lo/hi given, use them (global)."""
    if lo is None or hi is None:
        lo, hi = np.percentile(img2d, [p_lo, p_hi])
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            lo, hi = float(img2d.min()), float(img2d.max())
            if hi <= lo:
                return np.zeros_like(img2d, dtype=np.float32), lo, hi
    out = (img2d - lo) / (hi - lo) if hi > lo else np.zeros_like(img2d)
    return np.clip(out, 0.0, 1.0).astype(np.float32), lo, hi

# OLDER VERSION; keeep it here for the moment.
#def _norm01(img2d, p_lo=1.0, p_hi=99.5):
#    lo, hi = np.percentile(img2d, [p_lo, p_hi])
#    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
#        lo, hi = float(img2d.min()), float(img2d.max())
#        if hi <= lo:
#            return np.zeros_like(img2d, dtype=np.float32)
#    out = (img2d - lo) / (hi - lo)
#    return np.clip(out, 0.0, 1.0).astype(np.float32)


def starlet_scales_mosaic(
    scales,
    residual=None,
    *,
    lambda_index=None,       # if set, visualize this λ-slice; else aggregate across λ
    aggregate="sum",         # 'sum' | 'mean' | 'max' (used when lambda_index is None)
    ncols=None,              # number of columns for grid; None => single row (back-compat)
    gap_px=12,               # pixel gap between tiles (both x and y)
    normalize="per",         # 'per' (each tile) or 'global' (shared lo/hi)
    p_lo=1.0, p_hi=99.5,     # percentiles for normalization
    original=None,           # pass the original cube to append as the final tile
    original_label="Original",
):
    """
    Build a single 2D mosaic (no subplots) showing each scale, optional residual,
    and optionally the ORIGINAL cube as the last tile. Returns:
        mosaic : 2D float32 image in [0,1]
        labels : list[str] for each tile
        anchors: list[dict] with placement info per tile:
                 {'row','col','y','x','h','w','label'}
    """
    # --- 1) Collect 2D tiles from (x,y,λ) arrays
    def pick_plane(arr):
        if lambda_index is not None:
            return np.asarray(arr[:, :, lambda_index])
        if aggregate == "sum":
            return np.asarray(arr.sum(axis=2))
        if aggregate == "mean":
            return np.asarray(arr.mean(axis=2))
        if aggregate == "max":
            return np.asarray(arr.max(axis=2))
        raise ValueError("aggregate must be one of {'sum','mean','max'}")

    tiles_raw, labels = [], []
    for j, a in enumerate(scales, start=1):
        tiles_raw.append(pick_plane(a))
        labels.append(f"Scale {j}")

    if residual is not None:
        tiles_raw.append(pick_plane(residual))
        labels.append("Residual")

    if original is not None:
        tiles_raw.append(pick_plane(original))
        labels.append(original_label)

    if not tiles_raw:
        raise ValueError("No tiles to mosaic (empty inputs).")

    # --- 2) Normalize (per-tile or global)
    if normalize not in ("per", "global"):
        raise ValueError("normalize must be 'per' or 'global'.")

    gl_lo = gl_hi = None
    if normalize == "global":
        flat = np.concatenate([t.ravel() for t in tiles_raw])
        gl_lo, gl_hi = np.percentile(flat, [p_lo, p_hi])
        if not np.isfinite(gl_lo) or not np.isfinite(gl_hi) or gl_hi <= gl_lo:
            gl_lo, gl_hi = float(flat.min()), float(flat.max())
            if gl_hi <= gl_lo:
                gl_lo, gl_hi = 0.0, 1.0

    tiles = []
    for t in tiles_raw:
        if normalize == "per":
            res = _norm01(t, p_lo=p_lo, p_hi=p_hi)
        else:
            res = _norm01(t, p_lo=p_lo, p_hi=p_hi, lo=gl_lo, hi=gl_hi)
    
        # Accept both signatures: img01 OR (img01, lo, hi)
        if isinstance(res, tuple):
            t01 = res[0]
        else:
            t01 = res
        tiles.append(t01)

    # --- 3) Grid geometry
    n_tiles = len(tiles)
    ncols_eff = n_tiles if (ncols is None or ncols <= 0 or ncols >= n_tiles) else int(ncols)
    nrows = (n_tiles + ncols_eff - 1) // ncols_eff

    H_tile = max(t.shape[0] for t in tiles)
    W_tile = max(t.shape[1] for t in tiles)
    GY = int(gap_px); GX = int(gap_px)
    H = nrows * H_tile + (nrows - 1) * GY
    W = ncols_eff * W_tile + (ncols_eff - 1) * GX

    mosaic = np.zeros((H, W), dtype=np.float32)
    anchors = []

    # --- 4) Paste tiles and record anchors
    for idx, t in enumerate(tiles):
        r = idx // ncols_eff
        c = idx % ncols_eff
        y_base = r * (H_tile + GY)
        x_base = c * (W_tile + GX)
        h, w = t.shape
        y0 = y_base + (H_tile - h) // 2
        x0 = x_base + (W_tile - w) // 2
        mosaic[y0:y0 + h, x0:x0 + w] = t
        anchors.append(dict(row=r, col=c, y=y0, x=x0, h=h, w=w, label=labels[idx]))

    return mosaic, labels, anchors


def starlet_scale_powers(scales, residual=None, *, lambda_index=None, kind="l2"):
    """
    Compute total power per scale (and optionally residual).
    kind: 'l2' -> sum(w^2), 'l1' -> sum(|w|), 'sum' -> sum(w).
    Returns (powers, residual_power).
    """
    def power(arr):
        s = arr[:, :, lambda_index] if lambda_index is not None else arr
        if kind == "l2":
            return float(np.sum(s * s))
        elif kind == "l1":
            return float(np.sum(np.abs(s)))
        elif kind == "sum":
            return float(np.sum(s))
        else:
            raise ValueError("kind must be one of {'l2','l1','sum'}")
    p = [power(a) for a in scales]
    pres = power(residual) if residual is not None else None
    return np.array(p, dtype=np.float64), pres


def read_fits_cube(
    path: str,
    *,
    hdu=None,
    dtype=np.float64,
    memmap=True,
    squeeze_singletons=True,
    slice_axes=None,
    assume_spectral_first=True,
):
    """
    Read a FITS data cube and return (cube_xyλ, info).

    - Returns a numpy array ordered as (x, y, λ).
    - Only reorders axes; no spatial regridding or interpolation.
    - Detects the spectral axis from CTYPE keywords when possible (WAVE/FREQ/VELO/etc).
      If ambiguous, assumes the first numpy axis is spectral (typical FITS order: (λ, y, x)).

    Parameters
    ----------
    path : str
        Path to the FITS file.
    hdu : int | str | None
        HDU index or name. If None, the first HDU with ndim >= 3 is chosen.
    dtype : numpy dtype, default float64
        Cast output data to this dtype.
    memmap : bool, default True
        Open with memory-mapping (astropy).
    squeeze_singletons : bool, default True
        Squeeze size-1 axes before further handling.
    slice_axes : dict[int,int] | None
        For >3D data, select indices along extra axes BEFORE reordering (e.g. {0:0} for Stokes I).
    assume_spectral_first : bool, default True
        If CTYPE detection fails, treat the first numpy axis as spectral; if False, raise.

    Returns
    -------
    cube_xyλ : np.ndarray
        Data cube with shape (nx, ny, nλ).
    info : dict
        {'header', 'hdu', 'cunit', 'ctype',
         'spectral_header_axis', 'spectral_numpy_axis',
         'spectral_values', 'spectral_unit',
         'axis_order_numpy_in', 'axis_order_numpy_out'}
    """
    try:
        from astropy.io import fits
    except Exception as e:
        raise ImportError("read_fits_cube requires astropy. Please `pip install astropy`.") from e

    hdul = fits.open(path, memmap=memmap)
    try:
        # Choose HDU
        if hdu is None:
            sel = None
            for ext in hdul:
                if getattr(ext, "data", None) is not None and hasattr(ext.data, "ndim") and ext.data.ndim >= 3:
                    sel = ext
                    break
            if sel is None:
                raise ValueError("No HDU with ndim >= 3 found in the FITS file.")
        else:
            sel = hdul[hdu]

        data = sel.data
        hdr = sel.header.copy()
        if data is None:
            raise ValueError("Selected HDU has no data.")

        arr = np.asarray(data, dtype=dtype)

        # Optional slicing for higher-dim cubes (e.g., Stokes/time)
        if slice_axes:
            for ax, idx in sorted(slice_axes.items()):
                arr = np.take(arr, indices=idx, axis=ax)

        # Squeeze if requested
        if squeeze_singletons:
            arr = np.squeeze(arr)

        # Ensure 3D
        if arr.ndim != 3:
            squeezed = np.squeeze(arr)
            if squeezed.ndim == 3:
                arr = squeezed
            else:
                raise ValueError(
                    f"Expected a 3D cube, got {arr.ndim}D after slicing/squeezing. "
                    "Use 'slice_axes' to select extra dimensions."
                )

        # Detect spectral axis from FITS header CTYPE*
        n_header_axes = int(hdr.get("NAXIS", arr.ndim))
        ctypes = [str(hdr.get(f"CTYPE{i}", "")).upper() for i in range(1, n_header_axes + 1)]
        cunits = [str(hdr.get(f"CUNIT{i}", "")).upper() for i in range(1, n_header_axes + 1)]

        spectral_tokens = ("WAVE", "AWAV", "FREQ", "VELO", "VRAD", "BETA", "ENER", "WAVN", "ZOPT", "FELO")
        header_spec_axes = [i for i, ct in enumerate(ctypes, start=1) if any(tok in ct for tok in spectral_tokens)]
        spectral_header_axis = header_spec_axes[0] if len(header_spec_axes) == 1 else None

        # Map header axis (1-based) -> numpy axis (0-based) for current array
        # numpy axis index = (arr.ndim - header_axis_index)
        spectral_numpy_axis = (arr.ndim - spectral_header_axis) if spectral_header_axis is not None else (0 if assume_spectral_first else None)
        if spectral_numpy_axis is None:
            raise ValueError(
                "Could not detect spectral axis from CTYPE, and assume_spectral_first=False. "
                "Set assume_spectral_first=True or specify 'slice_axes'."
            )

        # Build transpose order to (x, y, λ)
        axes_all = [0, 1, 2]
        spatial_axes = [ax for ax in axes_all if ax != spectral_numpy_axis]
        # Typical FITS is (λ, y, x); we want (x, y, λ)
        order_in = (spatial_axes[1], spatial_axes[0], spectral_numpy_axis) if len(spatial_axes) == 2 else (0, 1, 2)
        cube_xyl = np.transpose(arr, order_in).copy()

        # Attempt linear spectral coordinates (CRVAL/CDELT/CRPIX); else None
        spectral_unit = None
        spectral_values = None
        if spectral_header_axis is not None:
            j = spectral_header_axis
            try:
                crval = float(hdr.get(f"CRVAL{j}"))
                cdelt = float(hdr.get(f"CDELT{j}"))
                crpix = float(hdr.get(f"CRPIX{j}"))
                spectral_unit = hdr.get(f"CUNIT{j}")
                nλ = cube_xyl.shape[2]
                pix = np.arange(nλ, dtype=np.float64) + 1.0  # FITS is 1-based
                spectral_values = (pix - crpix) * cdelt + crval
            except Exception:
                spectral_values = None

        info = dict(
            header=hdr,
            hdu=getattr(sel, "name", None),
            cunit={i+1: cunits[i] for i in range(min(len(cunits), arr.ndim))},
            ctype={i+1: ctypes[i] for i in range(min(len(ctypes), arr.ndim))},
            spectral_header_axis=spectral_header_axis,
            spectral_numpy_axis=spectral_numpy_axis,
            spectral_values=spectral_values,
            spectral_unit=spectral_unit,
            axis_order_numpy_in=order_in,
            axis_order_numpy_out=(0, 1, 2),
        )
        return cube_xyl, info
    finally:
        hdul.close()


def _saliency_from_scale(scale_xyl: np.ndarray, mode: str = "l2") -> np.ndarray:
    """Build a 2D saliency map from a (x,y,λ) scale cube."""
    if mode == "l2":
        return np.sqrt(np.sum(scale_xyl * scale_xyl, axis=2))
    elif mode == "l1":
        return np.sum(np.abs(scale_xyl), axis=2)
    elif mode == "sum":
        return np.sum(scale_xyl, axis=2)
    else:
        raise ValueError("mode must be one of {'l2','l1','sum'}")

def _spectra_matrix(cube_xyl: np.ndarray, *, normalize: str = "unit") -> np.ndarray:
    """Return spectra as (N, L) from (x,y,λ) and optionally normalize per pixel."""
    nx, ny, nl = cube_xyl.shape
    M = cube_xyl.reshape(nx * ny, nl).astype(np.float64, copy=False)
    eps = 1e-12
    if normalize == "unit":
        med = np.median(M, axis=1, keepdims=True)
        M = M - med
        norms = np.linalg.norm(M, axis=1, keepdims=True) + eps
        M = M / norms
    elif normalize == "zscore":
        mu = M.mean(axis=1, keepdims=True)
        sd = M.std(axis=1, keepdims=True) + eps
        M = (M - mu) / sd
        # L2 normalize as well so cosine distance is meaningful
        norms = np.linalg.norm(M, axis=1, keepdims=True) + eps
        M = M / norms
    elif normalize in (None, "none"):
        pass
    else:
        raise ValueError("normalize must be 'unit', 'zscore', or None")
    return M


def segment_starlet_scale_region_growing(
    cube_xyl: np.ndarray,
    scale_xyl: np.ndarray,
    *,
    spectra_source: str = "original",   # 'original' | 'detail'
    saliency_mode: str = "l2",
    seed_quantile: float = 0.90,        # start from top-Q saliency pixels
    connectivity: int = 4,              # 4 or 8
    cos_sim_min: float = 0.92,          # spectral similarity threshold (0..1)
    min_region: int = 25,               # drop or merge tiny blobs
    normalize: str = "unit",            # per-pixel spectral normalization
    merge_cos_min: float = 0.96,        # optional post-merge threshold
    return_debug: bool = False,
):
    """
    Region-growing segmentation at a given starlet scale.

    - Spatial contiguity: growth only to 4/8-neighbors in (x,y).
    - Spectral similarity: growth only if cosine(sim( pixel, region_centroid )) >= cos_sim_min.
    - Seeds: highest-saliency pixels first (from the chosen scale).
    - Spectra source: 'original' (recommended) or 'detail' (use the scale cube itself).

    Returns
    -------
    labels : (nx, ny) int32 array, 0=background (below seed threshold), regions start at 1
    info   : dict with debugging fields if return_debug=True
    """
    nx, ny, nl = cube_xyl.shape
    assert scale_xyl.shape == cube_xyl.shape, "scale_xyl must match cube shape"

    # 1) saliency & mask
    sal = _saliency_from_scale(scale_xyl, mode=saliency_mode)
    thr = np.quantile(sal, seed_quantile) if 0.0 <= seed_quantile < 1.0 else np.min(sal) - 1.0
    mask = sal >= thr

    # 2) spectra
    if spectra_source == "original":
        M = _spectra_matrix(cube_xyl, normalize=normalize)    # (N, L)
    elif spectra_source == "detail":
        M = _spectra_matrix(scale_xyl, normalize=normalize)
    else:
        raise ValueError("spectra_source must be 'original' or 'detail'")
    N, L = M.shape
    flat_idx = np.arange(N, dtype=np.int32)

    # 3) helpers
    def to_flat(i, j): return i * ny + j
    def to_ij(f): return (f // ny, f % ny)

    if connectivity == 8:
        neighbors = [(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(-1,1),(1,-1),(1,1)]
    elif connectivity == 4:
        neighbors = [(-1,0),(1,0),(0,-1),(0,1)]
    else:
        raise ValueError("connectivity must be 4 or 8")

    # 4) seed order (descending saliency)
    order = np.argsort(-sal, axis=None)  # flat indices sorted by saliency
    labels = np.zeros((nx, ny), dtype=np.int32)
    visited = np.zeros((nx, ny), dtype=bool)
    label_id = 0

    cos_min = float(cos_sim_min)

    # 5) region growing
    for f in order:
        i, j = divmod(f, ny)
        if labels[i, j] != 0 or not mask[i, j]:
            continue
        label_id += 1

        # region centroid (start at seed spectrum)
        centroid = M[f].copy()
        q = deque([f])
        labels[i, j] = label_id
        visited[i, j] = True
        n_pix = 1

        while q:
            pf = q.popleft()
            pi, pj = divmod(pf, ny)

            # explore neighbors
            for di, dj in neighbors:
                ni, nj = pi + di, pj + dj
                if ni < 0 or nj < 0 or ni >= nx or nj >= ny:
                    continue
                if visited[ni, nj] or not mask[ni, nj]:
                    continue

                nf = to_flat(ni, nj)
                # cosine similarity with current centroid (vectors already unit-norm)
                cos_sim = float(np.dot(M[nf], centroid))
                if cos_sim >= cos_min:
                    labels[ni, nj] = label_id
                    visited[ni, nj] = True
                    q.append(nf)
                    # incremental centroid update (keep normalized for stability)
                    n_pix += 1
                    centroid = centroid + (M[nf] - centroid) / n_pix
                    # re-normalize to unit length to keep dot-product as cosine
                    norm = np.linalg.norm(centroid) + 1e-12
                    centroid = centroid / norm
                else:
                    visited[ni, nj] = True  # mark checked to avoid revisiting

    # 6) small-region cleanup: relabel to 0 (or could reassign to best neighbor)
    if min_region and min_region > 1:
        sizes = np.bincount(labels.ravel())
        small = set(np.nonzero((sizes > 0) & (sizes < min_region))[0].tolist())
        if small:
            labels[np.isin(labels, list(small))] = 0

    # 7) optional merging of touching regions with very similar centroids
    if merge_cos_min is not None and merge_cos_min > cos_min and label_id > 1:
        # compute centroids
        centroids = {}
        for lab in range(1, label_id + 1):
            mask_lab = (labels == lab).ravel()
            if not np.any(mask_lab):
                continue
            v = M[mask_lab].mean(axis=0)
            v /= (np.linalg.norm(v) + 1e-12)
            centroids[lab] = v

        # adjacency pairs (touching labels)
        pairs = set()
        for i in range(nx):
            for j in range(ny):
                lab = labels[i, j]
                if lab == 0:
                    continue
                if i + 1 < nx and labels[i+1, j] != lab and labels[i+1, j] != 0:
                    pairs.add(tuple(sorted((lab, labels[i+1, j]))))
                if j + 1 < ny and labels[i, j+1] != lab and labels[i, j+1] != 0:
                    pairs.add(tuple(sorted((lab, labels[i, j+1]))))

        # simple union-find
        parent = {k: k for k in range(1, label_id + 1)}
        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a
        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb: parent[rb] = ra

        for a, b in pairs:
            if a in centroids and b in centroids:
                if float(np.dot(centroids[a], centroids[b])) >= float(merge_cos_min):
                    union(a, b)

        # apply unions
        for i in range(nx):
            for j in range(ny):
                lab = labels[i, j]
                if lab != 0:
                    labels[i, j] = find(lab)

        # compress labels to 1..K
        labs = np.unique(labels)
        labs = labs[labs != 0]
        remap = {lab: k for k, lab in enumerate(labs, start=1)}
        for lab, k in remap.items():
            labels[labels == lab] = k

    if return_debug:
        return labels, dict(saliency=sal, mask=mask, seed_threshold=thr)
    return labels

    
def compose_scale_and_seg_per_scale(scale_xyl, labels_xy, *,
                          lambda_index=None,   # or None to aggregate
                          aggregate="sum",     # 'sum' | 'mean' | 'max'
                          gap_px=16,
                          left_title="Scale",
                          right_title="Segmentation"):
    # Collapse (x,y,λ) -> (x,y)
    if lambda_index is not None:
        scale2d = scale_xyl[:, :, lambda_index]
    else:
        if   aggregate == "sum":  scale2d = scale_xyl.sum(axis=2)
        elif aggregate == "mean": scale2d = scale_xyl.mean(axis=2)
        elif aggregate == "max":  scale2d = scale_xyl.max(axis=2)
        else: raise ValueError("aggregate must be one of {'sum','mean','max'}")

    A = _norm01(np.asarray(scale2d))[0]
    B = _norm01(np.asarray(labels_xy).astype(np.float64))[0]  # just for display

    H = max(A.shape[0], B.shape[0])
    W = A.shape[1] + gap_px + B.shape[1]
    mosaic = np.zeros((H, W), dtype=np.float32)

    y0A = (H - A.shape[0]) // 2; x0A = 0
    y0B = (H - B.shape[0]) // 2; x0B = A.shape[1] + gap_px
    mosaic[y0A:y0A+A.shape[0], x0A:x0A+A.shape[1]] = A
    mosaic[y0B:y0B+B.shape[0], x0B:x0B+B.shape[1]] = B

    anchors = [
        {"x": x0A + A.shape[1]/2, "y": y0A - 6, "text": left_title},
        {"x": x0B + B.shape[1]/2, "y": y0B - 6, "text": right_title},
    ]
    return mosaic, anchors


def compose_scale_and_seg(
    scale_xyl, labels_xy,
    *, lambda_index=None, aggregate="sum", gap_px=16,
    left_title="Scale", right_title="Segmentation",
    overlay_boundaries=True, boundary_linewidth=1.2, boundary_alpha=1.0,
    seg_cmap_name="tab20"   # qualitative colormap
):
    """
    Single-axes composition (no subplots):
      - Left: scale image (2D), with optional polygonal boundaries (contours)
      - Right: segmentation (labels, ints) with a qualitative colormap
    Returns (fig, ax), mosaic size (H,W), and anchor dicts.
    """
    # Collapse (x,y,λ) -> (x,y)
    if lambda_index is not None:
        scale2d = scale_xyl[:, :, lambda_index]
    else:
        if   aggregate == "sum":  scale2d = scale_xyl.sum(axis=2)
        elif aggregate == "mean": scale2d = scale_xyl.mean(axis=2)
        elif aggregate == "max":  scale2d = scale_xyl.max(axis=2)
        else: raise ValueError("aggregate must be one of {'sum','mean','max'}")

    A = _norm01(np.asarray(scale2d))[0]
    B = np.asarray(labels_xy).astype(int)

    # Geometry
    H = max(A.shape[0], B.shape[0])
    W = A.shape[1] + gap_px + B.shape[1]
    y0A = (H - A.shape[0]) // 2; x0A = 0
    y0B = (H - B.shape[0]) // 2; x0B = A.shape[1] + gap_px

    # Figure (one axes)
    fig = plt.figure(dpi=150, figsize=(10, 4))
    ax = plt.gca()

    # Draw left tile by pasting into a 2D mosaic and showing once
    mosaic = np.zeros((H, W), dtype=np.float32)
    mosaic[y0A:y0A+A.shape[0], x0A:x0A+A.shape[1]] = A
    ax.imshow(mosaic, vmin=0.0, vmax=1.0)

    # Right tile: segmentation with a QUALITATIVE colormap (discrete)
    nlabels = int(B.max()) + 1
    if nlabels < 2: nlabels = 2
    seg_cmap = cm.get_cmap(seg_cmap_name, nlabels)  # e.g., 'tab20'
    ax.imshow(
        B,
        cmap=seg_cmap,
        vmin=-0.5, vmax=nlabels - 0.5,  # quantized bins per integer label
        interpolation="nearest",
        origin="upper",
        extent=(x0B, x0B + B.shape[1], y0B + B.shape[0], y0B)
    )

    # Optional polygonal boundaries over the LEFT tile (marching-squares via contour)
    if overlay_boundaries:
        Lmax = int(B.max())
        for lab in range(1, Lmax + 1):  # skip background 0
            ax.contour(
                (B == lab).astype(float),
                levels=[0.5],
                linewidths=boundary_linewidth,
                alpha=boundary_alpha,
                origin="upper",
                extent=(x0A, x0A + A.shape[1], y0A + A.shape[0], y0A),
                color="white"
            )

    # Titles
    ax.set_xticks([]); ax.set_yticks([])
    ax.text(x0A + A.shape[1]/2, y0A - 6, left_title,  ha="center", va="bottom", fontsize=10)
    ax.text(x0B + B.shape[1]/2, y0B - 6, right_title, ha="center", va="bottom", fontsize=10)

    plt.tight_layout(pad=0.5)
    return fig, ax, (H, W), dict(left=(x0A, y0A, A.shape[1], A.shape[0]), right=(x0B, y0B, B.shape[1], B.shape[0]))

        
def merge_scale_segmentations(
    seg_maps,                   # list of (nx, ny) int arrays, one per scale, 0=background
    *,
    scale_cubes=None,           # optional list of (nx, ny, λ) detail cubes w_j (same length as seg_maps)
    saliency_maps=None,         # optional list of (nx, ny) floats; overrides scale_cubes if provided
    saliency_mode="l2",         # when computing saliency from scale_cubes
    method="argmax",            # 'argmax' (per-pixel max-saliency) or 'first-wins' (hierarchical)
    priority="fine_first",      # tie/order: 'fine_first' or 'coarse_first'
    normalize="per",            # saliency normalization: 'per'|'global'|'none'
    relabel_compact=False,      # if True, remap final labels to 1..K (still unique globally)
):
    """
    Merge per-scale segmentations into one global label map.

    Returns
    -------
    labels_final : (nx, ny) int array, 0=background
    info : dict with offsets, chosen_scale, etc.
    """
    # --- validate & shapes
    n_scales = len(seg_maps)
    if n_scales == 0:
        raise ValueError("seg_maps is empty.")
    nx, ny = seg_maps[0].shape
    for L in seg_maps:
        if L.shape != (nx, ny):
            raise ValueError("All seg_maps must share the same (nx, ny) shape.")

    # --- compute (or take) saliency maps S_j
    S_list = None
    if saliency_maps is not None:
        if len(saliency_maps) != n_scales:
            raise ValueError("saliency_maps must match seg_maps length.")
        S_list = [np.asarray(S, dtype=np.float64) for S in saliency_maps]
    elif scale_cubes is not None:
        if len(scale_cubes) != n_scales:
            raise ValueError("scale_cubes must match seg_maps length.")
        S_list = [ _saliency_from_scale(np.asarray(W, dtype=np.float64), saliency_mode)
                   for W in scale_cubes ]
    else:
        # fallback: use binary masks as saliency
        S_list = [ (np.asarray(L) > 0).astype(np.float64) for L in seg_maps ]

    # --- optional normalization of S_j
    if normalize not in ("per", "global", "none"):
        raise ValueError("normalize must be 'per'|'global'|'none'")
    if normalize == "per":
        S_list = [
            (S / (np.percentile(S, 99.5) + 1e-12)) for S in S_list
        ]
    elif normalize == "global":
        stack = np.stack(S_list, axis=0)
        denom = np.percentile(stack, 99.5)
        S_list = [ S / (denom + 1e-12) for S in S_list ]
    # else 'none' keeps as-is

    # --- per-scale unique offsets so labels never collide
    offsets = []
    off = 0
    for L in seg_maps:
        offsets.append(off)
        off += int(np.max(L))

    # --- stacking helpers
    masks = [ (np.asarray(L, dtype=np.int32) > 0) for L in seg_maps ]
    # zero out saliency where that scale has no label
    S_masked = [ S * M for S, M in zip(S_list, masks) ]

    if method == "argmax":
        # compose per-pixel weights with deterministic tie-break
        weights = np.stack(S_masked, axis=0)  # (J, nx, ny)
        # tiny deltas for tie-breaking by priority
        if priority == "fine_first":   # j=0 (finest) wins ties
            deltas = np.linspace(1e-12, 0.0, n_scales, dtype=np.float64)[:,None,None]
        elif priority == "coarse_first":
            deltas = np.linspace(0.0, 1e-12, n_scales, dtype=np.float64)[:,None,None]
        else:
            raise ValueError("priority must be 'fine_first' or 'coarse_first'")
        weights = weights + deltas

        # argmax over scales
        j_star = np.argmax(weights, axis=0)            # (nx, ny) chosen scale index
        w_star = weights[j_star, np.arange(nx)[:,None], np.arange(ny)[None,:]]  # selected weight
        # background where no scale had a label (all zeros)
        any_mask = np.any(np.stack(masks, axis=0), axis=0)

        labels_final = np.zeros((nx, ny), dtype=np.int32)
        for j, Lj in enumerate(seg_maps):
            pick = (j_star == j) & any_mask
            # assign offset+label at chosen pixels
            labels_final[pick] = offsets[j] + np.asarray(Lj, dtype=np.int32)[pick]

    elif method == "first-wins":
        # fill in pixels from a priority-ordered list of scales, never overwriting
        if priority == "fine_first":
            order = list(range(n_scales))                # 0,1,2... (finest to coarsest)
        elif priority == "coarse_first":
            order = list(reversed(range(n_scales)))      # J-1,...,0
        else:
            raise ValueError("priority must be 'fine_first' or 'coarse_first'")

        labels_final = np.zeros((nx, ny), dtype=np.int32)
        for j in order:
            Lj = np.asarray(seg_maps[j], dtype=np.int32)
            write = (labels_final == 0) & (Lj > 0)
            labels_final[write] = offsets[j] + Lj[write]

        j_star = None  # not used in this mode
    else:
        raise ValueError("method must be 'argmax' or 'first-wins'")

    # --- optional compaction of ids to 1..K (still unique globally)
    mapping = None
    if relabel_compact:
        L = labels_final
        ids = np.unique(L)
        ids = ids[ids != 0]
        mapping = {old: new for new, old in enumerate(ids, start=1)}
        labels_compact = np.zeros_like(L)
        for old, new in mapping.items():
            labels_compact[L == old] = new
        labels_final = labels_compact

    info = dict(
        offsets=offsets,
        method=method,
        priority=priority,
        normalize=normalize,
        chosen_scale=j_star,      # None if method='first-wins'
        id_mapping=mapping,
    )
    return labels_final, info


def starlet_segment_cube(
    cube_xyl: np.ndarray,
    n_scales: int,
    *,
    # --- starlet params ---
    positive: bool = True,
    k_sigma: float | None = None,
    # --- per-scale segmentation params ---
    # Either a single dict applied to all scales, or a list[dict] with length n_scales
    segmentation_params: dict | list[dict] | None = None,
    # sensible defaults (used if segmentation_params is None or missing keys):
    spectra_source: str = "original",   # 'original' | 'detail'
    saliency_mode: str = "l2",
    seed_quantile: float = 0.90,
    connectivity: int = 4,
    cos_sim_min: float = 0.92,
    min_region: int = 25,
    normalize_spectra: str = "unit",    # 'unit' | 'zscore' | None
    merge_cos_min: float = 0.96,
    # --- fusion params ---
    fusion_method: str = "argmax",      # 'argmax' | 'first-wins'
    fusion_priority: str = "fine_first",# 'fine_first' | 'coarse_first'
    fusion_saliency_mode: str = "l2",   # how to compute per-scale saliency if needed
    fusion_normalize: str = "per",      # 'per' | 'global' | 'none'
    relabel_compact: bool = False,      # remap final ids to 1..K (still unique)
    # --- outputs ---
    return_intermediates: bool = True,
):
    """
    Full pipeline: starlet -> per-scale segmentation -> global fusion.

    Parameters
    ----------
    cube_xyl : np.ndarray
        Input cube, shape (x, y, λ).
    n_scales : int
        Number of starlet scales.
    positive, k_sigma :
        Starlet settings (positivity, optional soft-threshold for sparsity).
    segmentation_params :
        Dict of args for `segment_starlet_scale_region_growing`, applied to all scales,
        or a list of dicts of length n_scales for per-scale overrides.
        Recognized keys include: 'spectra_source','saliency_mode','seed_quantile',
        'connectivity','cos_sim_min','min_region','normalize','merge_cos_min'.
    spectra_source, saliency_mode, seed_quantile, connectivity, cos_sim_min,
    min_region, normalize_spectra, merge_cos_min :
        Default values used if not provided in segmentation_params.
    fusion_method, fusion_priority, fusion_saliency_mode, fusion_normalize, relabel_compact :
        Parameters forwarded to `merge_scale_segmentations`.

    Returns
    -------
    labels_final : (nx, ny) int32
        Unique global segmentation map (0 = background).
    info : dict (if return_intermediates=True)
        {
          'scales': list[(nx,ny,λ)] detail cubes,
          'residual': (nx,ny,λ),
          'per_scale_labels': list[(nx,ny)] of ints,
          'fusion_info': dict from merge_scale_segmentations,
          'starlet': {'positive':..., 'k_sigma':..., 'n_scales':...},
          'segmentation_params_resolved': list[dict] actually used per scale,
        }
    """
    # 1) Starlet decomposition
    scales, residual = starlet_decompose_cube(
        cube_xyl, n_scales, positive=positive, k_sigma=k_sigma, return_residual=True
    )

    # 2) Resolve per-scale segmentation parameters
    def _default_seg_params():
        return dict(
            spectra_source=spectra_source,
            saliency_mode=saliency_mode,
            seed_quantile=seed_quantile,
            connectivity=connectivity,
            cos_sim_min=cos_sim_min,
            min_region=min_region,
            normalize=normalize_spectra,
            merge_cos_min=merge_cos_min,
            return_debug=False,
        )

    if segmentation_params is None:
        seg_params_list = [_default_seg_params() for _ in range(n_scales)]
    elif isinstance(segmentation_params, dict):
        base = _default_seg_params()
        base.update(segmentation_params)
        seg_params_list = [base for _ in range(n_scales)]
    elif isinstance(segmentation_params, (list, tuple)):
        if len(segmentation_params) != n_scales:
            raise ValueError("segmentation_params list must have length n_scales.")
        seg_params_list = []
        for sp in segmentation_params:
            base = _default_seg_params()
            if sp is not None:
                base.update(sp)
            seg_params_list.append(base)
    else:
        raise ValueError("segmentation_params must be None, dict, or list[dict].")

    # 3) Per-scale segmentations
    per_scale_labels: list[np.ndarray] = []
    for j, (wj, sp) in enumerate(zip(scales, seg_params_list)):
        labels_j = segment_starlet_scale_region_growing(
            cube_xyl, wj,
            spectra_source=sp['spectra_source'],
            saliency_mode=sp['saliency_mode'],
            seed_quantile=sp['seed_quantile'],
            connectivity=sp['connectivity'],
            cos_sim_min=sp['cos_sim_min'],
            min_region=sp['min_region'],
            normalize=sp['normalize'],
            merge_cos_min=sp['merge_cos_min'],
            return_debug=False,
        )
        per_scale_labels.append(labels_j.astype(np.int32, copy=False))

    # 4) Fuse all scales into a single global map
    labels_final, fusion_info = merge_scale_segmentations(
        per_scale_labels,
        scale_cubes=scales,               # let fusion compute per-pixel saliency from w_j
        saliency_mode=fusion_saliency_mode,
        method=fusion_method,
        priority=fusion_priority,
        normalize=fusion_normalize,
        relabel_compact=relabel_compact,
    )

    if not return_intermediates:
        return labels_final

    info = dict(
        scales=scales,
        residual=residual,
        per_scale_labels=per_scale_labels,
        fusion_info=fusion_info,
        starlet=dict(positive=positive, k_sigma=k_sigma, n_scales=n_scales),
        segmentation_params_resolved=seg_params_list,
    )
    return labels_final, info


def refine_labels_by_sed_adjacency(
    labels: np.ndarray,
    cube_xyl: np.ndarray,
    *,
    connectivity: int = 8,      # 4 or 8
    min_size: int | None = 20,  # only attempt to merge regions with size <= min_size (None => consider all)
    sim_threshold: float = 0.97,# cosine similarity threshold to merge (0..1)
    max_passes: int = 3,        # run multiple passes until no merges or limit reached
    prefer_larger: bool = True, # merge smaller -> larger neighbor (ties broken by higher similarity)
    normalize: str = "unit",    # 'unit' => per-pixel continuum removal (median) + L2 norm
    compact_ids: bool = True,   # re-pack labels to 1..K at the end
    return_info: bool = True,
):
    """
    Refine a global segmentation map by merging spectrally-indistinguishable
    touching regions using original SEDs.

    Parameters
    ----------
    labels : (nx, ny) int
        Global label map (0 = background).
    cube_xyl : (nx, ny, λ) float
        Original data cube with per-pixel spectra.
    connectivity : 4|8
        Pixel connectivity to define region adjacency.
    min_size : int or None
        Only consider merging regions with size <= min_size. Use None to consider all.
    sim_threshold : float
        Cosine similarity threshold between region mean SEDs to merge.
    max_passes : int
        Maximum refinement passes (recompute adjacency/centroids each pass).
    prefer_larger : bool
        If True, merge into the larger touching neighbor (most similar among them).
    normalize : {'unit', None}
        'unit': per-pixel continuum removal (subtract median over λ) and L2 normalize.
        None: use raw spectra (not recommended if scales differ).
    compact_ids : bool
        If True, remap final labels to 1..K.
    return_info : bool
        If True, also returns a dict with stats.

    Returns
    -------
    labels_refined : (nx, ny) int
    info : dict (if return_info)
        {'passes': int, 'merges': int, 'sim_threshold': float, 'connectivity': int}
    """
    lbl = np.asarray(labels).copy()
    assert lbl.ndim == 2 and cube_xyl.ndim == 3 and lbl.shape == cube_xyl.shape[:2]
    nx, ny, nl = cube_xyl.shape

    # ---------- helpers ----------
    def _build_adjacency(L: np.ndarray, conn: int = 8):
        """Return neighbors dict: label -> set(neighbor_labels), excluding 0."""
        H, W = L.shape
        nbrs = {int(k): set() for k in np.unique(L) if k != 0}
        # horizontal/vertical
        a = L[:, :-1]; b = L[:, 1:]
        mask = (a != b)
        if np.any(mask):
            A = a[mask]; B = b[mask]
            for u, v in zip(A, B):
                if u != 0 and v != 0 and u != v:
                    nbrs[int(u)].add(int(v)); nbrs[int(v)].add(int(u))
        a = L[:-1, :]; b = L[1:, :]
        mask = (a != b)
        if np.any(mask):
            A = a[mask]; B = b[mask]
            for u, v in zip(A, B):
                if u != 0 and v != 0 and u != v:
                    nbrs[int(u)].add(int(v)); nbrs[int(v)].add(int(u))
        if conn == 8:
            a = L[:-1, :-1]; b = L[1:, 1:]
            mask = (a != b)
            if np.any(mask):
                A = a[mask]; B = b[mask]
                for u, v in zip(A, B):
                    if u != 0 and v != 0 and u != v:
                        nbrs[int(u)].add(int(v)); nbrs[int(v)].add(int(u))
            a = L[:-1, 1:]; b = L[1:, :-1]
            mask = (a != b)
            if np.any(mask):
                A = a[mask]; B = b[mask]
                for u, v in zip(A, B):
                    if u != 0 and v != 0 and u != v:
                        nbrs[int(u)].add(int(v)); nbrs[int(v)].add(int(u))
        return nbrs

    def _normalize_pixels(X: np.ndarray, mode: str = "unit"):
        """X: (N, L) -> normalized (N, L)."""
        if mode == "unit":
            med = np.median(X, axis=1, keepdims=True)
            Xc = X - med
            norms = np.linalg.norm(Xc, axis=1, keepdims=True) + 1e-12
            return Xc / norms
        elif mode in (None, "none"):
            return X
        else:
            raise ValueError("normalize must be 'unit' or None")

    def _region_centroids(L: np.ndarray, cube: np.ndarray, mode: str = "unit"):
        """Return centroids dict: label -> (L-dim unit vector), and sizes dict."""
        H, W, C = cube.shape
        lab = L.ravel()
        mask = lab > 0
        uniq = np.unique(lab[mask])
        # map labels -> 0..K-1 compact indices
        label_to_idx = {int(k): i for i, k in enumerate(uniq)}
        idx = np.array([label_to_idx[int(k)] for k in lab[mask]], dtype=np.int32)
        X = cube.reshape(-1, C)[mask].astype(np.float64, copy=False)
        Xn = _normalize_pixels(X, mode=mode)
        # accumulate normalized spectra
        K = len(uniq)
        sums = np.zeros((K, C), dtype=np.float64)
        np.add.at(sums, idx, Xn)
        counts = np.bincount(idx, minlength=K).astype(np.int64)
        # mean normalized spectra -> re-normalize to unit vectors
        with np.errstate(invalid="ignore", divide="ignore"):
            means = sums / counts[:, None]
        norms = np.linalg.norm(means, axis=1, keepdims=True) + 1e-12
        means /= norms
        centroids = {int(k): means[i] for i, k in enumerate(uniq)}
        sizes = {int(k): int(counts[i]) for i, k in enumerate(uniq)}
        return centroids, sizes

    def _relabel(L: np.ndarray, mapping: dict[int, int]) -> np.ndarray:
        """Apply label mapping; supports chained mappings."""
        if not mapping:
            return L
        ids = np.unique(L)
        max_id = int(ids.max(initial=0))
        lut = np.arange(max_id + 1, dtype=np.int32)
        for src, dst in mapping.items():
            if src <= max_id:
                lut[src] = dst
        # resolve chains: lut = lut[lut] repeatedly
        for _ in range(3):
            lut = lut[lut]
        out = lut[L]
        return out

    def _compact(L: np.ndarray) -> np.ndarray:
        ids = np.unique(L)
        ids = ids[ids != 0]
        remap = {int(old): i+1 for i, old in enumerate(ids)}
        out = L.copy()
        for old, new in remap.items():
            out[L == old] = new
        return out

    # ---------- iterative refinement ----------
    total_merges = 0
    passes_done = 0
    for _pass in range(max_passes):
        centroids, sizes = _region_centroids(lbl, cube_xyl, mode=normalize)
        nbrs = _build_adjacency(lbl, connectivity)

        # candidate order: small regions first (to simulate "fine" features)
        region_ids = sorted(nbrs.keys(), key=lambda k: sizes.get(k, np.inf))
        if min_size is not None:
            region_ids = [k for k in region_ids if sizes.get(k, 0) <= min_size]

        merges = {}
        for r in region_ids:
            if r in merges:  # already scheduled for merge this pass
                continue
            neigh = list(nbrs.get(r, ()))
            if not neigh:
                continue
            # pick neighbor with highest cosine similarity
            cr = centroids.get(r)
            if cr is None:
                continue
            sims = []
            for n in neigh:
                cn = centroids.get(n)
                if cn is None: 
                    continue
                sims.append((float(np.dot(cr, cn)), n))
            if not sims:
                continue
            sims.sort(reverse=True, key=lambda x: x[0])
            best_sim, best_n = sims[0]
            if best_sim < sim_threshold:
                continue
            # choose merge direction
            if prefer_larger and sizes.get(best_n, 0) >= sizes.get(r, 0):
                src, dst = r, best_n
            else:
                # if not prefer_larger, merge r -> most similar neighbor
                src, dst = r, best_n
            merges[src] = dst

        if not merges:
            break

        # apply merges
        lbl = _relabel(lbl, merges)
        # optional quick compress to remove holes in ids (helps next pass performance)
        lbl = _compact(lbl)
        total_merges += len(merges)
        passes_done += 1

    if compact_ids:
        lbl = _compact(lbl)

    if return_info:
        return lbl.astype(np.int32, copy=False), {
            "passes": passes_done,
            "merges": total_merges,
            "sim_threshold": sim_threshold,
            "connectivity": connectivity,
            "min_size": min_size,
        }
    return lbl.astype(np.int32, copy=False)
