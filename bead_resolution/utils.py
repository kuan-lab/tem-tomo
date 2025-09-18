import webknossos as wk

import tempfile
import numpy as np
import webknossos as wk

import uuid
from pathlib import Path


from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
from collections.abc import Iterable
from typing import List, Sequence, Tuple, Union, Any

from pathlib import Path
from typing import Dict, Tuple
import numpy as np

from pathlib import Path
import re
import yaml
import pandas as pd
import numpy as np

def get_annotation_points(annotation):
    """
    Extract annotation points (as Vec3Int) from all trees in the skeleton.

    Returns:
        List of Vec3Int(x, y, z)
    """
    skeleton = annotation.skeleton
    node_coords = []
    for tree in skeleton.trees:
        for _, node_attrs in tree.nodes.items():
            node_coords.append(node_attrs['position'])  # Vec3Int(x, y, z)
    return node_coords

def compute_bounding_boxes_flex(
    annotation_points: List[wk.Vec3Int],
    box_size: int,
    dataset_shape: Tuple[int, int, int],   # (X, Y, Z)
) -> List[wk.BoundingBox]:
    """
    Return same-sized BoundingBoxes centred on the given points, shifting them
    inwards (not shrinking) if they would exceed the dataset limits.

    Parameters
    ----------
    annotation_points : list[wk.Vec3Int]
        Centre voxels for the boxes.
    box_size : int
        Cube side length in voxels. Must be ≤ dataset extent in every axis.
    dataset_shape : tuple[int, int, int]
        Full dataset dimensions (x, y, z).

    Returns
    -------
    list[wk.BoundingBox]
    """
    half = box_size // 2
    max_idx = [d - 1 for d in dataset_shape]  # last valid index
    boxes = []

    if any(box_size > dim for dim in dataset_shape):
        raise ValueError(
            f"box_size ({box_size}) larger than dataset extent {dataset_shape}"
        )

    for p in annotation_points:
        # initial top-left and bottom-right
        tl = [p.x - half, p.y - half, p.z - half]
        br = [p.x + half, p.y + half, p.z + half]

        # shift per axis so the box lies inside [0, max_idx]
        for axis in range(3):
            if tl[axis] < 0:
                # shift forward
                shift = -tl[axis]
                tl[axis] += shift
                br[axis] += shift
            elif br[axis] > max_idx[axis]:
                # shift backward
                shift = br[axis] - max_idx[axis]
                tl[axis] -= shift
                br[axis] -= shift
            # after shifting, tl ≥ 0 and br ≤ max_idx

        boxes.append(wk.BoundingBox(tl, [box_size] * 3))

    return boxes

def save_subvolumes_to_npz(subvolumes, annotation_points, output_path):
    """
    Save subvolumes and annotation points to a .npz file.

    Parameters:
        subvolumes: list of np.ndarray (Z, Y, X)
        annotation_points: list of Vec3Int
        output_path: path to save .npz
    """
    arr_stack = np.stack(subvolumes)  # Shape: (N, Z, Y, X)
    coords = np.array([[pt.x, pt.y, pt.z] for pt in annotation_points])
    np.savez_compressed(output_path, subvolumes=arr_stack, points=coords)
    print(f"Saved {len(subvolumes)} subvolumes to {output_path}")

def _safe_name(s: str) -> str:
    """Mirror the downloader's filename sanitization."""
    return "".join(c if c.isalnum() or c in ("_", "-", ".") else "_" for c in s)

def get_bead_views(
    tomogram: str,
    recon_label: str,
    bead_index: int,
    root_dir: Path | str = "data",
) -> Dict[str, np.ndarray]:
    """
    Return orthogonal mid-plane views (XY, YZ, ZX) for a single bead's ROI.

    Parameters
    ----------
    tomogram : str
        The base tomogram name (folder created by the downloader).
    recon_label : str
        The reconstruction label (the downloader saved as <safe_name(recon_label)>.npz).
        This can be the dataset name or URL you used; it'll be sanitized the same way.
    bead_index : int
        Index of the bead within this tomogram/reconstruction (0-based).
    root_dir : Path | str, optional
        Root directory where the NPZ outputs were written (default: "data").

    Returns
    -------
    dict[str, np.ndarray]
        A dict with keys 'XY', 'YZ', 'ZX'.
        Shapes:
            - XY: (Y, X)
            - YZ: (Y, Z)
            - ZX: (Z, X)
        Dtypes match the stored subvolume dtype.

    Notes
    -----
    - Subvolumes are stored as (Z, Y, X).
    - Slices are taken at the *center* plane along each axis.
      If you want bead-centered slicing (robust to edge-shifted ROIs), we can add
      a variant that computes the bead's local index inside the ROI.
    """
    root_dir = Path(root_dir)
    npz_path = root_dir / _safe_name(tomogram) / f"{_safe_name(recon_label)}.npz"
    if not npz_path.exists():
        raise FileNotFoundError(f"NPZ not found: {npz_path}")

    data = np.load(npz_path, allow_pickle=True)
    if "subvolumes" not in data or "points" not in data:
        raise KeyError(f"{npz_path} missing expected keys 'subvolumes' and 'points'.")

    subvols = data["subvolumes"]  # typically object array of (Z,Y,X) ndarrays
    points = data["points"]       # global voxel coords (Z,Y,X) for each bead

    n = len(subvols)
    if not (0 <= bead_index < n):
        raise IndexError(f"bead_index {bead_index} out of range [0, {n-1}] for {npz_path.name}")

    vol = subvols[bead_index]
    if vol.ndim != 3:
        raise ValueError(f"Expected 3D subvolume (Z,Y,X), got shape {vol.shape}.")

    z, y, x = vol.shape
    zc, yc, xc = z // 2, y // 2, x // 2

    # XY at center Z
    xy = vol[zc, :, :]            # (Y, X)

    # YZ at center X (transpose to (Y, Z))
    yz = vol[:, :, xc].transpose(1, 0)  # (Z, Y) -> (Y, Z)

    # ZX at center Y (already (Z, X))
    zx = vol[:, yc, :]            # (Z, X)

    return {"XY": xy, "YZ": yz, "ZX": zx}
    
def plot_bead_views(
    tomogram: str,
    recon_label: str,
    bead_index: int,
    root_dir: str | Path = "data",
    *,
    clip_percentile: tuple[float, float] | None = (1, 99),
    cmap: str = "gray",
    title: str | None = None,
    save_path: str | Path | None = None,
):
    """
    One-liner: load a bead's ROI and plot XY / YZ / ZX mid-plane views side-by-side.

    clip_percentile: (low, high) intensity percentiles for contrast (None = raw).
    cmap: matplotlib colormap name.
    title: optional figure title.
    save_path: if given, saves the figure (PNG) instead of only showing it.

    Returns (fig, axes).
    """
    import matplotlib.pyplot as plt
    from pathlib import Path
    import numpy as np

    views = get_bead_views(tomogram, recon_label, bead_index, root_dir=root_dir)
    XY, YZ, ZX = views["XY"], views["YZ"], views["ZX"]

    def _vrange(img):
        if clip_percentile is None:
            return None
        lo, hi = np.percentile(img, clip_percentile)
        return (float(lo), float(hi))

    fig, axes = plt.subplots(1, 3, figsize=(9, 3), constrained_layout=True)
    v_xy = _vrange(XY); v_yz = _vrange(YZ); v_zx = _vrange(ZX)

    axes[0].imshow(XY, cmap=cmap, origin="lower", vmin=None if v_xy is None else v_xy[0], vmax=None if v_xy is None else v_xy[1])
    axes[0].set_title("XY")
    axes[1].imshow(YZ, cmap=cmap, origin="lower", vmin=None if v_yz is None else v_yz[0], vmax=None if v_yz is None else v_yz[1])
    axes[1].set_title("YZ")
    axes[2].imshow(ZX, cmap=cmap, origin="lower", vmin=None if v_zx is None else v_zx[0], vmax=None if v_zx is None else v_zx[1])
    axes[2].set_title("ZX")

    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([])

    if title is None:
        title = f"{tomogram} | {recon_label} | bead {bead_index}"
    fig.suptitle(title, fontsize=10)

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=200, bbox_inches="tight")

    return fig, axes

def extract_center_profiles_old(vol: np.ndarray, *, linewidth: int = 1
                            ) -> dict[str, np.ndarray]:
    """Return 1-D X/Y/Z profiles through the geometric centre of *vol*.
    `thickness` = number of pixels to average *perpendicularly* to each line."""
    z_sz, y_sz, x_sz = vol.shape
    zc, yc, xc = z_sz // 2, y_sz // 2, x_sz // 2

    if linewidth <= 1:           # single-voxel lines s
        return {"Z": vol[:,  yc,     xc],
                "Y": vol[zc, :,      xc],
                "X": vol[zc,  yc,    :]}
    half = linewidth // 2
    slc = lambda c, sz: slice(max(c - half, 0), min(c + half + 1, sz))
    return {"Z": vol[:,                slc(yc, y_sz), slc(xc, x_sz)].mean((1, 2)),
            "Y": vol[slc(zc, z_sz),    :,              slc(xc, x_sz)].mean((0, 2)),
            "X": vol[slc(zc, z_sz),    slc(yc, y_sz), :].mean((0, 1))}

def gaussian(x, A, x0, sigma, offset):
    """Simple 1-D Gaussian."""
    return A * np.exp(-((x - x0) ** 2) / (2 * sigma ** 2)) + offset

def fit_gaussian_and_compute_fwhm(profile: np.ndarray, *,
                                  bead_diam_nm: float = 0.0,
                                  pixel_size_nm: float | None = None,
                                  title: str = "",
                                  plot_fit: bool = False
                                  ) -> tuple[float, list[float], float]:
    """Fit a 1‑D profile with a Gaussian and compute (optionally deconvolved) FWHM.

    Parameters
    ----------
    profile
        1‑D or (linewidth × length) 2‑D intensity array.
    bead_diam_nm
        *D* — Physical diameter of the bead used for deconvolution **(nm)**.
        If ≤0, the raw measured FWHM is returned.
    pixel_size_nm
        Pixel size in nm / px. Required if *bead_diam_nm* > 0.
    plot_fit
        Plot data + Gaussian fit.
    title
        Figure title when *plot_fit* is True.

    Returns
    -------
    fwhm_effective_px : float
        FWHM after deconvolution if *bead_diam_nm* > 0, else the measured FWHM.
        Returned **in pixels** to keep downstream code unchanged.
    popt : list[float]
        Gaussian parameters ``[A, x0, sigma, offset]``.
    r_squared : float
        Coefficient of determination for the fit.
    """
    # ­­­Ensure NumPy array of floats for safe calculations
    profile_1d = np.asarray(profile, dtype=float)

    x = np.arange(profile_1d.size)
    x_dense = np.linspace(float(x.min()), float(x.max()), num=int((x.max() - x.min()) * 10) + 1)

    # Initial parameter guesses (per user‑supplied scheme)
    A_guess = float(np.max(profile_1d) - np.min(profile_1d))
    x0_guess = float(np.argmin(profile_1d))
    sigma_guess = 3.0  # pixels
    offset_guess = float(np.max(profile_1d))
    p0 = [A_guess, x0_guess, sigma_guess, offset_guess]

    try:
        popt, _ = curve_fit(gaussian, x, profile_1d, p0=p0)
        _, x0, sigma, _ = popt
        fwhm_meas_px = 2.355 * abs(sigma)

        # Default effective FWHM is the measured one
        fwhm_eff_px: float = fwhm_meas_px

        # Optional deconvolution
        if bead_diam_nm > 0:
            if pixel_size_nm is None or pixel_size_nm <= 0:
                raise ValueError("pixel_size_nm must be >0 when deconvolving with bead size")
            fwhm_psf_nm = estimate_psf_fwhm_from_top_hat(
                fwhm_meas_px, bead_diam_nm, pixel_size_nm
            )
            fwhm_eff_px = float(fwhm_psf_nm / pixel_size_nm)

        # Goodness of fit (R²)
        fitted = gaussian(x, *popt)
        ss_res = float(np.sum((profile_1d - fitted) ** 2))
        ss_tot = float(np.sum((profile_1d - profile_1d.mean()) ** 2))
        r_squared = 1 - ss_res / ss_tot if ss_tot else np.nan

        if plot_fit:
            plt.plot(x, profile_1d, "b.", label="Data")
            plt.plot(
                x_dense,
                gaussian(x_dense, *popt),
                "r-",
                label=(
                    f"Gaussian Fit\nFWHM = {fwhm_eff_px:.2f} px\nR² = {r_squared:.3f}"
                ),
            )
            plt.title(title)
            plt.xlabel("Position (px)")
            plt.ylabel("Intensity")
            #plt.legend() 
            plt.legend(bbox_to_anchor=(2, 0.5), loc='center right', borderaxespad=0.)
            plt.grid(True)
            plt.show()
            
        return fwhm_eff_px, popt, r_squared

    except Exception as e:
        print(f"⚠️ Gaussian fit failed: {e}")
        return np.nan, [np.nan] * 4, np.nan

def _trim_profile_nan_edges(profile: np.ndarray, *, min_len: int = 7) -> np.ndarray:
    """Keep the contiguous finite block from first to last finite sample."""
    y = np.asarray(profile, dtype=float)
    finite = np.isfinite(y)
    if not finite.any():
        return y[:0]
    i0 = int(np.argmax(finite))                          # first True
    i1 = len(y) - int(np.argmax(finite[::-1]))          # one past last True
    seg = y[i0:i1]
    return seg if seg.size >= min_len else y[:0]


def spot_check_fit_minimal(
    tomogram: str,
    recon_label: str,
    bead_index: int,
    root_dir: str | Path = "data",
    *,
    axes_to_show: Iterable[str] = ("X", "Z"),
    linewidth: int = 1,
    bead_diam_nm: float = 0.0,          # usually 0 (no deconv)
    pixel_size_nm: float | None = None, # needed only if bead_diam_nm > 0
) -> Dict[str, Dict[str, Any]]:
    """
    NaN-aware spot-check:
      - loads ROI,
      - extracts 1-D center profiles with your extract_center_profiles_old,
      - trims each profile to the finite run around its minimum (handles NaN padding),
      - calls your fit_gaussian_and_compute_fwhm (which plots per-axis),
      - returns per-axis (fwhm_px, popt, r2).

    Note: if you use linewidth > 1, the old extractor averages with np.mean and can
    still produce NaNs. Prefer linewidth=1 with NaN-padded ROIs, or we can add a
    nan-aware extractor later.
    """
    import matplotlib.pyplot as plt

    npz_path = Path(root_dir) / _safe_name(tomogram) / f"{_safe_name(recon_label)}.npz"
    data = np.load(npz_path, allow_pickle=True)
    subvols = data["subvolumes"]
    if not (0 <= bead_index < len(subvols)):
        raise IndexError(f"bead_index {bead_index} out of range [0, {len(subvols)-1}] for {npz_path.name}")

    vol = subvols[bead_index]
    if vol.ndim != 3:
        raise ValueError(f"Expected 3D (Z,Y,X), got {vol.shape}")

    # centerline profiles (Z,Y,X) -> dict {"X","Y","Z"}
    profiles = extract_center_profiles_old(vol, linewidth=linewidth)

    results: Dict[str, Dict[str, Any]] = {}
    for axis in axes_to_show:
        if axis not in ("X", "Y", "Z"):
            continue

        raw = profiles[axis]
        cleaned = _trim_profile_nan_edges(raw, min_len=7)

        title = f"{tomogram} | {recon_label} | bead {bead_index} | axis {axis}"
        if cleaned.size == 0:
            # not enough finite samples to fit
            print(f"⚠️  Not enough finite samples to fit ({axis}); skipping.")
            results[axis] = {"fwhm_px": np.nan, "popt": [np.nan]*4, "r2": np.nan}
            continue

        fwhm_px, popt, r2 = fit_gaussian_and_compute_fwhm(
            cleaned,
            bead_diam_nm=bead_diam_nm,
            pixel_size_nm=pixel_size_nm if bead_diam_nm > 0 else None,
            title=title + f" (n={cleaned.size})",
            plot_fit=True,
        )
        results[axis] = {"fwhm_px": fwhm_px, "popt": popt, "r2": r2}

    return results

def generate_skip_bead_indices(
    r2x_list: Sequence[Sequence[float]],
    r2z_list: Sequence[Sequence[float]],
    *,
    r2_cutoff: float = 0.9,
) -> List[int]:
    """Return a **single 1‑D list of bead indices to skip across *all* datasets.*

    We assume every dataset shares the same bead ordering (i.e. bead index 5
    refers to the same physical bead in every dataset).  A bead index appears
    in the skip list if **any** dataset has an X‑ or Z‑axis Gaussian fit with
    R² < `r2_cutoff` (or non‑finite) for that bead.  Duplicate indices are
    removed automatically.
    """
    if not r2x_list:
        return []

    # Longest dataset defines the maximum bead index
    n_beads = max(len(ds) for ds in r2x_list)
    skip: set[int] = set()

    for bead_idx in range(n_beads):
        for rx_ds, rz_ds in zip(r2x_list, r2z_list):
            #if bead_idx >= len(rx_ds):
            #    continue  # this dataset has fewer beads; ignore
            rx = rx_ds[bead_idx]
            rz = rz_ds[bead_idx]
            if (not np.isfinite(rx)) or (not np.isfinite(rz)) or rx < r2_cutoff or rz < r2_cutoff:
                skip.add(bead_idx)
                break  # no need to check other datasets for this bead
    
    return sorted(skip)

def filter_by_index_excluding(nested_lists, skip_indices):
    """Keep only elements whose index is NOT in skip_indices, for every sublist."""
    skip = set(skip_indices)
    return [[v for i, v in enumerate(sub) if i not in skip] for sub in nested_lists]


def parse_recon_label(recon_label: str) -> tuple[float, float]:
    """
    Extract (num_proj, max_angle) from a recon label, robust to prefixes and URLs.

    Matches the **last** occurrence of "<num>[_-]*lim<num>" so it works for:
      - "11_lim15_tomo10a_16bit"          -> (11, 15)
      - "Mouse..._tomo11a_11_lim15_16bit" -> (11, 15)
      - full URLs ".../Mouse..._11_lim15_16bit#view..." -> (11, 15)

    Returns (nan, nan) if not found.
    """
    # keep only the last path segment; strip fragments/query
    base = recon_label.rsplit("/", 1)[-1].split("#", 1)[0].split("?", 1)[0]

    # find the LAST "<digits>[optional _ or -]lim<digits>"
    pat = re.compile(r"(?P<num_proj>\d+)\s*[_-]*\s*lim\s*(?P<max_ang>\d+)", re.IGNORECASE)
    match = None
    for m in pat.finditer(base):
        match = m
    if match:
        return float(match.group("num_proj")), float(match.group("max_ang"))
    return float("nan"), float("nan")
# --- tiny parser for recon label -> (num_proj, max_angle) ---
def parse_recon_label_old(recon_label: str) -> tuple[float, float]:
    """
    Extract (num_proj, max_angle) from labels like:
      '11_lim15_tomo10a_16bit' or '..._41_lim60_tomo10a_16bit'
    Works even if recon_label is a full URL (uses last path segment).
    Returns (np.nan, np.nan) if not found.
    """
    base = recon_label.split("/")[-1]
    m = re.search(r"(\d+)_lim(\d+)", base)
    if m:
        return float(m.group(1)), float(m.group(2))
    # fallback: any number before 'lim' + limXX
    m = re.search(r"(\d+)[^\d]+lim(\d+)", base)
    if m:
        return float(m.group(1)), float(m.group(2))
    return float("nan"), float("nan")

# --- loader for a recon's saved NPZ (produced by the downloader) ---
def _load_recon_npz(root_dir: Path | str, tomogram: str, recon_label: str):
    p = Path(root_dir) / _safe_name(tomogram) / f"{_safe_name(recon_label)}.npz"
    if not p.exists():
        raise FileNotFoundError(f"Missing NPZ: {p}")
    d = np.load(p, allow_pickle=True)
    return d["subvolumes"], d.get("points", None)

# ---------------- MAIN: build a tidy DataFrame of bead fits ----------------
def fits_dataframe_from_config(
    config_path: Path | str,
    *,
    root_dir: Path | str = "data",
    linewidth: int = 1,
    r2_cutoff: float = 0.6,
    bead_diam_nm: float = 0.0,           # usually 0 (no deconvolution)
) -> pd.DataFrame:
    """
    Iterate all tomograms in YAML, fit X/Y/Z for every bead in every recon,
    compute per-tomogram skip indices from R² (X & Z), and return a tidy DataFrame
    with one row per (kept) bead fit. NaN-padded ROIs are handled by trimming each
    1-D profile to the finite run around its minimum before fitting.
    """
    cfg = yaml.safe_load(Path(config_path).read_text())
    rows: list[dict] = []

    for entry in cfg.get("tomograms", []):
        tomogram = entry["tomogram"]
        thickness = float(entry.get("thickness", float("nan")))
        pxl_size = float(entry.get("pxl_size", float("nan")))
        recon_labels = list(entry["recon_labels"])

        all_fwhms_x, all_fwhms_y, all_fwhms_z = [], [], []
        all_r2s_x,   all_r2s_y,   all_r2s_z   = [], [], []
        beads_per_recon: list[int] = []

        for recon in recon_labels:
            subvols, _ = _load_recon_npz(root_dir, tomogram, recon)
            beads_per_recon.append(len(subvols))

            fwhm_x, fwhm_y, fwhm_z = [], [], []
            r2_x,   r2_y,   r2_z   = [], [], []

            for vol in subvols:
                # Get geometric-center profiles (may include NaNs at the ends)
                profs = extract_center_profiles_old(vol, linewidth=linewidth)

                # NaN-aware trimming around the (finite) minimum for each axis
                px = _trim_profile_nan_edges(profs["X"], min_len=7)
                py = _trim_profile_nan_edges(profs["Y"], min_len=7)
                pz = _trim_profile_nan_edges(profs["Z"], min_len=7)

                # Fit each axis if we have enough finite samples; else record NaNs
                if px.size:
                    fx, _, rx = fit_gaussian_and_compute_fwhm(
                        px,
                        bead_diam_nm=bead_diam_nm,
                        pixel_size_nm=pxl_size if bead_diam_nm > 0 else None,
                        title="",
                        plot_fit=False,
                    )
                else:
                    fx, rx = np.nan, np.nan

                if py.size:
                    fy, _, ry = fit_gaussian_and_compute_fwhm(
                        py,
                        bead_diam_nm=bead_diam_nm,
                        pixel_size_nm=pxl_size if bead_diam_nm > 0 else None,
                        title="",
                        plot_fit=False,
                    )
                else:
                    fy, ry = np.nan, np.nan

                if pz.size:
                    fz, _, rz = fit_gaussian_and_compute_fwhm(
                        pz,
                        bead_diam_nm=bead_diam_nm,
                        pixel_size_nm=pxl_size if bead_diam_nm > 0 else None,
                        title="",
                        plot_fit=False,
                    )
                else:
                    fz, rz = np.nan, np.nan

                fwhm_x.append(fx); fwhm_y.append(fy); fwhm_z.append(fz)
                r2_x.append(rx);   r2_y.append(ry);   r2_z.append(rz)

            all_fwhms_x.append(fwhm_x); all_fwhms_y.append(fwhm_y); all_fwhms_z.append(fwhm_z)
            all_r2s_x.append(r2_x);     all_r2s_y.append(r2_y);     all_r2s_z.append(r2_z)

        # --- per-tomogram counts & skip set (unchanged) ---
        n_beads_total = int(beads_per_recon[0]) if beads_per_recon else 0
        if any(n != n_beads_total for n in beads_per_recon):
            print(f"⚠️ bead count mismatch across recons for {tomogram}: {beads_per_recon} "
                  f"(using first: {n_beads_total})")

        skip_indices = generate_skip_bead_indices(all_r2s_x, all_r2s_z, r2_cutoff=r2_cutoff)
        kept_indices = [i for i in range(n_beads_total) if i not in set(skip_indices)]
        n_beads_pass = len(kept_indices)

        pct = (100.0 * n_beads_pass / n_beads_total) if n_beads_total else 0.0
        print(f"[{tomogram}] beads: total={n_beads_total}, pass (R²≥{r2_cutoff})={n_beads_pass} ({pct:.1f}%)")

        # --- emit rows for every kept bead in every recon (unchanged except bead_id) ---
        for r_idx, recon in enumerate(recon_labels):
            num_proj, max_angle = parse_recon_label(recon)
            fxs = all_fwhms_x[r_idx]; fys = all_fwhms_y[r_idx]; fzs = all_fwhms_z[r_idx]

            for bead_id in kept_indices:
                if bead_id >= len(fxs) or bead_id >= len(fys) or bead_id >= len(fzs):
                    continue
                fx, fy, fz = fxs[bead_id], fys[bead_id], fzs[bead_id]
                rows.append({
                    "tomogram": tomogram,
                    "thickness": thickness,
                    "pxl_size": pxl_size,
                    "recon_label": recon,
                    "num_proj": num_proj,
                    "max_angle": max_angle,
                    "bead_id": bead_id,              # 0-based, stable within tomogram
                    "n_beads_total": n_beads_total,  # same for all rows of this tomogram
                    "n_beads_pass": n_beads_pass,    # same for all rows of this tomogram
                    "fwhm_x": fx,
                    "fwhm_y": fy,
                    "fwhm_z": fz,
                    "res_nm_x": fx * pxl_size if np.isfinite(fx) and np.isfinite(pxl_size) else np.nan,
                    "res_nm_y": fy * pxl_size if np.isfinite(fy) and np.isfinite(pxl_size) else np.nan,
                    "res_nm_z": fz * pxl_size if np.isfinite(fz) and np.isfinite(pxl_size) else np.nan,
                })

    return pd.DataFrame(rows)


def plot_fwhm_z_summary_line_by_thickness(
    df,
    *,
    units: str = "nm",                 # "nm" or "px"
    max_angs: list[int] | None = None, # x-axis order; if None, inferred
    theta_labels: list[str] | None = None,
    agg: str = "median",               # "median" or "mean"
    figsize: tuple[int, int] = (4, 2),
):
    """
    Single-axis plot: FWHM_Z vs tilt range, one line per thickness.

    df must have columns: ['thickness','pxl_size','max_angle','fwhm_z'].
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    if max_angs is None:
        max_angs = sorted(df["max_angle"].dropna().unique().tolist())
    if theta_labels is not None and len(theta_labels) != len(max_angs):
        raise ValueError("theta_labels must have same length as max_angs")

    # choose y column / label
    dfx = df.copy()
    if units == "nm":
        dfx["fwhm_z_nm"] = dfx["fwhm_z"] * dfx["pxl_size"]
        ycol, ylab = "fwhm_z_nm", "Resolution (nm)"
    elif units == "px":
        ycol, ylab = "fwhm_z", "Resolution (pixels)"
    else:
        raise ValueError("units must be 'nm' or 'px'")

    # aggregate across beads (and tomograms) for each (thickness, angle)
    agg_fun = {"median": np.nanmedian, "mean": np.nanmean}[agg]
    g = (dfx.groupby(["thickness", "max_angle"], as_index=False)[ycol]
            .agg(value=agg_fun)
            .rename(columns={"value": ycol}))

    # enforce x order
    g = g[g["max_angle"].isin(max_angs)].copy()
    g["max_angle"] = pd.Categorical(g["max_angle"], categories=max_angs, ordered=True)
    g.sort_values(["thickness", "max_angle"], inplace=True)

    # plot
    sns.set_theme(rc={"lines.linewidth": 0.9}, style="ticks")
    fig, ax = plt.subplots(1, 1, figsize=figsize, constrained_layout=True)

    for t, sub in g.groupby("thickness", sort=True):
        ax.plot(sub["max_angle"].cat.codes, sub[ycol], marker="o", label=f"{t:.0f} nm")

    xticks = list(range(len(max_angs)))
    ax.set_xticks(xticks)
    if theta_labels is None:
        theta_labels = [f"±{int(a)}°" for a in max_angs]
    ax.set_xticklabels(theta_labels)

    ax.set_xlabel("Tilt Range")
    ax.set_ylabel(ylab)
    ax.legend(title="Thickness", fontsize=8)
    sns.despine(fig)

    return fig, ax
