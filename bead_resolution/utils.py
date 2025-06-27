import webknossos as wk

import tempfile
import numpy as np
import webknossos as wk

import uuid
from pathlib import Path


from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
from collections.abc import Iterable
from typing import List, Sequence, Tuple, Union

# -----------------------------------
# Step 1: Extract annotation points
# -----------------------------------

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


# -----------------------------------
# Step 2: Compute bounding boxes
# -----------------------------------
import webknossos as wk
from typing import List, Tuple

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


def compute_bounding_boxes(annotation_points, box_size):
    """
    Compute BoundingBoxes centered at each annotation point.

    Parameters:
        annotation_points: list of Vec3Int
        box_size: int (e.g., 15)

    Returns:
        List of wk.BoundingBox objects
    """
    half = box_size // 2
    bounding_boxes = []

    for point in annotation_points:
        x, y, z = point.x, point.y, point.z  # Vec3Int from webknossos
        topleft = [x - half, y - half, z - half]
        size = [box_size, box_size, box_size]
        bounding_boxes.append(wk.BoundingBox(topleft, size))

    return bounding_boxes


# -----------------------------------
# Step 3: Download individual subvolumes
# -----------------------------------

def download_subvolumes_individually(
    dataset_name,
    bounding_boxes,
    api_token,
    organization_id,
    webknossos_url="https://webknossos.org",
    mag_level=1
):
    subvolumes = []

    with wk.webknossos_context(url=webknossos_url, token=api_token):
        remote_dataset = wk.Dataset.open_remote(dataset_name, organization_id=organization_id)
        available_layers = list(remote_dataset.layers.keys())
        if not available_layers:
            raise ValueError("No layers found in dataset.")
        layer_name = available_layers[0]
        print(f"Using layer: {layer_name}")

        for i, box in enumerate(bounding_boxes):
            #print('Downloading box %i' % i)
            try:
                # Use a truly unique temp directory
                temp_path = Path(tempfile.gettempdir()) / f"wk_dl_box_{i}_{uuid.uuid4().hex}"

                dataset = wk.Dataset.download(
                    dataset_name_or_url=dataset_name,
                    organization_id=organization_id,
                    bbox=box,
                    layers=[layer_name],
                    mags=[wk.Mag(mag_level)],
                    webknossos_url=webknossos_url,
                    path=str(temp_path)  # Must be string path
                )

                data = dataset.get_layer(layer_name).get_mag(wk.Mag(mag_level)).read()
                data = data.squeeze(0)         # (1, Y, X, Z) → (Y, X, Z)
                data = data.transpose(2, 0, 1) # → (Z, Y, X)
                subvolumes.append(data)

                # Optionally clean up the folder
                # import shutil; shutil.rmtree(temp_path)

            except Exception as e:
                print(f"⚠️ Failed to download box {i}: {e}")

    return subvolumes




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


def load_subvolumes_from_npz(path):
    """
    Load subvolumes and annotation points from .npz.

    Returns:
        (list of subvolumes, list of Vec3Int-like point arrays)
    """
    npz = np.load(path)
    subvolumes = list(npz['subvolumes'])  # Shape: (N, Z, Y, X)
    points = [wk.Vec3Int(int(x), int(y), int(z)) for x, y, z in npz['points']]
    return subvolumes, points


def plot_bead_triplet_montage(subvolumes, title="Bead Z/Y/X Montage", cmap='gray'):
    """
    Plot a montage where each bead gets 3 views: Z, Y, and X center slices.

    Parameters:
        subvolumes: list of np.ndarray of shape (Z, Y, X)
        title: Title of the whole figure
    """
    n = len(subvolumes)
    cols = 3  # Z, Y, X
    rows = n

    fig, axes = plt.subplots(rows, cols, figsize=(cols * 3, rows * 3))
    if rows == 1:
        axes = np.expand_dims(axes, 0)  # ensure 2D array even if one bead

    for i, vol in enumerate(subvolumes):
        z, y, x = vol.shape
        z_slice = vol[z // 2, :, :]
        y_slice = vol[:, y // 2, :]
        x_slice = vol[:, :, x // 2]

        axes[i][0].imshow(z_slice, cmap=cmap)
        axes[i][0].set_title(f"Bead {i} - Z")
        axes[i][1].imshow(y_slice, cmap=cmap)
        axes[i][1].set_title(f"Bead {i} - Y")
        axes[i][2].imshow(x_slice, cmap=cmap)
        axes[i][2].set_title(f"Bead {i} - X")

        for ax in axes[i]:
            ax.axis('off')

    plt.suptitle(title, fontsize=16)
    plt.tight_layout()
    plt.show()


def gaussian(x, A, x0, sigma, offset):
    """Simple 1-D Gaussian."""
    return A * np.exp(-((x - x0) ** 2) / (2 * sigma ** 2)) + offset

def _moving_average_edge_safe(arr: np.ndarray, window: int) -> np.ndarray:
    """Return a centered moving‑average that never pads with zeros.

    The window is symmetric; at the boundaries the effective window shrinks
    so only *real* pixels contribute to the average (no artificial zeros).
    """
    if window <= 1:
        return arr.astype(float)

    # Calculate padding sizes so the convolution output length matches input.
    pad_left = window // 2
    pad_right = window - 1 - pad_left

    # Pad using edge values rather than zeros to avoid underestimating edges.
    padded = np.pad(arr, (pad_left, pad_right), mode="edge")
    kernel = np.ones(window, dtype=float) / float(window)
    return np.convolve(padded, kernel, mode="valid")


def fit_gaussian_and_compute_fwhm(
    profile: np.ndarray,
    *,
    linewidth: int = 1,
    bead_diam_nm: float = 0.0,
    pixel_size_nm: float | None = None,
    plot_fit: bool = False,
    title: str = "",
) -> Tuple[float, List[float], float]:
    """Fit a 1‑D profile with a Gaussian and compute (optionally deconvolved) FWHM.

    Parameters
    ----------
    profile
        1‑D or (linewidth × length) 2‑D intensity array.
    linewidth
        Number of adjacent rows/columns to average for a thicker line profile.
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
    profile = np.asarray(profile, dtype=float)

    # Collapse to 1‑D according to requested linewidth
    if profile.ndim == 2:
        if linewidth > profile.shape[0]:
            raise ValueError("linewidth exceeds available profile thickness")
        profile_1d = profile[:linewidth, :].mean(axis=0)
    else:
        profile_1d = (
            _moving_average_edge_safe(profile, int(linewidth))
            if linewidth > 1
            else profile
        )

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
            plt.legend() 
            plt.grid(True)
            plt.show()

        return fwhm_eff_px, popt, r_squared

    except Exception as e:
        print(f"⚠️ Gaussian fit failed: {e}")
        return np.nan, [np.nan] * 4, np.nan
    
def fit_gaussian_and_compute_fwhm_old2(profile, linewidth: int = 1, plot_fit: bool = False, title: str = ""):
    """
    Fit a 1‑D profile with a Gaussian and compute FWHM.

    Parameters
    ----------
    profile : np.ndarray
        Either a 1‑D array (intensity values) or a 2‑D array where the first
        dimension corresponds to line thickness (e.g., multiple rows sampled
        perpendicular to the profile direction).
    linewidth : int, optional
        Thickness of the line profile. A value of 1 leaves the profile as‑is.
        Values >1 average over that many pixels *without* introducing zeros at
        the edges (edge pixels are averaged only with existing neighbors).
    plot_fit : bool, optional
        Whether to plot the data and Gaussian fit.
    title : str, optional
        Plot title, used only when *plot_fit* is True.

    Returns
    -------
    fwhm : float
        Full‑Width‑at‑Half‑Maximum in pixels.
    popt : list[float]
        Fitted Gaussian parameters (A, x0, sigma, offset).
    r_squared : float
        Coefficient of determination describing goodness of fit.
    """
    # Ensure NumPy array of floats for safe calculations
    profile = np.asarray(profile, dtype=float)

    # Collapse to 1‑D according to requested linewidth
    if profile.ndim == 2:
        if linewidth > profile.shape[0]:
            raise ValueError("linewidth exceeds available profile thickness")
        profile_1d = profile[:linewidth, :].mean(axis=0)
    else:  # 1‑D input
        profile_1d = _moving_average_edge_safe(profile, int(linewidth)) if linewidth > 1 else profile

    x = np.arange(profile_1d.size)
    x_dense = np.linspace(x.min(), x.max(), num=int((x.max() - x.min()) * 10) + 1)

    # Initial guesses
    A_guess = np.max(profile_1d) - np.min(profile_1d)
    x0_guess = np.argmin(profile_1d)
    sigma_guess = 3
    offset_guess = np.max(profile_1d)
    p0 = [A_guess, x0_guess, sigma_guess, offset_guess]

    try:
        popt, _ = curve_fit(gaussian, x, profile_1d, p0=p0)
        _, x0, sigma, _ = popt
        fwhm = 2.355 * abs(sigma)  # Convert sigma to FWHM

        # Goodness of fit (R²)
        fitted = gaussian(x, *popt)
        ss_res = np.sum((profile_1d - fitted) ** 2)
        ss_tot = np.sum((profile_1d - profile_1d.mean()) ** 2)
        r_squared = 1 - ss_res / ss_tot if ss_tot else np.nan

        if plot_fit:
            plt.plot(x, profile_1d, "b.", label="Data")
            plt.plot(x_dense, gaussian(x_dense, *popt), "r-",
                     label=f"Gaussian Fit\nFWHM = {fwhm:.2f} px\nR² = {r_squared:.3f}")
            plt.title(title)
            plt.xlabel("Position (px)")
            plt.ylabel("Intensity")
            plt.legend()
            plt.grid(True)
            plt.show()

        return fwhm, popt, r_squared

    except Exception as e:
        print(f"⚠️ Gaussian fit failed: {e}")
        return None, None, None



def fit_gaussian_and_compute_fwhm_old(profile, plot_fit=False, title=""):
    """
    Fit a 1D profile with a Gaussian and compute FWHM.

    Parameters:
        profile: 1D np.array of intensities
        plot_fit: bool, whether to plot data and fit
        title: str, plot title if plotting

    Returns:
        fwhm: Full Width Half Maximum
        popt: fitted Gaussian parameters (A, x0, sigma, offset)
    """
    x = np.arange(len(profile))
    x2 = np.arange(np.min(x), np.max(x), .1)

    # Initial guesses
    A_guess = np.max(profile) - np.min(profile)
    x0_guess = np.argmin(profile)
    sigma_guess = 3 #len(profile) / 5
    offset_guess = np.max(profile)
    p0 = [A_guess, x0_guess, sigma_guess, offset_guess]

    try:
        popt, _ = curve_fit(gaussian, x, profile, p0=p0)
        _, x0, sigma, _ = popt
        fwhm = 2.355 * abs(sigma)  # FWHM formula

        if plot_fit:
            plt.plot(x, profile, 'b.', label='Data')
            plt.plot(x2, gaussian(x2, *popt), 'r-', label=f'Gaussian Fit\nFWHM={fwhm:.2f} px')
            plt.title(title)
            plt.xlabel('Position (px)')
            plt.ylabel('Intensity')
            plt.legend()
            plt.grid(True)
            plt.show()

        return fwhm, popt

    except Exception as e:
        print(f"⚠️ Gaussian fit failed: {e}")
        return None, None


def compute_fwhms_from_subvolumes(
    subvolumes: Iterable[np.ndarray],
    *,
    linewidth: int = 1,
    bead_diam_nm: float = 0.0,
    pixel_size_nm: float | None = None,
    r2_cutoff: float = 0.0,
) -> Tuple[List[float], List[float], List[float], List[float], List[float], List[float]]:
    """Compute FWHM along X/Y/Z for each bead, **optionally rejecting poor fits**.

    A bead (sub‑volume) is *kept* only if its X‑ and Z‑axis Gaussian fits both
    have R² ≥ `r2_cutoff`.  Y‑axis R² is ignored when deciding inclusion, but the
    corresponding Y FWHM is returned for beads that pass the X+Z criterion.

    Parameters
    ----------
    subvolumes : Iterable[np.ndarray]
        Each volume has shape (Z, Y, X).
    linewidth : int, optional
        Profile‑averaging thickness.
    bead_diam_nm : float, optional
        Bead diameter used for deconvolution; ≤0 disables deconvolution.
    pixel_size_nm : float, optional
        Pixel size (nm/px); required if deconvolving.
    r2_cutoff : float, optional
        Minimum R² required **simultaneously** for X and Z fits. Defaults to 0
        (no filtering).

    Returns
    -------
    fwhms_x, fwhms_y, fwhms_z, r2s_x, r2s_y, r2s_z : list[float]
        Lists contain only the beads that passed the R² criterion.
    """
    fwhms_x: List[float] = []
    fwhms_y: List[float] = []
    fwhms_z: List[float] = []
    r2s_x: List[float] = []
    r2s_y: List[float] = []
    r2s_z: List[float] = []

    for vol in subvolumes:
        z, y, x = vol.shape
        profile_z = vol[:, y // 2, x // 2]
        profile_y = vol[z // 2, :, x // 2]
        profile_x = vol[z // 2, y // 2, :]

        fwhm_z, _, r2_z = fit_gaussian_and_compute_fwhm(
            profile_z,
            linewidth=linewidth,
            bead_diam_nm=bead_diam_nm,
            pixel_size_nm=pixel_size_nm,
        )
        fwhm_x, _, r2_x = fit_gaussian_and_compute_fwhm(
            profile_x,
            linewidth=linewidth,
            bead_diam_nm=bead_diam_nm,
            pixel_size_nm=pixel_size_nm,
        )
        # Y‑axis (computed even though it doesn't influence inclusion)
        fwhm_y, _, r2_y = fit_gaussian_and_compute_fwhm(
            profile_y,
            linewidth=linewidth,
            bead_diam_nm=bead_diam_nm,
            pixel_size_nm=pixel_size_nm,
        )

        keep = (
            (np.isfinite(r2_x) and np.isfinite(r2_z))
            and (r2_x >= r2_cutoff and r2_z >= r2_cutoff)
        )
        if keep:
            fwhms_x.append(fwhm_x)
            fwhms_y.append(fwhm_y)
            fwhms_z.append(fwhm_z)
            r2s_x.append(r2_x)
            r2s_y.append(r2_y)
            r2s_z.append(r2_z)

    return fwhms_x, fwhms_y, fwhms_z, r2s_x, r2s_y, r2s_z

def compute_fwhms_from_subvolumes_old(
    subvolumes: Iterable[np.ndarray],
    *,
    linewidth: int = 1,
    bead_diam_nm: float = 0.0,
    pixel_size_nm: float | None = None,
) -> Tuple[List[float], List[float], List[float], List[float], List[float], List[float]]:
    """Measure (and optionally deconvolve) FWHM along X/Y/Z for each sub‑volume.

    Parameters
    ----------
    subvolumes
        Sequence of volumes with shape ``(Z, Y, X)``.
    linewidth
        Thickness for the line‑profile average.
    bead_diam_nm
        Bead diameter *D* (nm). If ≤0, no deconvolution.
    pixel_size_nm
        Pixel size (nm/px). Required if *bead_diam_nm* > 0.

    Returns
    -------
    (fwhms_x, fwhms_y, fwhms_z, r2s_x, r2s_y, r2s_z)
        Six lists matching *subvolumes* order.
    """
    fwhms_x: List[float] = []
    fwhms_y: List[float] = []
    fwhms_z: List[float] = []
    r2s_x: List[float] = []
    r2s_y: List[float] = []
    r2s_z: List[float] = []

    for vol in subvolumes:
        z, y, x = vol.shape

        profile_z = vol[:, y // 2, x // 2]
        profile_y = vol[z // 2, :, x // 2]
        profile_x = vol[z // 2, y // 2, :]

        fwhm_z, _, r2_z = fit_gaussian_and_compute_fwhm(
            profile_z,
            linewidth=linewidth,
            bead_diam_nm=bead_diam_nm,
            pixel_size_nm=pixel_size_nm,
        )
        fwhm_y, _, r2_y = fit_gaussian_and_compute_fwhm(
            profile_y,
            linewidth=linewidth,
            bead_diam_nm=bead_diam_nm,
            pixel_size_nm=pixel_size_nm,
        )
        fwhm_x, _, r2_x = fit_gaussian_and_compute_fwhm(
            profile_x,
            linewidth=linewidth,
            bead_diam_nm=bead_diam_nm,
            pixel_size_nm=pixel_size_nm,
        )

        fwhms_z.append(fwhm_z)
        fwhms_y.append(fwhm_y)
        fwhms_x.append(fwhm_x)

        r2s_z.append(r2_z)
        r2s_y.append(r2_y)
        r2s_x.append(r2_x)

    return fwhms_x, fwhms_y, fwhms_z, r2s_x, r2s_y, r2s_z


def compute_fwhms_from_subvolumes_old (subvolumes, linewidth = 1):
    fwhms_x, fwhms_y, fwhms_z = [], [], []
    r2s_x, r2s_y, r2s_z = [], [], []
    for vol in subvolumes:
        z, y, x = vol.shape

        profile_z = vol[:, y // 2, x // 2]
        profile_y = vol[z // 2, :, x // 2]
        profile_x = vol[z // 2, y // 2, :]

        fwhm_z, popt_z, r2_z = fit_gaussian_and_compute_fwhm(profile_z, linewidth = linewidth)
        fwhm_y, popt_y, r2_y = fit_gaussian_and_compute_fwhm(profile_y, linewidth = linewidth)
        fwhm_x, popt_x, r2_x = fit_gaussian_and_compute_fwhm(profile_x, linewidth = linewidth)

        fwhms_z.append(fwhm_z)
        fwhms_y.append(fwhm_y)
        fwhms_x.append(fwhm_x)

        r2s_z.append(r2s_z)
        r2s_y.append(r2s_y)
        r2s_x.append(r2s_x)

    return fwhms_x, fwhms_y, fwhms_z, r2s_x, r2s_y, r2s_z



def plot_fwhm_summary(all_fwhms_x, all_fwhms_y, all_fwhms_z, dataset_names, group_labels):
    """
    Plot mean ± std of FWHM for X/Y/Z axes, grouped by manual categories.

    Parameters:
        all_fwhms_x/y/z: list of lists, shape = [n_datasets][n_beads]
        dataset_names: list of dataset name strings, len = n_datasets
        group_labels: dict mapping group name -> list of dataset name substrings (used to match group)
    """
    groups = {}
    for label, substrings in group_labels.items():
        group_indices = [i for i, name in enumerate(dataset_names)
                         if any(sub in name for sub in substrings)]
        groups[label] = group_indices

    means = {'X': [], 'Y': [], 'Z': []}
    stds = {'X': [], 'Y': [], 'Z': []}
    labels = []

    for group_name, indices in groups.items():
        # Flatten and ignore NaNs
        x_vals = np.concatenate([np.array(all_fwhms_x[i]) for i in indices])
        y_vals = np.concatenate([np.array(all_fwhms_y[i]) for i in indices])
        z_vals = np.concatenate([np.array(all_fwhms_z[i]) for i in indices])

        means['X'].append(np.nanmean(x_vals))
        stds['X'].append(np.nanstd(x_vals))

        means['Y'].append(np.nanmean(y_vals))
        stds['Y'].append(np.nanstd(y_vals))

        means['Z'].append(np.nanmean(z_vals))
        stds['Z'].append(np.nanstd(z_vals))

        labels.append(group_name)

    x = np.arange(len(labels))
    width = 0.25

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(x - width, means['X'], width, yerr=stds['X'], label='X', capsize=6)
    ax.bar(x,         means['Y'], width, yerr=stds['Y'], label='Y', capsize=6)
    ax.bar(x + width, means['Z'], width, yerr=stds['Z'], label='Z', capsize=6)

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel('FWHM (px)')
    ax.set_title('Grouped FWHM Mean ± Std by Axis')
    ax.legend()
    ax.grid(True, axis='y')
    plt.tight_layout()
    plt.show()


def report_fwhm_outliers(all_fwhms_x, all_fwhms_y, all_fwhms_z, threshold=3.0):
    """
    Detect outliers and return lists of (dataset_index, bead_index, value) for each axis.
    """
    def detect_outliers(fwhm_lists):
        all_vals = np.concatenate([np.array(lst) for lst in fwhm_lists])
        mean = np.nanmean(all_vals)
        std = np.nanstd(all_vals)
        lower, upper = mean - threshold * std, mean + threshold * std

        outliers = [
            (i, j, val)
            for i, lst in enumerate(fwhm_lists)
            for j, val in enumerate(lst)
            if not np.isnan(val) and (val < lower or val > upper)
        ]
        return outliers, mean, std

    out_x, mean_x, std_x = detect_outliers(all_fwhms_x)
    out_y, mean_y, std_y = detect_outliers(all_fwhms_y)
    out_z, mean_z, std_z = detect_outliers(all_fwhms_z)

    print(f"X-axis: mean={mean_x:.2f}, std={std_x:.2f}, outliers={len(out_x)}")
    print(f"Y-axis: mean={mean_y:.2f}, std={std_y:.2f}, outliers={len(out_y)}")
    print(f"Z-axis: mean={mean_z:.2f}, std={std_z:.2f}, outliers={len(out_z)}")

    return out_x, out_y, out_z

def remove_joint_outliers(all_fwhms_x, all_fwhms_y, all_fwhms_z, threshold=100.0):
    """
    For each dataset, remove entries (same bead index) where any axis exceeds `threshold` or is invalid.
    """
    new_x, new_y, new_z = [], [], []

    for fx, fy, fz in zip(all_fwhms_x, all_fwhms_y, all_fwhms_z):
        fx_new, fy_new, fz_new = [], [], []

        for x, y, z in zip(fx, fy, fz):
            if all([
                np.isfinite(x), np.isfinite(y), np.isfinite(z),
                x <= threshold, y <= threshold, z <= threshold,
                x > 0, y > 0, z > 0
            ]):
                fx_new.append(x)
                fy_new.append(y)
                fz_new.append(z)

        new_x.append(fx_new)
        new_y.append(fy_new)
        new_z.append(fz_new)

    return new_x, new_y, new_z


def plot_fwhm_outliers(subvolumes_list, all_fwhms_axis, axis_label, outlier_tuples):
    """
    Plot bead volumes for detected FWHM outliers along one axis.
    """
    for i_dataset, i_bead, value in outlier_tuples:
        vol = subvolumes_list[i_dataset][i_bead]
        z, y, x = vol.shape

        if axis_label == 'X':
            profile = vol[z//2, y//2, :]
        elif axis_label == 'Y':
            profile = vol[z//2, :, x//2]
        elif axis_label == 'Z':
            profile = vol[:, y//2, x//2]
        else:
            continue

        x_vals = np.arange(len(profile))
        plt.figure()
        plt.plot(x_vals, profile, 'b.')
        plt.title(f'Outlier: Dataset {i_dataset}, Bead {i_bead}, Axis {axis_label}, FWHM={value:.2f}')
        plt.xlabel('Position (px)')
        plt.ylabel('Intensity')
        plt.grid(True)
        plt.show()

import matplotlib.pyplot as plt

def plot_fwhm_outliers_image(subvolumes_list, outlier_tuples, axis_label):
    """
    Visualize outlier bead volumes as 3 central slices (XY, XZ, YZ).
    axis_label is only used for labeling the plot, not filtering.
    """
    for i_dataset, i_bead, fwhm_val in outlier_tuples:
        vol = subvolumes_list[i_dataset][i_bead]
        z, y, x = vol.shape

        slice_xy = vol[z // 2, :, :]
        slice_xz = vol[:, y // 2, :]
        slice_yz = vol[:, :, x // 2]

        fig, axs = plt.subplots(1, 3, figsize=(10, 3))
        axs[0].imshow(slice_xy, cmap='gray')
        axs[0].set_title('XY plane')
        axs[1].imshow(slice_xz, cmap='gray')
        axs[1].set_title('XZ plane')
        axs[2].imshow(slice_yz, cmap='gray')
        axs[2].set_title('YZ plane')

        for ax in axs:
            ax.axis('off')

        fig.suptitle(f'Dataset {i_dataset}, Bead {i_bead}, {axis_label} FWHM = {fwhm_val:.2f}', fontsize=12)
        plt.tight_layout()
        plt.show()


def plot_fwhm_scatter(all_fwhms_x, all_fwhms_y, all_fwhms_z):
    """
    Scatter plot of all individual FWHM values by dataset and axis.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    n_datasets = len(all_fwhms_x)
    fig, ax = plt.subplots(figsize=(10, 5))

    for i in range(n_datasets):
        # Add horizontal jitter per axis for visibility
        jitter = 0.15
        x_vals = np.full(len(all_fwhms_x[i]), i - jitter)
        y_vals = np.full(len(all_fwhms_y[i]), i)
        z_vals = np.full(len(all_fwhms_z[i]), i + jitter)

        ax.scatter(x_vals, all_fwhms_x[i], color='blue', label='X' if i == 0 else "", alpha=0.6)
        ax.scatter(y_vals, all_fwhms_y[i], color='green', label='Y' if i == 0 else "", alpha=0.6)
        ax.scatter(z_vals, all_fwhms_z[i], color='red', label='Z' if i == 0 else "", alpha=0.6)

    ax.set_xticks(range(n_datasets))
    ax.set_xticklabels([f'DS{i}' for i in range(n_datasets)], rotation=45)
    ax.set_ylabel('FWHM (px)')
    ax.set_title('Raw FWHM Scatter per Dataset and Axis')
    ax.grid(True, axis='y')
    ax.legend()
    plt.tight_layout()
    plt.show()


def plot_fwhm_summary_by_dataset(all_fwhms_x, all_fwhms_y, all_fwhms_z, dataset_names):
    """
    Plot mean ± std of FWHM per axis, one bar per dataset (no grouping).
    """
    means_x = [np.nanmean(fwhms) for fwhms in all_fwhms_x]
    means_y = [np.nanmean(fwhms) for fwhms in all_fwhms_y]
    means_z = [np.nanmean(fwhms) for fwhms in all_fwhms_z]

    stds_x = [np.nanstd(fwhms) for fwhms in all_fwhms_x]
    stds_y = [np.nanstd(fwhms) for fwhms in all_fwhms_y]
    stds_z = [np.nanstd(fwhms) for fwhms in all_fwhms_z]

    x = np.arange(len(dataset_names))
    width = 0.25

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(x - width, means_x, width, yerr=stds_x, label='X', capsize=6)
    ax.bar(x,         means_y, width, yerr=stds_y, label='Y', capsize=6)
    ax.bar(x + width, means_z, width, yerr=stds_z, label='Z', capsize=6)

    for i,xi in enumerate(x):
        ax.scatter(np.repeat(xi - width,len(all_fwhms_x[i])), all_fwhms_x[i])
        ax.scatter(np.repeat(xi,len(all_fwhms_y[i])), all_fwhms_y[i])
        ax.scatter(np.repeat(xi + width,len(all_fwhms_z[i])), all_fwhms_z[i])
              
    ax.set_xticks(x)
    ax.set_xticklabels(dataset_names, rotation=45, ha='right')
    ax.set_ylabel('FWHM (px)')
    ax.set_title('FWHM Mean ± Std per Dataset and Axis')
    ax.legend()
    ax.grid(True, axis='y')
    plt.tight_layout()
    #plt.show()

def estimate_psf_fwhm_from_top_hat(fwhm_meas_px: float,
                                   bead_diam_nm: float,
                                   pixel_size_nm: float) -> float:
    """
    Estimate PSF FWHM (nm) when the object is a top-hat bead.

    Parameters
    ----------
    fwhm_meas_px : float
        FWHM of the fitted Gaussian (in pixels).
    bead_diam_nm : float
        Physical bead diameter (nm), e.g. 5.0.
    pixel_size_nm : float
        Detector sampling in nm/px.

    Returns
    -------
    float
        Estimated PSF FWHM in nm.
    """
    fwhm_meas_px = np.asarray(fwhm_meas_px, dtype=float)

    fwhm_meas_nm = fwhm_meas_px * pixel_size_nm
    # Convert to sigma
    sigma_meas = fwhm_meas_nm / 2.355

    # Variance of a 1-D top-hat of width D
    var_bead = (bead_diam_nm ** 2) / 12.0
    if sigma_meas**2 <= var_bead:
        raise ValueError("Measured profile is narrower than the bead itself.")
    sigma_psf = np.sqrt(sigma_meas**2 - var_bead)
    return 2.355 * sigma_psf   # back to FWHM

# -----------------------------------------------------------------------------
# Skip‑list helpers (based on R² quality)
# -----------------------------------------------------------------------------

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
            if bead_idx >= len(rx_ds):
                continue  # this dataset has fewer beads; ignore
            rx = rx_ds[bead_idx]
            rz = rz_ds[bead_idx]
            if (not np.isfinite(rx)) or (not np.isfinite(rz)) or rx < r2_cutoff or rz < r2_cutoff:
                skip.add(bead_idx)
                break  # no need to check other datasets for this bead

    return sorted(skip)

# -----------------------------------------------------------------------------
# Plotting helpers
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Profile extraction (user‑annotated geometric centre)
# -----------------------------------------------------------------------------

def _extract_center_profiles(vol: np.ndarray, *, linewidth: int = 1) -> dict[str, np.ndarray]:
    """Return X, Y, Z line‑profiles that pass through the **geometric centre**
    of *vol* (assumes the sub‑volume is already centred on the bead by the
    user's annotation).

    If *linewidth* > 1, average symmetrically about the centre point; at the
    edges the window shrinks as needed.
    """
    z_sz, y_sz, x_sz = vol.shape
    zc, yc, xc = z_sz // 2, y_sz // 2, x_sz // 2

    if linewidth <= 1:
        return {
            "Z": vol[:, yc, xc],
            "Y": vol[zc, :, xc],
            "X": vol[zc, yc, :],
        }

    half = linewidth // 2
    slc = lambda c, sz: slice(max(c - half, 0), min(c + half + 1, sz))
    return {
        "Z": vol[:, slc(yc, y_sz), slc(xc, x_sz)].mean(axis=(1, 2)),
        "Y": vol[slc(zc, z_sz), :, slc(xc, x_sz)].mean(axis=(0, 2)),
        "X": vol[slc(zc, z_sz), slc(yc, y_sz), :].mean(axis=(0, 1)),
    }


# -----------------------------------------------------------------------------
# Grid‑plot helper (rows = beads, cols = X/Y/Z)
# -----------------------------------------------------------------------------

def plot_dataset_fit_grid(
    dataset_idx: int,
    subvolumes_list: Sequence[Sequence[np.ndarray]],
    *,
    linewidth: int = 1,
    bead_diam_nm: float = 0.0,
    pixel_size_nm: float | None = None,
    skip_indices_global: Sequence[int] | None = None,
):
    """Visualise Gaussian fits for every bead in *dataset_idx*.

    Parameters
    ----------
    dataset_idx : int
        Index of the dataset to plot (0‑based).
    subvolumes_list : list[list[np.ndarray]]
        Master list: ``[[vols_ds0], [vols_ds1], …]``.
    skip_indices_global : list[int] | None, optional
        Global bead indices to omit (e.g. from :pyfunc:`generate_skip_bead_indices`).
    """
    if dataset_idx < 0 or dataset_idx >= len(subvolumes_list):
        raise IndexError("dataset_idx out of range")

    vols = subvolumes_list[dataset_idx]
    skip_set = set(skip_indices_global or [])
    keep = [i for i in range(len(vols)) if i not in skip_set]
    if not keep:
        print("All beads skipped for this dataset.")
        return

    axes_labels = ("X", "Y", "Z")
    fig, axs = plt.subplots(len(keep), 3, figsize=(9, 2.5 * len(keep)), squeeze=False)

    for row, bead_idx in enumerate(keep):
        vol = vols[bead_idx]
        profiles = _extract_center_profiles(vol, linewidth=linewidth)
        for col, axis in enumerate(axes_labels):
            ax = axs[row, col]
            profile = profiles[axis]
            x = np.arange(len(profile))
            fwhm, popt, r2 = fit_gaussian_and_compute_fwhm(
                profile,
                linewidth=1,
                bead_diam_nm=bead_diam_nm,
                pixel_size_nm=pixel_size_nm,
            )
            ax.plot(x, profile, "b.")
            if np.isfinite(popt[0]):
                x_dense = np.linspace(0, len(profile) - 1, 10 * len(profile))
                ax.plot(x_dense, gaussian(x_dense, *popt), "r-")
            ax.set_title(f"Bead {bead_idx} {axis} | F={fwhm:.1f}px R²={r2:.2f}")
            ax.set_xlabel("px")
            if col == 0:
                ax.set_ylabel("intensity")
            ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.show()