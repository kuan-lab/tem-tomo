import webknossos as wk

import tempfile
import numpy as np
import webknossos as wk

import uuid
from pathlib import Path


from scipy.optimize import curve_fit

import matplotlib.pyplot as plt

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
    return A * np.exp(-(x - x0)**2 / (2 * sigma**2)) + offset

def fit_gaussian_and_compute_fwhm(profile, plot_fit=False, title=""):
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
            plt.plot(x, gaussian(x, *popt), 'r-', label=f'Gaussian Fit\nFWHM={fwhm:.2f} px')
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



def compute_fwhms_from_subvolumes(subvolumes):
    fwhms_x, fwhms_y, fwhms_z = [], [], []

    for vol in subvolumes:
        z, y, x = vol.shape

        profile_z = vol[:, y // 2, x // 2]
        profile_y = vol[z // 2, :, x // 2]
        profile_x = vol[z // 2, y // 2, :]

        fwhm_z, _ = fit_gaussian_and_compute_fwhm(profile_z)
        fwhm_y, _ = fit_gaussian_and_compute_fwhm(profile_y)
        fwhm_x, _ = fit_gaussian_and_compute_fwhm(profile_x)

        fwhms_z.append(fwhm_z)
        fwhms_y.append(fwhm_y)
        fwhms_x.append(fwhm_x)

    return fwhms_x, fwhms_y, fwhms_z



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

    ax.set_xticks(x)
    ax.set_xticklabels(dataset_names, rotation=45, ha='right')
    ax.set_ylabel('FWHM (px)')
    ax.set_title('FWHM Mean ± Std per Dataset and Axis')
    ax.legend()
    ax.grid(True, axis='y')
    plt.tight_layout()
    plt.show()
