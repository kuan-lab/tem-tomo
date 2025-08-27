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
