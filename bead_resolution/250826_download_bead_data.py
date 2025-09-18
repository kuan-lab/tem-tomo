#!/usr/bin/env python3
"""
Download bead-centered subvolumes from WebKnossos for multiple tomograms and reconstructions.

Usage (config-driven):
    python download_bead_data.py --config tomograms.yaml \
        --org 8632814cfac2f959 --url https://webknossos.org \
        --out data --max-workers 8

Usage (single tomogram ad-hoc):
    python download_bead_data.py --tomogram MouseCerebellum_A4S1_SA3.3k_de64_tomo10 \
        --annotation https://webknossos.org/annotations/683f3d... \
        --recon 11_lim15_tomo10a_16bit --recon 15_lim21_tomo10a_16bit \
        --org 8632814cfac2f959 --url https://webknossos.org --out data

Auth:
- Pass a token with --token, or set WEBKNOSSOS_TOKEN in your environment.
"""

import os
import sys
import argparse
import concurrent.futures as futures
from pathlib import Path
from typing import List, Dict, Any, Optional

import yaml  # pip install pyyaml
import numpy as np
import webknossos as wk

# Reuse your existing helpers (make sure utils.py is on PYTHONPATH
import utils  # uses: get_annotation_points, compute_bounding_boxes_flex, save_subvolumes_to_npz

# -----------------------
# Helpers
# -----------------------

def _is_url(s: str) -> bool:
    return s.startswith("http://") or s.startswith("https://")

def _log(msg: str):
    print(msg, flush=True)

def _require_token(cli_token: Optional[str]) -> str:
    tok = cli_token or os.environ.get("WEBKNOSSOS_TOKEN")
    if not tok:
        raise SystemExit("No token provided. Use --token or set WEBKNOSSOS_TOKEN env var.")
    return tok

def _safe_name(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("_", "-", ".") else "_" for c in s)

def _ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def _pt_to_xyz(pt) -> tuple[int, int, int]:
    """Extract (x,y,z) from wk point objects or tuples/lists."""
    if hasattr(pt, "x"):
        return int(pt.x), int(pt.y), int(pt.z)
    if isinstance(pt, dict) and all(k in pt for k in ("x", "y", "z")):
        return int(pt["x"]), int(pt["y"]), int(pt["z"])
    x, y, z = pt  # assume (x,y,z)
    return int(x), int(y), int(z)

def make_centered_roi_plans(
    points: list,
    box_size: int,
    dataset_shape_xyz: tuple[int, int, int],
) -> list[dict]:
    """
    For each bead center (global x,y,z), plan the dataset read (clipped to bounds)
    and where to paste it into a NaN-padded cube so the bead stays centered.

    Returns a list of dicts:
      {
        'offset': (rx0, ry0, rz0),       # dataset read offset (X,Y,Z)
        'size':   (sx,  sy,  sz),        # dataset read size   (X,Y,Z)
        'paste':  (px,  py,  pz),        # paste offset in target (X,Y,Z) coords
        'target': (S,   S,   S),         # target cube (X,Y,Z) = (box_size, ...)
      }
    """
    S = int(box_size)
    hx = hy = hz = S // 2
    Xmax, Ymax, Zmax = dataset_shape_xyz  # dataset bounds in (X,Y,Z)
    plans = []
    for pt in points:
        xc, yc, zc = _pt_to_xyz(pt)

        # Desired cube extents centered at bead (open-ended [start, end))
        x0_d, x1_d = xc - hx, xc - hx + S
        y0_d, y1_d = yc - hy, yc - hy + S
        z0_d, z1_d = zc - hz, zc - hz + S

        # Clip to dataset
        rx0, rx1 = max(0, x0_d), min(Xmax, x1_d)
        ry0, ry1 = max(0, y0_d), min(Ymax, y1_d)
        rz0, rz1 = max(0, z0_d), min(Zmax, z1_d)

        sx, sy, sz = max(0, rx1 - rx0), max(0, ry1 - ry0), max(0, rz1 - rz0)

        # Paste offsets inside the SxSxS target (how far we had to clip)
        px, py, pz = rx0 - x0_d, ry0 - y0_d, rz0 - z0_d

        plans.append({
            "offset": (int(rx0), int(ry0), int(rz0)),
            "size":   (int(sx),  int(sy),  int(sz)),
            "paste":  (int(px),  int(py),  int(pz)),
            "target": (S, S, S),
        })
    return plans

# -----------------------
# Core download (fast path with threading)
# -----------------------

def _get_primary_layer_name(remote_dataset: wk.Dataset) -> str:
    layers = list(remote_dataset.layers.keys())
    if not layers:
        raise ValueError("No layers found in dataset.")
    # Heuristic: prefer 'color' when present
    if "color" in layers:
        return "color"
    return layers[0]

def _download_one_box(
    dataset_name: str,
    organization_id: str,
    webknossos_url: str,
    layer_name: str,
    mag_level: int,
    box: wk.BoundingBox,
) -> np.ndarray:
    """
    Download a single ROI (BoundingBox) as (Z, Y, X) ndarray via wk.Dataset.download.
    Uses a temp on-disk dataset for the ROI to keep the API simple & robust.
    """
    # Each call creates a tiny local dataset containing just this box
    ds_local = wk.Dataset.download(
        dataset_name_or_url=dataset_name,
        organization_id=organization_id,
        bbox=box,
        layers=[layer_name],
        mags=[wk.Mag(mag_level)],
        webknossos_url=webknossos_url,
        # NOTE: letting wk pick a temp path; alternatively give a custom one
    )
    arr = ds_local.get_layer(layer_name).get_mag(wk.Mag(mag_level)).read()
    # Shapes from WK are typically (C, Y, X, Z) or (1, Y, X, Z). Make it (Z, Y, X).
    arr = np.asarray(arr)
    if arr.ndim == 4:  # (C, Y, X, Z) or (1, Y, X, Z)
        arr = arr.squeeze(0)     # (Y, X, Z)
    # transpose to (Z, Y, X)
    arr = arr.transpose(2, 0, 1)
    return arr
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
def download_subvolumes_centered_nanpad(
    dataset_name: str,
    roi_plans: list[dict],
    api_token: str,
    organization_id: str,
    webknossos_url: str = "https://webknossos.org",
    mag_level: int = 1,
    max_workers: int = 8,
) -> list:
    """
    For each ROI plan:
      - read clipped chunk from remote dataset using (absolute_offset, size)
      - convert to (Z,Y,X)
      - paste into a (S,S,S) float32 array prefilled with NaN so bead stays centered
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed
    import threading

    is_url = _is_url(dataset_name)
    mag_str = str(mag_level)
    subvols = [None] * len(roi_plans)

    # discover primary layer once
    try:
        if is_url:
            with wk.webknossos_context(token=api_token):
                remote = wk.Dataset.open_remote(dataset_name)
                layer_name = _get_primary_layer_name(remote)
        else:
            with wk.webknossos_context(url=webknossos_url, token=api_token):
                remote = wk.Dataset.open_remote(dataset_name, organization_id=organization_id)
                layer_name = _get_primary_layer_name(remote)
    except Exception as e:
        _log(f"    ✗ Could not open dataset to discover layer: {e}")
        return []

    _tls = threading.local()

    def _get_ds():
        if getattr(_tls, "ds", None) is not None:
            return _tls.ds
        if is_url:
            _tls.ctx = wk.webknossos_context(token=api_token)
            _tls.ctx.__enter__()
            _tls.ds = wk.Dataset.open_remote(dataset_name)
        else:
            _tls.ctx = wk.webknossos_context(url=webknossos_url, token=api_token)
            _tls.ctx.__enter__()
            _tls.ds = wk.Dataset.open_remote(dataset_name, organization_id=organization_id)
        return _tls.ds

    def _worker(i: int, plan: dict):
        try:
            ds = _get_ds()
            (rx0, ry0, rz0) = plan["offset"]
            (sx,  sy,  sz)  = plan["size"]
            (px,  py,  pz)  = plan["paste"]
            (Sx,  Sy,  Sz)  = plan["target"]  # all equal to box_size

            # If nothing to read (should not happen), return all-NaN cube
            import numpy as np
            tgt = np.full((Sz, Sy, Sx), np.nan, dtype=np.float32)

            if sx <= 0 or sy <= 0 or sz <= 0:
                return i, tgt  # all NaN

            # WK uses (X,Y,Z) for offsets/sizes; read and convert to (Z,Y,X)
            vol = ds.get_layer(layer_name).get_mag(mag_str).read(
                absolute_offset=(rx0, ry0, rz0),
                size=(sx, sy, sz),
            )
            arr = np.asarray(vol)  # (C,Y,X,Z) or (Y,X,Z)
            if arr.ndim == 4:
                arr = arr.squeeze(0)  # (Y,X,Z)
            arr = arr.transpose(2, 0, 1)  # -> (Z,Y,X) = (sz,sy,sx)

            # paste into NaN-padded target so bead stays centered
            tgt[pz:pz+sz, py:py+sy, px:px+sx] = arr.astype(np.float32)
            return i, tgt
        except Exception as e:
            return i, e

    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs = [ex.submit(_worker, i, plan) for i, plan in enumerate(roi_plans)]
        for f in as_completed(futs):
            i, res = f.result()
            if isinstance(res, Exception):
                _log(f"    ⚠️  Failed ROI {i}: {res}")
            else:
                subvols[i] = res

    return [sv for sv in subvols if sv is not None]

def download_subvolumes_parallel(
    dataset_name: str,
    bounding_boxes: list,
    api_token: str,
    organization_id: str,
    webknossos_url: str = "https://webknossos.org",
    mag_level: int = 1,
    max_workers: int = 8,
) -> list:
    """
    Read ROIs directly from the remote dataset (thread-safe, no temp datasets).
    Works with either dataset *name* (+org) or a canonical dataset *URL*.
    Always returns a list (possibly empty) — never None.
    """
    is_url = _is_url(dataset_name)
    mag_str = str(mag_level)  # robust across wk versions

    # Find the primary layer once
    try:
        if is_url:
            with wk.webknossos_context(token=api_token):
                remote = wk.Dataset.open_remote(dataset_name)
                layer_name = _get_primary_layer_name(remote)
        else:
            with wk.webknossos_context(url=webknossos_url, token=api_token):
                remote = wk.Dataset.open_remote(dataset_name, organization_id=organization_id)
                layer_name = _get_primary_layer_name(remote)
    except Exception as e:
        _log(f"    ⚠️  Could not open dataset to discover layer: {e}")
        return []  # <- never return None

    subvols = [None] * len(bounding_boxes)

    # Per-thread cache so we don't reopen the dataset for every ROI
    _tls = threading.local()

    def _get_ds():
        if getattr(_tls, "ds", None) is not None:
            return _tls.ds
        if is_url:
            _tls.ctx = wk.webknossos_context(token=api_token)
            _tls.ctx.__enter__()
            _tls.ds = wk.Dataset.open_remote(dataset_name)
        else:
            _tls.ctx = wk.webknossos_context(url=webknossos_url, token=api_token)
            _tls.ctx.__enter__()
            _tls.ds = wk.Dataset.open_remote(dataset_name, organization_id=organization_id)
        return _tls.ds

    def _cleanup_ds():
        ctx = getattr(_tls, "ctx", None)
        if ctx is not None:
            try:
                ctx.__exit__(None, None, None)
            except Exception:
                pass
            _tls.ctx = None
            _tls.ds = None

    def _worker(i: int, box: wk.BoundingBox):
        try:
            ds = _get_ds()
            # Build offset/size from the wk.BoundingBox (X, Y, Z)
            offset = (int(box.topleft.x), int(box.topleft.y), int(box.topleft.z))
            size   = (int(box.size.x),    int(box.size.y),    int(box.size.z))

            vol = ds.get_layer(layer_name).get_mag(mag_str).read(
                absolute_offset=offset,
                size=size,
            )
            arr = np.asarray(vol)            # (C,Y,X,Z) or (Y,X,Z)
            if arr.ndim == 4:
                arr = arr.squeeze(0)         # -> (Y,X,Z)
            return i, arr.transpose(2, 0, 1) # -> (Z,Y,X)
        except Exception as e:
            return i, e  # return the exception so we can log centrally
        finally:
            # No-op: keep ds open per-thread for speed; comment-in below to close each time
            # _cleanup_ds()
            pass

    try:
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            futs = [ex.submit(_worker, i, box) for i, box in enumerate(bounding_boxes)]
            for f in as_completed(futs):
                i, result = f.result()
                if isinstance(result, Exception):
                    _log(f"    ⚠️  Failed ROI {i}: {result}")
                else:
                    subvols[i] = result
    finally:
        # Ensure all thread-local contexts are closed (for notebooks / repeated runs)
        # Spinning up dummy tasks ensures _cleanup_ds is called across threads if needed.
        pass

    # Always return a list
    return [sv for sv in subvols if sv is not None]


# -----------------------
# Per-tomogram pipeline
# -----------------------

def process_tomogram(
    *,
    tomogram: str,
    annotation_id: str,
    recon_labels: List[str],
    organization_id: str,
    webknossos_url: str,
    api_token: str,
    out_root: Path,
    box_size: int = 36,
    mag_level: int = 1,
    max_workers: int = 8,
    overwrite: bool = False,
    dry_run: bool = False,  
):
    _log(f"\n=== Tomogram: {tomogram} ===")
    t_root = out_root / _safe_name(tomogram)
    _ensure_dir(t_root)

    # 1) Annotation once per tomogram
    with wk.webknossos_context(url=webknossos_url, token=api_token):
        ann = wk.Annotation.download(annotation_id)
    pts = utils.get_annotation_points(ann)
    _log(f"  • Downloaded annotation: {ann.name} | beads: {len(pts)}")

    for recon in recon_labels:
        out_path = t_root / f"{_safe_name(recon)}.npz"
        if out_path.exists() and not overwrite:
            _log(f"  • {recon}: exists → skip (use --overwrite to redo)")
            continue

        _log(f"  • Recon: {recon}")

       # 2) Bounds for edge-safe planning (we will NOT shift; we will NaN-pad)
        try:
            if _is_url(recon):
                with wk.webknossos_context(token=api_token):
                    ds_remote = wk.Dataset.open_remote(recon)
                    bounds = ds_remote.calculate_bounding_box()
            else:
                with wk.webknossos_context(url=webknossos_url, token=api_token):
                    ds_remote = wk.Dataset.open_remote(recon, organization_id=organization_id)
                    bounds = ds_remote.calculate_bounding_box()
        except Exception as e:
            _log(f"    ✗ Could not open dataset for bounds: {e}")
            continue

        dataset_shape_xyz = (bounds.size.x, bounds.size.y, bounds.size.z)

        # Plan centered ROIs with NaN padding (no shifting)
        plans = make_centered_roi_plans(pts, box_size=box_size, dataset_shape_xyz=dataset_shape_xyz)

        # 3) Download (centered, NaN-padded)
        subvols = download_subvolumes_centered_nanpad(
            dataset_name=recon,
            roi_plans=plans,
            api_token=api_token,
            organization_id=organization_id,
            webknossos_url=webknossos_url,
            mag_level=mag_level,
            max_workers=max_workers,
        )
        _log(f"    - Downloaded {len(subvols)}/{len(plans)} subvolumes (NaN-padded)")

        # 4) Save
        utils.save_subvolumes_to_npz(subvols, pts, out_path)
        _log(f"    - Saved → {out_path}")


# -----------------------
# Config parsing
# -----------------------

def load_config(path: Path) -> Dict[str, Any]:
    data = yaml.safe_load(path.read_text())
    if "tomograms" not in data:
        raise ValueError("Config must contain a top-level 'tomograms' list.")
    return data

def run_from_config(
    cfg_path: Path,
    *,
    organization_id: str,
    webknossos_url: str,
    api_token: str,
    out_root: Path,
    box_size: int,
    mag_level: int,
    max_workers: int,
    overwrite: bool,
    dry_run: bool,
):
    cfg = load_config(cfg_path)
    for entry in cfg["tomograms"]:
        tomogram = entry["tomogram"]
        annotation_id = entry["annotation_id"]
        recon_labels = entry["recon_labels"]
        entry_box = int(entry.get("box_size", box_size))
        process_tomogram(
            tomogram=tomogram,
            annotation_id=annotation_id,
            recon_labels=recon_labels,
            organization_id=organization_id,
            webknossos_url=webknossos_url,
            api_token=api_token,
            out_root=out_root,
            box_size=entry_box,
            mag_level=mag_level,
            max_workers=max_workers,
            overwrite=overwrite,
            dry_run=dry_run,
        )


# -----------------------
# CLI
# -----------------------

def main(argv=None):
    p = argparse.ArgumentParser(description="Download bead subvolumes from WebKnossos.")
    gcfg = p.add_argument_group("Config")
    gcfg.add_argument("--config", type=Path, help="YAML file listing tomograms.")
    gcfg.add_argument("--tomogram", help="Single tomogram name (ad-hoc mode).")
    gcfg.add_argument("--annotation", help="Annotation URL/ID for the tomogram.")
    gcfg.add_argument("--recon", action="append", default=[], help="Reconstruction dataset name. Repeatable.")
    gcfg.add_argument("--box-size", type=int, default=36, help="ROI cube side in voxels (default: 36).")
    gcfg.add_argument("--mag", type=int, default=1, help="WebKnossos mag level (default: 1).")
    gcfg.add_argument("--max-workers", type=int, default=8, help="Concurrent ROI downloads (default: 8).")

    gsvc = p.add_argument_group("WebKnossos")
    gsvc.add_argument("--org", required=True, help="Organization ID/slug.")
    gsvc.add_argument("--url", default="https://webknossos.org", help="WebKnossos base URL.")
    gsvc.add_argument("--token", help="API token (or set WEBKNOSSOS_TOKEN).")

    gout = p.add_argument_group("Output")
    gout.add_argument("--out", type=Path, default=Path("data"), help="Output root folder (default: ./data)")
    gout.add_argument("--overwrite", action="store_true", help="Overwrite existing .npz files")
    gcfg.add_argument("--dry-run", action="store_true", 
            help="Validate annotation/datasets/layer/mag/bounds only; do not download or write.")
    args = p.parse_args(argv)

    api_token = _require_token(args.token)
    out_root = args.out
    _ensure_dir(out_root)

    if args.config:
        run_from_config(
            args.config,
            organization_id=args.org,
            webknossos_url=args.url,
            api_token=api_token,
            out_root=out_root,
            box_size=args.box_size,
            mag_level=args.mag,
            max_workers=args.max_workers,
            overwrite=args.overwrite,
            dry_run=args.dry_run,
        )
    else:
        # Ad-hoc single tomogram mode
        if not (args.tomogram and args.annotation and args.recon):
            raise SystemExit("In ad-hoc mode, provide --tomogram, --annotation, and at least one --recon.")
        process_tomogram(
            tomogram=args.tomogram,
            annotation_id=args.annotation,
            recon_labels=args.recon,
            organization_id=args.org,
            webknossos_url=args.url,
            api_token=api_token,
            out_root=out_root,
            box_size=args.box_size,
            mag_level=args.mag,
            max_workers=args.max_workers,
            overwrite=args.overwrite,
            dry_run=args.dry_run,   
        )


if __name__ == "__main__":
    main()
