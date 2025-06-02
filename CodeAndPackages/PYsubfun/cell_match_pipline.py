"""cell_match_pipeline.py
===================================================
Python re‑implementation of Lu Zhang's MATLAB pipeline for identifying and
tracking the **same neurons** across multiple calcium‑imaging sessions
processed with Suite2p.

⚠️ 2025‑05‑15 — **SimpleITK demons was replaced by a thin wrapper that calls
MATLAB's `imregdemons` via the MATLAB Engine API**.  This guarantees identical
behaviour to the original MATLAB scripts and removes the VectorFloat pixel‑type
mismatch issues seen with SimpleITK.

Install once (from an *elevated* CMD or PowerShell where MATLAB is installed):
```bash
cd "C:\Program Files\MATLAB\R2024a\extern\engines\python"
python -m pip install .
```
Then the rest of this pipeline works entirely in Python while delegating
non‑rigid registration to MATLAB.
"""

from __future__ import annotations

from pathlib import Path
from dataclasses import dataclass
from typing import List, Tuple, Dict, Sequence

import numpy as np
import pandas as pd
from skimage import exposure
from skimage.measure import label, regionprops
from skimage.morphology import dilation, square
from scipy.stats import pearsonr
from scipy.ndimage import map_coordinates
import matlab.engine                # MATLAB Engine API

###############################################################################
# 0. MATLAB engine handle (singleton) -----------------------------------------
###############################################################################

_engine = None

def _get_matlab_engine():
    global _engine
    if _engine is None:
        running = matlab.engine.find_matlab()
        _engine = matlab.engine.connect_matlab(running[0]) if running else matlab.engine.start_matlab()
    return _engine

###############################################################################
# 1. ROI‑neighbourhood extraction --------------------------------------------
###############################################################################

def roi2neighbour(mean_img: np.ndarray,
                  cell_id_map: np.ndarray,
                  cell_ids: Sequence[int],
                  radius: int = 20,
                  intensity_normalise: bool = True
                  ) -> Tuple[np.ndarray, List[np.ndarray]]:
    """Extract a square *radius* neighbourhood around each ROI."""
    pad = radius
    mean_ext = np.pad(mean_img, pad, mode="constant", constant_values=np.nan)
    id_ext   = np.pad(cell_id_map, pad, mode="constant", constant_values=np.nan)

    neighs, bnds = [], []
    for cid in cell_ids:
        ys, xs = np.where(cell_id_map == cid)
        if not len(ys):
            continue
        cy, cx = int(np.median(ys))+pad, int(np.median(xs))+pad
        sl_y = slice(cy-pad, cy+pad+1)
        sl_x = slice(cx-pad, cx+pad+1)
        patch = mean_ext[sl_y, sl_x].copy()
        if intensity_normalise:
            patch = exposure.rescale_intensity(patch, in_range="image")
        neighs.append(patch)

        roi_mask = (id_ext[sl_y, sl_x] == cid)
        roi_mask = dilation(roi_mask, square(3))
        labeled  = label(roi_mask)
        coords   = np.vstack([prop.coords for prop in regionprops(labeled)])
        bnds.append(coords)

    return np.stack(neighs, axis=-1).astype(np.float32), bnds

###############################################################################
# 2. Session container --------------------------------------------------------
###############################################################################

@dataclass
class SessionData:
    mean_img: np.ndarray   # (H, W, nPlanes) float32
    cell_map: np.ndarray   # (H, W, nPlanes) int32 label image
    iscell:   np.ndarray   # N×2 from combined/iscell.npy
    stat:     List[dict]   # stat list

    @property
    def n_good(self):
        return int(self.iscell[:,0].sum())

###############################################################################
# 3. Demons via MATLAB --------------------------------------------------------
###############################################################################

def demons_register(moving: np.ndarray,
                    fixed:  np.ndarray,
                    n_iter: int = 500) -> Tuple[np.ndarray, np.ndarray]:
    """Call MATLAB `imregdemons` and return (disp, moving_registered).
    disp is (H, W, 2) float64, MATLAB order (rowΔ, colΔ)."""
    eng = _get_matlab_engine()

    mv = matlab.double(moving.astype(float).tolist())
    fx = matlab.double(fixed.astype(float).tolist())
    disp_mat, mov_reg = eng.imregdemons(mv, fx, int(n_iter), nargout=2)

    disp = np.asarray(disp_mat, dtype=np.float64)    # (H, W, 2)
    moving_reg = np.asarray(mov_reg, dtype=np.float32)
    return disp, moving_reg

###############################################################################
# 4. Utility: warp label image with displacement ------------------------------
###############################################################################

def warp_labels_nn(label_img: np.ndarray, disp: np.ndarray) -> np.ndarray:
    """Nearest‑neighbour warp of an integer label image using disp (Δy,Δx)."""
    h, w = label_img.shape
    grid_y, grid_x = np.mgrid[0:h, 0:w]
    new_y = grid_y + disp[..., 0]
    new_x = grid_x + disp[..., 1]
    warped = map_coordinates(label_img, [new_y.ravel(), new_x.ravel()], order=0, mode='nearest')
    return warped.reshape(label_img.shape).astype(label_img.dtype)

###############################################################################
# 5. Pairwise matching --------------------------------------------------------
###############################################################################

def match_sessions(sess_fixed: SessionData,
                   sess_moving: SessionData,
                   overlap_th: float = 0.5,
                   radius: int = 20) -> Dict:
    disp, mov_reg = demons_register(sess_moving.mean_img, sess_fixed.mean_img)

    moving_map_reg = warp_labels_nn(sess_moving.cell_map, disp)

    mapping, overlap_pix = [], []
    for mid in np.unique(moving_map_reg):
        if mid == 0:
            continue
        mask_mov = moving_map_reg == mid
        overlap_ids = sess_fixed.cell_map[mask_mov]
        overlap_ids = overlap_ids[overlap_ids > 0]
        if not overlap_ids.size:
            continue
        fid, counts = np.unique(overlap_ids, return_counts=True)
        idx = np.argmax(counts)
        best_fix = int(fid[idx]); best_cnt = counts[idx]
        overlap = best_cnt / (sess_fixed.cell_map == best_fix).sum()
        if overlap >= overlap_th:
            mapping.append((best_fix, int(mid)))
            overlap_pix.append(overlap)

    mapping = np.asarray(mapping, dtype=int)
    overlap_pix = np.asarray(overlap_pix, dtype=float)

    if mapping.size:
        fix_neigh, _ = roi2neighbour(sess_fixed.mean_img, sess_fixed.cell_map, mapping[:,0], radius)
        mov_neigh, _ = roi2neighbour(sess_moving.mean_img, sess_moving.cell_map, mapping[:,1], radius)
        corr = np.array([pearsonr(fix_neigh[...,i].ravel(), mov_neigh[...,i].ravel())[0]
                         for i in range(fix_neigh.shape[-1])])
    else:
        corr = np.empty(0)

    return {
        'mapping': mapping,
        'overlap_ratio': overlap_pix,
        'corr_match': corr,
        'disp': disp,                # H×W×2 float64
        'moving_img_reg': mov_reg,
    }

###############################################################################
# 6. Plane‑wise helpers & loading --------------------------------------------
###############################################################################

def load_session_from_suite2p(s2p_root: str | Path) -> SessionData:
    s2p_root = Path(s2p_root)
    plane_dirs = sorted(s2p_root.glob("plane*"), key=lambda p: int(p.name.replace("plane", "")))
    if not plane_dirs:
        raise FileNotFoundError(f"No plane* folders in {s2p_root}")
    mean, lbl = [], []
    for p in plane_dirs:
        ops = np.load(p/'ops.npy', allow_pickle=True).item()
        mean.append(ops['meanImgE'].astype(np.float32))
        m = np.zeros_like(ops['meanImgE'], dtype=np.int32)
        for i,s in enumerate(np.load(p/'stat.npy', allow_pickle=True), start=1):
            m[s['ypix'], s['xpix']] = i
        lbl.append(m)
    mean_stack = np.stack(mean, axis=-1)
    map_stack  = np.stack(lbl, axis=-1)
    comb = s2p_root/'combined'
    return SessionData(mean_stack, map_stack,
                       np.load(comb/'iscell.npy'),
                       np.load(comb/'stat.npy', allow_pickle=True).tolist())

