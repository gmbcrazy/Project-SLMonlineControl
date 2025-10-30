#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import suite2p
import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt
import scipy.io
import pandas as pd
import re
import os
import glob
from typing import Union, List
import shutil
from natsort import natsorted
from multiprocessing import Pool
from joblib import Parallel, delayed
from pathlib import Path
import gc
gc.collect()
import importlib
import sys
from pathlib import Path
# Add the directory containing your script to the Python path
from suite2p import run_s2p
import FastBin_Suite2p as FBS
import CalBin2p as CB
# Iterate over X and Y values and search for matching files
#PointFile=[None]*len(PointI_values)
import networkx as nx
from typing import List, Sequence, Optional, Union, Any, Dict
import tifffile


def slm_target_match_cell(SLM3D, NeuroPos3D, DistTh):
    """
    Match each SLM3D point to the closest neuron in NeuroPos3D on the same plane (z-value).
    
    Parameters:
        SLM3D (ndarray): N x 3 array of SLM target positions.
        NeuroPos3D (ndarray): M x 3 array of neuron positions.
        DistTh (float): Distance threshold. Matches beyond this are invalid.
    
    Returns:
        SLMtarget (ndarray): Indices of matched neurons, or -1 if too far.
        SLMtargetCellDist (ndarray): Distances to matched neurons.
    """
    num_points = SLM3D.shape[0]
    SLMtarget = np.full(num_points, -1, dtype=int)
    SLMtargetCellDist = np.full(num_points, np.nan)

    for i in range(num_points):
        z_match = NeuroPos3D[:, 2] == SLM3D[i, 2]
        plane_indices = np.where(z_match)[0]

        if plane_indices.size == 0:
            continue

        cell_xy = NeuroPos3D[plane_indices][:, :2]
        mp_xy = SLM3D[i, :2]
        dists = np.linalg.norm(cell_xy - mp_xy, axis=1)

        min_idx = np.argmin(dists)
        min_dist = dists[min_idx]

        if min_dist <= DistTh:
            SLMtarget[i] = plane_indices[min_idx]
        SLMtargetCellDist[i] = min_dist

    return SLMtarget, SLMtargetCellDist



def remove_duplicate_matches(iscell1, iscell2, SLMtargetCellDist):
    """
    Ensure only one SLM point matches to each neuron — keep the closest only.
    """
    match_dict = {}  # session2 neuron idx → (session1 SLM idx, distance)

    for s1, s2 in zip(iscell1, iscell2):  # SLM index, neuron index
        dist = SLMtargetCellDist[s1]      # use SLM index for distance
        if s2 not in match_dict or dist < match_dict[s2][1]:
            match_dict[s2] = (s1, dist)

    filtered_iscell2 = np.array(sorted(match_dict.keys()))
    filtered_iscell1 = np.array([match_dict[s2][0] for s2 in filtered_iscell2])

    return filtered_iscell1, filtered_iscell2




def load_slm_data(ProcessFolder):
    # Load SLMFunGroup.mat
    slm_fun_info = scipy.io.loadmat(os.path.join(ProcessFolder, "SLMFunGroup.mat"))
    pos3d_fun_group = slm_fun_info.get('FinalPos3D')
    fun_score = slm_fun_info.get('FinalFunScore')
    
    # Filter based on non-NaN in FunScore
    valid_fun_rows = ~np.isnan(fun_score[:, 1]) & ~np.isnan(fun_score[:, 2])
    filtered_pos3d = pos3d_fun_group[valid_fun_rows]
    filtered_fun_score = fun_score[valid_fun_rows]

    # Load SLMIncludedIndFromIscell.mat
    slm_test_info = scipy.io.loadmat(os.path.join(ProcessFolder, "SLMIncludedIndFromIscell.mat"))
    slm_test_pos3d = slm_test_info.get('Pos3Dneed')
    slm_test_tbl = slm_fun_info.get('SLMTable')  # Note: Still from SLMFunGroup.mat

    # Filter based on non-NaN in test table
    valid_test_rows = ~np.isnan(slm_test_tbl[:, 1])
    filtered_test_pos3d = slm_test_pos3d[valid_test_rows]
    filtered_test_tbl = slm_test_tbl[valid_test_rows]

    return filtered_pos3d, filtered_fun_score, filtered_test_pos3d, filtered_test_tbl


# run_s2p_by_range.py

# def run_s2p_by_range(root, ops_in, id_range):
#     """
#     root      : full path to the parent 'Data' folder
#     ops_in    : dict with Suite2p parameters you care about (e.g. {'tau': 1.5})
#     id_range  : (low, high) tuple of integers; keeps folders whose final
#                 3-digit suffix is between low and high inclusive.
#     """
#     low, high = id_range

#     # -------- pick sub-folders --------
#     pattern   = re.compile(r'^TSeries-.*-(\d{3})$')
#     subfolders = []
#     for p in Path(root).iterdir():
#         if not p.is_dir():                 # skip files
#             continue
#         m = pattern.match(p.name)
#         if m and low <= int(m.group(1)) <= high:
#             subfolders.append(p.name)

#     if not subfolders:
#         raise ValueError(f'No TSeries folders with suffix in [{low}, {high}] found')

#     # -------- build ops + db --------
#     ops = dict(look_one_level_down=True)   # <- needed when you give subfolders
#     ops.update(ops_in)                     # user-specified overrides

#     db  = dict(data_path=[str(root)],
#                subfolders=subfolders)

#     # -------- run Suite2p --------
#     run_s2p(ops=ops, db=db, server=False)


# --------------------------------------------------------------
#  extract_suite2p.py
#  Author:  (feel free to add your name)
# --------------------------------------------------------------
"""
Extract Suite2p calcium‑imaging results (Python version).

Parameters
----------
confSet : dict‑like
    Required keys
    ├── 'save_path0' : str | pathlib.Path
    ├── 'scan_Z'     : (nPlanes,) array‑like   – microscope focus depth (µm)
    └── 'ETL'        : (nPlanes,) array‑like   – ETL offset (µm)

Returns
-------
Pos3D        : (nCells, 3)  np.ndarray   – XYZ for accepted cells (iscell==1)
Pos3DRaw     : (nRoi , 3)  np.ndarray   – XYZ for *all* ROIs
CaData       : dict                         – combined data, mirroring MATLAB struct
CaDataPlane  : list[dict]                   – per‑plane data dictionaries
stat         : list[dict]                   – stat for accepted cells only
"""

# ----------------------------------------
def _load_npy(fname):
    """Numpy loader with pickling enabled for dict / list objects."""
    return np.load(fname, allow_pickle=True)

# ----------------------------------------

#def _load_npy(fname):
#    return np.load(fname, allow_pickle=True)

def _as_float_array(x):
    """
    Accept list/tuple/ndarray OR a single space/comma‑separated string.
    Returns: np.ndarray(dtype=float)
    """
    if isinstance(x, str):
        # split on space or comma
        parts = [p for p in re.split(r'[,\s]+', x.strip()) if p]
        return np.asarray(parts, dtype=float)
    else:
        return np.asarray(x, dtype=float)


def extract_suite2pNoPlane(s2p_dir):
    s2p_dir  = Path(s2p_dir)
    #etls    = _as_float_array(confSet['ETL'])
    #scan_Zs = _as_float_array(confSet['scan_Z'])
    # -------- 1.  load COMBINED data --------
    CaData = {
        'F'     : _load_npy(s2p_dir / 'F.npy'),
        'Fneu'  : _load_npy(s2p_dir / 'Fneu.npy'),
        'spks'  : _load_npy(s2p_dir / 'spks.npy'),
        'iscell': _load_npy(s2p_dir / 'iscell.npy'),
        'stat'  : _load_npy(s2p_dir / 'stat.npy').tolist(),
        'ops'   : _load_npy(s2p_dir / 'ops.npy' ).item(),
    }

    nPlanes = int(CaData['ops']['nplanes'])
    return CaData

def extract_suite2p(s2p_dir, confSet):
    s2p_dir  = Path(s2p_dir)
    comb_dir = s2p_dir / 'combined'
    etls    = _as_float_array(confSet['ETL'])
    scan_Zs = _as_float_array(confSet['scan_Z'])
    # -------- 1.  load COMBINED data --------
    CaData = {
        'F'     : _load_npy(comb_dir / 'F.npy'),
        'Fneu'  : _load_npy(comb_dir / 'Fneu.npy'),
        'spks'  : _load_npy(comb_dir / 'spks.npy'),
        'iscell': _load_npy(comb_dir / 'iscell.npy'),
        'stat'  : _load_npy(comb_dir / 'stat.npy').tolist(),
        'ops'   : _load_npy(comb_dir / 'ops.npy' ).item(),
    }

    nPlanes = int(CaData['ops']['nplanes'])
    print(nPlanes)
    # -------- 2.  per‑plane loading --------
    plane_dirs = sorted(s2p_dir.glob('plane*'),
                        key=lambda p: int(p.name.replace('plane', '')))
    print(plane_dirs)
    print(nPlanes)
    if len(plane_dirs) != nPlanes:
        raise RuntimeError(f'Plane folder count ({len(plane_dirs)}) '
                           f'does not match ops[nplanes] ({nPlanes})')

    CaDataPlane = []
    for p_dir in plane_dirs:
        plane_dict = {
            'F'     : _load_npy(p_dir / 'F.npy'),
            'Fneu'  : _load_npy(p_dir / 'Fneu.npy'),
            'spks'  : _load_npy(p_dir / 'spks.npy'),
            'iscell': _load_npy(p_dir / 'iscell.npy'),
            'stat'  : _load_npy(p_dir / 'stat.npy').tolist(),
            'ops'   : _load_npy(p_dir / 'ops.npy').item(),
        }
        CaDataPlane.append(plane_dict)

    # store per‑plane mean images in the combined dict
    CaData['PlaneMeanImg']  = np.stack([p['ops']['meanImg']  for p in CaDataPlane])
    CaData['PlaneMeanImgE'] = np.stack([p['ops']['meanImgE'] for p in CaDataPlane])

    # -------- 3.  concatenate stat + plane IDs --------
    stat_all     = []
    plane_id_all = []
    for idx_pl, p in enumerate(CaDataPlane, start=1):
        stat_all.extend(p['stat'])
        plane_id_all.extend([idx_pl] * len(p['stat']))
    plane_id_all = np.asarray(plane_id_all, dtype=int)

    xy_all = np.vstack([s['med'][::-1] for s in stat_all])  # (x, y)

    z_layer = etls + scan_Zs
    z_all   = z_layer[plane_id_all - 1]

    Pos3DRaw = np.column_stack((xy_all, z_all))

    # -------- 4.  accepted cells only --------
    good_mask = CaData['iscell'][:, 0].astype(bool)
    stat      = [s for s, keep in zip(stat_all, good_mask) if keep]

    CaData['statCell']        = stat
    CaData['CellPlaneID']     = plane_id_all[good_mask]
    CaData['CellPlaneIDRaw']  = plane_id_all

    Pos3D = Pos3DRaw[good_mask]

    return Pos3D, Pos3DRaw, CaData, CaDataPlane, stat


def subfolders_and_pairs(folder, KeyWords=None):
    """
    Parameters
    ----------
    folder : str | pathlib.Path
        The parent directory containing sub‑folders whose names start with the
        pattern 'sub<number>-<number>...'.

    KeyWords : str | list[str], optional
        One keyword **or** a list of keywords. Only folder names that contain at
        least one of these keywords (case‑insensitive substring match) will be
        included. If *None*, all matching folders are returned.

    Returns
    -------
    names : list[str]
        Sorted list of sub‑folder names (not full paths).

    pairs : np.ndarray, shape (n, 2), dtype=int
        Two‑column array with the numbers extracted from each name.
    """
    from pathlib import Path
    import re, numpy as np

    folder = Path(folder)
    if not folder.is_dir():
        raise NotADirectoryError(folder)

    # Normalize KeyWords to a list
    if KeyWords is None:
        KeyWords = []
    elif isinstance(KeyWords, str):
        KeyWords = [KeyWords]

    # Pre‑lower keywords for faster comparison
    KeyWords = [kw.lower() for kw in KeyWords]

    # regex:  sub   <firstDigits>  -  <secondDigits>  anythingElse
    pat = re.compile(r"^sub(\d+)-(\d+)", re.IGNORECASE)

    names, nums = [], []
    for p in sorted(folder.iterdir()):
        if p.is_dir():
            m = pat.match(p.name)
            if m:
                if not KeyWords or any(kw in p.name.lower() for kw in KeyWords):
                    names.append(p.name)
                    nums.append([int(m.group(1)), int(m.group(2))])

    pairs = np.asarray(nums, dtype=int).reshape(-1, 2)
    return names, pairs




def backup_iscell_npy(root: Union[str, Path], overwrite: bool = False) -> List[Path]:
    """
    Recursively copy every `iscell.npy` under `root` to `iscellRaw.npy`
    in the same directory.

    Parameters
    ----------
    root : str | Path   # path to your top-level folder
    overwrite : bool    # replace existing backups if True

    Returns
    -------
    List[Path]          # list of created `iscellRaw.npy` paths
    """
    root = Path(root)
    copied: List[Path] = []

    for iscell_path in root.rglob("iscell.npy"):
        raw_path = iscell_path.with_name("iscellRaw.npy")

        if raw_path.exists() and not overwrite:
            continue  # keep the existing backup

        shutil.copy2(iscell_path, raw_path)  # preserves timestamps on most OSes
        copied.append(raw_path)

    return copied



def _backup_local_iscell_npy(folder: Path, *, overwrite: bool = False) -> None:
    """Copy <folder>/iscell.npy → <folder>/iscellRaw.npy (if it exists)."""
    src = folder / "iscell.npy"
    dst = folder / "iscellRaw.npy"
    if src.exists() and (overwrite or not dst.exists()):
        shutil.copy2(src, dst)          # metadata-preserving copy

def update_iscell(
    save_folder: Union[str, Path],
    iscell_update: np.ndarray,
    nplanes: int,
    ROIcolors: [np.ndarray] = None,
):

    """
    Update the 'iscell.npy' flags in a Suite2p combined folder and every plane.

    Parameters
    ----------
    save_folder : str or Path
        Root folder that contains 'combined' plus sub-folders like 'plane0', 'plane1', etc.
    iscell_update : np.ndarray
        1D boolean array or integer indices into the combined units array.
        Values where this is True (or indices listed) will be flagged as real cells.


    Returns
    -------
    tuple : (iscell_combined, stat_combined, unit_plane)
        Updated iscell array, stat list, and plane ID mapping for all units.
    """
    save_folder = Path(save_folder)
    combined_path = save_folder / "combined"

    # Load data from combined folder
    iscell_combined = np.load(combined_path / "iscell.npy")
    stat_combined   = np.load(combined_path / "stat.npy", allow_pickle=True)
    ops_combined    = np.load(combined_path / "ops.npy",  allow_pickle=True)

    # Update iscell flags
    #print(iscell_combined.shape)
    #print(iscell_update.shape)
    iscell_combined[:, 0] = 0
    iscell_combined[iscell_update, 0] = 1
    np.save(combined_path / "iscell.npy", iscell_combined)

    print(
        f"{len(stat_combined)} units total — "
        f"{int(iscell_combined[:, 0].sum())} flagged as cells in combined data"
    )

    # Optional: save ROI colors
    if ROIcolors is not None:
        assert len(ROIcolors) == len(stat_combined), "ROIcolors must match number of ROIs"
        np.save(combined_path / "roi_colors.npy", ROIcolors)
        print("Saved custom ROI color map.")

    # Extract unit-plane mapping
    unit_plane = np.array([u["iplane"] for u in stat_combined], dtype=int)
    #print(ops_combined["nplanes"])


    #nplanes    = int(ops_combined["nplanes"])

    # Update each plane's iscell.npy and stat
    for plane_idx in range(nplanes):
        plane_path = save_folder / f"plane{plane_idx}"
        plane_mask = unit_plane == plane_idx
        np.save(plane_path / "iscell.npy", iscell_combined[plane_mask])

        stat_plane = np.load(plane_path / "stat.npy", allow_pickle=True)
        stat_combined[plane_mask] = stat_plane

    # Save metadata for traceability
    np.savez(
        combined_path / "statUpdate.npz",
        stat      = stat_combined,
        UnitPlane = unit_plane,
        iscell    = iscell_combined,
    )

    return iscell_combined, stat_combined, unit_plane



def ConsistentCell(matched_pairs):
    n_sessions = len(matched_pairs)+1          # total recording sessions
    G = nx.Graph()

    # ➊ add one node for every (session, cell_index) we ever saw
    # ➋ add an edge for every confirmed match between consecutive sessions
    for s, pairs in enumerate(matched_pairs):   # matched_pairs[0] links s=0↔1, etc.
        for c1, c2 in pairs:
            G.add_edge((s,   c1),               # node in session s
                        (s+1, c2))              # node in session s+1
    
    # ➌ connected components whose size == n_sessions have one node per session
    stable_cells = []
    for comp in nx.connected_components(G):
        if len(comp) == n_sessions:
            ordered = sorted(comp, key=lambda x: x[0])   # sort by session index
            cell_chain = [idx for (_, idx) in ordered]   # just the cell IDs
            stable_cells.append(cell_chain)

    print(f"{len(stable_cells)} cells found in ALL {n_sessions} sessions")
    return np.array(stable_cells, dtype=int)
    # stable_cells[k][s]  ⇒ index of that neuron in session s


def ConsistentCellColor(stable_cells, cellNum, seed=None):
    if seed is not None:
        np.random.seed(seed)

    N_sessions = len(cellNum)
    N_cells = len(stable_cells)  # number of stable/matched cells

    # Generate one random color per stable cell (normalized to 0–1)
    common_colors = np.random.permutation(N_cells) / N_cells*0.6+0.2

    # Initialize per-session color arrays
    colors1 = [np.zeros(n) for n in cellNum]       # where consistent values will go

    for color_val, cell_indices in zip(common_colors, stable_cells):
        for session_idx, cell_idx in enumerate(cell_indices):
            colors1[session_idx][cell_idx] = color_val

    return colors1  # list of 1D arrays, one per session





Array2D = np.ndarray  # alias for readability (expects int16 Ly×Lx)

def _load_refs_from_folder(folder: Path) -> List[Array2D]:
    """Collect **refImg** from each `plane*/ops.npy` under *folder*.

    Each plane folder must already contain a valid `ops['refImg']` (no fallback
    to `meanImg`).
    """
    refs: List[Array2D] = []
    plane_dirs = sorted(
        p for p in folder.iterdir() if p.is_dir() and p.name.startswith("plane")
    )
    if not plane_dirs:
        raise FileNotFoundError(f"No plane* sub‑folders found in {folder}")

    for pd in plane_dirs:
        ops_path = pd / "ops.npy"
        if not ops_path.exists():
            raise FileNotFoundError(f"Missing ops.npy in {pd}")
        ops_plane = np.load(ops_path, allow_pickle=True).item()
        if "refImg" not in ops_plane:
            raise KeyError(f"'refImg' missing in {ops_path}; cannot continue")
        refs.append(ops_plane["refImg"].astype(np.int16))
    return refs


def run_s2p_by_range(
    root: Path,
    ops_in: dict,
    id_range: Sequence[int],
    refImg_init: Optional[Union[Array2D, List[Array2D], str]] = None,
):
    """Run Suite2p on TSeries sub‑folders within a numeric suffix range.

    Parameters
    ----------
    root : Path
        Parent directory that contains *TSeries‑xxx‑NNN* sub‑folders.
    ops_in : dict
        User‑supplied overrides for Suite2p ``default_ops()``.
    id_range : tuple(int, int)
        Inclusive numeric range `(low, high)` applied to the **3‑digit suffix**
        in each folder name.  Example folder: ``TSeries‑20250521‑002`` → suffix
        `002`.
    refImg_init : None | 2‑D ndarray | list[2‑D] | str/Path
        * ``None`` – let Suite2p build its own reference(s).
        * 2‑D array – use same ref for all planes.
        * List of 2‑D arrays – one ref per plane.
        * Path – load refs from existing ``plane*/ops.npy`` in that folder.
    """

    low, high = id_range
    pattern = re.compile(r"^TSeries-.*-(\d{3})$")

    # ------------------------------------------------------------------
    # 1. Discover matching sub‑folders
    # ------------------------------------------------------------------
    subfolders: List[str] = []
    for p in Path(root).iterdir():
        if not p.is_dir():
            continue
        m = pattern.match(p.name)
        if m and low <= int(m.group(1)) <= high:
            subfolders.append(p.name)

    if not subfolders:
        raise ValueError(f"No TSeries folders with suffix in [{low}, {high}] found")

    # ------------------------------------------------------------------
    # 2. Assemble ops & db
    # ------------------------------------------------------------------
    ops = dict(look_one_level_down=True)
    ops.update(ops_in)

    db = dict(data_path=[str(root)], subfolders=subfolders)

    # ------------------------------------------------------------------
    # 3. Inject a forced reference (optional)
    # ------------------------------------------------------------------
    if refImg_init is not None:
        
        print('Load forced refImg')
        if isinstance(refImg_init, (str, Path)):
            ops["refImg"] = _load_refs_from_folder(Path(refImg_init))
        elif isinstance(refImg_init, list):
            ops["refImg"] = [np.asarray(r, dtype=np.int16) for r in refImg_init]
        elif isinstance(refImg_init, np.ndarray):
            ops["refImg"] = refImg_init.astype(np.int16)
        else:
            raise TypeError("refImg_init must be None, ndarray, list of ndarrays, or str/Path")

        ops["force_refImg"] = True  # ensure Suite2p uses the supplied reference
        # ---- reporter ------------------------------------------------
        if isinstance(ops["refImg"], list):
            shapes = [r.shape for r in ops["refImg"]]
            print(f"[run_s2p_by_range] Forced refImg list loaded for {len(shapes)} planes: {shapes}")
        else:
            print(f"[run_s2p_by_range] Forced refImg loaded with shape {ops['refImg'].shape}")
    else: 
        print('No forced refImg')
    # ------------------------------------------------------------------
    # 4. Launch Suite2p
    # ------------------------------------------------------------------
    run_s2p(ops=ops, db=db)