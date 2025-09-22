#!/usr/bin/env python
# coding: utf-8

# In[ ]:
#import suite2p
import numpy as np
#from scipy.io import loadmat
import matplotlib.pyplot as plt
import scipy.io
import pandas as pd
import re
import os
import glob
from natsort import natsorted
from multiprocessing import Pool
from joblib import Parallel, delayed
from pathlib import Path
import gc
gc.collect()
import importlib
import sys
# Add the directory containing your script to the Python path

from scipy.signal import welch
from scipy.stats import zscore
from suite2p.extraction import dcnv as s2p_dcnv
import math
from typing import Iterable, Optional, Tuple, List, Any,Sequence
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
#import FastBin_Suite2p as FBS
#import CalBin2p as CB
#import SLM_Offline as slmO


sys.path.append(r'C:\Users\zhangl33\Projects\Project-SLMonlineControl\CodeAndPackages\PYsubfun')

def cell_to_str(x):
    # unwrap MATLAB cell entries that come in as 0-d/1-d numpy arrays
    if isinstance(x, np.ndarray):
        # most loadmat cell entries are 0-d or 1-d arrays with a single string
        try:
            return str(x.item())
        except Exception:
            # fallback: join if someone stored more than one string
            return "\\".join(map(str, x.ravel().tolist()))
    return str(x)


def is_valid_windows_path(path: str) -> (bool, str):
    # Windows forbidden characters
    forbidden = r'[<>:"/\\|?*]'
    # Reserved device names (case-insensitive)
    reserved = {
        "CON","PRN","AUX","NUL",
        *(f"COM{i}" for i in range(1,10)),
        *(f"LPT{i}" for i in range(1,10))
    }

    # Check total path length
    if len(path) >= 260:  # old Win32 limit
        return False, f"path too long ({len(path)} chars)"

    # Check each component
    for part in path.split("\\"):
        if not part:
            continue
        if re.search(forbidden, part):
            return False, f"forbidden character in '{part}'"
        if part.upper() in reserved:
            return False, f"reserved name '{part}'"

    return True, "ok"




def flatten_csvnames(csv_field: Any) -> List[str]:
    """Return a flat list of CSV paths from a MATLAB-style field.

    Parameters
    ----------
    csv_field : Any
        The value corresponding to the MATLAB struct field (e.g., FOV['CSVname']).
        This may be a nested object/ndarray of dtype=object with 1xN cell strings.

    Returns
    -------
    List[str]
        A Python list of file paths (strings).
    """
    arr = np.array(csv_field, dtype=object)
    flat = arr.ravel()

    paths: List[str] = []
    for item in flat:
        # Unwrap MATLAB cell arrays that often appear as length-1 lists/arrays
        if isinstance(item, (list, tuple, np.ndarray)) and len(item) == 1:
            item = item[0]
        paths.append(str(item))

    return paths





def extract_metrics_from_stat(statRaw, iscell, recording_name):
    """Extract selected metrics from Suite2p statRaw list into a DataFrame."""
    df = pd.DataFrame({
        "recording": [recording_name] * len(statRaw),
        "iscell": iscell.astype(int),
        "footprint": [s.get('footprint') for s in statRaw],
        "compact": [s.get('compact') for s in statRaw],
        "aspect_ratio": [s.get('aspect_ratio') for s in statRaw],
        "npix": [s.get('npix') for s in statRaw],
        "skew": [s.get('skew') for s in statRaw],
        "mrs": [s.get('mrs') for s in statRaw],
        "mrs0": [s.get('mrs0') for s in statRaw],
        "solidity": [s.get('solidity') for s in statRaw],
        "skew": [s.get('skew') for s in statRaw],
        "npix_soma": [s.get('npix_soma') for s in statRaw],    
        "skew": [s.get('skew') for s in statRaw],
        "npix": [s.get('npix') for s in statRaw],
        "npix_norm": [s.get('npix_norm') for s in statRaw],
        "npix_norm_no_crop": [s.get('npix_norm_no_crop') for s in statRaw],
        "radius": [s.get('radius') for s in statRaw],
        "std": [s.get('std') for s in statRaw],
        # add other metrics as needed
    })
    return df


# -------- PSD noise (same as before) --------
def _ensure_2d(x):
    x = np.asarray(x)
    if x.ndim == 1:
        x = x[None, :]
    return x

def psd_noise_sigma(dff, fs, fmin=None, fmax=None, nperseg=None):
    X = _ensure_2d(dff)
    nyq = fs / 2.0

    lo = 0.7 * nyq if fmin is None else float(fmin)
    hi = 0.95 * nyq if fmax is None else float(fmax)
    lo = max(0.0, lo)
    hi = min(hi, nyq * 0.98)

    # compute once on the first trace to know frequency grid length
    f0, _ = welch(X[0], fs=fs, nperseg=nperseg, detrend='constant',
                  return_onesided=True, scaling='density')
    n_traces, n_freq = X.shape[0], f0.size

    sigmas  = np.zeros(n_traces, dtype=float)
    sigmasZ = np.zeros(n_traces, dtype=float)

    # store PSDs (density) if you want them
    psdZ_density = np.zeros((n_traces, n_freq), dtype=float)

    # store percent-normalized PSDs (sum to 100%)
    psd_pct   = np.zeros((n_traces, n_freq), dtype=float)
    psdZ_pct  = np.zeros((n_traces, n_freq), dtype=float)

    # frequency bin width (Welch uses uniform spacing)
    df = f0[1] - f0[0]

    for i in range(n_traces):
        f, Pxx = welch(X[i], fs=fs, nperseg=nperseg, detrend='constant',
                       return_onesided=True, scaling='density')
        _, PxxZ = welch(zscore(X[i]), fs=fs, nperseg=nperseg, detrend='constant',
                        return_onesided=True, scaling='density')

        # keep the z-PSD (density) if needed by caller
        psdZ_density[i, :] = PxxZ

        # ---- percent-normalization (convert density -> power per bin first) ----
        power_bins   = Pxx  * df
        power_bins_Z = PxxZ * df
        total_power   = power_bins.sum()
        total_power_Z = power_bins_Z.sum()
        # avoid divide-by-zero if a trace is all zeros
        if total_power   > 0: psd_pct[i, :]  = power_bins   / total_power
        if total_power_Z > 0: psdZ_pct[i, :] = power_bins_Z / total_power_Z
        # ------------------------------------------------------------------------

        # your band mask
        m = (f >= lo) & (f <= hi)
        if not np.any(m):
            k = max(1, f.size // 10)   # fallback: top 10% freqs
            m = np.zeros_like(f, dtype=bool); m[-k:] = True

        # your original flat-PSD approximation for variance
        P0  = np.median(Pxx[m])
        P0Z = np.median(PxxZ[m])
        var  = P0  * nyq
        varz = P0Z * nyq
        sigmas[i]  = np.sqrt(max(0.0, var))
        sigmasZ[i] = np.sqrt(max(0.0, varz))

    # returns: sigmas, sigmas on zscored signal, band info, and percent-normalized z-PSD
    return sigmas, sigmasZ, (lo, hi, nyq), psdZ_pct, f


# -------- Event parsing + robust amplitude --------
def _event_groups(mask_true):
    """yield contiguous index groups where mask_true is True"""
    idx = np.flatnonzero(mask_true)
    if idx.size == 0:
        return []
    cuts = np.where(np.diff(idx) > 1)[0] + 1
    return np.split(idx, cuts)

def robust_event_amplitude(spks, sigmas,  eps=1e-8, amp_method="max",
                           min_event_snr=2.5,  # drop events with amp < 2.5 * sigma
                           amp_percentile=0.8, # take 80th-percentile across retained events
                           min_event_frames=1):
    """
    spks: (N,T) Suite2p deconvolved signal
    sigmas: (N,) noise SD from PSD of dF
    Returns A_evt per ROI using noise-thresholded high-quantile peak (or sum) amplitudes.
    """
    S = np.asarray(spks, float)
    N, T = S.shape
    A_evt = np.zeros(N, float)

    for i in range(N):
        s = S[i]
        groups = _event_groups(s > eps)
        if not groups:
            A_evt[i] = 0.0
            continue

        # per-event amplitude
        if amp_method == "max":
            amps = [s[g].max() for g in groups if len(g) >= min_event_frames]
        elif amp_method == "sum":
            amps = [s[g].sum() for g in groups if len(g) >= min_event_frames]
        elif amp_method == "energy":  # optional: L2 energy; penalizes many small frames
            amps = [np.sqrt((s[g]**2).sum()) for g in groups if len(g) >= min_event_frames]
        else:
            raise ValueError("amp_method must be 'max', 'sum', or 'energy'")

        # noise-aware censoring: keep only “real” events
        thr = float(min_event_snr) * sigmas[i]
        amps = [a for a in amps if a >= thr]
        #print(amps)
        #print(thr)


        if len(amps) == 0:
            A_evt[i] = 0.0
        else:
            q = float(amp_percentile)
            q = min(max(q, 0.0), 1.0)
            A_evt[i] = np.quantile(amps, q)

    return A_evt

# -------- Main wrapper (uses Suite2p exactly like run_s2p.py) --------
def compute_snr_with_suite2p_preprocess(
    CaData,
    nperseg=512,
    fband=None,                 # e.g., (2.5, 3.3) for fs≈6.9; None => auto near Nyquist
    amp_method="sum",           # 'max' resists many-small-event cells; 'sum' available too
    min_event_snr=10,          # events must exceed this many sigmas to count
    amp_percentile=0.8,         # robust high-quantile across events
    eps=1e-8,
    recompute_spks=False
):
    F    = np.asarray(CaData['F'], float)
    Fneu = np.asarray(CaData['Fneu'], float)
    spks = np.asarray(CaData['spks'], float)
    ops  = CaData['ops']

    fs        = float(ops.get('fs', 10.0))
    print(fs)
    neucoeff  = float(ops.get('neucoeff', 0.7))
    print(neucoeff)
    baseline  = str(ops.get('baseline', 'maximin'))
    win_base  = float(ops.get('win_baseline', 60.0))
    print(win_base)
    sig_base  = float(ops.get('sig_baseline', 10.0))
    print(sig_base)
    prct_base = float(ops.get('prctile_baseline', 8.0))

    # Suite2p-like preprocess (exact calls)
    dF = F.copy() - neucoeff * Fneu
    dF = s2p_dcnv.preprocess(F=dF, baseline=baseline, win_baseline=win_base,
                             sig_baseline=sig_base, fs=fs, prctile_baseline=prct_base)

    if recompute_spks:
        spks = s2p_dcnv.oasis(F=dF, batch_size=int(ops.get('batch_size', 500)),
                              tau=float(ops.get('tau', 1.0)), fs=fs)

    # PSD noise sigma
    if fband is None:
        sigma_psd, sigmasZ, band, psdPerc, fre = psd_noise_sigma(dF, fs, nperseg=nperseg)
    else:
        sigma_psd, sigmasZ, band, psdPerc, fre = psd_noise_sigma(dF, fs, fmin=fband[0], fmax=fband[1], nperseg=nperseg)

    # Noise-thresholded, high-quantile event amplitude
    A_evt = robust_event_amplitude(spks, sigma_psd, eps=eps, amp_method=amp_method,
                                   min_event_snr=min_event_snr, amp_percentile=amp_percentile)

    snr = A_evt / np.maximum(sigma_psd, 1e-12)
    return snr, sigma_psd, A_evt, psdPerc, fre

def median_spks_per_roi_eventwise(spks, eps=1e-8, method='sum'):
    """
    Group contiguous positive frames as events; per-event amplitude = sum or max over those frames.
    Returns median amplitude across events, per ROI.
    """
    S = np.asarray(spks, dtype=float)
    n, T = S.shape
    out = np.zeros(n, dtype=float)
    for i in range(n):
        s = S[i]
        act = s > eps
        if not np.any(act):
            out[i] = 0.0
            continue
        idx = np.flatnonzero(act)
        # split events where there's a gap >1 frame
        cuts = np.where(np.diff(idx) > 1)[0] + 1
        groups = np.split(idx, cuts)
        if method == 'sum':
            amps = [s[g].sum() for g in groups]
        elif method == 'max':
            amps = [s[g].max() for g in groups]
        else:
            raise ValueError("method must be 'sum' or 'max'")
        out[i] = float(np.median(amps)) if len(amps) else 0.0
    return out


def save_metric_distributions_pdf(
    metrics_table: pd.DataFrame,
    iscell: Iterable,
    pdf_path: str,
    *,
    exclude_columns: Optional[Iterable[str]] = None,
    grid: Tuple[int, int] = (3, 4),
    plots_per_page: Optional[int] = None,
    kde: bool = True,
    alpha: float = 0.5,
    palette = {0: "orange", 1: "blue"},
    probability: bool = False,
) -> str:
    """
    Create a multi-page PDF of histogram distributions for all numeric columns in a
    DataFrame, colored by an independent iscell vector (0/1).

    Parameters
    ----------
    metrics_table : pd.DataFrame
        Input table containing numeric metrics.
    iscell : Iterable
        A required vector of 0/1 values (same length as metrics_table) used for hue coloring.
    pdf_path : str
        Destination path for the PDF. Parent directories are created if needed.
    exclude_columns : Iterable[str] or None, default None
        Columns to exclude from plotting.
    grid : (int, int), default (3, 4)
        Rows, columns for each page's subplot grid.
    plots_per_page : int or None, default None
        Number of plots per page. If None, computed as grid[0] * grid[1].
    kde : bool, default True
        Whether to overlay a kernel density estimate.
    alpha : float, default 0.5
        Bar transparency for histograms.
    palette : mapping or sequence, default {0: 'orange', 1: 'blue'}
        Colors for the hue levels.
    probability : bool, default False
        If True, y-axis will show probabilities; otherwise counts.

    Returns
    -------
    str
        Absolute path to the saved PDF.
    """
    if plots_per_page is None:
        plots_per_page = grid[0] * grid[1]

    # Ensure output directory exists
    os.makedirs(os.path.dirname(os.path.abspath(pdf_path)), exist_ok=True)

    # Attach iscell vector as a column to the table copy
    df = metrics_table.copy()
    df["iscell"] = iscell

    # Build list of numeric metric columns, excluding any extras
    exclude = {"iscell"}
    if exclude_columns:
        exclude.update(exclude_columns)

    metric_columns = (
        df
        .select_dtypes(include=[float, int])
        .columns
        .difference(exclude, sort=False)
    )

    if len(metric_columns) == 0:
        raise ValueError("No numeric columns found to plot after exclusions.")

    # Decide figure size based on grid
    fig_w = 5 * grid[1] / 2  # heuristic for compact pages
    fig_h = 5 * grid[0] / 2

    # Configure seaborn `stat` argument if requested
    sns_histplot_kwargs = {}
    if probability:
        try:
            sns.histplot(data=df, x=metric_columns[0], hue="iscell", stat="probability")
            plt.close()
            sns_histplot_kwargs["stat"] = "probability"
        except TypeError:
            sns_histplot_kwargs = {}

    with PdfPages(pdf_path) as pdf:
        num_pages = math.ceil(len(metric_columns) / plots_per_page)

        for page in range(num_pages):
            start_idx = page * plots_per_page
            end_idx = min(start_idx + plots_per_page, len(metric_columns))
            metrics_subset = metric_columns[start_idx:end_idx]

            fig, axes = plt.subplots(*grid, figsize=(fig_w, fig_h))
            axes = axes.flatten()

            for idx, col in enumerate(metrics_subset):
                ax = axes[idx]
                sns.histplot(
                    data=df,
                    x=col,
                    hue="iscell",
                    hue_order=sorted(df["iscell"].dropna().unique().tolist()),
                    common_norm=False,
                    kde=kde,
                    palette=palette,
                    alpha=alpha,
                    ax=ax,
                    **sns_histplot_kwargs,
                )
                ax.set_title(str(col))
                ylabel = "Probability" if probability else "Count"
                ax.set_ylabel(ylabel)
                ax.set_xlabel("")

            # Remove unused axes
            for j in range(len(metrics_subset), len(axes)):
                fig.delaxes(axes[j])

            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    return os.path.abspath(pdf_path)


def low_high_ratio(psd, fre, low_band, high_band):
    """
    Compute Low/High frequency band ratio for each ROI.

    Parameters
    ----------
    psd : ndarray, shape (n_rois, n_freqs)
        PSD values (already normalized or raw).
    fre : ndarray, shape (n_freqs,)
        Frequency bins corresponding to PSD columns.
    low_band : tuple (f_low, f_high)
        Frequency range for low band (Hz).
    high_band : tuple (f_low, f_high)
        Frequency range for high band (Hz).

    Returns
    -------
    ratios : ndarray, shape (n_rois,)
        Low/High ratio for each ROI.
    """
    psd = np.asarray(psd)
    fre = np.asarray(fre)

    # Masks for frequency ranges
    low_mask = (fre >= low_band[0]) & (fre <= low_band[1])
    high_mask = (fre >= high_band[0]) & (fre <= high_band[1])

    # Integrate (sum) power in each band
    low_power = psd[:, low_mask].mean(axis=1)
    high_power = psd[:, high_mask].mean(axis=1)

    # Avoid divide-by-zero
    ratios = low_power / (high_power + 1e-12)

    return ratios


def _auto_grid(n: int) -> Tuple[int, int]:
    """Choose a near-square (rows, cols) grid big enough for n plots."""
    cols = int(np.ceil(np.sqrt(n)))
    rows = int(np.ceil(n / cols))
    return rows, cols


def save_pairwise_scatter_pdf(
    data: pd.DataFrame,
    metrics_subset: Sequence[str],
    pdf_path: str,
    *,
    hue_col: str = "iscell",
    alpha: float = 0.6,
    point_size: float = 20,
    palette = {0: "orange", 1: "blue"},
    max_plots_per_page: int = 16,
) -> str:
    """
    Plot all pairwise combinations of metrics (columns) colored by `hue_col`
    and save to a multi-page PDF.

    Parameters
    ----------
    data : pd.DataFrame
        Table containing the metrics and a categorical/binary `hue_col`.
    metrics_subset : Sequence[str]
        Column names to pairwise-scatter against each other.
    pdf_path : str
        Destination PDF file path.
    hue_col : str, default 'iscell'
        Column used to color points.
    alpha : float, default 0.6
        Marker alpha.
    point_size : float, default 20
        Marker size passed to seaborn.
    palette : mapping, default {0: 'orange', 1: 'blue'}
        Colors for hue levels.
    max_plots_per_page : int, default 16
        Maximum number of subplots per PDF page (auto-paginates if exceeded).

    Returns
    -------
    str
        Absolute path to the saved PDF.
    """
    # Build all unique pairs
    pairs: List[Tuple[str, str]] = list(itertools.combinations(list(metrics_subset), 2))
    if not pairs:
        raise ValueError("metrics_subset needs at least two metrics to form pairs.")

    # Clean data: only keep columns we need to avoid KeyErrors / NaNs explosions
    cols_needed = set([hue_col]) | set(sum(([x, y] for x, y in pairs), []))
    df = data.loc[:, list(cols_needed)].copy()

    with PdfPages(pdf_path) as pdf:
        # Paginate
        for start in range(0, len(pairs), max_plots_per_page):
            chunk = pairs[start : start + max_plots_per_page]
            n = len(chunk)
            rows, cols = _auto_grid(n)

            fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 5 * rows))
            axes = np.atleast_1d(axes).ravel()

            for i, (x_metric, y_metric) in enumerate(chunk):
                ax = axes[i]
                # Drop rows with NaNs in either metric or the hue
                local = df[[x_metric, y_metric, hue_col]].dropna()
                sns.scatterplot(
                    data=local,
                    x=x_metric,
                    y=y_metric,
                    hue=hue_col,
                    palette=palette,
                    alpha=alpha,
                    s=point_size,
                    ax=ax,
                )
                ax.set_title(f"{x_metric} vs {y_metric}")
                ax.set_xlabel(x_metric)
                ax.set_ylabel(y_metric)
                ax.grid(True, alpha=0.2)

            # Hide any unused axes on the last page
            for j in range(i + 1, len(axes)):
                axes[j].axis("off")

            # Single shared legend at top
            handles, labels = axes[0].get_legend_handles_labels()
            if handles:
                fig.legend(handles, labels, loc="upper center", ncol=min(4, len(labels)), title=hue_col)
                for ax in axes:
                    leg = ax.get_legend()
                    if leg:
                        leg.remove()

            plt.tight_layout(rect=(0, 0, 1, 0.95))
            pdf.savefig(fig)
            plt.close(fig)

    return str(pd.Series([pdf_path]).iloc[0])
