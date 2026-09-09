"""Data loading for the four measured spectra. Reflectance returned as a fraction in [0,~1].

Filename mapping is explicit in ``SPECTRA`` below.
Columns: '波数 (cm-1)', '反射率 (%)'. The SiC 15-degree spectrum contains values >100%
(a calibration artefact) -- see
``load_spectrum`` percent handling.
"""
import os
import numpy as np
import pandas as pd

# Resolve repo root from this file location (src/io_data.py -> repo root)
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = os.path.join(_REPO_ROOT, "data", "raw")

DATASETS = {
    "SiC_10deg": dict(file="SiC_10deg_reflectance.xlsx", material="SiC", angle_deg=10.0),
    "SiC_15deg": dict(file="SiC_15deg_reflectance.xlsx", material="SiC", angle_deg=15.0),
    "Si_10deg":  dict(file="Si_10deg_reflectance.xlsx", material="Si",  angle_deg=10.0),
    "Si_15deg":  dict(file="Si_15deg_reflectance.xlsx", material="Si",  angle_deg=15.0),
}


def load_spectrum(key, percent_to_fraction=True, clamp=False):
    """Load one spectrum by key (e.g. 'SiC_10deg'). Returns (sigma_cm1, R) sorted by sigma.

    percent_to_fraction: divide the reflectance column by 100 (data are in %).
    clamp: if True, clip R into [0, 1] (documents the >100% handling). Default False to preserve raw data;
           the calibration artefact is handled explicitly by analyses that need it.
    """
    meta = DATASETS[key]
    path = os.path.join(RAW_DIR, meta["file"])
    df = pd.read_excel(path)
    cols = list(df.columns)
    wcol = next((c for c in cols if ("波数" in str(c)) or ("cm" in str(c))), cols[0])
    rcol = next((c for c in cols if ("反射" in str(c)) or ("%" in str(c))), cols[1])
    sigma = df[wcol].astype(float).to_numpy()
    R = df[rcol].astype(float).to_numpy()
    m = np.isfinite(sigma) & np.isfinite(R)
    sigma, R = sigma[m], R[m]
    order = np.argsort(sigma)
    sigma, R = sigma[order], R[order]
    if percent_to_fraction:
        R = R / 100.0
    if clamp:
        R = np.clip(R, 0.0, 1.0)
    return sigma, R


def dataset_meta(key):
    """Return a copy of the metadata dict (file, material, angle_deg) for a dataset key."""
    return dict(DATASETS[key])
