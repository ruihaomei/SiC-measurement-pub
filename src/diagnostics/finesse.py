"""Finesse diagnostics — the heart of the Table 4 vs Table B.7 reconciliation (Study A).

Two distinct quantities must not be compared under one threshold:

1. **Lineshape finesse** ``F_line = FSR / FWHM`` measured from fringe peaks in a spectrum
   (Table 4, F ~ 1.6-1.8 "significant"). This is an *observable* of the actual reflectance.
2. **Ideal-cavity finesse** ``F_cav = pi * sqrt(Re) / (1 - Re)`` with Re an effective reflectance
   (Table B.7, F ~ 0.77-0.98 "weak"). This is a *theoretical* constant of an ideal Fabry-Perot.

They are NOT the same definition, so a single threshold cannot classify both -- this module exposes both
explicitly plus the rigorous round-trip statistic ``|q|`` so the calibration can rest on physics rather
than on conflating amplitude and intensity quantities.
"""
import numpy as np
from scipy.signal import find_peaks, peak_widths


# ---------------------------------------------------------------- lineshape finesse (observable)
def lineshape_finesse(sigma_cm1, R, prominence_factor=0.30, min_peaks=3):
    """Lineshape finesse F = mean(FSR)/mean(FWHM) from fringe peaks (matches the reference finesse_from_data routine).

    Tries both R and 1-R (peaks vs dips) and returns the richer one. Returns a dict with finesse, FSR,
    FWHM, n_peaks. NaN finesse if too few peaks.
    """
    sigma_cm1 = np.asarray(sigma_cm1, float)
    R = np.asarray(R, float)

    def _calc(y):
        prom = max(np.std(y) * prominence_factor, 0.01)
        peaks, _ = find_peaks(y, prominence=prom, distance=max(5, len(y) // 200))
        if len(peaks) < min_peaks:
            return dict(finesse=np.nan, FSR=np.nan, FWHM=np.nan, n_peaks=int(len(peaks)))
        sigma_pk = sigma_cm1[peaks]
        FSR = float(np.nanmean(np.diff(sigma_pk)))
        widths_pts = peak_widths(y, peaks, rel_height=0.5)[0]
        dsigma = (sigma_cm1[-1] - sigma_cm1[0]) / (len(sigma_cm1) - 1)
        FWHM = float(np.nanmean(widths_pts * dsigma))
        fin = FSR / FWHM if FWHM > 0 else np.nan
        return dict(finesse=fin, FSR=FSR, FWHM=FWHM, n_peaks=int(len(peaks)))

    a = _calc(R)
    b = _calc(1.0 - R)
    return a if a["n_peaks"] >= b["n_peaks"] else b


def finesse_to_reflectance(F_line):
    """Map lineshape finesse -> effective reflectance via the cavity relation used in code.

    Code uses Fcoeff = (2 F_line / pi)^2 then R = (sqrt(Fcoeff+1)-1)^2 / Fcoeff
    (the reference R-from-finesse routines). Returned R is an *effective intensity* reflectance.
    """
    F_line = np.asarray(F_line, float)
    Fcoeff = (2.0 * F_line / np.pi) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        R = (np.sqrt(Fcoeff + 1.0) - 1.0) ** 2 / Fcoeff
    return R


# ---------------------------------------------------------------- ideal-cavity finesse (theoretical)
def cavity_finesse_from_R(Re):
    """Ideal Fabry-Perot finesse F = pi*sqrt(Re)/(1-Re) from an effective *intensity* reflectance Re.

    This is the Table B.7 / Eq.(11) form. Re here must be an intensity reflectance in [0,1).
    """
    Re = np.asarray(Re, float)
    return np.pi * np.sqrt(np.clip(Re, 0.0, 0.999999)) / (1.0 - np.clip(Re, 0.0, 0.999999))


def reflectance_to_finesse(Re):
    """Alias of cavity_finesse_from_R for readability."""
    return cavity_finesse_from_R(Re)


def cavity_finesse_from_roundtrip_amplitude(absq):
    """Ideal-cavity proxy from a generalized round-trip amplitude magnitude.

    For unequal and/or absorbing interfaces, ``|q|`` is the sigma-resolved
    round-trip amplitude magnitude. Treating ``|q|`` as the equivalent ideal
    mirror intensity reflectance gives ``pi*sqrt(|q|)/(1-|q|)``. This remains a
    theoretical proxy, not the measured lineshape finesse.
    """
    return cavity_finesse_from_R(absq)


# ---------------------------------------------------------------- rigorous round-trip statistic |q|
def band_q_statistics(q, mask=None):
    """Band statistics of the complex round-trip factor q = r10 r12 e^{2i delta}.

    Returns median and 95th percentile of |q| over the (optionally masked) band. |q| is the rigorous,
    sigma-resolved, absorption-inclusive control parameter for two-beam truncation error.
    """
    aq = np.abs(np.asarray(q))
    if mask is not None:
        aq = aq[mask]
    aq = aq[np.isfinite(aq)]
    if aq.size == 0:
        return dict(median_absq=np.nan, p95_absq=np.nan, max_absq=np.nan, n=0)
    return dict(
        median_absq=float(np.median(aq)),
        p95_absq=float(np.percentile(aq, 95)),
        max_absq=float(np.max(aq)),
        n=int(aq.size),
    )


def two_beam_relative_error(absq):
    """First-order two-beam truncation error magnitude ~ |q|/(1-|q|) (geometric-series tail)."""
    a = np.asarray(absq, float)
    return a / np.clip(1.0 - a, 1e-9, None)
