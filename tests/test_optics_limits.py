"""Physics-limit and units tests for the shared optics layer.

These guard the scientific foundation every Phase-3 study builds on:
* two-beam -> Airy/TMM in the weak-feedback limit (|q| -> 0);
* Airy single-film == Abeles TMM for a 3-layer stack (same physics, two code paths);
* finesse definition consistency (the Table 4 vs B.7 distinction is real, not a bug);
* wavenumber/units sanity.
"""
import numpy as np
import pytest

from src.optics import (
    reflectance_two_beam, reflectance_airy, reflectance_abeles,
    roundtrip_factor_q, eps_drude_si,
)
from src.optics.dielectric import n_from_eps
from src.diagnostics import (
    cavity_finesse_from_R, cavity_finesse_from_roundtrip_amplitude,
    finesse_to_reflectance, lineshape_finesse,
)


SIGMA = np.linspace(400.0, 4000.0, 2000)


def _weak_film():
    # Low-contrast, lossless-ish film on a low-contrast substrate -> small |q|.
    eps_f = (1.5 ** 2) * np.ones_like(SIGMA, dtype=complex)
    eps_s = (1.7 ** 2) * np.ones_like(SIGMA, dtype=complex)
    return eps_f, eps_s


def _strong_film():
    # High index contrast at BOTH interfaces -> large |q| -> two-beam departs from Airy.
    # (A modest SiC-like 2.6/3.4 contrast gives median|q|~0.06, still weak; need a genuine high-Q cavity.)
    eps_f = (3.5 ** 2) * np.ones_like(SIGMA, dtype=complex)
    eps_s = (1.5 ** 2) * np.ones_like(SIGMA, dtype=complex)
    return eps_f, eps_s


def test_two_beam_approaches_airy_in_weak_limit():
    eps_f, eps_s = _weak_film()
    d_m = 7.0e-6
    R_tb = reflectance_two_beam(eps_f, eps_s, d_m, SIGMA, theta0_deg=10.0)
    R_airy = reflectance_airy(eps_f, eps_s, d_m, SIGMA, theta0_deg=10.0)
    # |q| is small here; the two should agree closely.
    n0 = np.ones_like(SIGMA, complex)
    q = roundtrip_factor_q(n0, n_from_eps(eps_f), n_from_eps(eps_s), SIGMA, d_m, 10.0, "s")
    assert np.median(np.abs(q)) < 0.15
    assert np.max(np.abs(R_tb - R_airy)) < 0.02


def test_two_beam_deviates_in_strong_limit():
    eps_f, eps_s = _strong_film()
    d_m = 7.0e-6
    R_tb = reflectance_two_beam(eps_f, eps_s, d_m, SIGMA, theta0_deg=10.0)
    R_airy = reflectance_airy(eps_f, eps_s, d_m, SIGMA, theta0_deg=10.0)
    # Larger feedback -> the truncation visibly departs from the exact Airy line shape.
    assert np.max(np.abs(R_tb - R_airy)) > 0.02


def test_airy_equals_abeles_single_film():
    """Airy single-film and Abeles TMM are the same physics; they must agree to tight tolerance."""
    eps_f, eps_s = _strong_film()
    d_um = 7.0
    R_airy = reflectance_airy(eps_f, eps_s, d_um * 1e-6, SIGMA, theta0_deg=12.0, pol="avg")
    n_list = [np.ones_like(SIGMA, complex), n_from_eps(eps_f), n_from_eps(eps_s)]
    R_tmm = reflectance_abeles(SIGMA, 12.0, n_list, [d_um], pol="avg")
    assert np.max(np.abs(R_airy - R_tmm)) < 1e-9


def test_finesse_code_mappings_are_mutual_inverses():
    """The two finesse<->reflectance mappings used in the project code are the SAME physics.

    cavity_finesse_from_R: F = pi*sqrt(Re)/(1-Re)   (Eq. 11 / Table B.7)
    finesse_to_reflectance: Re = (sqrt(Fcoeff+1)-1)^2/Fcoeff with Fcoeff=(2F/pi)^2 (code mapping)
    These are exact mutual inverses (since (2F/pi)^2 = 4Re/(1-Re)^2 = coefficient of finesse).

    KEY FINDING (Study A, P0-B): therefore the Table 4 vs Table B.7 contradiction is NOT a formula bug.
    It arises because Table 4 uses *measured* lineshape finesse (FSR/FWHM from broadened, noisy spectra)
    while Table B.7 uses an *idealized theoretical* effective reflectance. Different inputs, same formula.
    """
    for Re in [0.045, 0.085, 0.161, 0.18, 0.30]:
        F_cav = cavity_finesse_from_R(Re)
        Re_back = finesse_to_reflectance(F_cav)
        assert np.isclose(Re, Re_back, atol=1e-6), (Re, Re_back)


def test_historical_threshold_mapping_is_approximately_correct():
    """F=0.7/1.0/1.5 -> Re ~ 0.045/0.085/0.161, as stated in the revision (sanity, not endorsement).

    This confirms the *arithmetic* of the mapping; Study A separately addresses whether an ideal-cavity
    Re derived from amplitude coefficients is the right physical control parameter (vs band |q|).
    """
    for F_target, Re_claim in [(0.7, 0.045), (1.0, 0.085), (1.5, 0.161)]:
        Re = finesse_to_reflectance(np.array([F_target]))[0]
        assert abs(Re - Re_claim) < 0.01, (F_target, Re, Re_claim)


def test_roundtrip_amplitude_proxy_is_named_and_not_squared():
    """Generalized |q| is already a round-trip amplitude magnitude."""
    absq = 0.2
    assert np.isclose(cavity_finesse_from_roundtrip_amplitude(absq), cavity_finesse_from_R(absq))
    assert not np.isclose(cavity_finesse_from_roundtrip_amplitude(absq), cavity_finesse_from_R(absq**2))


def test_units_wavenumber_positive_and_reflectance_bounded():
    eps_f, eps_s = _weak_film()
    R = reflectance_airy(eps_f, eps_s, 7e-6, SIGMA, theta0_deg=10.0)
    assert np.all(SIGMA > 0)
    assert np.all(R >= -1e-9) and np.all(R <= 1.0 + 1e-9)


def test_lineshape_finesse_on_synthetic_fringes():
    """A clean cosine fringe pattern yields a finite, positive lineshape finesse."""
    sigma = np.linspace(2000, 3200, 2000)
    R = 0.2 + 0.05 * np.cos(2 * np.pi * sigma / 40.0)  # FSR ~ 40 cm^-1
    out = lineshape_finesse(sigma, R)
    assert np.isfinite(out["finesse"]) and out["finesse"] > 0
    assert out["n_peaks"] >= 3
