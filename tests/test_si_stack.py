"""Tests for the Si 4-layer stack model and its high-frequency mechanisms."""
import numpy as np
import pytest

from src.optics.si_stack import si_reflectance
from src.optics.dielectric import (
    eps_drude_si,
    eps_si_total_empirical_hf_correction,
    si_extra_dispersion_real,
)
from src.reference_params import SI_STACK as S

SIGMA = np.linspace(700.0, 4000.0, 1200)
LOR = (S["S1"], S["w01"], S["g1"], S["S2"], S["w02"], S["g2"])


def _base(mech, **kw):
    return si_reflectance(SIGMA, 10.0, S["d_um"], S["t_ox_nm"], S["N_epi"], S["mu_epi"],
                          S["N_sub"], S["mu_sub"], S["eps_inf_e"], S["eps_inf_s"], LOR,
                          A_disp=S["A_disp"], mechanism=mech, **kw)


def test_reflectance_bounded_all_mechanisms():
    for mech in ["none", "empirical_real_correction", "empirical_imag_sensitivity",
                 "angle_avg", "thick_avg", "roughness"]:
        kw = {}
        if mech in {"empirical_real_correction", "empirical_imag_sensitivity"}:
            kw = dict(B_damp=0.5, B_width=600.0)
        elif mech == "roughness":
            kw = dict(roughness_nm=2.0)
        R = _base(mech, **kw)
        assert np.all(R >= -1e-6) and np.all(R <= 1.0 + 1e-6), mech
        assert np.all(np.isfinite(R)), mech


def test_empirical_real_correction_changes_high_sigma_reflectance_structure():
    """The faithful reference Re(epsilon) correction changes the high-sigma spectrum."""
    R0 = _base("none")
    Rd = _base("empirical_real_correction", B_damp=0.5, B_width=600.0)
    hi = SIGMA > 2000
    assert np.max(np.abs(Rd - R0)[hi]) > 1e-3


def test_zero_knob_mechanisms_reduce_to_baseline():
    """Averaging with zero spread / zero roughness must equal the baseline."""
    R0 = _base("none")
    assert np.allclose(_base("angle_avg", angle_sigma_deg=0.0), R0, atol=1e-9)
    assert np.allclose(_base("thick_avg", d_jitter_um=0.0), R0, atol=1e-9)
    assert np.allclose(_base("roughness", roughness_nm=0.0), R0, atol=1e-9)


def test_empirical_hf_placement_is_explicit_and_component_specific():
    base = (
        eps_drude_si(SIGMA, S["N_epi"], S["mu_epi"], S["eps_inf_e"])
        + si_extra_dispersion_real(SIGMA, S["A_disp"])
    )
    real_variant = eps_si_total_empirical_hf_correction(
        SIGMA, S["N_epi"], S["mu_epi"], S["eps_inf_e"], S["A_disp"], 0.5, 600.0,
        placement="real")
    imag_variant = eps_si_total_empirical_hf_correction(
        SIGMA, S["N_epi"], S["mu_epi"], S["eps_inf_e"], S["A_disp"], 0.5, 600.0,
        placement="imag")
    assert np.allclose(np.imag(real_variant), np.imag(base))
    assert np.max(np.abs(np.real(real_variant) - np.real(base))) > 0
    assert np.allclose(np.real(imag_variant), np.real(base))
    assert np.max(np.abs(np.imag(imag_variant) - np.imag(base))) > 0


def test_ambiguous_damping_mechanism_name_is_rejected():
    with pytest.raises(ValueError, match="unknown Si high-frequency mechanism"):
        _base("damping", B_damp=0.5, B_width=600.0)
