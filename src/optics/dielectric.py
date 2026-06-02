"""Dielectric-function models.

SiC  : MDF / L-STC (LO-TO phonon factor) + Drude free carriers  [from the reference SiC model's  eps_mdf]
Si   : Drude + optional empirical saturating-log high-frequency correction [from the reference Si model]
SiO2 : double Lorentz oscillator  [from the reference Si model]

All functions take wavenumber sigma [cm^-1] and return complex eps arrays.
"""
import numpy as np

from .constants import (
    TWOPI, C_CM_S, SIGMA_REF, E_CHARGE, EPS0, M_E,
    plasma_omega_from_N, gamma_p_from_mobility,
)


# --------------------------------------------------------------------------- SiC
def eps_mdf_sic(sigma_cm1, eps_inf, wT_cm, wL_cm, GT_cm, GL_cm, wp_rad_s, gamma_p_rad_s):
    """SiC modified-dielectric-function (L-STC LO-TO phonon factor + Drude).

    eps = eps_inf * [ (wL^2 - w^2 - i GL w)/(wT^2 - w^2 - i GT w)  -  wp^2/(w^2 + i gamma_p w) ]

    Phonon parameters (wT, wL, GT, GL) are passed in cm^-1 and converted to rad/s here, matching
    the reference SiC model's eps_mdf called via eps_from_params_advanced.
    """
    w = TWOPI * C_CM_S * np.asarray(sigma_cm1, dtype=np.float64)
    w = w.astype(np.complex128)
    wT = TWOPI * C_CM_S * wT_cm
    wL = TWOPI * C_CM_S * wL_cm
    GT = TWOPI * C_CM_S * GT_cm
    GL = TWOPI * C_CM_S * GL_cm
    phonon = (wL**2 - w**2 - 1j * GL * w) / (wT**2 - w**2 - 1j * GT * w)
    free_carrier = wp_rad_s**2 / (w**2 + 1j * gamma_p_rad_s * w)
    return eps_inf * (phonon - free_carrier)


def eps_mdf_sic_from_Nmu(sigma_cm1, eps_inf, wT_cm, wL_cm, GT_cm, GL_cm,
                         N_cm3, mu_cm2_Vs, m_conc_si, m_mob_si):
    """Convenience: build SiC eps from carrier density/mobility (as in the reference SiC dielectric builder)."""
    wp = plasma_omega_from_N(N_cm3, eps_inf, m_conc_si)
    gp = gamma_p_from_mobility(mu_cm2_Vs, m_mob_si)
    return eps_mdf_sic(sigma_cm1, eps_inf, wT_cm, wL_cm, GT_cm, GL_cm, wp, gp)


# --------------------------------------------------------------------------- Si
def eps_drude_si(sigma_cm1, N_cm3, mu_cm2Vs, eps_inf=11.7, mstar=0.26):
    """Silicon Drude dielectric (matches the reference Si model's drude_eps_Si)."""
    w = TWOPI * C_CM_S * np.asarray(sigma_cm1, dtype=np.float64)
    m_eff = mstar * M_E
    N_m3 = max(N_cm3, 1e10) * 1e6
    mu_m2 = max(mu_cm2Vs, 1e-6) * 1e-4
    wp2 = N_m3 * E_CHARGE * E_CHARGE / (EPS0 * m_eff)
    gamma = E_CHARGE / (m_eff * mu_m2)
    return eps_inf - wp2 / (w * (w + 1j * gamma))


def si_extra_dispersion_real(sigma_cm1, A_disp=0.10):
    """Low-frequency real-dispersion correction A_disp*(SIGMA_REF/sigma)^4 (reference Si model)."""
    return A_disp * (SIGMA_REF / np.maximum(sigma_cm1, 1.0)) ** 4


def sat_log_weight(x, cap=1.2, shape=1.0):
    """Saturating logarithmic weight w(x)=1-exp(-(ln(1+x)/cap)^shape), x>=0 (the reference Si model's _sat_log_weight)."""
    z = np.log1p(np.maximum(x, 0.0))
    return 1.0 - np.exp(-(z / max(cap, 1e-9)) ** max(shape, 1e-6))


def eps_si_total_empirical_hf_correction(
        sigma_cm1, N_cm3, mu_cm2Vs, eps_inf, A_disp, B_damp, B_width,
        sigma_thr=700.0, cap=1.2, shape=1.0, placement="real"):
    """Si epilayer eps with an empirical saturating-log high-frequency correction.

    ``the reference Si model's eps_Si_total`` computes ``-B_damp * w_sat`` as a real
    array and adds it directly to epsilon. The faithful reproduction is therefore
    ``placement="real"``. ``placement="imag"`` is an explicit sensitivity
    variant for auditing the previously documented, but not implemented,
    Im(epsilon) interpretation. Neither placement is a physical replacement.
    """
    base = eps_drude_si(sigma_cm1, N_cm3, mu_cm2Vs, eps_inf)
    low = si_extra_dispersion_real(sigma_cm1, A_disp)
    x = (np.asarray(sigma_cm1, float) - sigma_thr) / max(B_width, 1e-9)
    w_sat = sat_log_weight(x, cap=cap, shape=shape)
    high = -B_damp * w_sat
    if placement == "imag":
        high = 1j * high
    elif placement != "real":
        raise ValueError("placement must be 'real' or 'imag'")
    return base + low + high


# --------------------------------------------------------------------------- SiO2
def eps_lorentz_sio2(sigma_cm1, eps_inf=2.1, S1=0.3, w01=1070.0, g1=40.0,
                     S2=0.15, w02=800.0, g2=60.0, g_min=8.0):
    """Native-oxide SiO2 double-Lorentz dielectric (matches the reference Si model's lorentz_eps_SiO2)."""
    s = np.asarray(sigma_cm1, dtype=np.float64)
    g1e, g2e = max(g1, g_min), max(g2, g_min)
    t1 = S1 * w01**2 / ((w01**2 - s * s) - 1j * g1e * s)
    t2 = S2 * w02**2 / ((w02**2 - s * s) - 1j * g2e * s)
    return eps_inf + t1 + t2


def n_from_eps(eps):
    """Complex refractive index from complex permittivity (principal branch)."""
    return np.sqrt(np.asarray(eps, dtype=np.complex128) + 0j)
