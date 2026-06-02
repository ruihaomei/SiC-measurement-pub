"""Silicon 4-layer stack forward model (air / SiO2 / Si-epi / Si-substrate) with selectable corrections.

Mirrors the reference Si model's ideal_reflectance_with_averaging but exposes the high-frequency correction as a
*pluggable mechanism* so Study C can compare:
  - 'none'        : Drude epi + Drude sub, no high-frequency correction (baseline)
  - 'empirical_real_correction': faithful reference saturating-log subtraction from Re(eps_epi)
  - 'empirical_imag_sensitivity': optional subtraction from Im(eps_epi), sensitivity audit only
  - 'angle_avg'   : baseline + Gaussian averaging over incidence angle (objective/aperture spread)
  - 'thick_avg'   : baseline + Gaussian averaging over film thickness (spot non-uniformity)
  - 'roughness'   : baseline + Debye-Waller-like interface coherence attenuation (physical, KK-agnostic envelope)

All return reflectance over sigma (fraction). Thicknesses: d_um (epi), t_ox_nm (oxide). SiO2 double-Lorentz.
"""
import numpy as np

from .dielectric import (
    eps_drude_si, eps_lorentz_sio2, eps_si_total_empirical_hf_correction, n_from_eps,
)
from .forward import reflectance_abeles


def _gauss_nodes(mu, sig, n):
    if sig <= 0 or n <= 1:
        return np.array([mu], float), np.array([1.0], float)
    xs = np.linspace(mu - 3 * sig, mu + 3 * sig, n)
    w = np.exp(-0.5 * ((xs - mu) / sig) ** 2)
    w /= w.sum()
    return xs, w


def si_reflectance(sigma, angle_deg, d_um, t_ox_nm, N_e, mu_e, N_s, mu_s,
                   eps_inf_e, eps_inf_s, lorentz_params, A_disp,
                   mechanism="none", B_damp=0.0, B_width=600.0, sigma_thr=700.0,
                   angle_sigma_deg=0.40, d_jitter_um=0.03, roughness_nm=0.0,
                   n_angle=7, n_d=3, n_tox=3):
    """Si 4-layer reflectance with a selectable high-frequency mechanism.

    lorentz_params = (S1, w01, g1, S2, w02, g2) for the SiO2 oxide.
    """
    S1, w01, g1, S2, w02, g2 = lorentz_params
    eps_ox = eps_lorentz_sio2(sigma, 2.1, S1, w01, g1, S2, w02, g2)

    # the reference's implemented empirical term changes Re(eps). The Im(eps) case is an
    # explicitly separate sensitivity variant because earlier prose conflated them.
    if mechanism == "empirical_real_correction":
        eps_e = eps_si_total_empirical_hf_correction(
            sigma, N_e, mu_e, eps_inf_e, A_disp, B_damp, B_width, sigma_thr, placement="real")
    elif mechanism == "empirical_imag_sensitivity":
        eps_e = eps_si_total_empirical_hf_correction(
            sigma, N_e, mu_e, eps_inf_e, A_disp, B_damp, B_width, sigma_thr, placement="imag")
    else:
        # baseline epi = Drude + real low-frequency dispersion only
        from .dielectric import si_extra_dispersion_real
        eps_e = eps_drude_si(sigma, N_e, mu_e, eps_inf_e) + si_extra_dispersion_real(sigma, A_disp)
    eps_s = eps_drude_si(sigma, N_s, mu_s, eps_inf_s)

    n_air = np.ones_like(sigma, complex)
    n_list = [n_air, n_from_eps(eps_ox), n_from_eps(eps_e), n_from_eps(eps_s)]

    # averaging dimensions depend on mechanism
    if mechanism == "angle_avg":
        angs, wa = _gauss_nodes(angle_deg, angle_sigma_deg, n_angle)
    else:
        angs, wa = np.array([angle_deg]), np.array([1.0])
    if mechanism == "thick_avg":
        dvals, wd = _gauss_nodes(d_um, d_jitter_um, n_d)
    else:
        dvals, wd = np.array([d_um]), np.array([1.0])

    R = np.zeros_like(sigma, float); W = 0.0
    for th, wth in zip(angs, wa):
        for dd, wdi in zip(dvals, wd):
            d_list = [t_ox_nm * 1e-3, dd]
            Rc = reflectance_abeles(sigma, th, n_list, d_list, "avg")
            R += wth * wdi * Rc
            W += wth * wdi
    R = R / max(W, 1e-12)

    # roughness / coherence attenuation: Debye-Waller-like envelope on the *fringe contrast* at high sigma.
    # Physical surrogate: rms interface roughness scatters coherent light ~ exp(-(2 sigma_cm * 2pi * rough)^2/2)
    if mechanism == "roughness" and roughness_nm > 0:
        rough_cm = roughness_nm * 1e-7
        dw = np.exp(-0.5 * (2.0 * np.pi * 2.0 * sigma * rough_cm) ** 2)
        Rbar = np.mean(R)
        R = Rbar + dw * (R - Rbar)  # attenuate fringe contrast toward the mean, preserve baseline
    known = {
        "none", "empirical_real_correction", "empirical_imag_sensitivity",
        "angle_avg", "thick_avg", "roughness",
    }
    if mechanism not in known:
        raise ValueError(f"unknown Si high-frequency mechanism: {mechanism}")
    return R
