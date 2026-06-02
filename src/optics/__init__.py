"""Shared optical forward models for IR reflectance thickness metrology (MEAS-D-26-04372).

Grounded in the project's canonical code (the reference SiC model for SiC MDF+Airy; the reference Si model for
Si Drude+SiO2+Abeles TMM with saturating-log damping). This package provides a clean, tested,
reusable implementation so that every Phase-3 study (finesse, SiC uncertainty, damping, angular)
draws on identical physics.

Conventions
-----------
* Wavenumber ``sigma`` in cm^-1 (the data's native unit).
* Reflectance returned as a fraction in [0, ~1] (NOT percent).
* Complex refractive index ``n_tilde = n + i k``; complex dielectric ``eps = n_tilde**2``.
* Angles in degrees at the API boundary; radians internally.
* Thickness ``d`` in micrometres at the API boundary.
"""
from .constants import (
    E_CHARGE, EPS0, M_E, C_CM_S, TWOPI, SIGMA_REF,
    m_eff_4H_concentration, m_eff_4H_mobility,
)
from .dielectric import (
    eps_mdf_sic, eps_drude_si, eps_lorentz_sio2, eps_si_total_empirical_hf_correction,
    sat_log_weight, si_extra_dispersion_real,
)
from .fresnel import (
    fresnel_rij, cos_theta_complex, phase_thickness, roundtrip_factor_q,
)
from .forward import (
    reflectance_two_beam, reflectance_airy, reflectance_abeles,
    angles_in_stack,
)

__all__ = [
    "E_CHARGE", "EPS0", "M_E", "C_CM_S", "TWOPI", "SIGMA_REF",
    "m_eff_4H_concentration", "m_eff_4H_mobility",
    "eps_mdf_sic", "eps_drude_si", "eps_lorentz_sio2", "eps_si_total_empirical_hf_correction",
    "sat_log_weight", "si_extra_dispersion_real",
    "fresnel_rij", "cos_theta_complex", "phase_thickness", "roundtrip_factor_q",
    "reflectance_two_beam", "reflectance_airy", "reflectance_abeles",
    "angles_in_stack",
]
