"""Physical constants and SiC effective masses (CODATA / project values).

Values match the canonical scripts (the reference SiC model, the reference Si model) so reproductions are bit-comparable.
"""
import math

E_CHARGE = 1.602176634e-19      # C
EPS0 = 8.8541878128e-12         # F/m
M_E = 9.10938356e-31            # kg
C_CM_S = 2.99792458e10          # cm/s  (speed of light in cm/s)
TWOPI = 2.0 * math.pi
SIGMA_REF = 1000.0              # cm^-1, reference wavenumber for Si extra-dispersion term


def m_eff_4H_concentration(m_MK: float = 0.31, m_MG: float = 0.58) -> float:
    """Density-of-states effective mass for 4H-SiC (geometric mean), in kg."""
    return math.sqrt(m_MK * m_MG) * M_E


def m_eff_4H_mobility(m_MK: float = 0.31, m_MG: float = 0.58) -> float:
    """Conductivity (mobility) effective mass for 4H-SiC, in kg."""
    num = (m_MK ** -0.5 + m_MG ** -0.5)
    den = (m_MK ** -1.5 + m_MG ** -1.5)
    return (num / den) * M_E


def m_eff_6H_isotropic() -> float:
    """Isotropic effective mass for 6H-SiC, in kg."""
    return 0.35 * M_E


def wnum_to_omega(sigma_cm1):
    """Angular frequency omega [rad/s] from wavenumber sigma [cm^-1]."""
    return TWOPI * C_CM_S * sigma_cm1


def plasma_omega_from_N(N_cm3: float, eps_inf: float, m_eff_si: float) -> float:
    """Plasma angular frequency [rad/s] from carrier density N [cm^-3]."""
    N_m3 = max(N_cm3, 1e10) * 1e6
    return math.sqrt(N_m3 * E_CHARGE * E_CHARGE / (EPS0 * eps_inf * m_eff_si))


def gamma_p_from_mobility(mu_cm2_Vs: float, m_eff_mob_si: float) -> float:
    """Drude damping rate [rad/s] from mobility mu [cm^2/(V s)]."""
    mu_m2 = max(mu_cm2_Vs, 1e-6) * 1e-4
    return E_CHARGE / (m_eff_mob_si * mu_m2)
