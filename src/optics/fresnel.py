"""Fresnel coefficients, phase thickness, and the multi-beam round-trip factor q.

The round-trip factor ``q = r10 * r12 * exp(2 i beta1)`` is the quantity that rigorously governs the
two-beam truncation error (Study A): the Airy reflectance is the geometric series sum 1/(1-q-ish), and
truncating after the first internal beam leaves a relative error controlled by |q|/(1-|q|). q includes
interface reflection, absorption, and propagation phase, unlike an ideal-cavity finesse constant.
"""
import numpy as np


def cos_theta_complex(n_in, n_out, theta_in_rad):
    """cos(theta_out) via Snell's law with complex indices (matches the reference Snell-angle routine).

    n_in, n_out may be complex arrays; theta_in_rad is the (real) incidence angle in the incident medium.
    """
    sin_in = np.sin(theta_in_rad)
    ratio = (n_in / n_out) * sin_in
    return np.sqrt(1.0 - ratio * ratio + 0j)


def fresnel_rij(n_i, n_j, cos_i, cos_j, pol):
    """Fresnel amplitude reflection coefficient r_ij for 's' or 'p' polarisation.

    Sign/form convention follows the project's oblique-incidence reflectance reference implementation.
    """
    if pol == "s":
        return (n_i * cos_i - n_j * cos_j) / (n_i * cos_i + n_j * cos_j)
    elif pol == "p":
        return (n_j * cos_i - n_i * cos_j) / (n_j * cos_i + n_i * cos_j)
    raise ValueError("pol must be 's' or 'p'")


def phase_thickness(n1, d_m, sigma_cm1, cos1):
    """Single-pass phase thickness delta = k0 * n1 * d * cos(theta1), with k0 = 2*pi*sigma.

    sigma in cm^-1, d in metres -> use k0 = 2*pi/lambda = 2*pi*sigma*100 (1/m). Here we keep the exact
    convention used in the SiC reference model (lambda_m = 1/(sigma*100)), i.e. delta = (2*pi/lambda_m) * n1 * d_m * cos1.
    The round-trip phase is 2*delta.
    """
    lam_m = 1.0 / (np.asarray(sigma_cm1, float) * 100.0)
    k0 = 2.0 * np.pi / lam_m
    return k0 * n1 * d_m * cos1


def roundtrip_factor_q(n0, n1, n2, sigma_cm1, d_m, theta0_deg, pol):
    """Complex round-trip feedback factor q = r10 * r12 * exp(2 i delta) for a single film 0|1|2.

    Returns the complex array q(sigma). |q| is the rigorous local control parameter for two-beam
    truncation error. r10 = -r01 (amplitude), and we compute it directly from the interfaces.
    """
    theta0 = np.deg2rad(theta0_deg)
    cos0 = np.cos(theta0) + 0j
    cos1 = cos_theta_complex(n0, n1, theta0)
    cos2 = cos_theta_complex(n0, n2, theta0)
    # r10 is reflection from medium 1 into 0 = -r01
    r01 = fresnel_rij(n0, n1, cos0, cos1, pol)
    r10 = -r01
    r12 = fresnel_rij(n1, n2, cos1, cos2, pol)
    delta = phase_thickness(n1, d_m, sigma_cm1, cos1)
    return r10 * r12 * np.exp(2j * delta)
