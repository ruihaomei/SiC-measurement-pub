"""Forward reflectance models: two-beam, Airy (single film), Abeles TMM (multilayer).

Hierarchy (R2.1 framing):
* ``reflectance_abeles`` / ``reflectance_airy`` -- rigorous electromagnetic baseline (multiple beams).
* ``reflectance_two_beam`` -- the truncated (front + first internal beam) simplification whose adequacy
  is diagnosed. It is the SAME as Airy but with the denominator (1 + r01 r12 e^{2i delta}) set to 1.

Airy single-film unpolarised reflectance (the reference Airy implementation):
    r_pol = (r01 + r12 e^{2i delta}) / (1 + r01 r12 e^{2i delta}),   R = 0.5 (|r_s|^2 + |r_p|^2)
Two-beam: numerator only (denominator -> 1), i.e. r = r01 + t01 r12 t10 e^{2i delta} with the
thin-film identity t01 t10 = 1 - r01^2; we implement the standard truncation r ≈ r01 + (1-r01^2) r12 e^{2i delta}.
"""
import numpy as np

from .fresnel import cos_theta_complex, fresnel_rij, phase_thickness


# --------------------------------------------------------------- single film 0|1|2
def _film_amplitudes(n0, n1, n2, sigma_cm1, d_m, theta0_deg, pol):
    """Interface amplitude reflections (r01, r12) and round-trip phasor exp(2i*delta) for a 0|1|2 film."""
    theta0 = np.deg2rad(theta0_deg)
    cos0 = np.cos(theta0) + 0j
    cos1 = cos_theta_complex(n0, n1, theta0)
    cos2 = cos_theta_complex(n0, n2, theta0)
    r01 = fresnel_rij(n0, n1, cos0, cos1, pol)
    r12 = fresnel_rij(n1, n2, cos1, cos2, pol)
    delta = phase_thickness(n1, d_m, sigma_cm1, cos1)
    e2i = np.exp(2j * delta)
    return r01, r12, e2i


def reflectance_airy(eps_f, eps_s, d_m, sigma_cm1, theta0_deg=10.0, pol="avg", n0=1.0):
    """Exact single-film Airy reflectance (the rigorous baseline for a 3-layer 0|film|substrate stack)."""
    n0 = np.asarray(n0, dtype=np.complex128) * np.ones_like(np.asarray(sigma_cm1, float), dtype=np.complex128)
    n1 = np.sqrt(np.asarray(eps_f, np.complex128) + 0j)
    n2 = np.sqrt(np.asarray(eps_s, np.complex128) + 0j)

    def one(pol_):
        r01, r12, e2i = _film_amplitudes(n0, n1, n2, sigma_cm1, d_m, theta0_deg, pol_)
        r = (r01 + r12 * e2i) / (1.0 + r01 * r12 * e2i)
        return np.abs(r) ** 2

    if pol == "avg":
        return 0.5 * (one("s") + one("p"))
    return one(pol)


def reflectance_two_beam(eps_f, eps_s, d_m, sigma_cm1, theta0_deg=10.0, pol="avg", n0=1.0):
    """Two-beam (front + first internal beam) truncation of the Airy series.

    r ≈ r01 + (1 - r01^2) r12 e^{2i delta}   (denominator dropped). This is the simplified algorithm
    whose breakdown the finesse diagnostic is meant to detect.
    """
    n0 = np.asarray(n0, dtype=np.complex128) * np.ones_like(np.asarray(sigma_cm1, float), dtype=np.complex128)
    n1 = np.sqrt(np.asarray(eps_f, np.complex128) + 0j)
    n2 = np.sqrt(np.asarray(eps_s, np.complex128) + 0j)

    def one(pol_):
        r01, r12, e2i = _film_amplitudes(n0, n1, n2, sigma_cm1, d_m, theta0_deg, pol_)
        r = r01 + (1.0 - r01 * r01) * r12 * e2i
        return np.abs(r) ** 2

    if pol == "avg":
        return 0.5 * (one("s") + one("p"))
    return one(pol)


# --------------------------------------------------------------- multilayer Abeles TMM
def angles_in_stack(n_list, theta0_deg):
    """Snell angles through a stack given the list of (complex) indices (the reference stack-angle routine)."""
    th0 = np.deg2rad(theta0_deg)
    n0 = n_list[0]
    th = [th0 * np.ones_like(np.asarray(n0))]
    for j in range(1, len(n_list)):
        th.append(np.arcsin(n0 * np.sin(th0) / n_list[j]))
    return th


def _q_tilt(n, th, pol):
    """Tilted optical admittance for the Abeles matrix: n*cos(theta) for s-pol, n/cos(theta) for p-pol."""
    ct = np.cos(th)
    return n * ct if pol == "s" else n / np.cos(th)


def reflectance_abeles(sigma_cm1, theta0_deg, n_list, d_list_um, pol="avg"):
    """Abeles transfer-matrix reflectance for a stack (matches the reference Si model's reflectance_abeles).

    n_list: [n_incident, n_layer1, ..., n_substrate] (complex arrays over sigma).
    d_list_um: thicknesses of the *internal* layers (len = len(n_list) - 2), in micrometres.
    """
    sigma_cm1 = np.asarray(sigma_cm1, float)
    L = len(n_list) - 2
    th = angles_in_stack(n_list, theta0_deg)

    def calc(polx):
        q0 = _q_tilt(n_list[0], th[0], polx)
        qS = _q_tilt(n_list[-1], th[-1], polx)
        M11 = np.ones_like(sigma_cm1, complex); M12 = np.zeros_like(sigma_cm1, complex)
        M21 = np.zeros_like(sigma_cm1, complex); M22 = np.ones_like(sigma_cm1, complex)
        for j in range(1, L + 1):
            nj, thj = n_list[j], th[j]
            dj_cm = d_list_um[j - 1] * 1e-4              # micrometres -> centimetres (sigma is in cm^-1)
            # phase thickness delta = 2*pi * sigma[cm^-1] * n * cos(theta) * d[cm]  (dimensionless)
            delta = 2.0 * np.pi * sigma_cm1 * nj * np.cos(thj) * dj_cm
            c, s = np.cos(delta), np.sin(delta)
            qj = _q_tilt(nj, thj, polx)
            A11, A12, A21, A22 = c, 1j * s / qj, 1j * qj * s, c
            T11 = M11 * A11 + M12 * A21
            T12 = M11 * A12 + M12 * A22
            T21 = M21 * A11 + M22 * A21
            T22 = M21 * A12 + M22 * A22
            M11, M12, M21, M22 = T11, T12, T21, T22
        num = (q0 * M11 + q0 * qS * M12 - M21 - qS * M22)
        den = (q0 * M11 + q0 * qS * M12 + M21 + qS * M22)
        den = np.where(np.abs(den) < 1e-30, 1e-30 + 0j, den)
        r = num / den
        return np.real(r * np.conj(r))

    if pol == "avg":
        return 0.5 * (calc("s") + calc("p"))
    return calc(pol)
