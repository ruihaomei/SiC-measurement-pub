"""Thickness inversion under a chosen forward model.

For the finesse-calibration study we generate an *exact Airy* spectrum from a known optical response and
then invert it with the *two-beam* model, fitting thickness (+ a small linear dispersion / scale / offset
nuisance, mirroring the manuscript's Appendix-C protocol). The recovered-vs-true thickness bias is the
quantity calibrated against the finesse / |q| diagnostic.

These are deliberately lightweight, well-conditioned fits over a fixed optical response so that the bias
isolates *model-form* error (two-beam vs Airy), not optical-parameter mis-fit.
"""
import numpy as np
from scipy.optimize import least_squares

from ..optics.forward import reflectance_two_beam, reflectance_airy


def _make_eps_callables(eps_f, eps_s):
    eps_f = np.asarray(eps_f, complex)
    eps_s = np.asarray(eps_s, complex)
    return eps_f, eps_s


def fit_thickness_two_beam(sigma_cm1, R_target, eps_f, eps_s, d0_um,
                           theta0_deg=10.0, pol="avg", fit_nuisance=True,
                           d_bounds=(0.5, 30.0)):
    """Fit film thickness with the TWO-BEAM model to a target spectrum.

    Parameters fitted: d_um (thickness in micrometres) plus, if ``fit_nuisance``, a multiplicative
    ``scale`` and additive ``offset`` (instrument baseline) and a small uniform film-index shift ``dn``
    (a linear-dispersion nuisance). ``d_bounds`` is in micrometres and brackets physically plausible
    epilayer thicknesses; the scale/offset/dn bounds (set below) are deliberately tight so the fit
    isolates *model-form* error rather than absorbing it into nuisance parameters. Tolerances
    ``xtol/ftol=1e-12`` request a tight optimum on these well-conditioned synthetic problems.
    Returns dict with d_um, scale, offset, dn, mse, success.
    """
    eps_f, eps_s = _make_eps_callables(eps_f, eps_s)
    sigma = np.asarray(sigma_cm1, float)
    R_target = np.asarray(R_target, float)

    def model(p):
        d_um = p[0]
        scale = p[1] if fit_nuisance else 1.0
        offset = p[2] if fit_nuisance else 0.0
        dn = p[3] if fit_nuisance else 0.0
        eps_f_eff = (np.sqrt(eps_f) + dn) ** 2  # small uniform index shift (dispersion nuisance)
        R = reflectance_two_beam(eps_f_eff, eps_s, d_um * 1e-6, sigma, theta0_deg, pol)
        return scale * R + offset

    if fit_nuisance:
        p0 = np.array([d0_um, 1.0, 0.0, 0.0])
        # bounds: [d_um, scale, offset, dn] — d_um in micrometres (d_bounds); scale near 1,
        # offset near 0 (tight instrument baseline), dn a small +/-0.5 film-index dispersion shift.
        lb = np.array([d_bounds[0], 0.5, -0.2, -0.5])
        ub = np.array([d_bounds[1], 1.5, 0.2, 0.5])
    else:
        p0 = np.array([d0_um])
        lb = np.array([d_bounds[0]])
        ub = np.array([d_bounds[1]])

    res = least_squares(lambda p: model(p) - R_target, p0, bounds=(lb, ub),
                        max_nfev=4000, xtol=1e-12, ftol=1e-12)
    R_fit = model(res.x)
    mse = float(np.mean((R_fit - R_target) ** 2))
    return dict(
        d_um=float(res.x[0]),
        scale=float(res.x[1]) if fit_nuisance else 1.0,
        offset=float(res.x[2]) if fit_nuisance else 0.0,
        dn=float(res.x[3]) if fit_nuisance else 0.0,
        mse=mse, success=bool(res.success), cost=float(res.cost),
    )


def fit_thickness_airy(sigma_cm1, R_target, eps_f, eps_s, d0_um,
                       theta0_deg=10.0, pol="avg", fit_nuisance=True,
                       d_bounds=(0.5, 30.0)):
    """Same as fit_thickness_two_beam but with the exact AIRY model (consistency control)."""
    eps_f, eps_s = _make_eps_callables(eps_f, eps_s)
    sigma = np.asarray(sigma_cm1, float)
    R_target = np.asarray(R_target, float)

    def model(p):
        d_um = p[0]
        scale = p[1] if fit_nuisance else 1.0
        offset = p[2] if fit_nuisance else 0.0
        dn = p[3] if fit_nuisance else 0.0
        eps_f_eff = (np.sqrt(eps_f) + dn) ** 2
        R = reflectance_airy(eps_f_eff, eps_s, d_um * 1e-6, sigma, theta0_deg, pol)
        return scale * R + offset

    if fit_nuisance:
        p0 = np.array([d0_um, 1.0, 0.0, 0.0])
        # bounds: [d_um, scale, offset, dn] — d_um in micrometres (d_bounds); scale near 1,
        # offset near 0 (tight instrument baseline), dn a small +/-0.5 film-index dispersion shift.
        lb = np.array([d_bounds[0], 0.5, -0.2, -0.5])
        ub = np.array([d_bounds[1], 1.5, 0.2, 0.5])
    else:
        p0 = np.array([d0_um]); lb = np.array([d_bounds[0]]); ub = np.array([d_bounds[1]])

    res = least_squares(lambda p: model(p) - R_target, p0, bounds=(lb, ub),
                        max_nfev=4000, xtol=1e-12, ftol=1e-12)
    R_fit = model(res.x)
    return dict(d_um=float(res.x[0]), mse=float(np.mean((R_fit - R_target) ** 2)),
                success=bool(res.success))
