"""Jacobian-based conditional covariance, condition number, and parameter correlation matrix.

These are *conditional* on the measurement equation and estimator (R1.2 / R2.4 wording): a fitting
covariance, not a full GUM uncertainty budget. Mirrors the method in the reference Si model's ci_from_jac but exposed
as reusable, testable functions.
"""
import numpy as np


def jacobian_conditional_cov(jac, residuals, n_params=None):
    """Conditional covariance cov = s^2 (J^T J)^{-1}, s^2 = SSR/(m - p).

    jac: (m, p) residual Jacobian (d residual / d param). residuals: (m,) residual vector at the optimum.
    Returns (cov, sigma) where sigma = sqrt(diag(cov)) is the per-parameter conditional standard uncertainty.
    """
    J = np.asarray(jac, float)
    r = np.asarray(residuals, float)
    m, p = J.shape
    if n_params is None:
        n_params = p
    dof = max(m - n_params, 1)
    s2 = float(r @ r) / dof
    JTJ = J.T @ J
    try:
        cov = np.linalg.inv(JTJ) * s2
    except np.linalg.LinAlgError:
        # rank-deficient / singular normal matrix (a parameter is unidentified at the optimum):
        # use the Moore-Penrose pseudo-inverse so the identified directions still get a covariance.
        cov = np.linalg.pinv(JTJ, rcond=1e-12) * s2
    sigma = np.sqrt(np.maximum(np.diag(cov), 0.0))
    return cov, sigma


def condition_number(jac, weights=None):
    """Condition number kappa(J^T W J) (or kappa(J) if weights None). Large -> ill-conditioned/poorly identified."""
    J = np.asarray(jac, float)
    if weights is not None:
        w = np.sqrt(np.asarray(weights, float))
        J = J * w[:, None]
    JTJ = J.T @ J
    sv = np.linalg.svd(JTJ, compute_uv=False)
    sv = sv[sv > 0]
    if sv.size == 0:
        return np.inf
    return float(sv[0] / sv[-1])


def parameter_correlations(cov):
    """Full parameter correlation matrix from a covariance matrix."""
    cov = np.asarray(cov, float)
    d = np.sqrt(np.maximum(np.diag(cov), 1e-300))
    corr = cov / np.outer(d, d)
    return corr
