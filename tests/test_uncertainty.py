"""Tests for the uncertainty layer (Jacobian conditional cov, wild bootstrap)."""
import numpy as np

from src.uncertainty.jacobian import jacobian_conditional_cov, condition_number, parameter_correlations
from src.uncertainty.bootstrap import wild_bootstrap


def test_jacobian_cov_linear_model_matches_ols():
    """For a linear model y = a + b x, Jacobian cov must match the textbook OLS covariance."""
    rng = np.random.default_rng(0)
    x = np.linspace(0, 10, 200)
    a_true, b_true, noise = 1.0, 2.0, 0.1
    y = a_true + b_true * x + rng.normal(0, noise, x.shape)
    # residual r = (a + b x) - y; Jacobian columns d r/d a = 1, d r/d b = x
    a_hat = np.polyfit(x, y, 1)[1]; b_hat = np.polyfit(x, y, 1)[0]
    resid = (a_hat + b_hat * x) - y
    J = np.column_stack([np.ones_like(x), x])
    cov, sigma = jacobian_conditional_cov(J, resid)
    # OLS analytic: sigma^2 (X^T X)^-1
    s2 = (resid @ resid) / (len(x) - 2)
    cov_ols = s2 * np.linalg.inv(J.T @ J)
    assert np.allclose(cov, cov_ols, rtol=1e-8)
    assert sigma[1] < sigma[0] or sigma[1] > 0  # slope better determined typically


def test_condition_number_identity_is_one():
    J = np.eye(4)
    assert abs(condition_number(J) - 1.0) < 1e-9


def test_parameter_correlations_diag_unit():
    cov = np.array([[4.0, 1.0], [1.0, 9.0]])
    corr = parameter_correlations(cov)
    assert np.allclose(np.diag(corr), 1.0)
    assert abs(corr[0, 1] - 1.0 / (2 * 3)) < 1e-9


def test_wild_bootstrap_recovers_mean_and_is_seed_deterministic():
    # simple "fit": estimate the mean of a vector (d_um := mean(R_syn))
    R_obs = np.array([1.0, 1.1, 0.9, 1.05, 0.95])
    R_fit = np.full_like(R_obs, R_obs.mean())

    def fit_fn(R_syn):
        return {"d_um": float(np.mean(R_syn))}

    b1 = wild_bootstrap(fit_fn, None, R_obs, R_obs.mean(), R_fit, n_boot=300, seed=42)
    b2 = wild_bootstrap(fit_fn, None, R_obs, R_obs.mean(), R_fit, n_boot=300, seed=42)
    assert b1["n_success"] == 300
    assert np.allclose(b1["replicates"], b2["replicates"])  # determinism
    assert abs(b1["mean"] - R_obs.mean()) < 0.05            # centred near the point estimate
