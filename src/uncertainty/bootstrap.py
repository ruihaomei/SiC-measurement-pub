"""Seeded wild bootstrap for nonlinear least-squares thickness estimates (R1.2).

Wild bootstrap resamples residuals with random signs (Rademacher) preserving heteroskedasticity, refits,
and returns the distribution of the parameter(s) of interest. Seed is recorded for reproducibility; all
replicate values are returned so they can be saved to CSV.
"""
import numpy as np


def wild_bootstrap(fit_fn, sigma, R_obs, p_hat, R_fit, n_boot=1000, seed=20260601):
    """Generic wild bootstrap.

    fit_fn(R_synthetic) -> dict with at least key 'd_um' (or returns a float thickness). Must refit on a
    given synthetic reflectance vector and return the estimate(s).
    sigma, R_obs: data. p_hat: point estimate (for reference). R_fit: model at the optimum (for residuals).
    Returns dict(replicates=ndarray, mean, std, ci95=(lo,hi), seed, n_boot).
    """
    rng = np.random.default_rng(seed)
    resid = np.asarray(R_obs, float) - np.asarray(R_fit, float)
    Rf = np.asarray(R_fit, float)
    reps = []
    for _ in range(n_boot):
        signs = rng.choice([-1.0, 1.0], size=resid.shape)  # Rademacher
        R_syn = Rf + signs * resid
        try:
            out = fit_fn(R_syn)
            val = out["d_um"] if isinstance(out, dict) else float(out)
            if np.isfinite(val):
                reps.append(val)
            # non-finite refits are dropped; the kept count is reported as n_success below
        except Exception:
            continue  # a failed refit is skipped; failures = n_boot - n_success
    reps = np.asarray(reps, float)
    if reps.size == 0:
        return dict(replicates=reps, mean=np.nan, std=np.nan, ci95=(np.nan, np.nan),
                    seed=seed, n_boot=n_boot, n_success=0)
    lo, hi = np.percentile(reps, [2.5, 97.5])
    return dict(replicates=reps, mean=float(np.mean(reps)), std=float(np.std(reps, ddof=1)),
                ci95=(float(lo), float(hi)), seed=seed, n_boot=n_boot, n_success=int(reps.size))
