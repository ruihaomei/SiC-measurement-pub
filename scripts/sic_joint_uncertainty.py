#!/usr/bin/env python3
"""Study B — SiC joint multi-angle inversion + conditional uncertainty (addresses R1.2, R2.5).

Provides a single, formally-defensible reported SiC thickness with a conditional uncertainty interval,
WITHOUT any external-trueness claim. Three parts:

  (1) Independent-angle multi-beam (Airy) fits  -> retained as the cross-angle *invariance* diagnostic.
  (2) Joint multi-angle fit: a SHARED physical thickness d and shared optical-response parameters across
      10 deg and 15 deg, with per-angle scale nuisances -> ONE reported thickness.
  (3) Conditional uncertainty on the joint d:
        - Jacobian-based conditional standard uncertainty (cov = s^2 (J^T J)^-1),
        - seeded wild bootstrap 95% interval (replicates saved),
        - a small numerical perturbation budget (band choice, subsample/resolution) as model-choice
          sensitivity, clearly separated from the fitting interval.

Terminology: "conditional uncertainty" (conditional on the measurement equation + estimator); NOT a full
GUM budget, NOT a trueness claim.

Outputs: outputs/tables/sic_results.csv, outputs/tables/sic_uncertainty_budget.csv,
         outputs/bootstrap/sic_{10deg,15deg,joint}_replicates.csv,
         outputs/figures/sic_bootstrap.png, manifest rows.
Deterministic: bootstrap seed fixed and recorded.
"""
import os
import sys
import json
import csv
import datetime
import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.io_data import load_spectrum
from src.optics.forward import reflectance_airy
from src.optics.dielectric import eps_mdf_sic_from_Nmu
from src.optics.constants import m_eff_4H_concentration, m_eff_4H_mobility
from src.uncertainty.jacobian import jacobian_conditional_cov, condition_number
from src.uncertainty.bootstrap import wild_bootstrap
from src.reference_params import SIC_4H
from src.provenance import update_manifest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(REPO, "outputs", "tables")
BOOT = os.path.join(REPO, "outputs", "bootstrap")
FIG = os.path.join(REPO, "outputs", "figures")
SEED = 20260601
SUBSAMPLE = 3
N_BOOT = 1000
MC, MM = m_eff_4H_concentration(), m_eff_4H_mobility()
P = SIC_4H

# per-angle physical params: logNf, muf, logNs, mus, GTf, GLf, GTs, GLs  (8) + shared d
PNAMES_PHYS = ["log10_N_film", "mu_film", "log10_N_sub", "mu_sub", "GT_film", "GL_film", "GT_sub", "GL_sub"]


def _eps(sigma, logN, mu, GT, GL):
    return eps_mdf_sic_from_Nmu(sigma, P["eps_inf"], P["wT_cm"], P["wL_cm"], GT, GL, 10**logN, mu, MC, MM)


def _model_one(sigma, angle, d_um, phys, scale=1.0):
    logNf, muf, logNs, mus, GTf, GLf, GTs, GLs = phys
    ef = _eps(sigma, logNf, muf, GTf, GLf)
    es = _eps(sigma, logNs, mus, GTs, GLs)
    return scale * reflectance_airy(ef, es, d_um * 1e-6, sigma, angle, "avg")


def fit_independent(sigma, R, angle):
    def resid(x):
        d = x[0]; phys = x[1:9]
        return _model_one(sigma, angle, d, phys) - R
    x0 = np.array([P["d_um"], P["log10_N_film"], P["mu_film_cm2_Vs"], P["log10_N_sub"],
                   P["mu_sub_cm2_Vs"], P["GT_film_cm"], P["GL_film_cm"], P["GT_sub_cm"], P["GL_sub_cm"]])
    lb = np.array([5, 15, 10, 16, 10, 1, 1, 1, 1.])
    ub = np.array([10, 19.8, 800, 19.8, 500, 4, 15, 4, 15.])
    r = least_squares(resid, np.clip(x0, lb, ub), bounds=(lb, ub), loss="soft_l1", f_scale=0.01, max_nfev=4000)
    return r


def load_all(band=(400, 4000), sub=SUBSAMPLE):
    data = {}
    for key, ang in [("SiC_10deg", 10.0), ("SiC_15deg", 15.0)]:
        s, R = load_spectrum(key)
        m = (s >= band[0]) & (s <= band[1])
        data[ang] = (s[m][::sub], R[m][::sub])
    return data


def fit_joint(data):
    """Shared d + shared optical params across angles; per-angle scale. Returns (x, residual_fn, sizes)."""
    angles = sorted(data.keys())
    sigmas = [data[a][0] for a in angles]
    Rs = [data[a][1] for a in angles]
    # params: [d, 8 shared phys, scale_per_angle...]
    nphys = 8
    nsc = len(angles)

    def split(x):
        d = x[0]; phys = x[1:1 + nphys]; scales = x[1 + nphys:1 + nphys + nsc]
        return d, phys, scales

    def resid(x):
        d, phys, scales = split(x)
        parts = []
        for a, s, R, sc in zip(angles, sigmas, Rs, scales):
            parts.append(_model_one(s, a, d, phys, sc) - R)
        return np.concatenate(parts)

    x0 = np.array([P["d_um"], P["log10_N_film"], P["mu_film_cm2_Vs"], P["log10_N_sub"],
                   P["mu_sub_cm2_Vs"], P["GT_film_cm"], P["GL_film_cm"], P["GT_sub_cm"], P["GL_sub_cm"]]
                  + [1.0] * nsc)
    lb = np.array([5, 15, 10, 16, 10, 1, 1, 1, 1.] + [0.5] * nsc)
    ub = np.array([10, 19.8, 800, 19.8, 500, 4, 15, 4, 15.] + [1.5] * nsc)
    r = least_squares(resid, np.clip(x0, lb, ub), bounds=(lb, ub), loss="soft_l1", f_scale=0.01, max_nfev=6000)
    return r, resid, angles


def main():
    for d in (TAB, BOOT, FIG):
        os.makedirs(d, exist_ok=True)
    log = lambda m: print(f"[{datetime.datetime.now().isoformat(timespec='seconds')}] {m}", flush=True)

    data = load_all()

    # (1) independent-angle fits -> invariance diagnostic
    indep = {}
    indep_fits = {}
    for a in sorted(data):
        s, R = data[a]
        r = fit_independent(s, R, a)
        indep[a] = float(r.x[0])
        indep_fits[a] = r
        log(f"independent {a:.0f}deg: d={r.x[0]:.4f} um")
    d10, d15 = indep[10.0], indep[15.0]
    cross_angle_diff = abs(d10 - d15)

    # Saved independent-angle replicate files reproduce the per-angle intervals
    # separately from the joint shared-d result.
    independent_boot = {}
    for a in sorted(data):
        sigma_a, R_obs_a = data[a]
        fit_a = indep_fits[a]
        d_a = float(fit_a.x[0])
        R_fit_a = _model_one(sigma_a, a, d_a, fit_a.x[1:9])

        def refit_angle(R_syn, sigma=sigma_a, angle=a):
            rr = fit_independent(sigma, R_syn, angle)
            return {"d_um": float(rr.x[0])}

        boot_a = wild_bootstrap(refit_angle, sigma_a, R_obs_a, d_a, R_fit_a, n_boot=N_BOOT, seed=SEED)
        independent_boot[a] = boot_a
        path = os.path.join(BOOT, f"sic_{int(a)}deg_replicates.csv")
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["replicate_d_um"])
            for value in boot_a["replicates"]:
                w.writerow([value])
        log(f"independent {a:.0f}deg bootstrap ci95=({boot_a['ci95'][0]:.4f},{boot_a['ci95'][1]:.4f}) "
            f"n={boot_a['n_success']}")

    # (2) joint fit -> single reported thickness
    rj, resid_fn, angles = fit_joint(data)
    d_joint = float(rj.x[0])
    log(f"joint d={d_joint:.4f} um (shared across {angles})")

    # (3a) Jacobian conditional uncertainty on joint d
    J = rj.jac
    r_at_opt = rj.fun
    cov, sigma = jacobian_conditional_cov(J, r_at_opt)
    sigma_d = float(sigma[0])
    kappa = condition_number(J)
    log(f"Jacobian sigma_d={sigma_d:.4f} um, kappa(J^TJ)={kappa:.3e}")

    # (3b) seeded wild bootstrap on joint d
    R_fit_concat = np.concatenate([_model_one(data[a][0], a, d_joint, rj.x[1:9], rj.x[9 + i])
                                   for i, a in enumerate(angles)])
    R_obs_concat = np.concatenate([data[a][1] for a in angles])

    def refit(R_syn):
        # rebuild per-angle synthetic, refit joint, return d
        sizes = [len(data[a][0]) for a in angles]
        idx = np.cumsum([0] + sizes)
        syn = {a: (data[a][0], R_syn[idx[i]:idx[i + 1]]) for i, a in enumerate(angles)}
        rr, _, _ = fit_joint(syn)
        return {"d_um": float(rr.x[0])}

    boot = wild_bootstrap(refit, None, R_obs_concat, d_joint, R_fit_concat, n_boot=N_BOOT, seed=SEED)
    log(f"bootstrap d: mean={boot['mean']:.4f} ci95=({boot['ci95'][0]:.4f},{boot['ci95'][1]:.4f}) "
        f"n={boot['n_success']}")

    # save bootstrap replicates
    with open(os.path.join(BOOT, "sic_joint_replicates.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["replicate_d_um"])
        for v in boot["replicates"]:
            w.writerow([v])

    # (3c) model-choice perturbation budget (band + resolution/subsample)
    pert = []
    for band in [(400, 4000), (1500, 4000), (600, 4000)]:
        for sub in [2, 3, 5]:
            dd = fit_joint(load_all(band=band, sub=sub))[0].x[0]
            pert.append(float(dd))
    pert = np.array(pert)
    model_choice_spread = float(pert.max() - pert.min())
    log(f"model-choice spread (band x subsample): {model_choice_spread*1000:.1f} nm")

    # ---- write tables ----
    with open(os.path.join(TAB, "sic_results.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["quantity", "value_um", "note"])
        w.writerow(["d_independent_10deg", round(d10, 4), "invariance diagnostic"])
        w.writerow(["d_independent_15deg", round(d15, 4), "invariance diagnostic"])
        w.writerow(["bootstrap_10deg_ci95_lo", round(independent_boot[10.0]["ci95"][0], 4), "wild bootstrap, seed=%d" % SEED])
        w.writerow(["bootstrap_10deg_ci95_hi", round(independent_boot[10.0]["ci95"][1], 4), "wild bootstrap, seed=%d" % SEED])
        w.writerow(["bootstrap_15deg_ci95_lo", round(independent_boot[15.0]["ci95"][0], 4), "wild bootstrap, seed=%d" % SEED])
        w.writerow(["bootstrap_15deg_ci95_hi", round(independent_boot[15.0]["ci95"][1], 4), "wild bootstrap, seed=%d" % SEED])
        w.writerow(["cross_angle_diff", round(cross_angle_diff, 4), "|d10-d15| (consistency)"])
        w.writerow(["d_joint_reported", round(d_joint, 4), "SHARED-d joint fit (reported thickness)"])
        w.writerow(["jacobian_sigma_d", round(sigma_d, 4), "conditional std uncertainty"])
        w.writerow(["bootstrap_ci95_lo", round(boot["ci95"][0], 4), "wild bootstrap, seed=%d" % SEED])
        w.writerow(["bootstrap_ci95_hi", round(boot["ci95"][1], 4), "wild bootstrap, seed=%d" % SEED])

    with open(os.path.join(TAB, "sic_uncertainty_budget.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["component", "type", "value_um", "note"])
        w.writerow(["fitting (Jacobian conditional std)", "A/conditional", round(sigma_d, 4), "cov=s^2(JTJ)^-1"])
        w.writerow(["fitting (bootstrap half-width)", "A/conditional",
                    round((boot["ci95"][1] - boot["ci95"][0]) / 2, 4), "wild bootstrap 95%"])
        w.writerow(["model-choice (band+resolution)", "B/sensitivity", round(model_choice_spread, 4),
                    "range over band x subsample grid"])
        w.writerow(["cross-angle consistency (half-range)", "diagnostic", round(cross_angle_diff / 2, 4),
                    "NOT added; invariance check"])
        w.writerow(["optical-model-form / instrument / sample", "UNASSESSED", "NA",
                    "no independent reference; not a full GUM budget"])

    # ---- figure ----
    import matplotlib
    matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(boot["replicates"], bins=40, color="C0", alpha=0.8)
    ax.axvline(d_joint, c="C3", lw=1.5, label=f"joint d={d_joint:.3f} µm")
    ax.axvline(boot["ci95"][0], c="0.4", ls="--", lw=1)
    ax.axvline(boot["ci95"][1], c="0.4", ls="--", lw=1, label="95% CI")
    ax.set_xlabel("SiC epilayer thickness (µm)"); ax.set_ylabel("bootstrap count")
    ax.set_title("SiC joint-fit thickness: wild bootstrap (seed %d, n=%d)" % (SEED, boot["n_success"]))
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "sic_bootstrap.png"), dpi=200)
    fig.savefig(os.path.join(FIG, "sic_bootstrap.pdf")); plt.close(fig)

    update_manifest(REPO, [
        dict(claim_id="SIC_JOINT_D", description="SiC joint multi-angle reported thickness (shared d)",
             value=f"{d_joint:.4f} um", input_data="data/raw/附件1.xlsx, 附件2.xlsx",
             script="scripts/sic_joint_uncertainty.py", output_file="outputs/tables/sic_results.csv", seed=SEED),
        dict(claim_id="SIC_JOINT_SIGMA", description="SiC joint Jacobian conditional std uncertainty",
             value=f"{sigma_d:.4f} um", input_data="data/raw/附件1.xlsx, 附件2.xlsx",
             script="scripts/sic_joint_uncertainty.py", output_file="outputs/tables/sic_results.csv", seed=SEED),
        dict(claim_id="SIC_JOINT_BOOT_CI", description="SiC joint wild-bootstrap 95% CI",
             value=f"({boot['ci95'][0]:.4f},{boot['ci95'][1]:.4f}) um",
             input_data="data/raw/附件1.xlsx, 附件2.xlsx", script="scripts/sic_joint_uncertainty.py",
             output_file="outputs/bootstrap/sic_joint_replicates.csv", seed=SEED),
        dict(claim_id="SIC_10DEG_BOOT_CI", description="SiC independent 10deg wild-bootstrap 95% CI",
             value=f"({independent_boot[10.0]['ci95'][0]:.4f},{independent_boot[10.0]['ci95'][1]:.4f}) um",
             input_data="data/raw/附件1.xlsx", script="scripts/sic_joint_uncertainty.py",
             output_file="outputs/bootstrap/sic_10deg_replicates.csv", seed=SEED),
        dict(claim_id="SIC_15DEG_BOOT_CI", description="SiC independent 15deg wild-bootstrap 95% CI",
             value=f"({independent_boot[15.0]['ci95'][0]:.4f},{independent_boot[15.0]['ci95'][1]:.4f}) um",
             input_data="data/raw/附件2.xlsx", script="scripts/sic_joint_uncertainty.py",
             output_file="outputs/bootstrap/sic_15deg_replicates.csv", seed=SEED),
        dict(claim_id="SIC_CROSS_ANGLE_DIFF", description="SiC independent-angle |d10-d15|",
             value=f"{cross_angle_diff:.4f} um", input_data="data/raw/附件1.xlsx, 附件2.xlsx",
             script="scripts/sic_joint_uncertainty.py", output_file="outputs/tables/sic_results.csv", seed=SEED),
        dict(claim_id="SIC_MODEL_CHOICE_SPREAD", description="SiC band/resolution sensitivity",
             value=f"{model_choice_spread:.4f} um", input_data="data/raw/附件1.xlsx, 附件2.xlsx",
             script="scripts/sic_joint_uncertainty.py", output_file="outputs/tables/sic_uncertainty_budget.csv", seed=SEED),
    ])
    log("Study B complete.")


if __name__ == "__main__":
    main()
