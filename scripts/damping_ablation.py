#!/usr/bin/env python3
"""Study C — Si empirical high-frequency correction audit (R1.4, R2.3).

The canonical reference implementation subtracts a real-valued saturating-log term from Re(eps), although the
revision documents called it Im(eps) damping. We preserve the implemented real correction, add an explicit
imaginary-placement sensitivity case, and compare both against simpler reflectance-level alternatives.

  1. Fit the Si stack to the real 附件3/附件4 data under six high-frequency mechanisms on the SAME data/
     objective: none, empirical_real_correction (faithful reference), empirical_imag_sensitivity (audit only),
     angle_avg, thick_avg, and roughness (reflectance-level alternatives).
  2. For each: recovered epilayer thickness d, MSE, AICc/BIC, parameter count, held-out-band prediction
     (fit on 1500-3000, predict 3000-4000), Jacobian condition number, and the FULL parameter correlation
     matrix (not just 6 corr-with-d).
  3. Decision rule (locked): adopt a physical alternative for the *reported* model iff it is reproducible
     and >= comparable on held-out fit + thickness stability + identifiability + parsimony. Here the
     thickness-jitter fit reaches its upper bound, so it remains a sensitivity comparator. The conservative
     chosen reporting model is the no-empirical-correction baseline.

Outputs: outputs/tables/damping_ablation.csv, outputs/tables/damping_corr_matrix.csv,
         outputs/tables/damping_corr_matrix_imag_sensitivity.csv, outputs/tables/si_results.csv,
         outputs/figures/damping_residuals.png, outputs/figures/damping_corr_heatmap.png, manifest rows.
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
from src.optics.si_stack import si_reflectance
from src.uncertainty.jacobian import jacobian_conditional_cov, condition_number, parameter_correlations
from src.reference_params import SI_STACK as S
from src.provenance import update_manifest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(REPO, "outputs", "tables")
FIG = os.path.join(REPO, "outputs", "figures")

MECHANISMS = [
    "none",
    "empirical_real_correction",
    "empirical_imag_sensitivity",
    "angle_avg",
    "thick_avg",
    "roughness",
]
# free parameters common to all: d, t_ox, N_e, mu_e, N_s, mu_s + the mechanism's own knob
LORENTZ = (S["S1"], S["w01"], S["g1"], S["S2"], S["w02"], S["g2"])


def load_si(band=(700, 4000), sub=3):
    out = {}
    for key, ang in [("Si_10deg", 10.0), ("Si_15deg", 15.0)]:
        s, R = load_spectrum(key)
        m = (s >= band[0]) & (s <= band[1])
        out[ang] = (s[m][::sub], R[m][::sub])
    return out


def build_forward(mechanism):
    """Return (forward(params, sigma, angle), x0, lb, ub, pnames) for a mechanism."""
    # base params: d, t_ox, log10N_e, mu_e, log10N_s, mu_s
    base_x0 = [S["d_um"], S["t_ox_nm"], np.log10(S["N_epi"]), S["mu_epi"], np.log10(S["N_sub"]), S["mu_sub"]]
    base_lb = [S["d_um"] * 0.6, 0.0, 16.0, 50.0, 18.5, 10.0]
    base_ub = [S["d_um"] * 1.6, 12.0, 18.5, 800.0, 20.5, 400.0]
    base_names = ["d_um", "t_ox_nm", "log10N_e", "mu_e", "log10N_s", "mu_s"]
    A_disp = S["A_disp"]

    if mechanism == "none":
        knob0, knoblb, knobub, knobname = [], [], [], []
    elif mechanism in {"empirical_real_correction", "empirical_imag_sensitivity"}:
        knob0, knoblb, knobub, knobname = [S["B_damp"], S["B_width"]], [0.0, 100.0], [0.8, 2000.0], ["B_damp", "B_width"]
    elif mechanism == "angle_avg":
        knob0, knoblb, knobub, knobname = [0.40], [0.05], [2.0], ["angle_sigma_deg"]
    elif mechanism == "thick_avg":
        knob0, knoblb, knobub, knobname = [0.05], [0.0], [0.3], ["d_jitter_um"]
    elif mechanism == "roughness":
        knob0, knoblb, knobub, knobname = [2.0], [0.0], [20.0], ["roughness_nm"]
    else:
        raise ValueError(mechanism)

    x0 = np.array(base_x0 + knob0, float)
    lb = np.array(base_lb + knoblb, float)
    ub = np.array(base_ub + knobub, float)
    names = base_names + knobname

    def forward(params, sigma, angle):
        d, tox, lNe, mue, lNs, mus = params[:6]
        knob = params[6:]
        kw = dict(mechanism=mechanism, A_disp=A_disp)
        if mechanism in {"empirical_real_correction", "empirical_imag_sensitivity"}:
            kw.update(B_damp=knob[0], B_width=knob[1])
        elif mechanism == "angle_avg":
            kw.update(angle_sigma_deg=knob[0])
        elif mechanism == "thick_avg":
            kw.update(d_jitter_um=knob[0])
        elif mechanism == "roughness":
            kw.update(roughness_nm=knob[0])
        return si_reflectance(sigma, angle, d, tox, 10**lNe, mue, 10**lNs, mus,
                              S["eps_inf_e"], S["eps_inf_s"], LORENTZ, **kw)

    return forward, x0, lb, ub, names


def fit_mechanism(data, mechanism, band_fit=None):
    forward, x0, lb, ub, names = build_forward(mechanism)
    angles = sorted(data.keys())

    def resid(p):
        parts = []
        for a in angles:
            s, R = data[a]
            if band_fit is not None:
                mfit = (s >= band_fit[0]) & (s <= band_fit[1])
                parts.append(forward(p, s[mfit], a) - R[mfit])
            else:
                parts.append(forward(p, s, a) - R)
        return np.concatenate(parts)

    r = least_squares(resid, np.clip(x0, lb, ub), bounds=(lb, ub), loss="soft_l1", f_scale=0.03, max_nfev=4000)
    return r, forward, names


def metrics(data, forward, p):
    angles = sorted(data.keys())
    R_all, F_all = [], []
    for a in angles:
        s, R = data[a]
        R_all.append(R); F_all.append(forward(p, s, a))
    R_all = np.concatenate(R_all); F_all = np.concatenate(F_all)
    resid = F_all - R_all
    sse = float(resid @ resid); n = len(R_all); k = len(p)
    mse = sse / n
    # AICc/BIC under Gaussian assumption
    aic = n * np.log(sse / n) + 2 * k
    aicc = aic + (2 * k * (k + 1)) / max(n - k - 1, 1)
    bic = n * np.log(sse / n) + k * np.log(n)
    return dict(mse=mse, sse=sse, n=n, k=k, aicc=aicc, bic=bic)


def heldout_band(data, mechanism):
    """Fit on 1500-3000, predict 3000-4000; return held-out MSE."""
    r, forward, names = fit_mechanism(data, mechanism, band_fit=(1500, 3000))
    angles = sorted(data.keys())
    res = []
    for a in angles:
        s, R = data[a]
        m = (s >= 3000) & (s <= 4000)
        if m.sum() > 5:
            res.append(forward(r.x, s[m], a) - R[m])
    if not res:
        return np.nan
    res = np.concatenate(res)
    return float(np.mean(res ** 2))


def main():
    for d in (TAB, FIG):
        os.makedirs(d, exist_ok=True)
    log = lambda m: print(f"[{datetime.datetime.now().isoformat(timespec='seconds')}] {m}", flush=True)
    data = load_si()

    rows = []
    corr_by_empirical_variant = {}
    names_by_empirical_variant = {}
    fit_by_mechanism = {}
    resid_curves = {}
    for mech in MECHANISMS:
        r, forward, names = fit_mechanism(data, mech)
        d_um = float(r.x[0])
        met = metrics(data, forward, r.x)
        cov, sigma = jacobian_conditional_cov(r.jac, r.fun)
        kappa = condition_number(r.jac)
        ho = heldout_band(data, mech)
        fit_by_mechanism[mech] = dict(result=r, forward=forward, names=names, sigma=sigma, metrics=met, heldout=ho)
        knob_value = float(r.x[6]) if len(r.x) > 6 else None
        knob_at_bound = bool(
            len(r.x) > 6 and (
                np.isclose(r.x[6], build_forward(mech)[2][6], rtol=0, atol=1e-5) or
                np.isclose(r.x[6], build_forward(mech)[3][6], rtol=0, atol=1e-5)
            )
        )
        rows.append(dict(
            mechanism=mech, n_params=met["k"], d_um=round(d_um, 4),
            sigma_d_um=round(float(sigma[0]), 4), mse=round(met["mse"], 7),
            aicc=round(met["aicc"], 1), bic=round(met["bic"], 1),
            heldout_mse=round(ho, 7) if np.isfinite(ho) else None,
            kappa_JTJ=f"{kappa:.2e}",
            mechanism_knob=(names[6] if len(names) > 6 else ""),
            knob_value=(round(knob_value, 6) if knob_value is not None else ""),
            knob_at_bound=("yes" if knob_at_bound else "no"),
            correction_placement=(
                "Re(epsilon), faithful to the reference implementation"
                if mech == "empirical_real_correction"
                else "Im(epsilon), sensitivity only"
                if mech == "empirical_imag_sensitivity"
                else "not applicable"
            ),
        ))
        log(f"{mech}: d={d_um:.4f} k={met['k']} mse={met['mse']:.2e} aicc={met['aicc']:.1f} "
            f"heldout_mse={ho:.2e} kappa={kappa:.1e}")
        # residual curve at 10deg
        s10, R10 = data[10.0]
        resid_curves[mech] = (s10, forward(r.x, s10, 10.0) - R10)
        if mech in {"empirical_real_correction", "empirical_imag_sensitivity"}:
            corr_by_empirical_variant[mech] = parameter_correlations(cov)
            names_by_empirical_variant[mech] = names

    # ablation table
    with open(os.path.join(TAB, "damping_ablation.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["mechanism", "n_params", "d_um", "sigma_d_um", "mse",
                                          "aicc", "bic", "heldout_mse", "kappa_JTJ",
                                          "mechanism_knob", "knob_value", "knob_at_bound",
                                          "correction_placement"])
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # Full matrices for both empirical placements. The primary matrix is the
    # faithful Re(epsilon) implementation from the Si reference model.
    matrix_outputs = {
        "empirical_real_correction": "damping_corr_matrix.csv",
        "empirical_imag_sensitivity": "damping_corr_matrix_imag_sensitivity.csv",
    }
    for mechanism, filename in matrix_outputs.items():
        corr = corr_by_empirical_variant[mechanism]
        names = names_by_empirical_variant[mechanism]
        with open(os.path.join(TAB, filename), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow([""] + names)
            for i, nm in enumerate(names):
                w.writerow([nm] + [round(float(corr[i, j]), 3) for j in range(len(names))])

    # Thickness shifts when removing each empirical correction vs baseline.
    d_none = next(r["d_um"] for r in rows if r["mechanism"] == "none")
    d_real = next(r["d_um"] for r in rows if r["mechanism"] == "empirical_real_correction")
    d_imag = next(r["d_um"] for r in rows if r["mechanism"] == "empirical_imag_sensitivity")
    shift_real_nm = abs(d_none - d_real) * 1000
    shift_imag_nm = abs(d_none - d_imag) * 1000
    log(f"thickness shift removing faithful Re(epsilon) correction: {shift_real_nm:.1f} nm")
    log(f"thickness shift removing Im(epsilon) sensitivity correction: {shift_imag_nm:.1f} nm")

    # Conservative final Si table. The thickness-jitter comparator reaches its
    # fitted upper bound, so these spectra do not establish it as a replacement.
    chosen = fit_by_mechanism["none"]
    chosen_result = chosen["result"]
    sensitivity_d = np.array([float(r["d_um"]) for r in rows])
    with open(os.path.join(TAB, "si_results.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["quantity", "value", "unit", "note"])
        w.writerow(["chosen_model", "none", "", "baseline: no empirical HF correction; conservative conditional reporting model"])
        w.writerow(["d_chosen", round(float(chosen_result.x[0]), 4), "um", "conditional on baseline Si stack"])
        w.writerow(["jacobian_sigma_d", round(float(chosen["sigma"][0]), 4), "um", "conditional std uncertainty"])
        w.writerow(["oxide_thickness_fit", round(float(chosen_result.x[1]), 4), "nm", "nuisance estimate reaches upper bound; do not report as validated oxide thickness"])
        w.writerow(["candidate_d_range", round(float(sensitivity_d.max() - sensitivity_d.min()), 4), "um", "range across ablation candidates, including sensitivity-only empirical placements"])
        w.writerow(["faithful_real_correction_removal_shift", round(shift_real_nm, 1), "nm", "the reference implements subtraction from Re(epsilon)"])
        w.writerow(["imag_correction_removal_shift", round(shift_imag_nm, 1), "nm", "sensitivity-only Im(epsilon) placement"])
        w.writerow(["replacement_status", "not established", "", "reflectance-level alternatives remain sensitivity comparators; thickness-jitter knob reaches upper bound"])

    # figures
    import matplotlib
    matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 4.2))
    for mech, (s, res) in resid_curves.items():
        ax.plot(s, res, lw=0.8, alpha=0.8, label=mech)
    ax.axhline(0, c="k", lw=0.6, ls="--"); ax.set_xlabel("σ (cm⁻¹)"); ax.set_ylabel("residual (model−data), 10°")
    ax.set_title("Study C: Si residual morphology by high-frequency mechanism")
    ax.legend(frameon=False, fontsize=8, ncol=5)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "damping_residuals.png"), dpi=200)
    fig.savefig(os.path.join(FIG, "damping_residuals.pdf")); plt.close(fig)

    corr_for_real = corr_by_empirical_variant["empirical_real_correction"]
    names_for_real = names_by_empirical_variant["empirical_real_correction"]
    if corr_for_real is not None:
        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.imshow(corr_for_real, vmin=-1, vmax=1, cmap="RdBu_r")
        ax.set_xticks(range(len(names_for_real))); ax.set_xticklabels(names_for_real, rotation=90, fontsize=7)
        ax.set_yticks(range(len(names_for_real))); ax.set_yticklabels(names_for_real, fontsize=7)
        for i in range(len(names_for_real)):
            for j in range(len(names_for_real)):
                ax.text(j, i, f"{corr_for_real[i,j]:.2f}", ha="center", va="center", fontsize=6)
        fig.colorbar(im, ax=ax, shrink=0.8); ax.set_title("Si empirical Re(epsilon)-correction full correlation")
        fig.tight_layout(); fig.savefig(os.path.join(FIG, "damping_corr_heatmap.png"), dpi=200)
        fig.savefig(os.path.join(FIG, "damping_corr_heatmap.pdf")); plt.close(fig)

    update_manifest(REPO, [
        dict(claim_id="DAMPING_ABLATION",
             description="Si empirical HF-correction audit: faithful Re(epsilon), Im(epsilon) sensitivity, and reflectance-level alternatives",
             value=f"faithful Re(epsilon) correction removal shifts d by {shift_real_nm:.1f} nm",
             input_data="data/raw/附件3.xlsx, 附件4.xlsx", script="scripts/damping_ablation.py",
             output_file="outputs/tables/damping_ablation.csv", seed=None),
        dict(claim_id="DAMPING_CORR_MATRIX",
             description="Si faithful Re(epsilon)-correction FULL parameter correlation matrix",
             value="full matrix saved", input_data="data/raw/附件3.xlsx, 附件4.xlsx",
             script="scripts/damping_ablation.py", output_file="outputs/tables/damping_corr_matrix.csv", seed=None),
        dict(claim_id="DAMPING_CORR_MATRIX_IMAG_SENSITIVITY",
             description="Si Im(epsilon)-correction sensitivity FULL parameter correlation matrix",
             value="full sensitivity matrix saved", input_data="data/raw/附件3.xlsx, 附件4.xlsx",
             script="scripts/damping_ablation.py",
             output_file="outputs/tables/damping_corr_matrix_imag_sensitivity.csv", seed=None),
        dict(claim_id="SI_CHOSEN_MODEL",
             description="Conservative Si chosen-model result without empirical HF correction",
             value=f"{float(chosen_result.x[0]):.4f} um", input_data="data/raw/附件3.xlsx, 附件4.xlsx",
             script="scripts/damping_ablation.py", output_file="outputs/tables/si_results.csv", seed=None),
    ])
    log("Study C complete.")


if __name__ == "__main__":
    main()
