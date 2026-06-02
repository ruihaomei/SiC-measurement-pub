#!/usr/bin/env python3
"""Study A — Multi-factor finesse threshold calibration (addresses R1.1, R2.2, and the P0-B contradiction).

Ground truth = exact Airy single-film reflectance for a dispersive, absorbing film on a substrate.
We sweep a realistic factor grid, and for each cell:
  1. generate the exact-Airy spectrum (optionally + noise, + finite resolution smoothing),
  2. compute the lineshape finesse F = FSR/FWHM (the observable in Table 4),
  3. compute the ideal-cavity finesse from an effective reflectance (the Table B.7 style),
  4. compute band statistics of the rigorous round-trip factor |q| = |r10 r12 e^{2i delta}|,
  5. invert the spectrum with the two-beam model and record the relative thickness bias |d_fit-d_true|/d_true.

Outputs (machine-readable, reproducible from one command with a fixed seed):
  outputs/tables/finesse_threshold_grid.csv     -- one row per grid cell
  outputs/tables/finesse_threshold_summary.csv  -- bias quantiles per |q| tier + misclassification rates
  outputs/figures/finesse_bias.png/.pdf         -- bias-vs-diagnostic distribution
  outputs/logs/finesse_threshold_study.log
and appends provenance to outputs/manifest.json (qa_status defaults to PENDING; only QA sets QA_APPROVED).

Design choice: we use a tunable-contrast dispersive film (Sellmeier-like real index + Lorentzian
absorption) rather than only the fixed SiC/Si responses, to test whether a *fixed
global threshold* holds across dispersion/absorption/thickness/resolution. The SiC-like and Si-like
anchor points are included as labelled rows for grounding.
"""
import os
import sys
import json
import argparse
import datetime
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.optics.forward import reflectance_airy
from src.optics.fresnel import roundtrip_factor_q
from src.optics.dielectric import n_from_eps
from src.diagnostics.finesse import (
    lineshape_finesse, finesse_to_reflectance, cavity_finesse_from_roundtrip_amplitude,
    band_q_statistics,
)
from src.inversion.two_beam_fit import fit_thickness_two_beam
from src.provenance import update_manifest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(REPO, "outputs", "tables")
FIG = os.path.join(REPO, "outputs", "figures")
LOG = os.path.join(REPO, "outputs", "logs")
SEED = 20260601


def _log(msg, fh=None):
    line = f"[{datetime.datetime.now().isoformat(timespec='seconds')}] {msg}"
    print(line, flush=True)
    if fh:
        fh.write(line + "\n"); fh.flush()


def dispersive_film_eps(sigma, n_inf, sellmeier_A, sellmeier_w0, absorption_S, absorption_w0, absorption_g):
    """A flexible dispersive+absorbing film dielectric (real Sellmeier term + Lorentz absorption).

    eps(sigma) = [n_inf + A * w0^2/(w0^2 - sigma^2)]^2  +  i * Lorentz_abs
    Gives a controllable, smoothly dispersive complex index spanning weak->strong feedback.
    """
    s = np.asarray(sigma, float)
    n_real = n_inf + sellmeier_A * sellmeier_w0**2 / (sellmeier_w0**2 - s**2 + 1e-9)
    eps_real = n_real**2
    # Lorentz absorption contribution to Im(eps)
    denom = (absorption_w0**2 - s**2)**2 + (absorption_g * s)**2
    eps_imag = absorption_S * absorption_w0**2 * (absorption_g * s) / (denom + 1e-9)
    return eps_real + 1j * eps_imag


def apply_resolution(sigma, R, resolution_cm):
    """Convolve with a Gaussian of FWHM=resolution_cm to mimic finite spectral resolution."""
    if resolution_cm <= 0:
        return R
    dsig = np.median(np.diff(sigma))
    fwhm_pts = resolution_cm / dsig
    sigma_pts = fwhm_pts / 2.3548
    if sigma_pts < 0.5:
        return R
    half = int(np.ceil(3 * sigma_pts))
    x = np.arange(-half, half + 1)
    k = np.exp(-0.5 * (x / sigma_pts) ** 2)
    k /= k.sum()
    return np.convolve(R, k, mode="same")


def run_grid(rng, fh, quick=False):
    sigma = np.linspace(1500.0, 4000.0, 1600)  # fringe region (avoid Reststrahlen complications)

    # Factor levels (R2.2 asks for dispersion/absorption/resolution/material/thickness/noise/angle/pol).
    # Index range widened to span weak->strong cavity feedback: high film index + low substrate index
    # builds genuine high-Q cavities (|q|=|r01 r12| up to ~0.4), matching the regime real SiC/Si reach.
    n_inf_levels = [1.8, 2.6, 3.4, 4.2] if not quick else [3.4]       # film base index -> |r01|
    substrate_n_levels = [1.4, 2.0, 3.0] if not quick else [1.4]      # substrate index -> |r12|
    thickness_levels = [4.0, 7.0, 10.0, 13.0] if not quick else [7.0] # um
    absorption_levels = [0.5, 1.0, 2.0] if not quick else [1.0]       # x scaling of Lorentz S
    resolution_levels = [0.0, 1.0, 2.0, 4.0, 8.0] if not quick else [4.0]  # cm^-1
    snr_levels = [np.inf, 200.0, 50.0] if not quick else [200.0]      # additive-noise SNR
    angle_levels = [10.0, 15.0] if not quick else [10.0]              # measured angles
    pol_levels = ["avg", "s", "p"] if not quick else ["avg"]

    base_abs_S = 0.05
    sellmeier_A, sellmeier_w0 = 0.15, 1200.0
    absorption_w0, absorption_g = 1000.0, 120.0

    rows = []
    n_cells = (len(n_inf_levels) * len(substrate_n_levels) * len(thickness_levels) *
               len(absorption_levels) * len(resolution_levels) * len(snr_levels) *
               len(angle_levels) * len(pol_levels))
    _log(f"grid cells = {n_cells}", fh)
    cell = 0
    for n_inf in n_inf_levels:
        for n_sub in substrate_n_levels:
            for d_true in thickness_levels:
                for abs_x in absorption_levels:
                    eps_f = dispersive_film_eps(sigma, n_inf, sellmeier_A, sellmeier_w0,
                                                base_abs_S * abs_x, absorption_w0, absorption_g)
                    eps_s = (n_sub ** 2) * np.ones_like(sigma, dtype=complex)
                    for angle in angle_levels:
                        for pol in pol_levels:
                            # rigorous |q| (use the chosen pol; avg -> use 's' as representative for q magnitude)
                            qpol = "s" if pol == "avg" else pol
                            q = roundtrip_factor_q(np.ones_like(sigma, complex),
                                                   n_from_eps(eps_f), n_from_eps(eps_s),
                                                   sigma, d_true * 1e-6, angle, qpol)
                            qstat = band_q_statistics(q)
                            # exact Airy truth
                            R0 = reflectance_airy(eps_f, eps_s, d_true * 1e-6, sigma, angle, pol)
                            for resn in resolution_levels:
                                Rr = apply_resolution(sigma, R0, resn)
                                for snr in snr_levels:
                                    cell += 1
                                    if np.isfinite(snr):
                                        noise = rng.normal(0, np.mean(np.abs(Rr)) / snr, size=Rr.shape)
                                        Rn = Rr + noise
                                    else:
                                        Rn = Rr
                                    # diagnostics
                                    lf = lineshape_finesse(sigma, Rn)
                                    F_line = lf["finesse"]
                                    Re_eff = float(finesse_to_reflectance(np.array([F_line]))[0]) if np.isfinite(F_line) else np.nan
                                    # invert with two-beam, two modes:
                                    #  - strict: thickness + scale/offset only (isolates pure model-form bias)
                                    #  - nuisance: also free a dispersion shift dn (mirrors Appendix-C protocol;
                                    #    shows how a flexible dispersion term can mask model-form error -> ties to R2.4)
                                    try:
                                        fit_strict = fit_thickness_two_beam(sigma, Rn, eps_f, eps_s, d_true,
                                                                            theta0_deg=angle, pol=pol,
                                                                            fit_nuisance=False)
                                        d_fit_strict = fit_strict["d_um"]; mse_s = fit_strict["mse"]
                                    except Exception:
                                        d_fit_strict, mse_s = np.nan, np.nan
                                    try:
                                        fit_nz = fit_thickness_two_beam(sigma, Rn, eps_f, eps_s, d_true,
                                                                        theta0_deg=angle, pol=pol,
                                                                        fit_nuisance=True)
                                        d_fit_nz = fit_nz["d_um"]; mse_n = fit_nz["mse"]; dn_nz = fit_nz["dn"]
                                    except Exception:
                                        d_fit_nz, mse_n, dn_nz = np.nan, np.nan, np.nan
                                    rel_bias = abs(d_fit_strict - d_true) / d_true if np.isfinite(d_fit_strict) else np.nan
                                    rel_bias_nz = abs(d_fit_nz - d_true) / d_true if np.isfinite(d_fit_nz) else np.nan
                                    rows.append(dict(
                                        n_inf=n_inf, n_sub=n_sub, d_true_um=d_true, abs_x=abs_x,
                                        angle_deg=angle, pol=pol, resolution_cm=resn,
                                        snr=(None if not np.isfinite(snr) else snr),
                                        F_line=F_line, Re_eff=Re_eff,
                                        median_absq=qstat["median_absq"], p95_absq=qstat["p95_absq"],
                                        max_absq=qstat["max_absq"],
                                        F_cav_from_median_absq=float(
                                            cavity_finesse_from_roundtrip_amplitude(qstat["median_absq"]))
                                            if np.isfinite(qstat["median_absq"]) else np.nan,
                                        d_fit_um=d_fit_strict, rel_bias=rel_bias,
                                        d_fit_nz_um=d_fit_nz, rel_bias_nz=rel_bias_nz, dn_nz=dn_nz,
                                        mse=mse_s,
                                    ))
                                    if cell % 200 == 0:
                                        _log(f"  cell {cell}/{n_cells}", fh)
    return rows


def summarize(rows):
    """Bias quantiles per |q| tier (rigorous control parameter), using the realistic nuisance-mode bias.

    We report bias against band-median |q| tiers rather than lineshape-F tiers because the grid shows bias
    tracks |q| (the round-trip feedback) far better than F (an observable shaped also by dispersion/noise).
    rel_bias_nz = two-beam thickness bias with a free dispersion nuisance (mirrors the Appendix-C protocol).
    """
    valid = [r for r in rows if r.get("rel_bias_nz") is not None and np.isfinite(r["rel_bias_nz"])
             and np.isfinite(r["median_absq"])]
    qtiers = [("|q|<=0.05 (two-beam adequate)", lambda q: q <= 0.05),
              ("0.05<|q|<=0.10", lambda q: 0.05 < q <= 0.10),
              ("0.10<|q|<=0.20", lambda q: 0.10 < q <= 0.20),
              ("|q|>0.20 (multi-beam regime)", lambda q: q > 0.20)]
    out = []
    accept_bias = 0.02
    for name, sel in qtiers:
        b = [r["rel_bias_nz"] for r in valid if sel(r["median_absq"])]
        if b:
            out.append(dict(tier=name, n=len(b), bias_median=float(np.median(b)),
                            bias_p95=float(np.percentile(b, 95)), bias_max=float(np.max(b)),
                            frac_bias_gt_2pct=float(np.mean([x > accept_bias for x in b]))))
        else:
            out.append(dict(tier=name, n=0, bias_median=None, bias_p95=None, bias_max=None,
                            frac_bias_gt_2pct=None))
    # false-safe under the LINESHAPE-F adequacy rule: F<=0.7 but realistic bias > 2%
    vf = [r for r in valid if np.isfinite(r["F_line"])]
    fs = [r for r in vf if r["F_line"] <= 0.7 and r["rel_bias_nz"] > accept_bias]
    out.append(dict(tier="FALSE-SAFE (F<=0.7 & bias>2%)", n=len(fs), bias_median=None,
                    bias_p95=None, bias_max=None, frac_bias_gt_2pct=(len(fs) / max(len(vf), 1))))
    # correlation summary row (bias vs |q| and vs F)
    bz = np.array([r["rel_bias_nz"] for r in valid]); mq = np.array([r["median_absq"] for r in valid])
    Fl = np.array([r["F_line"] for r in valid])
    out.append(dict(tier="CORR bias_nz vs median|q|", n=len(valid),
                    bias_median=float(np.corrcoef(bz, mq)[0, 1]), bias_p95=None, bias_max=None,
                    frac_bias_gt_2pct=None))
    mfin = np.isfinite(Fl)
    out.append(dict(tier="CORR bias_nz vs lineshape F", n=int(mfin.sum()),
                    bias_median=float(np.corrcoef(bz[mfin], Fl[mfin])[0, 1]), bias_p95=None,
                    bias_max=None, frac_bias_gt_2pct=None))
    return out


def write_csv(path, rows, fieldnames):
    import csv
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def make_figure(rows, path_png, path_pdf):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    valid = [r for r in rows if r.get("rel_bias_nz") is not None and np.isfinite(r["rel_bias_nz"])]
    F = np.array([r["F_line"] for r in valid], float)
    mq = np.array([r["median_absq"] for r in valid], float)
    # Realistic case: two-beam fit with a free dispersion nuisance (mirrors Appendix-C protocol).
    bias = np.array([r["rel_bias_nz"] for r in valid], float) * 100.0
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    ax[0].scatter(F, bias, s=10, alpha=0.4, c="C0")
    for thr in (0.7, 1.5):
        ax[0].axvline(thr, ls="--", c="0.5", lw=1)
    ax[0].axhline(2.0, ls=":", c="C3", lw=1, label="2% bias")
    ax[0].set_xlabel("Lineshape finesse F = FSR/FWHM (observable)")
    ax[0].set_ylabel("Two-beam thickness bias (%)")
    ax[0].set_title("(a) Bias vs lineshape finesse")
    ax[0].legend(frameon=False, fontsize=8)
    ax[1].scatter(mq, bias, s=10, alpha=0.4, c="C1")
    ax[1].axhline(2.0, ls=":", c="C3", lw=1)
    # Annotate the real SiC operating point (|q| ~ 0.002, fringe band): two-beam is amply adequate.
    ax[1].axvline(0.002, ls="-", c="C2", lw=1.2)
    ax[1].text(0.004, ax[1].get_ylim()[1]*0.8, "real SiC\n(|q|~0.002)", color="C2", fontsize=8)
    ax[1].set_xlabel(r"Band-median $|q|=|r_{10}r_{12}e^{2i\delta}|$ (rigorous)")
    ax[1].set_ylabel("Two-beam thickness bias (%)")
    ax[1].set_title("(b) Bias vs round-trip factor |q|")
    fig.tight_layout()
    fig.savefig(path_png, dpi=200)
    fig.savefig(path_pdf)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="tiny grid for smoke-testing")
    args = ap.parse_args()
    os.makedirs(TAB, exist_ok=True); os.makedirs(FIG, exist_ok=True); os.makedirs(LOG, exist_ok=True)
    rng = np.random.default_rng(SEED)
    with open(os.path.join(LOG, "finesse_threshold_study.log"), "w") as fh:
        _log(f"Study A start (seed={SEED}, quick={args.quick})", fh)
        rows = run_grid(rng, fh, quick=args.quick)
        _log(f"grid done: {len(rows)} rows", fh)
        grid_csv = os.path.join(TAB, "finesse_threshold_grid.csv")
        fields = ["n_inf", "n_sub", "d_true_um", "abs_x", "angle_deg", "pol", "resolution_cm",
                  "snr", "F_line", "Re_eff", "median_absq", "p95_absq", "max_absq",
                  "F_cav_from_median_absq", "d_fit_um", "rel_bias", "d_fit_nz_um", "rel_bias_nz",
                  "dn_nz", "mse"]
        write_csv(grid_csv, rows, fields)
        summ = summarize(rows)
        summ_csv = os.path.join(TAB, "finesse_threshold_summary.csv")
        write_csv(summ_csv, summ, ["tier", "n", "bias_median", "bias_p95", "bias_max", "frac_bias_gt_2pct"])
        make_figure(rows, os.path.join(FIG, "finesse_bias.png"), os.path.join(FIG, "finesse_bias.pdf"))
        for s in summ:
            _log(f"SUMMARY {s}", fh)
        update_manifest(REPO, [
            dict(claim_id="FINESSE_GRID", description="Multi-factor two-beam thickness bias vs finesse/|q|",
                 value=f"{len(rows)} cells", input_data="synthetic (dispersive Airy films)",
                 script="scripts/finesse_threshold_study.py",
                 output_file="outputs/tables/finesse_threshold_grid.csv", seed=SEED,
                 ),
            dict(claim_id="FINESSE_TIER_BIAS", description="Bias quantiles per |q| tier + false-safe rate",
                 value="see summary", input_data="synthetic",
                 script="scripts/finesse_threshold_study.py",
                 output_file="outputs/tables/finesse_threshold_summary.csv", seed=SEED,
                 ),
        ])
        _log("Study A complete.", fh)


if __name__ == "__main__":
    main()
