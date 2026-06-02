#!/usr/bin/env python3
"""Study D — angular-separation identifiability (virtual) for the SiC system (R1.3, R2.4).

Motivation: 10 deg and 15 deg are too close (5 deg) to decouple thickness from refractive-index
dispersion; the Jacobian may stay ill-conditioned. We quantify, as a *virtual/simulated design analysis*
(the available measured-data angles are 10/15 only), how three identifiability metrics improve with angular separation:

  - |Delta d|  : thickness discrepancy when an angle PAIR is fit under a deliberately misspecified model
                 (here: two-beam inversion of an Airy-truth spectrum) -- larger = more sensitive to model error;
  - kappa(J^T W J) : condition number of the thickness+dispersion sensitivity matrix for the angle pair;
  - corr(d, dispersion) : off-diagonal parameter correlation between thickness and a linear dispersion term.

Truth model: SiC MDF system (reference_params). We build synthetic Airy spectra at each angle, then form a
2-parameter local sensitivity model (thickness d, dispersion shift dn) for the *pair* and report metrics
vs separation, for s / p / unpolarised.

IMPORTANT LABELLING: every output row is SIMULATED/VIRTUAL. The 10/15 SiC row is tagged only as matching
the angle pair available in the measured dataset; it is not labelled as an experiment.

Outputs: outputs/tables/angular_identifiability.csv, outputs/figures/angular_identifiability.png, manifest.
"""
import os
import sys
import json
import csv
import datetime
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.io_data import load_spectrum  # noqa (not used directly; kept for parity)
from src.optics.forward import reflectance_airy, reflectance_two_beam
from src.optics.dielectric import eps_mdf_sic_from_Nmu, n_from_eps
from src.optics.constants import m_eff_4H_concentration, m_eff_4H_mobility
from src.uncertainty.jacobian import condition_number, parameter_correlations, jacobian_conditional_cov
from src.inversion.two_beam_fit import fit_thickness_two_beam
from src.reference_params import SIC_4H
from src.provenance import update_manifest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(REPO, "outputs", "tables")
FIG = os.path.join(REPO, "outputs", "figures")
MC, MM = m_eff_4H_concentration(), m_eff_4H_mobility()
P = SIC_4H
SIGMA = np.linspace(1500.0, 4000.0, 1400)   # fringe band for thickness/dispersion identifiability


def sic_eps(sigma, dn=0.0):
    ef = eps_mdf_sic_from_Nmu(sigma, P["eps_inf"], P["wT_cm"], P["wL_cm"], P["GT_film_cm"], P["GL_film_cm"],
                              10**P["log10_N_film"], P["mu_film_cm2_Vs"], MC, MM)
    es = eps_mdf_sic_from_Nmu(sigma, P["eps_inf"], P["wT_cm"], P["wL_cm"], P["GT_sub_cm"], P["GL_sub_cm"],
                              10**P["log10_N_sub"], P["mu_sub_cm2_Vs"], MC, MM)
    if dn != 0.0:
        ef = (np.sqrt(ef) + dn) ** 2
    return ef, es


def strong_eps(sigma, dn=0.0):
    """High-contrast synthetic cavity (film n~3.4 on substrate n~1.5) -> genuine model-form error,
    so the angular-separation test has something to detect (contrast case for R2.4)."""
    nf = 3.4 + 0.05 * (1000.0 / np.maximum(sigma, 1.0)) ** 2  # mild dispersion
    ef = (nf + dn) ** 2 + 0j
    es = (1.5 ** 2) * np.ones_like(sigma, complex)
    return ef, es


def pair_sensitivity(angle_a, angle_b, pol, d0=None, eps_fn=sic_eps):
    """Build the 2-param (d, dn) Jacobian for an angle pair using Airy truth; return kappa, corr(d,dn)."""
    if d0 is None:
        d0 = P["d_um"]

    def model_pair(d, dn):
        efd, esd = eps_fn(SIGMA, dn=dn)
        Ra = reflectance_airy(efd, esd, d * 1e-6, SIGMA, angle_a, pol)
        Rb = reflectance_airy(efd, esd, d * 1e-6, SIGMA, angle_b, pol)
        return np.concatenate([Ra, Rb])

    # finite-difference Jacobian wrt (d, dn)
    base = model_pair(d0, 0.0)
    hd, hn = 1e-3, 1e-4
    Jd = (model_pair(d0 + hd, 0.0) - base) / hd
    Jn = (model_pair(d0, hn) - base) / hn
    J = np.column_stack([Jd, Jn])
    kappa = condition_number(J)
    cov, _ = jacobian_conditional_cov(J, np.zeros(J.shape[0]) + 1e-6)  # residual ~0; structure only
    corr = parameter_correlations(cov)
    return kappa, float(abs(corr[0, 1]))


def pair_delta_d(angle_a, angle_b, pol, eps_fn=sic_eps, d0=None):
    """|Delta d| between the two angles when each is fit with a MISSPECIFIED two-beam model on Airy truth."""
    if d0 is None:
        d0 = P["d_um"]
    ef, es = eps_fn(SIGMA)
    Ra = reflectance_airy(ef, es, d0 * 1e-6, SIGMA, angle_a, pol)
    Rb = reflectance_airy(ef, es, d0 * 1e-6, SIGMA, angle_b, pol)
    fa = fit_thickness_two_beam(SIGMA, Ra, ef, es, d0, theta0_deg=angle_a, pol=pol, fit_nuisance=True)
    fb = fit_thickness_two_beam(SIGMA, Rb, ef, es, d0, theta0_deg=angle_b, pol=pol, fit_nuisance=True)
    return abs(fa["d_um"] - fb["d_um"]) * 1000.0  # nm


def main():
    for d in (TAB, FIG):
        os.makedirs(d, exist_ok=True)
    log = lambda m: print(f"[{datetime.datetime.now().isoformat(timespec='seconds')}] {m}", flush=True)

    base_angle = 10.0
    partner_angles = [15.0, 20.0, 30.0, 40.0]  # all simulated; only 10/15 matches measured-data angles
    systems = [("SiC_MDF", sic_eps, P["d_um"]), ("strong_contrast", strong_eps, 7.0)]
    rows = []
    for sysname, eps_fn, d0 in systems:
        for pol in ["avg", "s", "p"]:
            for b in partner_angles:
                sep = b - base_angle
                kappa, corr_d_dn = pair_sensitivity(base_angle, b, pol, d0=d0, eps_fn=eps_fn)
                dd = pair_delta_d(base_angle, b, pol, eps_fn=eps_fn, d0=d0)
                pair_matches_available_measured_angles = (b == 15.0 and sysname == "SiC_MDF")
                rows.append(dict(
                    system=sysname, pol=pol, angle_a=base_angle, angle_b=b, separation_deg=sep,
                    analysis_kind="simulated_virtual",
                    pair_matches_available_measured_angles=("yes" if pair_matches_available_measured_angles else "no"),
                    delta_d_nm=round(dd, 2), kappa_d_dn=f"{kappa:.3e}", corr_d_dispersion=round(corr_d_dn, 4),
                ))
            log(f"[{sysname}] pol={pol}: " + " ".join(
                f"{base_angle:.0f}/{b:.0f}:|dd|={r['delta_d_nm']}nm,k={r['kappa_d_dn']}"
                for b, r in zip(partner_angles, rows[-len(partner_angles):])))

    with open(os.path.join(TAB, "angular_identifiability.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["system", "pol", "angle_a", "angle_b", "separation_deg",
                                          "analysis_kind", "pair_matches_available_measured_angles",
                                          "delta_d_nm", "kappa_d_dn", "corr_d_dispersion"])
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # figure: metrics vs separation (unpolarised, strong-contrast system where there is error to detect)
    import matplotlib
    matplotlib.use("Agg"); import matplotlib.pyplot as plt
    avg = [r for r in rows if r["pol"] == "avg" and r["system"] == "strong_contrast"]
    seps = [r["separation_deg"] for r in avg]
    dds = [r["delta_d_nm"] for r in avg]
    corrs = [r["corr_d_dispersion"] for r in avg]
    kaps = [float(r["kappa_d_dn"]) for r in avg]
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    ax[0].plot(seps, dds, "o-"); ax[0].set_xlabel("angular separation (°)"); ax[0].set_ylabel("|Δd| (nm)")
    ax[0].set_title("(a) two-beam misfit |Δd| vs separation"); ax[0].axvline(5, ls=":", c="C2")
    ax[1].semilogy(seps, kaps, "s-", c="C1"); ax[1].set_xlabel("angular separation (°)")
    ax[1].set_ylabel(r"$\kappa(J^TWJ)$"); ax[1].set_title("(b) conditioning vs separation"); ax[1].axvline(5, ls=":", c="C2")
    ax[2].plot(seps, corrs, "^-", c="C3"); ax[2].set_xlabel("angular separation (°)")
    ax[2].set_ylabel("|corr(d, dispersion)|"); ax[2].set_title("(c) d–dispersion correlation"); ax[2].axvline(5, ls=":", c="C2")
    for a in ax:
        a.text(5.2, a.get_ylim()[0], "10/15° measured-data angles\\n(curve simulated)",
               color="C2", fontsize=7, rotation=90, va="bottom")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "angular_identifiability.png"), dpi=200)
    fig.savefig(os.path.join(FIG, "angular_identifiability.pdf")); plt.close(fig)

    update_manifest(REPO, [dict(
        claim_id="ANGULAR_IDENTIFIABILITY",
        description="Virtual angular-separation identifiability: |Δd|, kappa(JTWJ), corr(d,dispersion) vs separation",
        value="see CSV; all rows simulated/virtual; 10/15 is the only measured-data angle pair available",
        input_data="synthetic (SiC MDF Airy truth)",
        script="scripts/angular_identifiability_study.py",
        output_file="outputs/tables/angular_identifiability.csv", seed=None,
    )])
    log("Study D complete.")


if __name__ == "__main__":
    main()
