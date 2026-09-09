#!/usr/bin/env python3
"""P0 discrepancy investigation: where does the SiC 0.293 um (3.8%) two-beam->multi-beam shift come from?

We compare, on the real 10°/15° SiC data, a ladder of forward models that progressively add the
ingredients present in the original multi-beam fit (V5_3.py) but absent from the two-beam baseline:

  M0  two-beam, semi-infinite substrate (front epi interface only, truncated)        [baseline]
  M1  Airy,     semi-infinite substrate (front epi interface, full series)           [pure model-form]
  M2  Airy,     + finite-substrate backside reflection (air|film|substrate|air)      [extra interface]
  M3  M2        + angular/thickness Gaussian averaging (instrument/aperture/jitter)   [V5_3-like]

Each model re-fits the SAME free parameters (epi/sub doping, mobilities, phonon dampings, thickness) so
the thickness differences isolate the *ingredient* responsible. If M0->M1 is a few nm but M2/M3 produce
the larger shift, then the reported 0.293 um is a substrate-backside + averaging (model-choice) effect,
NOT the two-beam->multi-beam (epilayer interface) effect the manuscript attributes it to.

Outputs: outputs/tables/sic_discrepancy_ladder.csv + manifest row. Deterministic (DE seeded).
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
from src.optics.forward import reflectance_two_beam, reflectance_airy
from src.optics.dielectric import eps_mdf_sic_from_Nmu, n_from_eps
from src.optics.fresnel import cos_theta_complex, fresnel_rij
from src.optics.constants import m_eff_4H_concentration, m_eff_4H_mobility
from src.reference_params import SIC_4H
from src.provenance import update_manifest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(REPO, "outputs", "tables")
SUBSAMPLE = 3
SUB_THICK_UM = 350.0          # typical SiC wafer thickness (finite substrate backside)
SUB_SIGMA_UM = 5.0            # substrate thickness spread -> backside fringe washout
FILM_SIGMA_NM = 30.0          # film thickness non-uniformity within spot
ANGLE_SIGMA_DEG = 0.40        # aperture angular spread
N_AVG = 5                     # Gaussian nodes per averaged dimension

MC, MM = m_eff_4H_concentration(), m_eff_4H_mobility()
P = SIC_4H


def _eps(sigma, logN, mu, GT, GL):
    return eps_mdf_sic_from_Nmu(sigma, P["eps_inf"], P["wT_cm"], P["wL_cm"], GT, GL, 10**logN, mu, MC, MM)


def reflectance_airy_finite_sub(eps_f, eps_s, d_um, d_sub_um, sigma, theta0_deg, pol):
    """Airy with a finite substrate backside (air|film|substrate|air), recursive amplitudes (V5_3 form)."""
    n0 = np.ones_like(sigma, complex)
    n1 = n_from_eps(eps_f); n2 = n_from_eps(eps_s); n3 = np.ones_like(sigma, complex)
    th0 = np.deg2rad(theta0_deg)
    c0 = np.cos(th0) + 0j
    c1 = cos_theta_complex(n0, n1, th0)
    c2 = cos_theta_complex(n0, n2, th0)
    c3 = c0
    lam = 1.0 / (sigma * 100.0)
    d1 = d_um * 1e-6; d2 = d_sub_um * 1e-6
    delta1 = 2 * np.pi * n1 * d1 * c1 / lam
    delta2 = 2 * np.pi * n2 * d2 * c2 / lam

    def stack(pol_):
        r01 = fresnel_rij(n0, n1, c0, c1, pol_)
        r12 = fresnel_rij(n1, n2, c1, c2, pol_)
        r23 = fresnel_rij(n2, n3, c2, c3, pol_)
        e2 = np.exp(2j * delta2)
        r_eff_12 = (r12 + r23 * e2) / (1 + r12 * r23 * e2)
        e1 = np.exp(2j * delta1)
        r_tot = (r01 + r_eff_12 * e1) / (1 + r01 * r_eff_12 * e1)
        return np.abs(r_tot) ** 2

    return 0.5 * (stack("s") + stack("p"))


def _gauss(mu, sig, n):
    if sig <= 0 or n <= 1:
        return np.array([mu]), np.array([1.0])
    x = np.linspace(mu - 3 * sig, mu + 3 * sig, n)
    w = np.exp(-0.5 * ((x - mu) / sig) ** 2); w /= w.sum()
    return x, w


def reflectance_airy_finite_sub_averaged(eps_f, eps_s, d_um, sigma, theta0_deg, pol):
    dvals, wd = _gauss(d_um, FILM_SIGMA_NM * 1e-3, N_AVG)
    svals, ws = _gauss(SUB_THICK_UM, SUB_SIGMA_UM, N_AVG)
    avals, wa = _gauss(theta0_deg, ANGLE_SIGMA_DEG, N_AVG)
    R = np.zeros_like(sigma, float); W = 0.0
    for dd, wdi in zip(dvals, wd):
        for ss, wsi in zip(svals, ws):
            for aa, wai in zip(avals, wa):
                R += wdi * wsi * wai * reflectance_airy_finite_sub(eps_f, eps_s, dd, ss, sigma, aa, pol)
                W += wdi * wsi * wai
    return R / max(W, 1e-12)


def fit_model(sigma, R, angle, model_name):
    def forward(x):
        logNf, muf, d, logNs, mus, GTf, GLf, GTs, GLs = x
        ef = _eps(sigma, logNf, muf, GTf, GLf)
        es = _eps(sigma, logNs, mus, GTs, GLs)
        if model_name == "M0_two_beam":
            return reflectance_two_beam(ef, es, d * 1e-6, sigma, angle, "avg")
        if model_name == "M1_airy_semiinf":
            return reflectance_airy(ef, es, d * 1e-6, sigma, angle, "avg")
        if model_name == "M2_airy_finite_sub":
            return reflectance_airy_finite_sub(ef, es, d, SUB_THICK_UM, sigma, angle, "avg")
        if model_name == "M3_airy_finite_sub_avg":
            return reflectance_airy_finite_sub_averaged(ef, es, d, sigma, angle, "avg")
        raise ValueError(model_name)

    x0 = np.array([P["log10_N_film"], P["mu_film_cm2_Vs"], P["d_um"], P["log10_N_sub"],
                   P["mu_sub_cm2_Vs"], P["GT_film_cm"], P["GL_film_cm"], P["GT_sub_cm"], P["GL_sub_cm"]])
    lb = np.array([15, 10, 5, 16, 10, 1, 1, 1, 1], float)
    ub = np.array([19.8, 800, 10, 19.8, 500, 4, 15, 4, 15], float)
    r = least_squares(lambda x: forward(x) - R, np.clip(x0, lb, ub), bounds=(lb, ub),
                      loss="soft_l1", f_scale=0.01, max_nfev=3000)
    rms = float(np.sqrt(np.mean((forward(r.x) - R) ** 2)))
    return r.x[2], rms


def main():
    os.makedirs(TAB, exist_ok=True)
    rows = []
    for key, ang in [("SiC_10deg", 10.0), ("SiC_15deg", 15.0)]:
        sigma, R = load_spectrum(key)
        m = (sigma >= 400) & (sigma <= 4000)
        sigma, R = sigma[m][::SUBSAMPLE], R[m][::SUBSAMPLE]
        res = {}
        for model in ["M0_two_beam", "M1_airy_semiinf", "M2_airy_finite_sub", "M3_airy_finite_sub_avg"]:
            d, rms = fit_model(sigma, R, ang, model)
            res[model] = (d, rms)
        rows.append(dict(
            dataset=key, angle_deg=ang,
            d_M0_two_beam=round(res["M0_two_beam"][0], 4),
            d_M1_airy=round(res["M1_airy_semiinf"][0], 4),
            d_M2_finite_sub=round(res["M2_airy_finite_sub"][0], 4),
            d_M3_finite_sub_avg=round(res["M3_airy_finite_sub_avg"][0], 4),
            shift_M0_to_M1_nm=round((res["M1_airy_semiinf"][0] - res["M0_two_beam"][0]) * 1000, 1),
            shift_M0_to_M3_nm=round((res["M3_airy_finite_sub_avg"][0] - res["M0_two_beam"][0]) * 1000, 1),
            rms_M0=round(res["M0_two_beam"][1], 5), rms_M3=round(res["M3_airy_finite_sub_avg"][1], 5),
        ))
    with open(os.path.join(TAB, "sic_discrepancy_ladder.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"[{datetime.datetime.now().isoformat(timespec='seconds')}] SiC model-ladder discrepancy:")
    for r in rows:
        print("  ", r)

    update_manifest(REPO, [dict(
        claim_id="SIC_DISCREPANCY_LADDER",
        description="Model ladder isolating the origin of the SiC two-beam->multi-beam thickness shift",
        value="see CSV (M0->M1 pure model-form vs M0->M3 with finite-substrate+averaging)",
        input_data="data/raw/SiC_10deg_reflectance.xlsx, SiC_15deg_reflectance.xlsx",
        script="scripts/sic_discrepancy_investigation.py",
        output_file="outputs/tables/sic_discrepancy_ladder.csv",
        seed=None,
    )])


if __name__ == "__main__":
    main()
