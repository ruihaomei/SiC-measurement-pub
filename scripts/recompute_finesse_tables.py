#!/usr/bin/env python3
"""Recompute the real-data finesse diagnostics behind Table 4 and Appendix Table B.7.

``finesse_table4.csv`` stores measured and fitted-spectrum *lineshape* finesse
``F_line = mean(FSR) / mean(FWHM)`` on the explicit 1500--4000 cm^-1 fringe
band. ``finesse_tableB7.csv`` stores the distinct theoretical round-trip
diagnostic ``|q(sigma, theta)|`` and an equivalent ideal-cavity proxy
``pi*sqrt(median(|q|))/(1-median(|q|))``.

The theoretical proxy is not a measured lineshape quantity. The output marks
Appendix Table B.7 for drop-or-replacement with ``|q|`` statistics rather than
direct comparison to Table 4 thresholds.
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.damping_ablation import fit_mechanism, load_si
from scripts.sic_joint_uncertainty import _eps, _model_one, fit_independent, load_all
from src.diagnostics.finesse import (
    band_q_statistics,
    cavity_finesse_from_roundtrip_amplitude,
    lineshape_finesse,
)
from src.io_data import load_spectrum
from src.optics.dielectric import (
    eps_drude_si,
    eps_lorentz_sio2,
    n_from_eps,
    si_extra_dispersion_real,
)
from src.optics.fresnel import cos_theta_complex, fresnel_rij, phase_thickness, roundtrip_factor_q
from src.provenance import update_manifest
from src.reference_params import SI_STACK as S

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(REPO, "outputs", "tables")
FRINGE_BAND = (1500.0, 4000.0)
LORENTZ = (S["S1"], S["w01"], S["g1"], S["S2"], S["w02"], S["g2"])


def _band(sigma):
    return (sigma >= FRINGE_BAND[0]) & (sigma <= FRINGE_BAND[1])


def _finesse_row(dataset, material, angle, sigma, measured, fitted, model_source):
    mask = _band(sigma)
    measured_metrics = lineshape_finesse(sigma[mask], measured[mask])
    fitted_metrics = lineshape_finesse(sigma[mask], fitted[mask])
    return dict(
        dataset=dataset,
        material=material,
        angle_deg=angle,
        diagnostic_kind="lineshape_observable_FSR_over_FWHM",
        band_cm1="1500-4000",
        measured_F_line=round(float(measured_metrics["finesse"]), 6),
        measured_FSR_cm1=round(float(measured_metrics["FSR"]), 6),
        measured_FWHM_cm1=round(float(measured_metrics["FWHM"]), 6),
        measured_n_peaks=int(measured_metrics["n_peaks"]),
        model_F_line=round(float(fitted_metrics["finesse"]), 6),
        model_FSR_cm1=round(float(fitted_metrics["FSR"]), 6),
        model_FWHM_cm1=round(float(fitted_metrics["FWHM"]), 6),
        model_n_peaks=int(fitted_metrics["n_peaks"]),
        model_source=model_source,
    )


def _q_summary(dataset, material, angle, q, model_source):
    stats = band_q_statistics(q)
    return dict(
        dataset=dataset,
        material=material,
        angle_deg=angle,
        diagnostic_kind="theoretical_roundtrip_proxy_not_lineshape",
        band_cm1="1500-4000",
        q_definition="|r_upper*r_lower*exp(2i*delta_film)|",
        median_absq=round(stats["median_absq"], 8),
        p95_absq=round(stats["p95_absq"], 8),
        max_absq=round(stats["max_absq"], 8),
        ideal_cavity_proxy_F_from_median_absq=round(
            float(cavity_finesse_from_roundtrip_amplitude(stats["median_absq"])), 8),
        model_source=model_source,
        tableB7_action="report explicit q statistics; do not compare directly with measured F_line thresholds",
    )


def _sic_rows():
    lineshape_rows, q_rows = [], []
    fit_data = load_all()
    for dataset, angle in [("SiC_10deg", 10.0), ("SiC_15deg", 15.0)]:
        fit = fit_independent(*fit_data[angle], angle)
        sigma, measured = load_spectrum(dataset)
        mask = _band(sigma)
        sigma, measured = sigma[mask], measured[mask]
        d_um, phys = float(fit.x[0]), fit.x[1:9]
        fitted = _model_one(sigma, angle, d_um, phys)
        source = "independent-angle SiC MDF Airy fit"
        lineshape_rows.append(_finesse_row(dataset, "SiC", angle, sigma, measured, fitted, source))
        ef = _eps(sigma, phys[0], phys[1], phys[4], phys[5])
        es = _eps(sigma, phys[2], phys[3], phys[6], phys[7])
        n0 = np.ones_like(sigma, complex)
        qs = roundtrip_factor_q(n0, n_from_eps(ef), n_from_eps(es), sigma, d_um * 1e-6, angle, "s")
        qp = roundtrip_factor_q(n0, n_from_eps(ef), n_from_eps(es), sigma, d_um * 1e-6, angle, "p")
        q_rows.append(_q_summary(dataset, "SiC", angle, np.concatenate([qs, qp]), source))
    return lineshape_rows, q_rows


def _si_roundtrip_q(sigma, angle, params):
    """Epilayer round-trip factor including the thin oxide in the upper reflection."""
    d_um, tox_nm, logNe, mue, logNs, mus = params[:6]
    eps_ox = eps_lorentz_sio2(sigma, 2.1, *LORENTZ)
    eps_epi = (
        eps_drude_si(sigma, 10**logNe, mue, S["eps_inf_e"])
        + si_extra_dispersion_real(sigma, S["A_disp"])
    )
    eps_sub = eps_drude_si(sigma, 10**logNs, mus, S["eps_inf_s"])
    n_air = np.ones_like(sigma, complex)
    n_ox, n_epi, n_sub = map(n_from_eps, [eps_ox, eps_epi, eps_sub])
    theta = np.deg2rad(angle)
    c_air = np.cos(theta) + 0j
    c_ox = cos_theta_complex(n_air, n_ox, theta)
    c_epi = cos_theta_complex(n_air, n_epi, theta)
    c_sub = cos_theta_complex(n_air, n_sub, theta)

    def one(pol):
        r_epi_ox = fresnel_rij(n_epi, n_ox, c_epi, c_ox, pol)
        r_ox_air = fresnel_rij(n_ox, n_air, c_ox, c_air, pol)
        delta_ox = phase_thickness(n_ox, tox_nm * 1e-9, sigma, c_ox)
        e2_ox = np.exp(2j * delta_ox)
        r_upper = (r_epi_ox + r_ox_air * e2_ox) / (1.0 + r_epi_ox * r_ox_air * e2_ox)
        r_lower = fresnel_rij(n_epi, n_sub, c_epi, c_sub, pol)
        delta_epi = phase_thickness(n_epi, d_um * 1e-6, sigma, c_epi)
        return r_upper * r_lower * np.exp(2j * delta_epi)

    return np.concatenate([one("s"), one("p")])


def _si_rows():
    lineshape_rows, q_rows = [], []
    data = load_si()
    fit, forward, _ = fit_mechanism(data, "none")
    source = "joint-angle Si baseline Abeles fit without empirical HF correction"
    for dataset, angle in [("Si_10deg", 10.0), ("Si_15deg", 15.0)]:
        sigma, measured = load_spectrum(dataset)
        mask = _band(sigma)
        sigma, measured = sigma[mask], measured[mask]
        fitted = forward(fit.x, sigma, angle)
        lineshape_rows.append(_finesse_row(dataset, "Si", angle, sigma, measured, fitted, source))
        q_rows.append(_q_summary(dataset, "Si", angle, _si_roundtrip_q(sigma, angle, fit.x), source))
    return lineshape_rows, q_rows


def _write(path, rows):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    os.makedirs(TAB, exist_ok=True)
    sic_table4, sic_b7 = _sic_rows()
    si_table4, si_b7 = _si_rows()
    table4 = sic_table4 + si_table4
    tableb7 = sic_b7 + si_b7
    _write(os.path.join(TAB, "finesse_table4.csv"), table4)
    _write(os.path.join(TAB, "finesse_tableB7.csv"), tableb7)
    update_manifest(REPO, [
        dict(
            claim_id="FINESSE_TABLE4_REALDATA",
            description="Measured and fitted-spectrum lineshape finesse FSR/FWHM on the explicit fringe band",
            value="see CSV", input_data="data/raw/SiC_10deg_reflectance.xlsx, SiC_15deg_reflectance.xlsx, Si_10deg_reflectance.xlsx, Si_15deg_reflectance.xlsx",
            script="scripts/recompute_finesse_tables.py", output_file="outputs/tables/finesse_table4.csv",
            seed=None,
        ),
        dict(
            claim_id="FINESSE_TABLEB7_THEORY",
            description="Theoretical round-trip |q| statistics and distinct equivalent ideal-cavity proxy",
            value="report explicit q statistics; not directly comparable with measured F_line",
            script="scripts/recompute_finesse_tables.py", output_file="outputs/tables/finesse_tableB7.csv",
            seed=None,
        ),
    ])
    print("Table 4 lineshape and Table B.7 theoretical round-trip diagnostics recomputed.")


if __name__ == "__main__":
    main()
