# Provenance Summary

Human-readable companion to the machine-readable record `outputs/manifest.json`. Each reported quantity is
traced to its input data, generating script, stored output, and random seed. Every row was independently
reviewed (`qa_status = QA_APPROVED`) at release time. Re-running the script in `scripts/` regenerates the
value from the data in `data/raw/`.

The shipped machine-readable rows preserve producing commit ids from the private audited development
lineage. That private history is intentionally excluded from this sanitized orphan release. A reproduction
run records the public checkout commit (or `unavailable-exported-source` for an archive export) and resets
the regenerated producer rows to `PENDING`.

| Claim | Value | Input data | Script | Output file | Seed | QA |
|---|---|---|---|---|---|---|
| DAMPING_ABLATION | faithful Re(epsilon) correction removal shifts d by 45.3 nm | data/raw/Si_10deg_reflectance.xlsx, Si_15deg_reflectance.xlsx | `scripts/damping_ablation.py` | `outputs/tables/damping_ablation.csv` | — | QA_APPROVED |
| DAMPING_CORR_MATRIX | full matrix saved | data/raw/Si_10deg_reflectance.xlsx, Si_15deg_reflectance.xlsx | `scripts/damping_ablation.py` | `outputs/tables/damping_corr_matrix.csv` | — | QA_APPROVED |
| DAMPING_CORR_MATRIX_IMAG_SENSITIVITY | full sensitivity matrix saved | data/raw/Si_10deg_reflectance.xlsx, Si_15deg_reflectance.xlsx | `scripts/damping_ablation.py` | `outputs/tables/damping_corr_matrix_imag_sensitivity.csv` | — | QA_APPROVED |
| SI_CHOSEN_MODEL | 4.3048 um | data/raw/Si_10deg_reflectance.xlsx, Si_15deg_reflectance.xlsx | `scripts/damping_ablation.py` | `outputs/tables/si_results.csv` | — | QA_APPROVED |
| FINESSE_TABLE4_REALDATA | see CSV | data/raw/SiC_10deg_reflectance.xlsx, SiC_15deg_reflectance.xlsx, Si_10deg_reflectance.xlsx, Si_15deg_reflectance.xlsx | `scripts/recompute_finesse_tables.py` | `outputs/tables/finesse_table4.csv` | — | QA_APPROVED |
| FINESSE_TABLEB7_THEORY | report explicit q statistics; not directly comparable with measured F_line | data/raw/SiC_10deg_reflectance.xlsx, SiC_15deg_reflectance.xlsx, Si_10deg_reflectance.xlsx, Si_15deg_reflectance.xlsx | `scripts/recompute_finesse_tables.py` | `outputs/tables/finesse_tableB7.csv` | — | QA_APPROVED |
| ANGULAR_IDENTIFIABILITY | see CSV; all rows simulated/virtual; 10/15 is the only measured-data angle pair available | synthetic (SiC MDF Airy truth) | `scripts/angular_identifiability_study.py` | `outputs/tables/angular_identifiability.csv` | — | QA_APPROVED |
| SIC_TWO_BEAM_VS_AIRY_REAL | diff 1.41 nm (10deg), 3.21 nm (15deg) | data/raw/SiC_10deg_reflectance.xlsx, SiC_15deg_reflectance.xlsx | `scripts/sic_two_beam_vs_airy_realdata.py` | `outputs/tables/sic_two_beam_vs_airy_realdata.csv` | — | QA_APPROVED |
| SIC_DISCREPANCY_LADDER | see CSV (M0->M1 pure model-form vs M0->M3 with finite-substrate+averaging) | data/raw/SiC_10deg_reflectance.xlsx, SiC_15deg_reflectance.xlsx | `scripts/sic_discrepancy_investigation.py` | `outputs/tables/sic_discrepancy_ladder.csv` | — | QA_APPROVED |
| SIC_JOINT_D | 7.4051 um | data/raw/SiC_10deg_reflectance.xlsx, SiC_15deg_reflectance.xlsx | `scripts/sic_joint_uncertainty.py` | `outputs/tables/sic_results.csv` | 20260601 | QA_APPROVED |
| SIC_JOINT_SIGMA | 0.0068 um | data/raw/SiC_10deg_reflectance.xlsx, SiC_15deg_reflectance.xlsx | `scripts/sic_joint_uncertainty.py` | `outputs/tables/sic_results.csv` | 20260601 | QA_APPROVED |
| SIC_JOINT_BOOT_CI | (7.3898,7.4112) um | data/raw/SiC_10deg_reflectance.xlsx, SiC_15deg_reflectance.xlsx | `scripts/sic_joint_uncertainty.py` | `outputs/bootstrap/sic_joint_replicates.csv` | 20260601 | QA_APPROVED |
| SIC_10DEG_BOOT_CI | (7.4221,7.4571) um | data/raw/SiC_10deg_reflectance.xlsx | `scripts/sic_joint_uncertainty.py` | `outputs/bootstrap/sic_10deg_replicates.csv` | 20260601 | QA_APPROVED |
| SIC_15DEG_BOOT_CI | (7.3667,7.4153) um | data/raw/SiC_15deg_reflectance.xlsx | `scripts/sic_joint_uncertainty.py` | `outputs/bootstrap/sic_15deg_replicates.csv` | 20260601 | QA_APPROVED |
| SIC_CROSS_ANGLE_DIFF | 0.0533 um | data/raw/SiC_10deg_reflectance.xlsx, SiC_15deg_reflectance.xlsx | `scripts/sic_joint_uncertainty.py` | `outputs/tables/sic_results.csv` | 20260601 | QA_APPROVED |
| SIC_MODEL_CHOICE_SPREAD | 0.0596 um | data/raw/SiC_10deg_reflectance.xlsx, SiC_15deg_reflectance.xlsx | `scripts/sic_joint_uncertainty.py` | `outputs/tables/sic_uncertainty_budget.csv` | 20260601 | QA_APPROVED |
| FINESSE_GRID | 12960 cells | synthetic (dispersive Airy films) | `scripts/finesse_threshold_study.py` | `outputs/tables/finesse_threshold_grid.csv` | 20260601 | QA_APPROVED |
| FINESSE_TIER_BIAS | see summary | synthetic | `scripts/finesse_threshold_study.py` | `outputs/tables/finesse_threshold_summary.csv` | 20260601 | QA_APPROVED |

## How to verify a value

1. Find the claim's `script` and `output_file` above (full record in `outputs/manifest.json`).
2. Run the script from the repository root (see `REPRODUCIBILITY.md`).
3. Compare the regenerated `output_file` under `outputs/` with the released value.

Stochastic claims (non-empty `Seed`) use a fixed seed; reproduction is value-reproducible to the precision
reported in the paper (trailing-digit float differences ~1e-6 do not affect any reported digit).

**Conditional results.** All thicknesses are conditional on the adopted reflectance model and estimators;
no independent reference measurement was available, so no absolute-trueness claim is made. Wider-angle
(20°–40°) entries are simulated/virtual design studies, not measurements.
