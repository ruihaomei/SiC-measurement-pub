# Model-adequacy–guided infrared reflectance thickness metrology

Reproducibility package for the manuscript:

> **"Model-adequacy–guided infrared reflectance thickness metrology: diagnosing two-beam breakdown and
> robust multi-beam inversion with uncertainty"** — *Measurement* (Elsevier), manuscript **MEAS-D-26-04372**
> (under review).

This repository contains the code, measured data, numerical outputs, and tests for the reported model-
adequacy studies and headline thickness results. It is designed so that a new researcher can regenerate the
released tables and figures **without contacting the authors**.

> ⚠️ **Scope of claims.** All thickness results are *conditional on the adopted reflectance model and
> estimators*. No independent reference measurement (SEM/TEM/ellipsometry) was available, so the package
> reports **internal consistency, model adequacy, and conditional uncertainty** — not absolute trueness.
> Wider-angle (20°–40°) analyses are **fully simulated/virtual design studies**, never additional
> experiments. The only measured incidence angles are **10° and 15°**.

---

## 1. Scientific background

Infrared (IR) reflectance interferometry infers the thickness of an epitaxial film from the fringe pattern
in its specular reflectance spectrum. A common simplification is the **two-beam** (single round-trip)
model, whereas the rigorous electromagnetic baseline is the **Airy / transfer-matrix method (TMM)**, which
sums all internal reflections. When internal feedback is non-negligible, the two-beam simplification breaks
down and biases the extracted thickness.

This work casts thickness extraction as a **model-adequacy problem**: it diagnoses *when* the two-beam model
is adequate using the wavenumber-resolved round-trip feedback factor `q(σ,θ,pol) = r₁₀ r₁₂ exp(2iβ₁)` and the
measured lineshape finesse `F = FSR / FWHM`, then performs robust multi-beam inversion with conditional
uncertainty for two materials (SiC and Si epilayers).

## 2. Research question

1. **When** does the two-beam model give an adequate thickness, and how can adequacy be diagnosed from the
   spectrum itself (`|q|` tiers and finesse), rather than from engineering rules of thumb?
2. **How large** is the two-beam-vs-Airy thickness bias on real spectra, under a like-for-like comparison?
3. **What conditional uncertainty** can be reported for the multi-beam (Airy) thickness, comparably for SiC
   and Si, in the absence of an independent reference?
4. **How sensitive** is angular-invariance diagnosis to angular separation (a simulated design study)?

## 3. Repository structure

```
src/            Reproduction library (importable package)
  optics/       Fresnel, dielectric models, Si stack, forward Airy/two-beam, constants
  diagnostics/  Lineshape finesse F = FSR/FWHM
  inversion/    Two-beam / Airy fitting
  uncertainty/  Jacobian conditional uncertainty, wild bootstrap
  io_data.py    Loads the four spectra from data/raw/
  provenance.py Run-provenance + manifest helpers
  reference_params.py  SiC-like / Si-like "truth" systems for synthetic studies
scripts/        Reproduction entry points (see §6)
tests/          pytest suite (limits, units, finesse consistency, manifest/CSV checks, release hygiene)
data/raw/       Measured spectra 附件1–4.xlsx (SiC 10°/15°, Si 10°/15°)  [CC BY 4.0]
outputs/        Generated tables/, figures/, bootstrap/, logs/, and manifest.json
docs/provenance/ Human-readable provenance notes; outputs/manifest.json is the machine-readable record
```

A detailed map is in [`REPOSITORY_STRUCTURE.md`](REPOSITORY_STRUCTURE.md).

## 4. Installation

Requires **Python ≥ 3.10**.

```bash
git clone https://github.com/ruihaomei/SiC-measurement.git
cd SiC-measurement
python -m venv .venv && source .venv/bin/activate    # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Core dependencies: `numpy`, `scipy`, `pandas`, `matplotlib`, `openpyxl` (to read the `.xlsx` spectra),
`pytest`. No GPU or deep-learning framework is required. After dependencies are installed, the numerical
workflow runs offline.

## 5. Quick start

```bash
pytest tests/                              # limits, units, finesse, provenance, release hygiene
python scripts/reproduce_main_tables.py    # regenerate the main tables
python scripts/reproduce_main_figures.py   # regenerate the main figures
```

Outputs are written under `outputs/`. Every run records its random seed and provenance.
The test suite finishes in seconds. The full table and figure workflows perform seeded optimization and
bootstrap refits and can take substantial CPU time; runtime depends on the machine and numerical libraries.

## 6. Full reproduction workflow

Run from the repository root with the virtual environment active. Each script regenerates its outputs under
`outputs/` and records provenance into `outputs/manifest.json`.

| Step | Command | Reproduces |
|---|---|---|
| Tests | `pytest tests/` | Physical limits, unit handling, finesse consistency, manifest/CSV checks, release hygiene |
| Tables | `python scripts/reproduce_main_tables.py` | Main-text tables (finesse tiers, etc.) |
| Figures | `python scripts/reproduce_main_figures.py` | Main-text figures |
| Study A | `python scripts/finesse_threshold_study.py` | `|q|`-tier adequacy calibration grid + summary |
| Study B | `python scripts/sic_joint_uncertainty.py` | SiC joint Airy thickness + Jacobian/bootstrap conditional uncertainty |
| Study C | `python scripts/damping_ablation.py` | Si empirical high-frequency correction ablation (Re(ε) canonical / Im(ε) sensitivity) |
| Study D | `python scripts/angular_identifiability_study.py` | Simulated/virtual angular-separation identifiability |

Supporting scripts: `recompute_finesse_tables.py`, `sic_discrepancy_investigation.py`,
`sic_two_beam_vs_airy_realdata.py`.

**Determinism note.** Stochastic steps (bootstrap, differential evolution) are seeded (see
`outputs/manifest.json`). Reproduction is *value-reproducible* to the precision reported in the paper;
trailing-digit differences (~1e-6 relative) can arise from BLAS/threading and do not affect any reported
digit. The reproduction scripts overwrite files under `outputs/` and reset recomputed manifest rows to a
`PENDING` review status — this is expected. To preserve the shipped independently reviewed manifest, run
reproductions on a copy or restore with `git restore outputs/` afterward.

## 7. Data description

`data/raw/` holds four measured unpolarized specular IR reflectance spectra (Thermo Scientific Nicolet iS50
FTIR with a variable-angle specular reflectance accessory; gold-mirror reference; 4 cm⁻¹ resolution,
128 co-added scans):

| File | Material | Incidence angle |
|---|---|---|
| `附件1.xlsx` | SiC epilayer | 10° |
| `附件2.xlsx` | SiC epilayer | 15° |
| `附件3.xlsx` | Si epilayer | 10° |
| `附件4.xlsx` | Si epilayer | 15° |

Each file has two columns: wavenumber `波数 (cm-1)` and reflectance `反射率 (%)`. Reflectance is stored in
percent and converted to a fraction for analysis (`src/io_data.py`). The SiC 15° spectrum contains a few
values above 100% (a calibration/reference artifact); these are **retained, not silently clamped** — see
`src/io_data.py` and the manuscript's preprocessing statement. Non-scientific Excel source-path and author
metadata were removed with `scripts/sanitize_xlsx_metadata.py`; worksheet cell values were preserved.
Checksums for the released workbooks are recorded in `data/raw/SHA256SUMS`. Data are licensed **CC BY 4.0**.

## 8. Figure generation

`python scripts/reproduce_main_figures.py` regenerates the main figures into `outputs/figures/` (both
`.png` and `.pdf`). Individual studies (A–D) also emit their own diagnostic figures. Figures derived from
the simulated angular study are labelled as virtual design guidance.

## 9. Table generation

`python scripts/reproduce_main_tables.py` regenerates the main tables into `outputs/tables/` as CSV. The
`tests/test_manifest_consistency.py` checks selected headline manifest values against the stored CSVs and
bootstrap replicates. The article itself is not included in this repository.

## 10. Testing

```bash
pytest tests/
```

The suite covers: optical-limit sanity (`test_optics_limits.py`), Si-stack forward model
(`test_si_stack.py`), uncertainty machinery (`test_uncertainty.py`), manifest/CSV consistency
(`test_manifest_consistency.py`), and public-release hygiene (`test_public_release_hygiene.py`). All tests
run offline in seconds.

## 11. Limitations

- **Conditional, not absolute.** No independent reference (SEM/TEM/ellipsometry) was available; results are
  conditional on the adopted reflectance model and estimators. No trueness/accuracy claim is made.
- **Two measured angles only (10°, 15°).** All wider-angle (20°–40°) results are **fully simulated virtual
  design studies**, not experiments.
- **Empirical Si high-frequency correction** is reported as **sensitivity-only**; the canonical placement is
  in `Re(ε)`, and the reported Si thickness uses the no-correction baseline. No Kramers–Kronig-consistent
  reporting model is claimed.
- **Value-reproducible, not bit-identical** across platforms/BLAS (see §6).

## 12. Citation

If you use this code or data, please cite the manuscript. Machine-readable metadata is in
[`CITATION.cff`](CITATION.cff). The manuscript is currently under review; the citation will be updated with
the final volume/DOI upon publication.

## 13. License

- **Code:** MIT License — see [`LICENSE`](LICENSE).
- **Data and generated research outputs** (`data/raw/`, `outputs/`): Creative Commons Attribution 4.0
  International (CC BY 4.0).

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for how to report issues or propose changes, and
[`CHANGELOG.md`](CHANGELOG.md) for the release history.
