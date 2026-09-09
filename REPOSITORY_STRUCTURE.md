# Repository Structure

Map of the public reproducibility package.

```
.
├── README.md                 Start here (background, install, full reproduction)
├── LICENSE                   MIT (code) + CC BY 4.0 (data in data/raw/)
├── CITATION.cff              Machine-readable citation metadata
├── CONTRIBUTING.md           How to report issues / propose changes
├── REPRODUCIBILITY.md        Step-by-step reproduction + provenance guide
├── REPOSITORY_STRUCTURE.md   This file
├── CHANGELOG.md              Release history
├── requirements.txt          Python dependencies
├── requirements-validated.txt Exact direct-dependency versions used for release-candidate verification
├── pyproject.toml            Project metadata + pytest config
│
├── src/                      Importable reproduction library
│   ├── io_data.py            Load the four spectra from data/raw/ (percent→fraction)
│   ├── provenance.py         Run provenance + manifest helpers
│   ├── reference_params.py   SiC-like / Si-like "truth" systems for synthetic studies
│   ├── optics/               constants, fresnel, dielectric, forward (Airy/two-beam), si_stack
│   ├── diagnostics/          finesse.py  (lineshape finesse F = FSR/FWHM)
│   ├── inversion/            two_beam_fit.py  (two-beam / Airy fitting)
│   └── uncertainty/          jacobian.py, bootstrap.py  (conditional uncertainty)
│
├── scripts/                  Reproduction entry points
│   ├── reproduce_main_tables.py        Regenerate main tables
│   ├── reproduce_main_figures.py       Regenerate main figures
│   ├── finesse_threshold_study.py      Study A: |q|-tier adequacy calibration
│   ├── sic_joint_uncertainty.py        Study B: SiC joint Airy + conditional uncertainty
│   ├── damping_ablation.py             Study C: Si HF-correction ablation (Re(ε)/Im(ε))
│   ├── angular_identifiability_study.py Study D: simulated angular identifiability
│   ├── recompute_finesse_tables.py     Finesse Table 4 / B.7 recomputation
│   ├── sic_discrepancy_investigation.py Forward-model ladder on SiC data
│   ├── sic_two_beam_vs_airy_realdata.py Like-for-like two-beam vs Airy on real data
│   └── sanitize_xlsx_metadata.py       Remove non-scientific XLSX source-path metadata
│
├── tests/                    pytest suite
│   ├── test_optics_limits.py           Physical-limit + finesse sanity
│   ├── test_si_stack.py                Si stack forward model
│   ├── test_uncertainty.py             Jacobian / bootstrap machinery
│   ├── test_manifest_consistency.py    Selected headline manifest value == stored CSV
│   └── test_public_release_hygiene.py  Public XLSX metadata + export fallback
│
├── data/
│   └── raw/                  Measured spectra (CC BY 4.0)
│       ├── README.md          Data dictionary + sanitation note
│       ├── SHA256SUMS         Released-workbook checksums
│       ├── SiC_10deg_reflectance.xlsx        SiC 10°
│       ├── SiC_15deg_reflectance.xlsx        SiC 15°
│       ├── Si_10deg_reflectance.xlsx        Si 10°
│       └── Si_15deg_reflectance.xlsx        Si 15°
│
├── outputs/                  Generated artifacts (regenerable; see REPRODUCIBILITY.md)
│   ├── tables/               CSV tables (13 generated files)
│   ├── figures/              PNG + PDF figures (10 generated files)
│   ├── bootstrap/            Seeded bootstrap replicate files (3 generated files)
│   ├── logs/                 Run logs
│   └── manifest.json         Provenance: claim → file → value → script → seed → commit → qa_status
│
└── docs/
    └── provenance/           Human-readable provenance notes (machine record: outputs/manifest.json)
```

## Notes

- **`data/raw/` uses a lowercase path** matching `src/io_data.py` (`RAW_DIR = .../data/raw`); this is
  required on case-sensitive filesystems (Linux/GitHub Actions).
- **`outputs/` is regenerable.** Shipped outputs are the released, independently-reviewed versions; running
  the scripts overwrites them (see `REPRODUCIBILITY.md` §6).
- **`data/raw/` metadata is sanitized.** The worksheet values are unchanged; only optional Excel
  source-path metadata is removed by `scripts/sanitize_xlsx_metadata.py`.
