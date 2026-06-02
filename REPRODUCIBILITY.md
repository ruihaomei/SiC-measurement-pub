# Reproducibility Guide

This document tells an independent researcher exactly how to reproduce the quantitative results of
manuscript MEAS-D-26-04372 from a clean checkout, and how to verify provenance.

## 1. Environment

| Requirement | Value |
|---|---|
| Python | ≥ 3.10 |
| OS | Linux, macOS, or Windows (paths are case-correct for case-sensitive filesystems) |
| Network | Required only to install dependencies; numerical execution is offline afterward |
| Accelerators | None (CPU only) |
| External tools | None |

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

`requirements.txt` records minimum supported versions. `requirements-validated.txt` records the exact direct
dependency versions used for the release-candidate verification. For archival runs, record the fully resolved
environment with `pip freeze > environment-resolved.txt`.

## 2. One-command sanity check

```bash
pytest tests/
```

This verifies physical limits, unit handling, lineshape-finesse consistency, selected headline
manifest-to-CSV values, and public-release hygiene.

## 3. Full reproduction

```bash
python scripts/reproduce_main_tables.py
python scripts/reproduce_main_figures.py
```

Outputs land in `outputs/{tables,figures,bootstrap,logs}/`.
The two full workflows perform seeded optimization and bootstrap refits and can take substantial CPU time;
runtime depends on the machine and numerical libraries. The one-command test suite is the fast sanity check.

## 4. Provenance: how every number is traced

The machine-readable record is **`outputs/manifest.json`**. Each entry links:

```
claim_id → description → value → input_data → script → output_file → seed → commit → qa_status
```

- `input_data` points into `data/raw/`.
- `script` is the exact generating script in `scripts/`.
- `output_file` is the stored CSV/JSON/figure under `outputs/`.
- `seed` is the random seed for stochastic steps (e.g. `20260601` for the SiC bootstraps).
- `commit` identifies the checked-out source commit used to generate the row. In a GitHub ZIP or
  `git archive` export without `.git`, regenerated rows use the explicit `unavailable-exported-source`
  sentinel.
- `qa_status` records independent review status at release time.

The shipped `QA_APPROVED` rows preserve producing commit ids from the private audited development lineage.
That private history is intentionally excluded from this sanitized orphan release. Re-running a script
replaces its rows with the public checkout commit (or the export sentinel) and resets them to `PENDING`.

A human-readable summary is in [`docs/provenance/PROVENANCE.md`](docs/provenance/PROVENANCE.md).
Released raw-data checksums are in [`data/raw/SHA256SUMS`](data/raw/SHA256SUMS).

## 5. Determinism and expected variation

Stochastic procedures (wild bootstrap, differential evolution) are seeded. Reproduction is
**value-reproducible to the precision reported in the paper**. Trailing-digit differences (~1e-6 relative)
can arise from BLAS/threading order and floating-point non-associativity; they are far below any reported
digit and do not change any conclusion.

## 6. Important: outputs are regenerated in place

Running a reproduction script **overwrites** files under `outputs/`, including `outputs/manifest.json`, and
resets recomputed rows to a `PENDING` review status (a producer cannot self-approve its own output). This is
expected. To keep the shipped, independently-reviewed manifest intact:

```bash
git restore outputs/        # discard local regeneration, restore the released outputs
```

or run reproductions against a copy of the repository.

## 7. What is NOT reproduced here

- The manuscript DOCX/PDF (the article itself) is not in this repository.
- Wider-angle (20°–40°) results are **simulated** design studies, not measurements, and are labelled as such
  in code and outputs.
