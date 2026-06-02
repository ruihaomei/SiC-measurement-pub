# Contributing

Thank you for your interest in this reproducibility package for *Measurement* manuscript MEAS-D-26-04372.

This repository accompanies a manuscript under peer review. Its primary purpose is **transparent
reproducibility**, so the bar for changes that alter any reported scientific value is intentionally high.

## Reporting issues

Please open a GitHub issue for:

- a reproduction step that fails or is unclear;
- a discrepancy between a reported value and what the code produces;
- environment/installation problems (please include OS, Python version, and `pip freeze`);
- documentation errors.

Include the exact command, the full error output, and your platform.

## Proposing changes

1. Fork and create a feature branch.
2. Keep changes focused and well-described.
3. Run `pytest tests/` and ensure all tests pass.
4. If you touch any analysis code, explain the scientific effect and update `outputs/manifest.json`
   provenance if values change. **Do not hand-edit generated outputs to match a desired number.**
5. Open a pull request describing the motivation and the verification you ran.

## Scientific-integrity rules

- Never introduce a numerical value that is not produced by a script in this repository from the data in
  `data/raw/`.
- Keep the distinction explicit between **measured** (10°/15°) and **simulated/virtual** (wider-angle)
  results.
- Do not present conditional, model-dependent results as absolute trueness/accuracy.

## Code style

- Python ≥ 3.10, standard scientific stack (`numpy`/`scipy`/`pandas`/`matplotlib`).
- Prefer small, readable functions with docstrings; match the style of the surrounding code.
- No hard-coded absolute paths or machine-specific assumptions; resolve paths relative to the repository
  root (see `src/io_data.py`).
