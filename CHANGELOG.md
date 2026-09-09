# Changelog

All notable changes to this reproducibility package are documented here. This project adheres loosely to
[Keep a Changelog](https://keepachangelog.com/) and [Semantic Versioning](https://semver.org/).

## [v1.0.0] — 2026-09-04

First public reproducibility release accompanying *Measurement* manuscript MEAS-D-26-04372.

### Added
- Curated open-source package: `src/`, `scripts/`, `tests/`, `outputs/`, `data/raw/`.
- Publication-grade documentation: `README.md`, `LICENSE` (MIT + CC BY 4.0 for data), `CITATION.cff`,
  `CONTRIBUTING.md`, `REPRODUCIBILITY.md`, `REPOSITORY_STRUCTURE.md`, this `CHANGELOG.md`.
- Machine-readable provenance in `outputs/manifest.json` (claim → file → value → script → seed → commit →
  qa_status) and human-readable notes under `docs/provenance/`.
- Raw-data checksums, an exact validated direct-dependency snapshot, and a lightweight GitHub Actions test
  workflow.

### Reproducibility / portability
- Data directory normalized to lowercase `data/raw/` so `src/io_data.py` resolves on case-sensitive
  filesystems.
- Removed non-scientific Excel source-path metadata from the four measured spectra and added an automated
  release-hygiene test.
- Added export-safe provenance handling for archive downloads and `git archive` exports.
- Fixed Python package discovery so wheels include all `src` subpackages.

### Notes
- The public branch is built from a sanitized **orphan root**.
- Results are conditional on the adopted reflectance model; no absolute-trueness claim is made; wider-angle
  analyses are simulated/virtual.
