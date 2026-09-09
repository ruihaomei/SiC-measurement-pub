# Measured Spectra

This directory contains the four measured unpolarized specular IR reflectance spectra used by the
reproduction scripts.

| File | Material | Incidence angle | Columns |
|---|---|---|---|
| `SiC_10deg_reflectance.xlsx` | SiC epilayer | 10 degrees | `波数 (cm-1)`, `反射率 (%)` |
| `SiC_15deg_reflectance.xlsx` | SiC epilayer | 15 degrees | `波数 (cm-1)`, `反射率 (%)` |
| `Si_10deg_reflectance.xlsx` | Si epilayer | 10 degrees | `波数 (cm-1)`, `反射率 (%)` |
| `Si_15deg_reflectance.xlsx` | Si epilayer | 15 degrees | `波数 (cm-1)`, `反射率 (%)` |

The released workbooks are sanitized copies. `scripts/sanitize_xlsx_metadata.py` removes non-scientific
Excel source-path and author metadata while preserving worksheet cell values. Released-workbook SHA-256
checksums are in `SHA256SUMS`.

The data in this directory are licensed under CC BY 4.0; see the repository `LICENSE`.
