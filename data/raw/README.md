# Measured Spectra

This directory contains the four measured unpolarized specular IR reflectance spectra used by the
reproduction scripts.

| File | Material | Incidence angle | Columns |
|---|---|---|---|
| `附件1.xlsx` | SiC epilayer | 10 degrees | `波数 (cm-1)`, `反射率 (%)` |
| `附件2.xlsx` | SiC epilayer | 15 degrees | `波数 (cm-1)`, `反射率 (%)` |
| `附件3.xlsx` | Si epilayer | 10 degrees | `波数 (cm-1)`, `反射率 (%)` |
| `附件4.xlsx` | Si epilayer | 15 degrees | `波数 (cm-1)`, `反射率 (%)` |

The released workbooks are sanitized copies. `scripts/sanitize_xlsx_metadata.py` removes non-scientific
Excel source-path and author metadata while preserving worksheet cell values. Released-workbook SHA-256
checksums are in `SHA256SUMS`.

The data in this directory are licensed under CC BY 4.0; see the repository `LICENSE`.
