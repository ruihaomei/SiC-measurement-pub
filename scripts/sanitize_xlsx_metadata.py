#!/usr/bin/env python3
"""Remove non-scientific source-path and author metadata from XLSX workbooks.

The measured cell values are preserved byte-for-byte inside the worksheet XML.
Only optional workbook/core-properties metadata is removed.
"""
from __future__ import annotations

import argparse
import re
import tempfile
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


ABS_PATH_RE = re.compile(br"<x15ac:absPath\b[^>]*/>")
CORE_PROPERTY_RE = re.compile(
    br"<(?:dc:creator|cp:lastModifiedBy)>.*?</(?:dc:creator|cp:lastModifiedBy)>"
)


def sanitize_workbook(path: Path) -> bool:
    """Strip non-scientific metadata from one XLSX file in place."""
    with ZipFile(path, "r") as src:
        members = [(info, src.read(info.filename)) for info in src.infolist()]

    changed = False
    with tempfile.NamedTemporaryFile(suffix=".xlsx", delete=False, dir=path.parent) as fh:
        tmp = Path(fh.name)
    try:
        with ZipFile(tmp, "w", compression=ZIP_DEFLATED) as dst:
            for info, payload in members:
                if info.filename == "xl/workbook.xml":
                    cleaned = ABS_PATH_RE.sub(b"", payload)
                    changed |= cleaned != payload
                    payload = cleaned
                elif info.filename == "docProps/core.xml":
                    cleaned = CORE_PROPERTY_RE.sub(b"", payload)
                    changed |= cleaned != payload
                    payload = cleaned
                dst.writestr(info, payload)
        if changed:
            tmp.replace(path)
    finally:
        if tmp.exists():
            tmp.unlink()
    return changed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbooks", nargs="+", type=Path)
    args = parser.parse_args()
    for path in args.workbooks:
        print(f"{path}: {'sanitized' if sanitize_workbook(path) else 'already clean'}")


if __name__ == "__main__":
    main()
