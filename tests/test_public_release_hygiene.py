"""Public-release privacy and export-reproduction guards."""
import shutil
from pathlib import Path
from zipfile import ZipFile

import src.provenance as provenance
from scripts.sanitize_xlsx_metadata import sanitize_workbook


REPO = Path(__file__).resolve().parents[1]


def test_raw_workbooks_do_not_embed_source_paths_or_authors():
    for path in sorted((REPO / "data" / "raw").glob("*.xlsx")):
        with ZipFile(path) as workbook:
            workbook_xml = workbook.read("xl/workbook.xml")
            core_xml = workbook.read("docProps/core.xml")
        assert b"absPath" not in workbook_xml, path.name
        assert b"<dc:creator>" not in core_xml, path.name
        assert b"<cp:lastModifiedBy>" not in core_xml, path.name


def test_source_commit_has_export_fallback(monkeypatch):
    def fail(*args, **kwargs):
        raise provenance.subprocess.CalledProcessError(128, ["git"])

    monkeypatch.setattr(provenance.subprocess, "check_output", fail)
    assert provenance.source_commit(REPO) == "unavailable-exported-source"


def test_workbook_sanitizer_is_idempotent(tmp_path):
    source = REPO / "data" / "raw" / "附件1.xlsx"
    workbook = tmp_path / source.name
    shutil.copyfile(source, workbook)
    before = workbook.read_bytes()
    assert sanitize_workbook(workbook) is False
    assert workbook.read_bytes() == before
