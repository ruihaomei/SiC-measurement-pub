"""CSV-backed manifest reproducibility guards.

These tests validate that the provenance manifest is internally consistent and that practical scalar
manifest claims match the generated CSV output files.
"""
import csv
import json
import os
import re

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANIFEST = os.path.join(REPO, "outputs", "manifest.json")
ALLOWED_QA = {"PENDING", "UNDER_REVIEW", "REJECTED", "REQUIRES_REVISION", "QA_APPROVED"}


def _claims():
    with open(MANIFEST, encoding="utf-8") as f:
        return json.load(f)["claims"]


def _claim(claim_id):
    return next(c for c in _claims() if c["claim_id"] == claim_id)


def _rows(path):
    with open(os.path.join(REPO, path), newline="") as f:
        return list(csv.DictReader(f))


def _quantity_rows(path):
    return {row["quantity"]: row for row in _rows(path)}


def _single_um(value):
    match = re.fullmatch(r"(-?[0-9.]+) um", value)
    assert match, value
    return float(match.group(1))


def _ci_um(value):
    match = re.fullmatch(r"\((-?[0-9.]+),(-?[0-9.]+)\) um", value)
    assert match, value
    return tuple(float(v) for v in match.groups())


def test_manifest_parses_and_has_expected_claims():
    ids = {claim["claim_id"] for claim in _claims()}
    assert {
        "FINESSE_GRID", "FINESSE_TABLE4_REALDATA", "FINESSE_TABLEB7_THEORY",
        "SIC_JOINT_D", "SIC_10DEG_BOOT_CI", "SIC_15DEG_BOOT_CI",
        "DAMPING_ABLATION", "SI_CHOSEN_MODEL", "ANGULAR_IDENTIFIABILITY",
    } <= ids


def test_manifest_claims_have_producer_provenance_and_existing_outputs():
    for claim in _claims():
        assert claim.get("qa_status") in ALLOWED_QA, claim["claim_id"]
        assert claim.get("script"), f"{claim['claim_id']} has no generating script"
        commit = claim.get("commit", "")
        assert re.fullmatch(r"[0-9a-f]{40}", commit) or commit == "unavailable-exported-source", claim["claim_id"]
        output = claim.get("output_file")
        assert output, f"{claim['claim_id']} has no output_file"
        assert os.path.exists(os.path.join(REPO, output)), f"{claim['claim_id']} output missing: {output}"
        assert claim.get("manuscript_location") is None


def test_qa_approved_claims_have_review_notes():
    for claim in _claims():
        if claim.get("qa_status") == "QA_APPROVED":
            assert claim.get("qa_note"), f"{claim['claim_id']} approved without a QA note"


def test_grid_manifest_count_matches_csv():
    rows = _rows("outputs/tables/finesse_threshold_grid.csv")
    assert _claim("FINESSE_GRID")["value"] == f"{len(rows)} cells"


def test_sic_manifest_scalars_match_csv_and_replicates():
    quantities = _quantity_rows("outputs/tables/sic_results.csv")
    assert _single_um(_claim("SIC_JOINT_D")["value"]) == float(quantities["d_joint_reported"]["value_um"])
    assert _single_um(_claim("SIC_JOINT_SIGMA")["value"]) == float(quantities["jacobian_sigma_d"]["value_um"])
    assert _single_um(_claim("SIC_CROSS_ANGLE_DIFF")["value"]) == float(quantities["cross_angle_diff"]["value_um"])

    for claim_id, path in [
        ("SIC_JOINT_BOOT_CI", "outputs/bootstrap/sic_joint_replicates.csv"),
        ("SIC_10DEG_BOOT_CI", "outputs/bootstrap/sic_10deg_replicates.csv"),
        ("SIC_15DEG_BOOT_CI", "outputs/bootstrap/sic_15deg_replicates.csv"),
    ]:
        values = np.array([float(row["replicate_d_um"]) for row in _rows(path)])
        expected = tuple(round(float(v), 4) for v in np.percentile(values, [2.5, 97.5]))
        assert _ci_um(_claim(claim_id)["value"]) == expected


def test_damping_and_si_manifest_scalars_match_csv():
    ablation = {row["mechanism"]: row for row in _rows("outputs/tables/damping_ablation.csv")}
    none_d = float(ablation["none"]["d_um"])
    real_d = float(ablation["empirical_real_correction"]["d_um"])
    shift_nm = round(abs(none_d - real_d) * 1000, 1)
    assert _claim("DAMPING_ABLATION")["value"] == (
        f"faithful Re(epsilon) correction removal shifts d by {shift_nm:.1f} nm"
    )
    si = _quantity_rows("outputs/tables/si_results.csv")
    assert _single_um(_claim("SI_CHOSEN_MODEL")["value"]) == float(si["d_chosen"]["value"])


def test_finesse_realdata_outputs_preserve_definition_distinction():
    table4 = _rows("outputs/tables/finesse_table4.csv")
    tableb7 = _rows("outputs/tables/finesse_tableB7.csv")
    assert len(table4) == len(tableb7) == 4
    assert {row["diagnostic_kind"] for row in table4} == {"lineshape_observable_FSR_over_FWHM"}
    assert {row["diagnostic_kind"] for row in tableb7} == {"theoretical_roundtrip_proxy_not_lineshape"}
    assert all("report explicit q statistics" in row["tableB7_action"] for row in tableb7)


def test_angular_identifiability_rows_are_all_simulated():
    rows = _rows("outputs/tables/angular_identifiability.csv")
    assert rows
    assert {row["analysis_kind"] for row in rows} == {"simulated_virtual"}
    tagged = [row for row in rows if row["pair_matches_available_measured_angles"] == "yes"]
    assert tagged
    assert all(row["system"] == "SiC_MDF" and row["angle_a"] == "10.0" and row["angle_b"] == "15.0"
               for row in tagged)


def test_no_claim_uses_rejected_number_silently():
    for claim in _claims():
        note = (claim.get("qa_note") or "").lower()
        description = (claim.get("description") or "").lower()
        if "0.293" in description or "0.293" in str(claim.get("value", "")):
            assert "reject" in note or "non-like-for-like" in note or "not" in note
