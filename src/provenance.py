"""Producer-side manifest helpers.

Analysis scripts may only propose fresh ``PENDING`` rows. Independent QA owns every
subsequent status transition.
"""
import json
import os
import subprocess


def source_commit(repo):
    """Return the checked-out source commit, or an explicit export sentinel.

    GitHub ZIP downloads and ``git archive`` exports do not include ``.git``.
    Reproduction remains supported there, but the generated manifest must state
    that the source commit metadata is unavailable rather than failing.
    """
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=repo, text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        return "unavailable-exported-source"


def update_manifest(repo, entries):
    """Replace manifest rows by claim id and reset them to producer-side PENDING."""
    path = os.path.join(repo, "outputs", "manifest.json")
    with open(path, encoding="utf-8") as f:
        manifest = json.load(f)
    commit = source_commit(repo)
    claim_ids = {entry["claim_id"] for entry in entries}
    manifest["claims"] = [
        claim for claim in manifest["claims"]
        if claim.get("claim_id") not in claim_ids
    ]
    for entry in entries:
        row = dict(entry)
        row["commit"] = commit
        row["manuscript_location"] = None
        row["qa_status"] = "PENDING"
        row["qa_note"] = ""
        manifest["claims"].append(row)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)
