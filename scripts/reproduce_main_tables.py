#!/usr/bin/env python3
"""Reproduce all machine-readable tables behind the revised manuscript (one command).

Runs the Phase-3 study scripts in order and prints a manifest summary. Each study writes its own CSVs to
outputs/tables/ and bootstrap/, and appends provenance to outputs/manifest.json.
"""
import os
import sys
import subprocess

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = [
    "scripts/finesse_threshold_study.py",
    "scripts/recompute_finesse_tables.py",
    "scripts/sic_two_beam_vs_airy_realdata.py",
    "scripts/sic_discrepancy_investigation.py",
    "scripts/sic_joint_uncertainty.py",
    "scripts/damping_ablation.py",
    "scripts/angular_identifiability_study.py",
]


def main():
    for s in SCRIPTS:
        print(f"\n===== running {s} =====", flush=True)
        rc = subprocess.call([sys.executable, os.path.join(REPO, s)])
        if rc != 0:
            print(f"!! {s} exited with {rc}", flush=True)
            sys.exit(rc)
    print("\nAll tables reproduced. See outputs/tables/, outputs/bootstrap/, outputs/manifest.json.")


if __name__ == "__main__":
    main()
