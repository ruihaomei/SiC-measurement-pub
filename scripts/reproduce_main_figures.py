#!/usr/bin/env python3
"""Reproduce all manuscript figures (one command). Figures are produced as a side-effect of the study
scripts (run scripts/reproduce_main_tables.py first, or this script which re-runs the figure-bearing ones).

Figures produced:
  outputs/figures/finesse_bias.{png,pdf}              (Study A)
  outputs/figures/sic_bootstrap.{png,pdf}             (Study B)
  outputs/figures/damping_residuals.{png,pdf}         (Study C)
  outputs/figures/damping_corr_heatmap.{png,pdf}      (Study C)
  outputs/figures/angular_identifiability.{png,pdf}   (Study D)
"""
import os
import sys
import subprocess

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = [
    "scripts/finesse_threshold_study.py",
    "scripts/sic_joint_uncertainty.py",
    "scripts/damping_ablation.py",
    "scripts/angular_identifiability_study.py",
]


def main():
    for s in SCRIPTS:
        print(f"\n===== running {s} (figures) =====", flush=True)
        rc = subprocess.call([sys.executable, os.path.join(REPO, s)])
        if rc != 0:
            sys.exit(rc)
    figs = sorted(f for f in os.listdir(os.path.join(REPO, "outputs", "figures")) if f.endswith(".png"))
    print("\nFigures present:", figs)


if __name__ == "__main__":
    main()
