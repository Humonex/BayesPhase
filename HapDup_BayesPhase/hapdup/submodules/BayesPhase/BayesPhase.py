#!/usr/bin/env python3
"""Compatibility entry point for the HapDup-BayesPhase workflow.

The maintained BayesPhase command-line implementation lives at the repository
root as `BayesPhase_joint_phase.py`. This wrapper lets the HapDup integration
call a BayesPhase script from the historical `submodules/BayesPhase/` location
without duplicating the implementation.
"""

from pathlib import Path
import runpy
import sys


if __name__ == "__main__":
    repo_root = Path(__file__).resolve().parents[4]
    target = repo_root / "BayesPhase_joint_phase.py"
    if not target.exists():
        sys.exit(f"Cannot find BayesPhase entry point: {target}")
    runpy.run_path(str(target), run_name="__main__")
