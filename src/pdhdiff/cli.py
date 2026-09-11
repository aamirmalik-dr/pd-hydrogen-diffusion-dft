"""Command-line entry point: ``pdhdiff <command>``.

Commands:

* ``profile``   extract the energy profile from the calculations tree
* ``analyze``   fit the profile and evaluate rates and the diffusion constant
* ``electronic``  plot the DOS and band structure of bulk Pd
* ``path``      write the relaxed path trajectory and geometry figures
* ``diagnostics``  residual forces and convergence diagnostics of the relaxations
* ``inputs``    generate CP-PAW structure files for a constrained path
* ``check``     compare regenerated result files with the committed ones
* ``all``       run every step in order
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "scripts"


def _run(name: str, argv: list[str]) -> int:
    import runpy

    sys.argv = [name, *argv]
    runpy.run_path(str(SCRIPTS / name), run_name="__main__")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Dispatch to the scripts in ``scripts/``."""
    parser = argparse.ArgumentParser(prog="pdhdiff", description=__doc__)
    parser.add_argument(
        "command",
        choices=[
            "profile",
            "analyze",
            "electronic",
            "path",
            "diagnostics",
            "inputs",
            "check",
            "all",
        ],
    )
    parser.add_argument("rest", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    table = {
        "profile": "extract_profile.py",
        "analyze": "analyze_diffusion.py",
        "electronic": "plot_electronic_structure.py",
        "path": "make_path_figures.py",
        "diagnostics": "relaxation_diagnostics.py",
        "inputs": "build_path_inputs.py",
        "check": "check_reproduction.py",
        "all": "run_all.py",
    }
    if not (SCRIPTS / table[args.command]).exists():
        raise SystemExit(
            "the pdhdiff console script drives the scripts/ directory of a source checkout; "
            f"{SCRIPTS} was not found. Clone the repository and install it with pip install -e ."
        )
    return _run(table[args.command], args.rest)


if __name__ == "__main__":
    raise SystemExit(main())
