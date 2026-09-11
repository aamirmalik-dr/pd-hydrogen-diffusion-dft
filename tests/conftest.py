"""Shared fixtures: paths into the committed calculations tree."""

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def calc_root() -> Path:
    return ROOT / "calculations"


@pytest.fixture(scope="session")
def results_root() -> Path:
    return ROOT / "results"


@pytest.fixture(scope="session")
def path_root(calc_root) -> Path:
    return calc_root / "04_path_octa_to_tetra"
