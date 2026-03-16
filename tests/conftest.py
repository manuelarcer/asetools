"""Shared fixtures for asetools test suite."""

from pathlib import Path

import pytest

TESTS_DIR = Path(__file__).parent
DATA_DIR = TESTS_DIR / "data"


@pytest.fixture
def tests_dir():
    """Path to the tests directory."""
    return TESTS_DIR


@pytest.fixture
def data_dir():
    """Path to the tests/data directory."""
    return DATA_DIR


@pytest.fixture
def sample_slab():
    """Create a simple FCC(111) slab for testing."""
    from ase.build import fcc111
    from ase.constraints import FixAtoms

    slab = fcc111("Cu", size=(2, 2, 4), vacuum=10.0)
    # Fix bottom two layers
    z = slab.positions[:, 2]
    mask = z < z.min() + 2.5
    slab.set_constraint(FixAtoms(mask=mask))
    return slab


@pytest.fixture
def sample_molecule():
    """Create a simple H2O molecule."""
    from ase import Atoms

    return Atoms("H2O", positions=[[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]])
