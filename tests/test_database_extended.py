"""Extended tests for asetools.database.databases module.

Tests check_if_exists_in_db and db_to_pandas with real ASE database objects.
"""

import numpy as np
import pandas as pd
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.db import connect


@pytest.fixture
def db_with_atoms(tmp_path):
    """Create a temp ASE database with one entry."""
    db_path = tmp_path / "test.db"
    db = connect(str(db_path))

    atoms = Atoms("H2O", positions=[[0, 0, 0], [0, 0, 0.96], [0, 0.96, 0]],
                  cell=[5, 5, 5], pbc=True)
    forces = np.array([[0.1, 0.2, 0.3], [-0.1, -0.2, -0.3], [0.0, 0.0, 0.0]])
    calc = SinglePointCalculator(atoms, energy=-10.5, forces=forces,
                                 free_energy=-10.6, magmom=0.0)
    atoms.calc = calc
    db.write(atoms, name="water")
    return db, atoms, forces


@pytest.fixture
def empty_db(tmp_path):
    """Create an empty temp ASE database."""
    db_path = tmp_path / "empty.db"
    return connect(str(db_path))


class TestCheckIfExistsInDb:
    def test_finds_existing_atoms(self, db_with_atoms):
        from asetools.database.databases import check_if_exists_in_db
        db, atoms, _ = db_with_atoms
        found, idx = check_if_exists_in_db(db, atoms)
        assert found is True
        assert idx is not None

    def test_not_found_different_atoms(self, db_with_atoms):
        from asetools.database.databases import check_if_exists_in_db
        db, _, _ = db_with_atoms

        other = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]],
                      cell=[5, 5, 5], pbc=True)
        forces = np.array([[0.5, 0.5, 0.5], [-0.5, -0.5, -0.5]])
        calc = SinglePointCalculator(other, energy=-5.0, forces=forces)
        other.calc = calc

        found, idx = check_if_exists_in_db(db, other)
        assert found is False
        assert idx is None

    def test_empty_db(self, empty_db):
        from asetools.database.databases import check_if_exists_in_db

        atoms = Atoms("H", positions=[[0, 0, 0]], cell=[5, 5, 5], pbc=True)
        forces = np.array([[0.0, 0.0, 0.0]])
        calc = SinglePointCalculator(atoms, energy=-1.0, forces=forces)
        atoms.calc = calc

        found, idx = check_if_exists_in_db(empty_db, atoms)
        assert found is False
        assert idx is None


class TestDbToPandas:
    def test_returns_dataframe(self, db_with_atoms):
        from asetools.database.databases import db_to_pandas
        db, _, _ = db_with_atoms
        df = db_to_pandas(db)
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 1

    def test_default_columns(self, db_with_atoms):
        from asetools.database.databases import db_to_pandas
        db, _, _ = db_with_atoms
        df = db_to_pandas(db)
        for col in ["name", "id", "energy", "free_energy", "magmom"]:
            assert col in df.columns

    def test_custom_columns(self, db_with_atoms):
        from asetools.database.databases import db_to_pandas
        db, _, _ = db_with_atoms
        df = db_to_pandas(db, columns=["name", "id"])
        assert list(df.columns) == ["name", "id"]

    def test_values_correct(self, db_with_atoms):
        from asetools.database.databases import db_to_pandas
        db, _, _ = db_with_atoms
        df = db_to_pandas(db)
        assert df.iloc[0]["name"] == "water"
        assert abs(df.iloc[0]["energy"] - (-10.5)) < 0.01

    def test_multiple_entries(self, tmp_path):
        from asetools.database.databases import db_to_pandas
        db = connect(str(tmp_path / "multi.db"))

        for i, name in enumerate(["entry1", "entry2", "entry3"]):
            atoms = Atoms("H", positions=[[0, 0, float(i)]], cell=[5, 5, 5], pbc=True)
            calc = SinglePointCalculator(atoms, energy=-float(i), forces=np.zeros((1, 3)),
                                         free_energy=-float(i), magmom=0.0)
            atoms.calc = calc
            db.write(atoms, name=name)

        df = db_to_pandas(db)
        assert len(df) == 3
