"""Tests for CLI tool functions that can be tested without VASP files."""

import os
import shutil
import tempfile

import pytest


class TestRemoveSlashes:
    """Test remove_slashes CLI tool functions."""

    def setup_method(self):
        self.test_dir = tempfile.mkdtemp()

    def teardown_method(self):
        shutil.rmtree(self.test_dir)

    def test_remove_slashes_basic(self):
        from asetools.cli.remove_slashes import remove_slashes_from_file

        infile = os.path.join(self.test_dir, "input.txt")
        outfile = os.path.join(self.test_dir, "output.txt")
        with open(infile, "w") as f:
            f.write("path/to/some/file\n")

        remove_slashes_from_file(infile, outfile)

        with open(outfile) as f:
            assert f.read() == "pathtosomefile\n"

    def test_remove_slashes_no_slashes(self):
        from asetools.cli.remove_slashes import remove_slashes_from_file

        infile = os.path.join(self.test_dir, "input.txt")
        outfile = os.path.join(self.test_dir, "output.txt")
        with open(infile, "w") as f:
            f.write("no slashes here\n")

        remove_slashes_from_file(infile, outfile)

        with open(outfile) as f:
            assert f.read() == "no slashes here\n"

    def test_remove_slashes_inplace(self):
        from asetools.cli.remove_slashes import remove_slashes_from_file

        infile = os.path.join(self.test_dir, "input.txt")
        with open(infile, "w") as f:
            f.write("a/b/c\n")

        remove_slashes_from_file(infile)  # No output file → overwrite

        with open(infile) as f:
            assert f.read() == "abc\n"

    def test_remove_slashes_empty_file(self):
        from asetools.cli.remove_slashes import remove_slashes_from_file

        infile = os.path.join(self.test_dir, "empty.txt")
        outfile = os.path.join(self.test_dir, "output.txt")
        with open(infile, "w") as f:
            f.write("")

        remove_slashes_from_file(infile, outfile)

        with open(outfile) as f:
            assert f.read() == ""

    def test_remove_slashes_only_slashes(self):
        from asetools.cli.remove_slashes import remove_slashes_from_file

        infile = os.path.join(self.test_dir, "slashes.txt")
        outfile = os.path.join(self.test_dir, "output.txt")
        with open(infile, "w") as f:
            f.write("///\n")

        remove_slashes_from_file(infile, outfile)

        with open(outfile) as f:
            assert f.read() == "\n"


class TestReorderAtoms:
    """Test reorder_atoms CLI functions."""

    def test_reorder_atoms_default_alphabetical(self):
        from ase import Atoms

        from asetools.cli.reorder_atoms import reorder_atoms

        atoms = Atoms("OHH", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        reordered = reorder_atoms(atoms)
        symbols = reordered.get_chemical_symbols()
        assert symbols == ["H", "H", "O"]

    def test_reorder_atoms_custom_order(self):
        from ase import Atoms

        from asetools.cli.reorder_atoms import reorder_atoms

        atoms = Atoms("CuOH", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        reordered = reorder_atoms(atoms, element_order=["O", "H", "Cu"])
        symbols = reordered.get_chemical_symbols()
        assert symbols == ["O", "H", "Cu"]

    def test_reorder_atoms_z_order_top_bottom(self):
        from ase import Atoms

        from asetools.cli.reorder_atoms import reorder_atoms

        atoms = Atoms("HH", positions=[[0, 0, 1], [0, 0, 5]])
        reordered = reorder_atoms(atoms, z_order="top-bottom")
        # Highest z first
        assert reordered.positions[0][2] == pytest.approx(5.0)
        assert reordered.positions[1][2] == pytest.approx(1.0)

    def test_reorder_atoms_z_order_bottom_top(self):
        from ase import Atoms

        from asetools.cli.reorder_atoms import reorder_atoms

        atoms = Atoms("HH", positions=[[0, 0, 5], [0, 0, 1]])
        reordered = reorder_atoms(atoms, z_order="bottom-top")
        assert reordered.positions[0][2] == pytest.approx(1.0)
        assert reordered.positions[1][2] == pytest.approx(5.0)

    def test_reorder_atoms_preserves_cell(self):
        from ase import Atoms

        from asetools.cli.reorder_atoms import reorder_atoms

        atoms = Atoms("OH", positions=[[0, 0, 0], [1, 0, 0]], cell=[10, 10, 10], pbc=True)
        reordered = reorder_atoms(atoms)
        assert (reordered.cell == atoms.cell).all()
        assert all(reordered.pbc)

    def test_reorder_atoms_missing_element_in_order(self):
        from ase import Atoms

        from asetools.cli.reorder_atoms import reorder_atoms

        atoms = Atoms("CuOH", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        # Only specify Cu and O — H should be appended
        reordered = reorder_atoms(atoms, element_order=["Cu", "O"])
        symbols = reordered.get_chemical_symbols()
        assert symbols[0] == "Cu"
        assert symbols[1] == "O"
        assert symbols[2] == "H"

    def test_get_fixed_indices_empty(self):
        from ase import Atoms

        from asetools.cli.reorder_atoms import get_fixed_indices

        atoms = Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
        assert get_fixed_indices(atoms) == []

    def test_get_fixed_indices_with_constraints(self):
        from ase import Atoms
        from ase.constraints import FixAtoms

        from asetools.cli.reorder_atoms import get_fixed_indices

        atoms = Atoms("H3", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]])
        atoms.set_constraint(FixAtoms(indices=[0, 2]))
        indices = get_fixed_indices(atoms)
        assert list(indices) == [0, 2]


class TestCompareVaspParam:
    """Test comparevaspparam CLI functions."""

    def setup_method(self):
        self.test_dir = tempfile.mkdtemp()

    def teardown_method(self):
        shutil.rmtree(self.test_dir)

    def test_read_outcar_extracts_keywords(self):
        from asetools.cli.comparevaspparam import read_outcar

        outcar = os.path.join(self.test_dir, "OUTCAR")
        with open(outcar, "w") as f:
            f.write("   ENCUT  =  500.00 eV\n")
            f.write("   EDIFF  =  0.1E-05\n")
            f.write("   ISPIN  =      2\n")
            f.write("   random line without keywords\n")

        params = read_outcar(outcar)
        assert "ENCUT" in params
        assert "EDIFF" in params
        assert "ISPIN" in params
        assert len(params) == 3

    def test_read_outcar_empty_file(self):
        from asetools.cli.comparevaspparam import read_outcar

        outcar = os.path.join(self.test_dir, "OUTCAR")
        with open(outcar, "w") as f:
            f.write("")

        params = read_outcar(outcar)
        assert params == {}

    def test_compare_parameters_identical(self, capsys):
        from asetools.cli.comparevaspparam import compare_parameters

        p1 = {"ENCUT": "ENCUT = 500", "EDIFF": "EDIFF = 1E-5"}
        p2 = {"ENCUT": "ENCUT = 500", "EDIFF": "EDIFF = 1E-5"}
        compare_parameters(p1, p2)
        captured = capsys.readouterr()
        assert captured.out == ""

    def test_compare_parameters_different(self, capsys):
        from asetools.cli.comparevaspparam import compare_parameters

        p1 = {"ENCUT": "ENCUT = 500"}
        p2 = {"ENCUT": "ENCUT = 600"}
        compare_parameters(p1, p2)
        captured = capsys.readouterr()
        assert "Difference in ENCUT" in captured.out

    def test_compare_parameters_missing_key(self, capsys):
        from asetools.cli.comparevaspparam import compare_parameters

        p1 = {"ENCUT": "ENCUT = 500"}
        p2 = {"EDIFF": "EDIFF = 1E-5"}
        compare_parameters(p1, p2)
        captured = capsys.readouterr()
        assert "Difference" in captured.out


class TestVaspBackup:
    """Test vaspbackup CLI functions."""

    def setup_method(self):
        self.test_dir = tempfile.mkdtemp()
        self.orig_dir = os.getcwd()
        os.chdir(self.test_dir)

    def teardown_method(self):
        os.chdir(self.orig_dir)
        shutil.rmtree(self.test_dir)

    def test_compress_file(self):
        import gzip

        from asetools.cli.vaspbackup import compress_file

        testfile = os.path.join(self.test_dir, "testfile.txt")
        with open(testfile, "w") as f:
            f.write("test content " * 100)

        compress_file(testfile)

        assert not os.path.exists(testfile)
        assert os.path.exists(testfile + ".gz")

        with gzip.open(testfile + ".gz", "rt") as f:
            content = f.read()
        assert content == "test content " * 100

    def test_main_creates_backup(self):
        """Test that main creates backup folder and copies/moves files."""
        import sys

        # Create dummy VASP files
        with open("CONTCAR", "w") as f:
            f.write("dummy contcar\n")
        with open("INCAR", "w") as f:
            f.write("dummy incar\n")
        with open("POSCAR", "w") as f:
            f.write("dummy poscar\n")

        from asetools.cli.vaspbackup import main

        sys.argv = ["vaspbackup", "test_backup"]
        main()

        assert os.path.isdir("test_backup")
        assert os.path.exists("test_backup/CONTCAR")
        assert os.path.exists("test_backup/INCAR")

    def test_main_existing_backup_aborts(self, capsys):
        """Test that main doesn't overwrite existing backup."""
        import sys

        os.makedirs("existing_backup")

        from asetools.cli.vaspbackup import main

        sys.argv = ["vaspbackup", "existing_backup"]
        main()

        captured = capsys.readouterr()
        assert "already exists" in captured.out


class TestGibbsFE:
    """Test gibbsFE (thermochem) helper functions."""

    def test_extract_vib_info_real_mode(self):
        """Test vibration extraction with synthetic OUTCAR lines."""
        from asetools.cli.gibbsFE import extract_vib_info

        lines = [
            "   1 f  =   95.123456 THz   597.651 2PiTHz 3172.838 cm-1   393.299 meV\n",
            "             X         Y         Z           dx          dy          dz\n",
            "      0.000   0.000   0.000    0.100    0.200    0.300\n",
            "      1.000   0.000   0.000   -0.100   -0.200   -0.300\n",
            "\n",
        ]
        vib, vasp6 = extract_vib_info(lines)
        assert "1" in vib
        assert vib["1"]["freq"] == pytest.approx(3172.838)
        assert vib["1"]["e"] == pytest.approx(0.393299)
        assert not vasp6

    def test_extract_vib_info_imaginary(self):
        from asetools.cli.gibbsFE import extract_vib_info

        lines = [
            "   1 f/i=   10.123456 THz    63.601 2PiTHz  337.638 cm-1    41.856 meV\n",
            "             X         Y         Z           dx          dy          dz\n",
            "      0.000   0.000   0.000    0.100    0.200    0.300\n",
            "\n",
        ]
        vib, _ = extract_vib_info(lines)
        assert vib["1"]["freq"] < 0  # Imaginary → negative

    def test_compute_corrections(self):
        from asetools.cli.gibbsFE import compute_corrections

        vib = {
            "1": {"e": 0.4, "freq": 3000},
            "2": {"e": 0.2, "freq": 1500},
            "3": {"e": 0.1, "freq": 500},
        }
        zpe, cpT, S = compute_corrections(vib, T=298)
        assert zpe == pytest.approx(0.35)  # (0.4 + 0.2 + 0.1) / 2
        assert cpT >= 0
        assert S >= 0

    def test_parse_atom_indices_range(self):
        from asetools.cli.gibbsFE import parse_atom_indices

        indices = parse_atom_indices("0-4")
        assert indices == {0, 1, 2, 3, 4}

    def test_parse_atom_indices_comma(self):
        from asetools.cli.gibbsFE import parse_atom_indices

        indices = parse_atom_indices("0,2,5")
        assert indices == {0, 2, 5}

    def test_parse_atom_indices_mixed(self):
        from asetools.cli.gibbsFE import parse_atom_indices

        indices = parse_atom_indices("0-2,5,8-10")
        assert indices == {0, 1, 2, 5, 8, 9, 10}

    def test_parse_atom_indices_none(self):
        from asetools.cli.gibbsFE import parse_atom_indices

        assert parse_atom_indices(None) is None

    def test_calculate_mode_character_all_atoms(self):
        from asetools.cli.gibbsFE import calculate_mode_character

        eigenvectors = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        char = calculate_mode_character(eigenvectors, {0, 1, 2})
        assert char == pytest.approx(1.0)

    def test_calculate_mode_character_subset(self):
        from asetools.cli.gibbsFE import calculate_mode_character

        eigenvectors = [[1, 0, 0], [1, 0, 0]]
        char = calculate_mode_character(eigenvectors, {0})
        assert char == pytest.approx(0.5)

    def test_calculate_mode_character_no_selection(self):
        from asetools.cli.gibbsFE import calculate_mode_character

        eigenvectors = [[1, 0, 0]]
        char = calculate_mode_character(eigenvectors, set())
        assert char == pytest.approx(1.0)  # Empty set → include all

    def test_is_close(self):
        from asetools.cli.gibbsFE import is_close

        assert is_close(1.0, 1.0005) is True
        assert is_close(1.0, 1.01) is False

    def test_determineshift(self):
        from asetools.cli.gibbsFE import determineshift

        shift = [[1, 2, 3], [4, 5, 6]]
        result = determineshift(shift, k=2)
        assert result == [[2, 4, 6], [8, 10, 12]]

    def test_determineshift_zero(self):
        from asetools.cli.gibbsFE import determineshift

        shift = [[1, 2, 3]]
        result = determineshift(shift, k=0)
        assert result == [[0, 0, 0]]


class TestSummaryFolders:
    """Test summaryfolders CLI functions."""

    def test_get_pyatoms_steps(self):
        from asetools.cli.summaryfolders import get_pyatoms_steps

        test_dir = tempfile.mkdtemp()
        try:
            # Create step folders
            for i in range(3):
                os.makedirs(os.path.join(test_dir, f"step_{i:02d}"))
            os.makedirs(os.path.join(test_dir, "other_dir"))

            steps = get_pyatoms_steps(test_dir)
            assert len(steps) == 3
            assert "step_00" in steps[0]
            assert "step_02" in steps[2]
        finally:
            shutil.rmtree(test_dir)

    def test_get_pyatoms_steps_empty(self):
        from asetools.cli.summaryfolders import get_pyatoms_steps

        test_dir = tempfile.mkdtemp()
        try:
            steps = get_pyatoms_steps(test_dir)
            assert steps == []
        finally:
            shutil.rmtree(test_dir)
