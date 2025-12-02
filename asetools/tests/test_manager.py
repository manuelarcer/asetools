import os
import shutil
import tempfile
import pytest
from ase.io import read, write
from asetools.manager.manager import load_structure


class TestLoadStructure:
    """Test the load_structure() function with various scenarios."""

    def setup_method(self):
        """Create a temporary directory for each test."""
        self.test_dir = tempfile.mkdtemp()
        self.original_dir = os.getcwd()
        os.chdir(self.test_dir)

    def teardown_method(self):
        """Clean up temporary directory after each test."""
        os.chdir(self.original_dir)
        shutil.rmtree(self.test_dir)

    def test_load_from_outcar(self):
        """Test loading from OUTCAR when available (priority 1)."""
        # Copy test OUTCAR to temp directory
        src_outcar = os.path.join(self.original_dir, 'asetools/data/OUTCAR')
        shutil.copy(src_outcar, 'OUTCAR')

        # Load structure
        atoms = load_structure()

        # Verify it loaded from OUTCAR
        assert atoms is not None
        assert len(atoms) > 0

        # Verify by comparing to direct OUTCAR read
        expected = read('OUTCAR', format='vasp-out', index=-1)
        assert len(atoms) == len(expected)
        assert atoms.get_chemical_formula() == expected.get_chemical_formula()

    def test_load_from_outcar_vasp6(self):
        """Test loading from VASP 6 OUTCAR."""
        # Copy test OUTCAR_vasp6 to temp directory
        src_outcar = os.path.join(self.original_dir, 'asetools/data/OUTCAR_vasp6')
        shutil.copy(src_outcar, 'OUTCAR')

        # Load structure
        atoms = load_structure()

        # Verify it loaded successfully
        assert atoms is not None
        assert len(atoms) > 0

        # Verify by comparing to direct OUTCAR read
        expected = read('OUTCAR', format='vasp-out', index=-1)
        assert len(atoms) == len(expected)

    def test_fallback_to_contcar(self):
        """Test fallback to CONTCAR when OUTCAR doesn't exist (priority 2)."""
        # Copy test CONTCAR to temp directory (no OUTCAR)
        src_contcar = os.path.join(self.original_dir, 'asetools/data/CONTCAR_CuNiOOH_001')
        shutil.copy(src_contcar, 'CONTCAR')

        # Load structure
        atoms = load_structure()

        # Verify it loaded from CONTCAR
        assert atoms is not None
        assert len(atoms) > 0

        # Verify by comparing to direct CONTCAR read
        expected = read('CONTCAR', format='vasp')
        assert len(atoms) == len(expected)
        assert atoms.get_chemical_formula() == expected.get_chemical_formula()

    def test_fallback_to_initial_pattern(self):
        """Test fallback to initial pattern when neither OUTCAR nor CONTCAR exist (priority 3)."""
        # Copy test POSCAR to temp directory
        src_poscar = os.path.join(self.original_dir, 'asetools/data/POSCAR')
        shutil.copy(src_poscar, 'POSCAR')

        # Load structure
        atoms = load_structure('POSCAR')

        # Verify it loaded from POSCAR
        assert atoms is not None
        assert len(atoms) > 0

    def test_fallback_with_glob_pattern(self):
        """Test fallback with glob pattern for initial structure."""
        # Copy test file with .vasp extension
        src_file = os.path.join(self.original_dir, 'asetools/data/POSCAR')
        shutil.copy(src_file, 'structure.vasp')

        # Load structure with pattern
        atoms = load_structure('*.vasp')

        # Verify it loaded
        assert atoms is not None
        assert len(atoms) > 0

    def test_outcar_priority_over_contcar(self):
        """Test that OUTCAR takes priority when both OUTCAR and CONTCAR exist."""
        # Copy both OUTCAR and CONTCAR
        src_outcar = os.path.join(self.original_dir, 'asetools/data/OUTCAR')
        src_contcar = os.path.join(self.original_dir, 'asetools/data/CONTCAR_CuNiOOH_001')
        shutil.copy(src_outcar, 'OUTCAR')
        shutil.copy(src_contcar, 'CONTCAR')

        # Load structure
        atoms = load_structure()

        # Verify it loaded from OUTCAR (not CONTCAR)
        expected_outcar = read('OUTCAR', format='vasp-out', index=-1)
        expected_contcar = read('CONTCAR', format='vasp')

        # Should match OUTCAR
        assert len(atoms) == len(expected_outcar)
        assert atoms.get_chemical_formula() == expected_outcar.get_chemical_formula()

    def test_corrupted_outcar_fallback_to_contcar(self):
        """Test fallback to CONTCAR when OUTCAR is corrupted."""
        # Create a corrupted OUTCAR file
        with open('OUTCAR', 'w') as f:
            f.write('CORRUPTED DATA\n')

        # Copy valid CONTCAR
        src_contcar = os.path.join(self.original_dir, 'asetools/data/CONTCAR_CuNiOOH_001')
        shutil.copy(src_contcar, 'CONTCAR')

        # Load structure - should fall back to CONTCAR
        atoms = load_structure()

        # Verify it loaded from CONTCAR
        assert atoms is not None
        assert len(atoms) > 0
        expected = read('CONTCAR', format='vasp')
        assert len(atoms) == len(expected)

    def test_corrupted_outcar_and_contcar_fallback_to_poscar(self):
        """Test fallback chain when both OUTCAR and CONTCAR are corrupted."""
        # Create corrupted files
        with open('OUTCAR', 'w') as f:
            f.write('CORRUPTED DATA\n')
        with open('CONTCAR', 'w') as f:
            f.write('CORRUPTED DATA\n')

        # Copy valid POSCAR
        src_poscar = os.path.join(self.original_dir, 'asetools/data/POSCAR')
        shutil.copy(src_poscar, 'POSCAR')

        # Load structure - should fall back to POSCAR
        atoms = load_structure()

        # Verify it loaded from POSCAR
        assert atoms is not None
        assert len(atoms) > 0

    def test_no_files_exits_with_error(self):
        """Test that function exits when no structure files exist."""
        # No files in temp directory

        # Should exit with error
        with pytest.raises(SystemExit) as exc_info:
            load_structure()

        assert exc_info.value.code == 1

    def test_pattern_not_matching_exits_with_error(self):
        """Test that function exits when pattern doesn't match any files."""
        # Create a file that doesn't match the pattern
        src_file = os.path.join(self.original_dir, 'asetools/data/POSCAR')
        shutil.copy(src_file, 'structure.xyz')

        # Should exit with error when looking for POSCAR
        with pytest.raises(SystemExit) as exc_info:
            load_structure('POSCAR')

        assert exc_info.value.code == 1


def test_load_structure_from_data_directory():
    """Integration test: Load structure from actual test data."""
    from asetools.manager.manager import load_structure

    # Save current directory
    original_dir = os.getcwd()

    # Change to data directory
    data_dir = 'asetools/data/summary_test/0_1.2'
    if os.path.exists(data_dir):
        os.chdir(data_dir)

        # This directory should have both OUTCAR and CONTCAR
        atoms = load_structure()

        # Verify structure loaded
        assert atoms is not None
        assert len(atoms) > 0

        # Return to original directory
        os.chdir(original_dir)
    else:
        # Skip if test data not available
        os.chdir(original_dir)
        pytest.skip("Test data directory not found")
