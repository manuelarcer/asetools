
import pytest
import numpy as np

def test_doscar_module():

    from asetools.doscar_analysis import extract_dos, extract_pdos_perstate, extract_pdos_perorbital

    doscar = 'asetools/data/DOSCAR_fe3o4_Feoct_2x2.doscar'
    dos = extract_dos(doscar)
    #print(dos)
    assert dos
    perstate = extract_pdos_perstate(dos, [0,1], ['s_states', 'p_states', 'd_states'])
    #print(perstate)
    assert perstate
    perorbital = extract_pdos_perorbital(dos, [0,1], 'all-d')
    #print(perorbital)
    assert perorbital

def test_dos_class_basic():
    """Test basic DOS class functionality."""
    from asetools.doscar_analysis import DOS
    
    doscar = 'asetools/data/DOSCAR_fe3o4_Feoct_2x2.doscar'
    dos = DOS(doscar)
    
    # Test basic properties
    assert dos.has_partial_dos == True
    assert dos.natoms > 0
    assert dos.fermi_energy is not None
    assert len(dos.energy) > 0
    assert len(dos.dos_up) > 0
    assert len(dos.dos_down) > 0

def test_band_center_d_orbitals():
    """Test d-band center calculation."""
    from asetools.doscar_analysis import DOS
    
    doscar = 'asetools/data/DOSCAR_fe3o4_Feoct_2x2.doscar'
    dos = DOS(doscar)
    
    # Test d-band center for first atom
    d_center = dos.calculate_band_center([0], orbitals='all-d')
    assert isinstance(d_center, float)
    assert np.isfinite(d_center)
    
    # Test different d-orbital subsets
    t2g_center = dos.calculate_band_center([0], orbitals='t2g')
    eg_center = dos.calculate_band_center([0], orbitals='eg')
    assert isinstance(t2g_center, float)
    assert isinstance(eg_center, float)
    assert np.isfinite(t2g_center)
    assert np.isfinite(eg_center)

def test_band_center_p_orbitals():
    """Test p-band center calculation."""
    from asetools.doscar_analysis import DOS
    
    doscar = 'asetools/data/DOSCAR_fe3o4_Feoct_2x2.doscar'
    dos = DOS(doscar)
    
    # Test p-band center
    p_center = dos.calculate_band_center([0], orbitals='all-p')
    assert isinstance(p_center, float)
    assert np.isfinite(p_center)

def test_band_center_states():
    """Test band center calculation using states."""
    from asetools.doscar_analysis import DOS
    
    doscar = 'asetools/data/DOSCAR_fe3o4_Feoct_2x2.doscar'
    dos = DOS(doscar)
    
    # Test using states instead of orbitals
    d_center = dos.calculate_band_center([0], states=['d_states'])
    p_center = dos.calculate_band_center([0], states=['p_states'])
    
    assert isinstance(d_center, float)
    assert isinstance(p_center, float)
    assert np.isfinite(d_center)
    assert np.isfinite(p_center)

def test_band_center_spin_treatment():
    """Test different spin treatments for band center calculation."""
    from asetools.doscar_analysis import DOS
    
    doscar = 'asetools/data/DOSCAR_fe3o4_Feoct_2x2.doscar'
    dos = DOS(doscar)
    
    # Test different spin treatments
    combined = dos.calculate_band_center([0], orbitals='all-d', spin_treatment='combined')
    up_only = dos.calculate_band_center([0], orbitals='all-d', spin_treatment='up')
    down_only = dos.calculate_band_center([0], orbitals='all-d', spin_treatment='down')
    separate = dos.calculate_band_center([0], orbitals='all-d', spin_treatment='separate')
    
    assert isinstance(combined, float)
    assert isinstance(up_only, float)
    assert isinstance(down_only, float)
    assert isinstance(separate, dict)
    
    # Check separate result structure
    assert 'up' in separate
    assert 'down' in separate
    assert isinstance(separate['up'], float)
    assert isinstance(separate['down'], float)

def test_band_center_energy_range():
    """Test band center calculation with energy range filtering."""
    from asetools.doscar_analysis import DOS
    
    doscar = 'asetools/data/DOSCAR_fe3o4_Feoct_2x2.doscar'
    dos = DOS(doscar)
    
    # Test with energy range
    full_range = dos.calculate_band_center([0], orbitals='all-d')
    limited_range = dos.calculate_band_center([0], orbitals='all-d', energy_range=(-5, 2))
    
    assert isinstance(full_range, float)
    assert isinstance(limited_range, float)
    assert np.isfinite(full_range)
    assert np.isfinite(limited_range)
    
    # They should be different (assuming DOS extends beyond the range)
    assert full_range != limited_range

def test_band_center_multiple_atoms():
    """Test band center calculation for multiple atoms."""
    from asetools.doscar_analysis import DOS
    
    doscar = 'asetools/data/DOSCAR_fe3o4_Feoct_2x2.doscar'
    dos = DOS(doscar)
    
    # Test with multiple atoms
    single_atom = dos.calculate_band_center([0], orbitals='all-d')
    multiple_atoms = dos.calculate_band_center([0, 1], orbitals='all-d')
    
    assert isinstance(single_atom, float)
    assert isinstance(multiple_atoms, float)
    assert np.isfinite(single_atom)
    assert np.isfinite(multiple_atoms)

def test_band_center_error_handling():
    """Test error handling in band center calculations."""
    from asetools.doscar_analysis import DOS
    
    doscar = 'asetools/data/DOSCAR_fe3o4_Feoct_2x2.doscar'
    dos = DOS(doscar)
    
    # Test parameter validation
    with pytest.raises(ValueError, match="Cannot specify both states and orbitals"):
        dos.calculate_band_center([0], orbitals='all-d', states=['d_states'])
    
    with pytest.raises(ValueError, match="Must specify either states or orbitals"):
        dos.calculate_band_center([0])
    
    with pytest.raises(ValueError, match="spin_treatment must be"):
        dos.calculate_band_center([0], orbitals='all-d', spin_treatment='invalid')

def test_band_center_legacy_function():
    """Test the legacy function wrapper."""
    from asetools.doscar_analysis import calculate_band_center
    
    doscar = 'asetools/data/DOSCAR_fe3o4_Feoct_2x2.doscar'
    
    # Test legacy function
    d_center = calculate_band_center(doscar, [0], orbitals='all-d')
    
    assert isinstance(d_center, float)
    assert np.isfinite(d_center)
    
    # Test with different parameters
    separate_spins = calculate_band_center(doscar, [0], orbitals='all-d', spin_treatment='separate')
    assert isinstance(separate_spins, dict)
    assert 'up' in separate_spins
    assert 'down' in separate_spins