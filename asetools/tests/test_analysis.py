
def test_analysis_vasp6():
    from asetools.analysis import check_energy_and_maxforce, check_outcar_convergence

    outcar = 'asetools/data/OUTCAR_vasp6'
    
    conv, vasp = check_outcar_convergence(outcar, verbose=True)
    assert conv
    assert vasp

    energy, maxforce = check_energy_and_maxforce(outcar, magmom=False, verbose=False)
    print(vasp, energy, maxforce)
    assert energy
    assert maxforce

def test_analysis_vasp5():
    from asetools.analysis import check_energy_and_maxforce, check_outcar_convergence

    outcar = 'asetools/data/OUTCAR'
    
    conv, vasp = check_outcar_convergence(outcar, verbose=True)
    assert conv
    assert vasp

    energy, maxforce = check_energy_and_maxforce(outcar, magmom=False, verbose=False)
    print(vasp, energy, maxforce)
    assert energy
    assert maxforce