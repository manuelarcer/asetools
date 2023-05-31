
def test_analysis():
    from asetools.analysis import check_energy_and_maxforce, check_outcar_convergence

    outcar = 'asetools/data/OUTCAR'
    
    conv = check_outcar_convergence(outcar, verbose=True)
    assert conv

    energy, maxforce = check_energy_and_maxforce(outcar, magmom=False, verbose=False)
    print(energy, maxforce)
    assert energy
    assert maxforce