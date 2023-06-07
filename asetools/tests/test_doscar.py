
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