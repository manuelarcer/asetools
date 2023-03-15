# TODO
# File containing the functions to analyze OUTCARs and VASP runs

from ase.io import read
import pandas as pd
from vasprun import vasprun


def check_convergence(xml='vasprun.xml'):
    vasp = vasprun('vasprun.xml')
    
