#!/usr/bin/env python

from asetools.analysis import check_energy_and_maxforce, check_outcar_convergence
import sys

if len( sys.argv ) > 1:
    outcar = sys.argv[1]
else:
    outcar = 'OUTCAR'

converged = check_outcar_convergence(outcar, verbose=True)
energy, maxforce = check_energy_and_maxforce(outcar, magmom=False, verbose=False)

print('{} {} {}'.format('Converged', 'MaxForce', 'Energy'))
print(f'{converged:10s} {maxforce:10.3f} {energy:10.4f}')
