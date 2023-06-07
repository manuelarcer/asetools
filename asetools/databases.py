from asetools.analysis_copy2 import check_outcar_convergence
from ase.io import read
import pandas as pd
    
def check_if_exists_in_db(db, atoms):
    
    forces = atoms.get_forces(apply_constraint=False)
    indb = False ; index = None
    if len(db) > 0:
        for row in db.select():
            db_forces = row.forces
            if len(atoms) == row.natoms:
                if ( db_forces == forces ).all():
                    indb = True
                    index = row.id
    return indb, index

def add_config_to_db(db, outcar, idname=None, update=False):
    # db: database file, existent or new
    # outcar: path/to/outcar/OUTCAR
    # idname: Assigned name to the config
    # update: wheter to update an existing configuration within the DB
     
    converged = check_outcar_convergence(outcar)
    if converged:
        if idname == None:
            idname = outcar.split('/')[-2]
        atoms = read(outcar, format='vasp-out', index=-1)
        indb, index = check_if_exists_in_db(db, atoms)
        if not indb:
            print('New configuration, added to the DB')
            return db.write(atoms, name=idname)
        elif indb and update:
            print('Config ALREADY in DB, UPDATING...')
            return db.update(index, atoms=atoms, name=idname)
        else:
            return print('Config ALREADY in DB, skipped....')
        
def db_to_pandas(db, columns=['name', 'id', 'energy', 'free_energy', 'magmom']):
    # db: must be a loaded ASE-DB
    # columns; columns to extract from the DB to the DataFrame
    
    dic = {}
    for c in columns:
        dic[c] = []

    for i, row in enumerate(db.select(columns='all')):
        for c in columns:
            dic[c].append(row[c])        
    df = pd.DataFrame.from_dict(dic)
    return df
