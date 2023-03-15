from ase.io import read

def check_outcar_convergence(outcar):
    out = open(outcar, 'r')
    lines = out.readlines()
    ibrion = None
    nsw = None
    opt = False
    convergence = False
    for line in lines:
        if 'IBRION' in line:
            ibrion = int( line.split()[2] )
            if ibrion == 1 or ibrion == 2:
                opt = True
            else:
                print(f'IBRION --> {ibrion}, not an OPTIMIZATION job')
                opt = False
        elif 'NSW' in line:
            nsw = int( line.split()[2] )
            if nsw > 0:
                opt = True
        elif 'reached required accuracy - stopping' in line and opt:
            convergence = True

    print(f'IBRION --> {ibrion}, NSW --> {nsw}')    
    
    if ibrion == 1 or ibrion == 2:
        if nsw > 0 and convergence:
            print('Optimization Job --> CONVERGED')
        elif nsw > 0 and not convergence:
            print('Optimization Job --> *NOT* converged')
            
    return convergence
    
def check_if_exists_in_db(db, atoms):
    
    forces = atoms.get_forces(apply_constraint=False)
    indb = False ; index = None
    if len(db) > 0:
        for row in db.select():
            db_forces = row.forces
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