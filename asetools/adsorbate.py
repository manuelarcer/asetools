
import numpy as np
from ase import Atom
from ase.io import read, write
from ase.visualize import view


class SurfaceAnalyzer:
    def __init__(self, atoms):
        self.atoms = atoms
        self.surface_indices = self.find_surface_atoms()
        self.surface_neighbors = self.find_surface_neighbors()
        self.symbols = np.unique( self.atoms.get_chemical_symbols() )
    
    def find_surface_atoms(self):
        z_max = self.atoms.positions[:, 2].max()  # Maximum z-coordinate
        surface_indices = []
        for i, atom in enumerate(self.atoms):
            if atom.position[2] > z_max - 0.5:
                surface_indices.append(i)
        return surface_indices
    
    def find_surface_neighbors(self, thr=3.0):
        surface_indices = sorted(self.surface_indices)
        surface_neighbors = []
        for i in range(len(surface_indices)):
            for j in range(i+1, len(surface_indices)):
                d_ij = self.atoms.get_distance(surface_indices[i], surface_indices[j], mic=True)
                if d_ij < thr:
                    for k in range(j+1, len(surface_indices)):
                        d_ik = self.atoms.get_distance(surface_indices[i], surface_indices[k], mic=True)
                        d_jk = self.atoms.get_distance(surface_indices[j], surface_indices[k], mic=True)
                        if d_ik < thr and d_jk < thr:
                            surface_neighbors.append((surface_indices[i], surface_indices[j], surface_indices[k]))
        return surface_neighbors

    def midpoint_three_atoms(self, three):
        # threeat is in (i, j, k) form detailing the indexes of atoms
        v_ij = self.atoms.get_distance(three[0], three[1],mic=True, vector=True)
        v_ik = self.atoms.get_distance(three[0], three[2],mic=True, vector=True)
        n = np.cross(v_ij, v_ik)

        # Calculate midpoint of triangle
        midpoint = (v_ij + v_ik) / 2.0

        # In a equilaterus triangle mid point is at 2/3 of R
        # http://www.treenshop.com/Treenshop/ArticlesPages/FiguresOfInterest_Article/The%20Equilateral%20Triangle.htm 
        return self.atoms.get_positions()[three[0]] + midpoint * 2 / 3

    def add_adsorbate_to_mid(self, three, adsorbate='H', z_off=1., thr=4.4):
        # three is in (i, j, k) form detailing the indexes of atoms
        # adsorbate, The Symbol of the atom to add as adsorbate
        # z_off, the shift in Z for the adsorbate. for H use 1.0, bigger atoms require higher numbers
        # thr, This is the threshold to find neighbors around the adsorbate, (generate cluster)
        
        mid = self.midpoint_three_atoms(three)
        self.atoms.append(Atom(adsorbate, mid + [0, 0, z_off]))
        self.adsorbate = adsorbate                      # Store the Symbol of Ads
        self.adsorbateindex = len( self.atoms ) - 1     # Store index of adsorbate
        self.adsorbatecluster = self.get_cluster_around_adsorbate(thr=thr)
        self.adsneighdistances = self.adsorbate_to_cluster_neighbors(thr=thr)

    def get_cluster_around_adsorbate(self, thr=4.4):
        try:
            dist_mtx = self.atoms.get_distances(self.adsorbateindex, indices=range(len(self.atoms)), mic=True)
            cluster = self.atoms[ (dist_mtx < thr) ]
            return cluster
        except:
            print('Possibly no adsorbate defined. Use SurfaceAnalizer.add_adsorbate_to_mid')

    def adsorbate_to_cluster_neighbors(self, thr=4.4, prec=3):
        dist_dic = {}
        for s in self.symbols:
            dist_neigh = self.adsorbatecluster.get_distances( len(self.adsorbatecluster)-1, range( len(self.adsorbatecluster) ), mic=True)
            dist_neigh = dist_neigh[ (dist_neigh > 0.01) & \
                                    (np.array(self.adsorbatecluster.get_chemical_symbols()) == s)]
            dist_neigh = [round(elem, prec) for elem in dist_neigh]
            dist_dic[s] = sorted(dist_neigh)
        return dist_dic