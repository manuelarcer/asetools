import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Dict, Union, Optional, Tuple

class DOS:
    """Class for handling and analyzing VASP DOSCAR files.
    
    This class provides a more modular approach to DOS analysis,
    storing DOS data and providing methods for extraction and visualization.
    """
    
    def __init__(self, doscarfile: str):
        """Initialize DOS object from DOSCAR file.
        
        Args:
            doscarfile: Path to the DOSCAR file
        """
        self.doscarfile = doscarfile
        self._load_doscar()
    
    def _load_doscar(self):
        """Load and parse DOSCAR file."""
        with open(self.doscarfile, 'r') as doscar:
            lines = doscar.readlines()
        
        # Parse header information
        self.natoms = int(lines[0].split()[0])
        typedos = int(lines[0].split()[2])
        self.has_partial_dos = typedos == 1
        
        # Get Fermi energy
        if len(lines) > 5 and len(lines[5].split()) >= 4:
            self.fermi_energy = float(lines[5].split()[3])
        else:
            raise ValueError(f"Invalid DOSCAR format: cannot extract Fermi energy from line 5")
        
        # Initialize data structure
        self.data = {'energy': [], 'DOSup': [], 'DOSdown': []}
        
        # Initialize per-atom data if partial DOS is available
        if self.has_partial_dos:
            for i in range(self.natoms):
                self.data[f'at-{i}'] = {
                    'energy': [], 's+': [], 's-': [], 'py+': [], 'py-': [],
                    'pz+': [], 'pz-': [], 'px+': [], 'px-': [],
                    'dxy+': [], 'dxy-': [], 'dyz+': [], 'dyz-': [],
                    'dz2+': [], 'dz2-': [], 'dxz+': [], 'dxz-': [],
                    'dx2+': [], 'dx2-': []
                }
        
        # Parse DOS data
        self._parse_dos_data(lines)
        
        # Convert to numpy arrays
        self._convert_to_arrays()
    
    def _parse_dos_data(self, lines: List[str]):
        """Parse DOS data from DOSCAR lines."""
        repline = lines[5]
        goDOS = False
        gopDOS = False
        count = -1
        
        for i, line in enumerate(lines):
            if goDOS:
                self.data['energy'].append(float(line.split()[0]) - self.fermi_energy)
                self.data['DOSup'].append(float(line.split()[1]))
                self.data['DOSdown'].append(-float(line.split()[2]))
            
            if gopDOS and line != repline:
                atom_key = f'at-{count}'
                for k, key in enumerate(self.data[atom_key]):
                    if k == 0:
                        self.data[atom_key]['energy'].append(
                            float(line.split()[k]) - self.fermi_energy
                        )
                    else:
                        self.data[atom_key][key].append(float(line.split()[k]))
            
            if i == 5:
                goDOS = True
            
            if i > 5 and line == repline and self.has_partial_dos:
                goDOS = False
                gopDOS = True
                count += 1
    
    def _convert_to_arrays(self):
        """Convert lists to numpy arrays."""
        for key, val in self.data.items():
            if isinstance(val, list):
                self.data[key] = np.array(val)
            elif isinstance(val, dict):
                for subkey, subval in val.items():
                    self.data[key][subkey] = np.array(subval)
    
    @property
    def energy(self) -> np.ndarray:
        """Energy grid (relative to Fermi energy)."""
        return self.data['energy']
    
    @property
    def dos_up(self) -> np.ndarray:
        """Spin-up DOS."""
        return self.data['DOSup']
    
    @property
    def dos_down(self) -> np.ndarray:
        """Spin-down DOS."""
        return self.data['DOSdown']
    
    @property
    def total_dos(self) -> np.ndarray:
        """Total DOS (sum of spin-up and absolute spin-down)."""
        return self.dos_up + np.abs(self.dos_down)
    
    def get_pdos_by_states(self, atoms: List[int], states: List[str]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Extract projected DOS by orbital states.
        
        Args:
            atoms: List of atom indices
            states: List of states ('s_states', 'p_states', 'd_states')
            
        Returns:
            Tuple of (energy, dos_up, dos_down)
        """
        if not self.has_partial_dos:
            raise ValueError("No partial DOS data available")
        
        state_orbitals = {
            's_states': {'s+': 1, 's-': 2},
            'p_states': {'py+': 3, 'py-': 4, 'pz+': 5, 'pz-': 6, 'px+': 7, 'px-': 8},
            'd_states': {
                'dxy+': 9, 'dxy-': 10, 'dyz+': 11, 'dyz-': 12,
                'dz2+': 13, 'dz2-': 14, 'dxz+': 15, 'dxz-': 16,
                'dx2+': 17, 'dx2-': 18
            }
        }
        
        energy = self.data['at-0']['energy']
        sum_plus = np.zeros(len(energy))
        sum_minus = np.zeros(len(energy))
        
        for at in atoms:
            for state in states:
                for orbital in state_orbitals[state]:
                    if orbital[-1] == '+':
                        sum_plus += self.data[f'at-{at}'][orbital]
                    elif orbital[-1] == '-':
                        sum_minus -= self.data[f'at-{at}'][orbital]
        
        return energy, sum_plus, sum_minus
    
    def get_pdos_by_orbitals(self, atoms: List[int], orbitals: Union[List[str], str]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Extract projected DOS by specific orbitals.
        
        Args:
            atoms: List of atom indices
            orbitals: List of orbital names or shorthand ('all-s', 'all-p', 'all-d', 't2g', 'eg')
            
        Returns:
            Tuple of (energy, dos_up, dos_down)
        """
        if not self.has_partial_dos:
            raise ValueError("No partial DOS data available")
        
        # Handle shorthand orbital specifications
        if orbitals == 'all-s':
            orbitals = ['s+', 's-']
        elif orbitals == 'all-p':
            orbitals = ['py+', 'py-', 'pz+', 'pz-', 'px+', 'px-']
        elif orbitals == 'all-d':
            orbitals = ['dxy+', 'dxy-', 'dyz+', 'dyz-', 'dxz+', 'dxz-', 'dz2+', 'dz2-', 'dx2+', 'dx2-']
        elif orbitals == 'all':
            orbitals = [
                's+', 's-', 'py+', 'py-', 'pz+', 'pz-', 'px+', 'px-',
                'dxy+', 'dxy-', 'dyz+', 'dyz-', 'dxz+', 'dxz-', 'dz2+', 'dz2-', 'dx2+', 'dx2-'
            ]
        elif orbitals == 't2g':
            orbitals = ['dxy+', 'dxy-', 'dyz+', 'dyz-', 'dxz+', 'dxz-']
        elif orbitals == 'eg':
            orbitals = ['dz2+', 'dz2-', 'dx2+', 'dx2-']
        
        energy = self.data['at-0']['energy']
        sum_plus = np.zeros(len(energy))
        sum_minus = np.zeros(len(energy))
        
        for at in atoms:
            for orbital in orbitals:
                if orbital[-1] == '+':
                    sum_plus += self.data[f'at-{at}'][orbital]
                elif orbital[-1] == '-':
                    sum_minus -= self.data[f'at-{at}'][orbital]
        
        return energy, sum_plus, sum_minus
    
    def plot_total_dos(self, ax: Optional[plt.Axes] = None, **kwargs) -> plt.Axes:
        """Plot total DOS.
        
        Args:
            ax: Matplotlib axes object (optional)
            **kwargs: Additional arguments passed to plot
            
        Returns:
            Matplotlib axes object
        """
        if ax is None:
            _, ax = plt.subplots(figsize=(8, 6))
        
        ax.plot(self.energy, self.dos_up, label='Spin Up', **kwargs)
        ax.plot(self.energy, self.dos_down, label='Spin Down', **kwargs)
        ax.axhline(y=0, color='k', linestyle='-', alpha=0.3)
        ax.axvline(x=0, color='k', linestyle='--', alpha=0.5, label='Fermi Level')
        
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('DOS (states/eV)')
        ax.set_title('Total Density of States')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        return ax
    
    def plot_pdos_by_states(self, atoms: List[int], states: List[str], 
                           ax: Optional[plt.Axes] = None, **kwargs) -> plt.Axes:
        """Plot projected DOS by orbital states.
        
        Args:
            atoms: List of atom indices
            states: List of states ('s_states', 'p_states', 'd_states')
            ax: Matplotlib axes object (optional)
            **kwargs: Additional arguments passed to plot
            
        Returns:
            Matplotlib axes object
        """
        if ax is None:
            _, ax = plt.subplots(figsize=(8, 6))
        
        energy, dos_up, dos_down = self.get_pdos_by_states(atoms, states)
        
        label = f"Atoms {atoms}, States {states}"
        ax.plot(energy, dos_up, label=f'{label} (Up)', **kwargs)
        ax.plot(energy, dos_down, label=f'{label} (Down)', **kwargs)
        ax.axhline(y=0, color='k', linestyle='-', alpha=0.3)
        ax.axvline(x=0, color='k', linestyle='--', alpha=0.5, label='Fermi Level')
        
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('PDOS (states/eV)')
        ax.set_title('Projected Density of States')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        return ax
    
    def plot_pdos_by_orbitals(self, atoms: List[int], orbitals: Union[List[str], str], 
                             ax: Optional[plt.Axes] = None, **kwargs) -> plt.Axes:
        """Plot projected DOS by specific orbitals.
        
        Args:
            atoms: List of atom indices
            orbitals: List of orbital names or shorthand
            ax: Matplotlib axes object (optional)
            **kwargs: Additional arguments passed to plot
            
        Returns:
            Matplotlib axes object
        """
        if ax is None:
            _, ax = plt.subplots(figsize=(8, 6))
        
        energy, dos_up, dos_down = self.get_pdos_by_orbitals(atoms, orbitals)
        
        label = f"Atoms {atoms}, Orbitals {orbitals}"
        ax.plot(energy, dos_up, label=f'{label} (Up)', **kwargs)
        ax.plot(energy, dos_down, label=f'{label} (Down)', **kwargs)
        ax.axhline(y=0, color='k', linestyle='-', alpha=0.3)
        ax.axvline(x=0, color='k', linestyle='--', alpha=0.5, label='Fermi Level')
        
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('PDOS (states/eV)')
        ax.set_title('Projected Density of States')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        return ax
    
    def to_dataframe(self) -> pd.DataFrame:
        """Convert DOS data to pandas DataFrame.
        
        Returns:
            DataFrame with energy and DOS data
        """
        df_data = {
            'energy': self.energy,
            'dos_up': self.dos_up,
            'dos_down': self.dos_down,
            'total_dos': self.total_dos
        }
        
        return pd.DataFrame(df_data)


# Backward compatibility functions
def extract_dos(doscarfile: str) -> Dict:
    """Extract DOS data from DOSCAR file (legacy function).
    
    Args:
        doscarfile: Path to DOSCAR file
        
    Returns:
        Dictionary with DOS data
    """
    dos = DOS(doscarfile)
    return dos.data

def extract_fermi_e(doscarfile):
    try:
        with open(doscarfile, 'r') as file:
            lines = file.readlines()
            fermie = float( lines[5].split()[3] )
        return fermie
    except FileNotFoundError:
        print(f"'{doscarfile}' not found.")
    except IndexError:
        print("Unexpected file format. Couldn't extract Fermi energy.")
    except ValueError:
        print("Error converting Fermi energy to float.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
    return None


def extract_pdos_perstate(data: Dict, atoms: List[int], states: List[str]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract projected DOS by orbital states (legacy function).
    
    Args:
        data: DOS data dictionary from extract_dos()
        atoms: List of atom indices
        states: List of states ('s_states', 'p_states', 'd_states')
        
    Returns:
        Tuple of (energy, dos_up, dos_down)
    """
    dicstates = {
        's_states': {'s+': 1, 's-': 2},
        'p_states': {'py+': 3, 'py-': 4, 'pz+': 5, 'pz-': 6, 'px+': 7, 'px-': 8},
        'd_states': {
            'dxy+': 9, 'dxy-': 10, 'dyz+': 11, 'dyz-': 12,
            'dz2+': 13, 'dz2-': 14, 'dxz+': 15, 'dxz-': 16,
            'dx2+': 17, 'dx2-': 18
        }
    }
    
    sum_plus = np.zeros(len(data['at-0']['py+']))
    sum_minus = np.zeros(len(data['at-0']['py+']))
    e = data['at-0']['energy']
    
    for at in atoms:
        for sss in states:
            for key in dicstates[sss]:
                if key[-1] == '+':
                    sum_plus += data['at-'+str(at)][key]
                elif key[-1] == '-':
                    sum_minus -= data['at-'+str(at)][key]
    
    return e, sum_plus, sum_minus

def extract_pdos_perorbital(data: Dict, atoms: List[int], orbitals: Union[List[str], str]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract projected DOS by specific orbitals (legacy function).
    
    Args:
        data: DOS data dictionary from extract_dos()
        atoms: List of atom indices
        orbitals: List of orbital names or shorthand
        
    Returns:
        Tuple of (energy, dos_up, dos_down)
    """
    if orbitals == 'all-s':
        orbitals = ['s+', 's-']
    elif orbitals == 'all-p':
        orbitals = ['py+', 'py-', 'pz+', 'pz-', 'px+', 'px-']
    elif orbitals == 'all-d':
        orbitals = ['dxy+', 'dxy-', 'dyz+', 'dyz-', 'dxz+', 'dxz-', 'dz2+', 'dz2-', 'dx2+', 'dx2-']
    elif orbitals == 'all':
        orbitals = [
            's+', 's-', 'py+', 'py-', 'pz+', 'pz-', 'px+', 'px-',
            'dxy+', 'dxy-', 'dyz+', 'dyz-', 'dxz+', 'dxz-', 'dz2+', 'dz2-', 'dx2+', 'dx2-'
        ]
    elif orbitals == 't2g':
        orbitals = ['dxy+', 'dxy-', 'dyz+', 'dyz-', 'dxz+', 'dxz-']
    elif orbitals == 'eg':
        orbitals = ['dz2+', 'dz2-', 'dx2+', 'dx2-']
    
    sum_plus = np.zeros(len(data['at-0']['py+']))
    sum_minus = np.zeros(len(data['at-0']['py+']))
    e = data['at-0']['energy']
    
    for at in atoms:
        for sss in orbitals:
            if sss[-1] == '+':
                sum_plus += data['at-'+str(at)][sss]
            elif sss[-1] == '-':
                sum_minus -= data['at-'+str(at)][sss]
    
    return e, sum_plus, sum_minus