#!/usr/bin/env python
"""
Bader Charge Analysis Tool

This script calculates Bader charges for atoms in a VASP calculation by:
1. Running Bader analysis if ACF.dat doesn't exist
2. Extracting ZVAL (valence electrons) from OUTCAR
3. Calculating final charges as: Charge = ZVAL - Bader_charge
4. Outputting results in a formatted table

Requirements:
- VASP output files: CONTCAR, OUTCAR, AECCAR0, AECCAR2
- Bader analysis tools: chgsum.pl, bader executable
"""

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
from ase import Atoms
from ase.io import read


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class BaderChargeAnalyzer:
    """
    A class to handle Bader charge analysis for VASP calculations.
    
    This class provides methods to:
    - Extract valence electron information from OUTCAR
    - Run Bader analysis if needed
    - Calculate final atomic charges
    - Export results in various formats
    """
    
    def __init__(self, structure_file: str = 'CONTCAR', 
                 outcar_file: str = 'OUTCAR',
                 acf_file: str = 'ACF.dat'):
        """
        Initialize the Bader charge analyzer.
        
        Args:
            structure_file: Path to structure file (CONTCAR/POSCAR)
            outcar_file: Path to OUTCAR file
            acf_file: Path to ACF.dat file (Bader output)
        """
        self.structure_file = Path(structure_file)
        self.outcar_file = Path(outcar_file)
        self.acf_file = Path(acf_file)
        
        # Initialize data containers
        self.atoms: Optional[Atoms] = None
        self.zval_dict: Dict[str, float] = {}
        self.bader_charges: np.ndarray = np.array([])
        self.final_charges: np.ndarray = np.array([])
        
    def load_structure(self) -> None:
        """Load atomic structure from file."""
        try:
            self.atoms = read(str(self.structure_file), format='vasp')
            logger.info(f"Loaded structure with {len(self.atoms)} atoms from {self.structure_file}")
        except Exception as e:
            logger.error(f"Failed to load structure from {self.structure_file}: {e}")
            raise
    
    def extract_zval_from_outcar(self) -> Dict[str, float]:
        """
        Extract valence electron information (ZVAL) from OUTCAR file.
        
        Returns:
            Dictionary mapping element symbols to their valence electron counts
            
        Raises:
            FileNotFoundError: If OUTCAR file doesn't exist
            ValueError: If ZVAL information cannot be extracted
        """
        if not self.outcar_file.exists():
            raise FileNotFoundError(f"OUTCAR file not found: {self.outcar_file}")
        
        symbols = self.atoms.get_chemical_symbols()
        unique_elements = list(dict.fromkeys(symbols))  # Preserve order
        n_unique = len(unique_elements)
        
        elements = []
        zval_values = []
        
        try:
            with open(self.outcar_file, 'r') as f:
                lines = f.readlines()
                
            for i, line in enumerate(lines):
                # Extract element types from POTCAR information
                if 'POTCAR:' in line and len(elements) < n_unique:
                    element = line.split()[2].split('_')[0]
                    elements.append(element)
                
                # Extract valence electron counts
                if 'Ionic Valenz' in line:
                    valence_line = lines[i + 1]
                    values = valence_line.split()[2:]
                    zval_values = [float(val) for val in values]
                    break
            
            if len(elements) != len(zval_values):
                raise ValueError(f"Mismatch between elements ({len(elements)}) and ZVAL values ({len(zval_values)})")
            
            self.zval_dict = dict(zip(elements, zval_values))
            logger.info(f"Extracted ZVAL for elements: {self.zval_dict}")
            return self.zval_dict
            
        except Exception as e:
            logger.error(f"Failed to extract ZVAL from OUTCAR: {e}")
            raise
    
    def run_bader_analysis(self) -> None:
        """
        Run Bader charge analysis if ACF.dat doesn't exist.
        
        This method executes the following commands:
        1. chgsum.pl AECCAR0 AECCAR2  (creates CHGCAR_sum)
        2. bader -vac off CHGCAR -ref CHGCAR_sum  (creates ACF.dat)
        
        Raises:
            RuntimeError: If Bader analysis fails
        """
        if self.acf_file.exists():
            logger.info(f"ACF.dat already exists, skipping Bader analysis")
            return
        
        logger.info("ACF.dat not found, running Bader charge analysis...")
        
        # Check required files exist
        required_files = ['AECCAR0', 'AECCAR2', 'CHGCAR']
        missing_files = [f for f in required_files if not Path(f).exists()]
        if missing_files:
            raise FileNotFoundError(f"Required files missing for Bader analysis: {missing_files}")
        
        try:
            # Generate charge density sum
            logger.info("Running chgsum.pl to create CHGCAR_sum...")
            ret1 = os.system('chgsum.pl AECCAR0 AECCAR2')
            if ret1 != 0:
                raise RuntimeError("chgsum.pl failed")
            
            # Run Bader analysis
            logger.info("Running Bader analysis...")
            ret2 = os.system('bader -vac off CHGCAR -ref CHGCAR_sum')
            if ret2 != 0:
                raise RuntimeError("Bader analysis failed")
            
            if not self.acf_file.exists():
                raise RuntimeError("Bader analysis completed but ACF.dat not created")
                
            logger.info("Bader analysis completed successfully")
            
        except Exception as e:
            logger.error(f"Bader analysis failed: {e}")
            raise
    
    def extract_bader_charges(self) -> np.ndarray:
        """
        Extract Bader charges from ACF.dat file.
        
        Returns:
            Array of Bader charges for each atom
            
        Raises:
            FileNotFoundError: If ACF.dat doesn't exist
            ValueError: If charge data cannot be extracted
        """
        if not self.acf_file.exists():
            raise FileNotFoundError(f"ACF.dat file not found: {self.acf_file}")
        
        try:
            with open(self.acf_file, 'r') as f:
                lines = f.readlines()
            
            n_atoms = len(self.atoms)
            charges = np.zeros(n_atoms)
            
            # Extract charges from lines 2 to n_atoms+2 (0-indexed)
            for i, line in enumerate(lines[2:2+n_atoms]):
                charges[i] = float(line.split()[4])
            
            self.bader_charges = charges
            logger.info(f"Extracted Bader charges for {len(charges)} atoms")
            return charges
            
        except Exception as e:
            logger.error(f"Failed to extract Bader charges: {e}")
            raise
    
    def calculate_final_charges(self) -> np.ndarray:
        """
        Calculate final atomic charges as: Charge = ZVAL - Bader_charge
        
        Returns:
            Array of final charges for each atom
        """
        if len(self.bader_charges) == 0:
            raise ValueError("Bader charges not loaded. Run extract_bader_charges() first.")
        
        symbols = self.atoms.get_chemical_symbols()
        final_charges = np.zeros(len(symbols))
        
        for i, symbol in enumerate(symbols):
            if symbol not in self.zval_dict:
                raise ValueError(f"ZVAL not found for element: {symbol}")
            final_charges[i] = self.zval_dict[symbol] - self.bader_charges[i]
        
        self.final_charges = final_charges
        logger.info("Calculated final charges for all atoms")
        return final_charges
    
    def write_results(self, output_file: str = 'bader_charges.txt', 
                     format_type: str = 'table') -> None:
        """
        Write charge analysis results to file.
        
        Args:
            output_file: Output filename
            format_type: Output format ('table', 'csv', 'json')
        """
        if len(self.final_charges) == 0:
            raise ValueError("Final charges not calculated. Run calculate_final_charges() first.")
        
        symbols = self.atoms.get_chemical_symbols()
        positions = self.atoms.get_positions()
        
        output_path = Path(output_file)
        
        try:
            if format_type == 'table':
                self._write_table_format(output_path, symbols, positions)
            elif format_type == 'csv':
                self._write_csv_format(output_path, symbols, positions)
            elif format_type == 'json':
                self._write_json_format(output_path, symbols, positions)
            else:
                raise ValueError(f"Unsupported format: {format_type}")
                
            logger.info(f"Results written to {output_path}")
            
        except Exception as e:
            logger.error(f"Failed to write results: {e}")
            raise
    
    def _write_table_format(self, output_path: Path, symbols: List[str], 
                           positions: np.ndarray) -> None:
        """Write results in formatted table format."""
        with open(output_path, 'w') as f:
            f.write(f"{'#':<3} {'Symbol':<8} {'X':<12} {'Y':<12} {'Z':<12} "
                   f"{'Bader_Chg':<12} {'Final_Chg':<12}\n")
            f.write("-" * 75 + "\n")
            
            for i, symbol in enumerate(symbols):
                f.write(f"{i:<3} {symbol:<8} {positions[i,0]:<12.6f} "
                       f"{positions[i,1]:<12.6f} {positions[i,2]:<12.6f} "
                       f"{self.bader_charges[i]:<12.6f} {self.final_charges[i]:<12.6f}\n")
    
    def _write_csv_format(self, output_path: Path, symbols: List[str], 
                         positions: np.ndarray) -> None:
        """Write results in CSV format."""
        import csv
        
        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Index', 'Symbol', 'X', 'Y', 'Z', 'Bader_Charge', 'Final_Charge'])
            
            for i, symbol in enumerate(symbols):
                writer.writerow([i, symbol, positions[i,0], positions[i,1], 
                               positions[i,2], self.bader_charges[i], self.final_charges[i]])
    
    def _write_json_format(self, output_path: Path, symbols: List[str], 
                          positions: np.ndarray) -> None:
        """Write results in JSON format."""
        import json
        
        results = {
            'zval_dict': self.zval_dict,
            'atoms': []
        }
        
        for i, symbol in enumerate(symbols):
            atom_data = {
                'index': i,
                'symbol': symbol,
                'position': positions[i].tolist(),
                'bader_charge': float(self.bader_charges[i]),
                'final_charge': float(self.final_charges[i])
            }
            results['atoms'].append(atom_data)
        
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
    
    def print_summary(self) -> None:
        """Print a summary of the charge analysis."""
        if len(self.final_charges) == 0:
            logger.warning("No charge data available for summary")
            return
        
        symbols = self.atoms.get_chemical_symbols()
        unique_symbols = sorted(set(symbols))
        
        print("\n" + "="*60)
        print("BADER CHARGE ANALYSIS SUMMARY")
        print("="*60)
        print(f"Total atoms: {len(self.atoms)}")
        print(f"Elements: {', '.join(unique_symbols)}")
        print(f"ZVAL values: {self.zval_dict}")
        print()
        
        # Statistics by element
        for element in unique_symbols:
            element_indices = [i for i, s in enumerate(symbols) if s == element]
            element_charges = self.final_charges[element_indices]
            
            print(f"{element:>4}: {len(element_indices):3d} atoms, "
                  f"charge range: {element_charges.min():8.4f} to {element_charges.max():8.4f}, "
                  f"average: {element_charges.mean():8.4f}")
        
        print(f"\nTotal charge: {self.final_charges.sum():8.4f}")
        print("="*60)


def main():
    """Main entry point for the command-line interface."""
    parser = argparse.ArgumentParser(
        description='Calculate Bader charges for VASP calculations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                           # Use default files
  %(prog)s -s POSCAR -o OUTCAR      # Specify input files
  %(prog)s --format csv -o charges.csv  # Output as CSV
  %(prog)s --skip-bader              # Skip Bader analysis (ACF.dat exists)
        """
    )
    
    parser.add_argument('-s', '--structure', default='CONTCAR',
                       help='Structure file (default: CONTCAR)')
    parser.add_argument('-o', '--outcar', default='OUTCAR',
                       help='OUTCAR file (default: OUTCAR)')
    parser.add_argument('-a', '--acf', default='ACF.dat',
                       help='ACF.dat file (default: ACF.dat)')
    parser.add_argument('--output', default='bader_charges.txt',
                       help='Output file (default: bader_charges.txt)')
    parser.add_argument('--format', choices=['table', 'csv', 'json'], default='table',
                       help='Output format (default: table)')
    parser.add_argument('--skip-bader', action='store_true',
                       help='Skip Bader analysis (assume ACF.dat exists)')
    parser.add_argument('--quiet', '-q', action='store_true',
                       help='Suppress output messages')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Configure logging level
    if args.quiet:
        logging.getLogger().setLevel(logging.ERROR)
    elif args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        # Initialize analyzer
        analyzer = BaderChargeAnalyzer(
            structure_file=args.structure,
            outcar_file=args.outcar,
            acf_file=args.acf
        )
        
        # Load structure
        analyzer.load_structure()
        
        # Extract ZVAL information
        analyzer.extract_zval_from_outcar()
        
        # Run Bader analysis if needed
        if not args.skip_bader:
            analyzer.run_bader_analysis()
        
        # Extract Bader charges
        analyzer.extract_bader_charges()
        
        # Calculate final charges
        analyzer.calculate_final_charges()
        
        # Write results
        analyzer.write_results(output_file=args.output, format_type=args.format)
        
        # Print summary
        if not args.quiet:
            analyzer.print_summary()
        
        logger.info("Bader charge analysis completed successfully!")
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()