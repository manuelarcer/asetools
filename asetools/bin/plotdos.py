#!/usr/bin/env python

"""
Plot DOS/PDOS from VASP DOSCAR files

This tool provides a convenient command-line interface for visualizing
Density of States data from VASP calculations, including total DOS,
partial DOS (PDOS), and band center calculations.
"""

import argparse
import sys
import os
import matplotlib.pyplot as plt
from asetools.doscar_analysis import DOS


def parse_atoms(atom_str):
    """
    Parse atom specification string into list of integers.

    Supports:
        - Single atoms: "0" -> [0]
        - Comma-separated: "0,2,5" -> [0, 2, 5]
        - Ranges: "0-5" -> [0, 1, 2, 3, 4, 5]
        - Mixed: "0,2-4,7" -> [0, 2, 3, 4, 7]

    Args:
        atom_str: String specification of atoms

    Returns:
        List of atom indices
    """
    if not atom_str:
        return None

    atoms = []
    for part in atom_str.split(','):
        part = part.strip()
        if '-' in part:
            start, end = map(int, part.split('-'))
            atoms.extend(range(start, end + 1))
        else:
            atoms.append(int(part))

    return sorted(list(set(atoms)))  # Remove duplicates and sort


def parse_range(range_str):
    """
    Parse range string like "-5,5" into tuple of floats.

    Args:
        range_str: String like "-5,5" or "-10.5,10.5"

    Returns:
        Tuple of (min, max) floats
    """
    if not range_str:
        return None

    parts = range_str.split(',')
    if len(parts) != 2:
        raise ValueError(f"Range must be two values separated by comma: {range_str}")

    return (float(parts[0]), float(parts[1]))


def parse_colors(color_str):
    """
    Parse color specification string.

    Args:
        color_str: Comma-separated color names/codes

    Returns:
        List of color strings
    """
    if not color_str:
        return None

    return [c.strip() for c in color_str.split(',')]


def main():
    parser = argparse.ArgumentParser(
        description='Plot DOS/PDOS from VASP DOSCAR file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Total DOS
  plotdos

  # d-orbitals for atom 0
  plotdos --atoms 0 --orbitals d

  # Multiple atoms and orbitals
  plotdos --atoms 0,2,5 --orbitals d
  plotdos --atoms 0-5 --orbitals d

  # Crystal field orbitals
  plotdos --atoms 0 --orbitals t2g
  plotdos --atoms 0 --orbitals eg

  # With customization
  plotdos --atoms 0 --orbitals d --xlim -5,5 --ylim -10,10
  plotdos --atoms 0 --orbitals d --output dos.png --dpi 600

  # Calculate d-band center
  plotdos --atoms 0 --orbitals d --band-center
  plotdos --atoms 0 --orbitals d --band-center --energy-range -10,0

  # Spin-resolved band center
  plotdos --atoms 0 --orbitals d --band-center --spin-treatment separate

  # Multiple orbitals (plots first one, calculates all)
  plotdos --atoms 0 --orbitals t2g,eg --band-center

  # Custom colors
  plotdos --atoms 0,1 --orbitals d --colors red,blue

Orbital options:
  s, p, d              Individual orbital types
  all-s, all-p, all-d  All orbitals of given type
  t2g                  t2g orbitals (dxy, dyz, dxz)
  eg                   eg orbitals (dz2, dx2-y2)

Spin treatment options (for --band-center):
  combined   Combine both spins (default)
  separate   Return both spin-up and spin-down
  up         Only spin-up
  down       Only spin-down
        '''
    )

    # Input/output
    parser.add_argument(
        '--doscar',
        default='DOSCAR',
        help='Path to DOSCAR file (default: DOSCAR)'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output file (default: show plot interactively). Format determined by extension.'
    )

    # What to plot
    parser.add_argument(
        '--atoms',
        help='Atom indices: "0", "0,2,5", or "0-5" (default: plot total DOS)'
    )
    parser.add_argument(
        '--orbitals',
        help='Orbitals: s, p, d, t2g, eg, all-s, all-p, all-d (comma-separated for multiple)'
    )
    parser.add_argument(
        '--total',
        action='store_true',
        help='Plot total DOS (overrides --atoms/--orbitals)'
    )

    # Appearance
    parser.add_argument(
        '--xlim',
        help='Energy range in eV, e.g., "-5,5" (default: full data range)'
    )
    parser.add_argument(
        '--ylim',
        help='DOS range, e.g., "-10,10" (default: auto)'
    )
    parser.add_argument(
        '--title',
        help='Plot title (default: auto-generated)'
    )
    parser.add_argument(
        '--dpi',
        type=int,
        default=300,
        help='Resolution for saved figure (default: 300)'
    )
    parser.add_argument(
        '--figsize',
        default='8,6',
        help='Figure size in inches as "width,height" (default: "8,6")'
    )
    parser.add_argument(
        '--linewidth',
        type=float,
        default=1.5,
        help='Line width (default: 1.5)'
    )
    parser.add_argument(
        '--colors',
        help='Custom colors (comma-separated), e.g., "red,blue,green"'
    )
    parser.add_argument(
        '--same-color-spins',
        action='store_true',
        help='Use same color for spin-up and spin-down'
    )
    parser.add_argument(
        '--no-legend',
        action='store_true',
        help='Disable legend'
    )

    # Analysis
    parser.add_argument(
        '--band-center',
        action='store_true',
        help='Calculate and print band center (requires --atoms and --orbitals)'
    )
    parser.add_argument(
        '--energy-range',
        help='Energy range for band center calculation, e.g., "-10,0"'
    )
    parser.add_argument(
        '--spin-treatment',
        choices=['combined', 'separate', 'up', 'down'],
        default='combined',
        help='Spin treatment for band center: combined, separate, up, or down (default: combined)'
    )

    args = parser.parse_args()

    # Check if DOSCAR exists
    if not os.path.exists(args.doscar):
        print(f'ERROR: DOSCAR file not found: {args.doscar}')
        print('Please specify the correct path using --doscar')
        sys.exit(1)

    # Load DOS
    try:
        dos = DOS(args.doscar)
        print(f'Loaded DOSCAR: {args.doscar}')
        print(f'  Atoms: {dos.natoms}')
        print(f'  Fermi energy: {dos.fermi_energy:.3f} eV')
        print(f'  Has partial DOS: {dos.has_partial_dos}')
        print(f'  Energy range: {dos.energy.min():.2f} to {dos.energy.max():.2f} eV')
    except Exception as e:
        print(f'ERROR: Failed to load DOSCAR: {e}')
        sys.exit(1)

    # Parse atoms
    atoms = parse_atoms(args.atoms) if args.atoms else None

    # Validate atom indices
    if atoms:
        invalid_atoms = [a for a in atoms if a < 0 or a >= dos.natoms]
        if invalid_atoms:
            print(f'ERROR: Invalid atom indices: {invalid_atoms}')
            print(f'Valid range: 0 to {dos.natoms - 1}')
            sys.exit(1)

    # Parse orbitals
    orbitals = [o.strip() for o in args.orbitals.split(',')] if args.orbitals else None

    # Check if PDOS is required but not available
    if (atoms or orbitals) and not args.total:
        if not dos.has_partial_dos:
            print('ERROR: Partial DOS requested but not available in DOSCAR')
            print('Please run VASP with LORBIT >= 10 to generate partial DOS')
            sys.exit(1)

    # Validate requirements for PDOS
    if not args.total and (atoms or orbitals):
        if not atoms:
            print('ERROR: --atoms required when plotting partial DOS')
            sys.exit(1)
        if not orbitals:
            print('ERROR: --orbitals required when plotting partial DOS')
            sys.exit(1)

    # Parse figure size
    figsize = tuple(map(float, args.figsize.split(',')))

    # Parse colors
    colors = parse_colors(args.colors)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot
    if args.total or (not atoms and not orbitals):
        # Total DOS
        dos.plot_total_dos(ax=ax)
        plot_title = args.title if args.title else 'Total DOS'
    else:
        # Partial DOS
        try:
            # Use first orbital for plotting
            primary_orbital = orbitals[0]

            # Plot using multi-atom plot for better handling
            dos.plot_multi_atom_pdos(
                atoms=atoms,
                orbitals=primary_orbital,
                ax=ax,
                same_color_spins=args.same_color_spins,
                colors=colors,
                linewidth=args.linewidth
            )

            # Generate title
            if args.title:
                plot_title = args.title
            else:
                atom_str = f"atom {atoms[0]}" if len(atoms) == 1 else f"atoms {','.join(map(str, atoms))}"
                plot_title = f'{primary_orbital.upper()} PDOS - {atom_str}'

        except Exception as e:
            print(f'ERROR: Failed to plot PDOS: {e}')
            sys.exit(1)

    # Set title
    ax.set_title(plot_title)

    # Set axis limits
    if args.xlim:
        xlim = parse_range(args.xlim)
        ax.set_xlim(xlim)
    else:
        # Use full data range by default
        ax.set_xlim(dos.energy.min(), dos.energy.max())

    if args.ylim:
        ylim = parse_range(args.ylim)
        ax.set_ylim(ylim)

    # Handle legend
    if args.no_legend:
        legend = ax.get_legend()
        if legend:
            legend.remove()

    # Tight layout
    plt.tight_layout()

    # Calculate band center if requested
    if args.band_center:
        if not atoms or not orbitals:
            print('\nWARNING: --band-center requires --atoms and --orbitals')
        else:
            energy_range = parse_range(args.energy_range) if args.energy_range else None

            print('\n' + '='*60)
            print('BAND CENTER ANALYSIS')
            print('='*60)

            atom_str = ','.join(map(str, atoms))
            print(f'Atoms: {atom_str}')
            if energy_range:
                print(f'Energy range: {energy_range[0]:.2f} to {energy_range[1]:.2f} eV')
            else:
                print(f'Energy range: full ({dos.energy.min():.2f} to {dos.energy.max():.2f} eV)')
            print(f'Spin treatment: {args.spin_treatment}')
            print()

            # Calculate for all specified orbitals
            for orbital in orbitals:
                try:
                    result = dos.calculate_band_center(
                        atoms=atoms,
                        orbitals=orbital,
                        energy_range=energy_range,
                        spin_treatment=args.spin_treatment
                    )

                    if args.spin_treatment == 'separate':
                        print(f'{orbital}-band center:')
                        print(f'  Spin-up:   {result["up"]:>8.4f} eV')
                        print(f'  Spin-down: {result["down"]:>8.4f} eV')
                        print(f'  Splitting: {result["up"] - result["down"]:>8.4f} eV')
                    else:
                        print(f'{orbital}-band center: {result:>8.4f} eV')

                except Exception as e:
                    print(f'ERROR calculating {orbital}-band center: {e}')

            print('='*60 + '\n')

    # Save or show
    if args.output:
        try:
            # Create output directory if needed
            output_dir = os.path.dirname(args.output)
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir)

            plt.savefig(args.output, dpi=args.dpi, bbox_inches='tight')
            print(f'\nSaved plot to: {args.output}')
        except Exception as e:
            print(f'ERROR: Failed to save figure: {e}')
            sys.exit(1)
    else:
        print('\nShowing interactive plot...')
        plt.show()


if __name__ == "__main__":
    main()
