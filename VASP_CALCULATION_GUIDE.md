# VASP/DFT Calculation Guide

Personal reference for computational materials science calculations and lessons learned.

## Purpose

This document serves as a personal knowledge base for VASP calculations, collecting tips, best practices, and solutions discovered through experience. It will be expanded over time as new techniques and insights are gained.

The focus is on practical, actionable advice that has proven useful in real calculations, particularly for challenging systems and convergence issues.

## Spin-Polarized Calculations with Magnetic Materials

Description: Guidelines and tips for handling magnetic systems in VASP, particularly for challenging convergence cases and DFT+U calculations.

### Key Challenges
- Achieving convergence to desired magnetic moments in atoms
- Handling strongly correlated systems with DFT+U
- Dealing with magnetic frustration and multiple magnetic states
- Convergence issues in complex magnetic structures

### Tips and Best Practices

#### Two-Step Convergence Approach
For difficult magnetic systems, especially with DFT+U, use a staged approach:

1. **First step**: Non-spin-polarized single-point energy calculation
   ```bash
   # INCAR settings for step 1
   NSW = 0          # Single point calculation
   NELM = 800       # High number of electronic steps
   ISPIN = 1        # Non-spin-polarized
   LCHARG = .TRUE.  # Save charge density to CHGCAR
   ```

2. **Second step**: Spin-polarized calculation with initial magnetic moment
   ```bash
   # INCAR settings for step 2
   ICHARG = 1       # Read charge density from CHGCAR
   ISPIN = 2        # Spin-polarized
   MAGMOM = ...     # Set initial magnetic moments
   # Include DFT+U parameters if needed
   LDAU = .TRUE.
   LDAUTYPE = 2
   LDAUL = ...
   LDAUU = ...
   ```

**Why this works**: The first step provides a good starting charge density without the complexity of magnetic degrees of freedom, making the second step more likely to converge to the desired magnetic state.

#### Convergence Parameter Tuning
When standard SCF convergence fails, adjust magnetic mixing parameters as recommended in the [VASP wiki](https://www.vasp.at/wiki/index.php/AMIX_MAG):

```bash
# Conservative mixing for difficult magnetic systems
AMIX = 0.2        # Linear mixing parameter
BMIX = 0.0001     # Cutoff wave vector for Kerker mixing
AMIX_MAG = 0.8    # Linear mixing parameter for magnetization
BMIX_MAG = 0.0001 # Cutoff wave vector for magnetic mixing
```

**When to use**: 
- Oscillating magnetic moments during SCF
- Convergence failures in strongly correlated systems
- Systems with competing magnetic states

#### Initial Magnetic Moment Specification
Proper MAGMOM initialization is crucial for convergence:

```bash
# Example for a system with Fe (4 μB) and O (0 μB) atoms
MAGMOM = 4*4.0 8*0.0  # 4 Fe atoms, 8 O atoms

# For mixed valence systems, specify per atom
MAGMOM = 4.5 3.5 0.0 0.0  # Different Fe oxidation states
```

**Tips**:
- Start with reasonable magnetic moments based on expected oxidation states
- For transition metals, use typical values (Fe: 4-5 μB, Co: 3 μB, Ni: 2 μB)
- Set non-magnetic atoms (O, C, etc.) to 0.0
- For unknown systems, try different initial values if convergence fails

### Common Issues and Solutions

| Problem | Solution |
|---------|----------|
| Oscillating magnetic moments | Reduce AMIX and AMIX_MAG |
| Wrong magnetic ground state | Try different MAGMOM initialization |
| SCF not converging | Increase NELM, use two-step approach |
| Unexpected magnetic moments | Check for charge transfer, verify DFT+U parameters |

### Additional Notes

- Always check final magnetic moments in OUTCAR to verify they match expectations
- For antiferromagnetic systems, initialize with alternating spin directions
- Consider using NUPDOWN for systems with known net magnetic moment
- Monitor convergence in OSZICAR for both energy and magnetic moments

---

*This section will be expanded with additional calculation types and techniques as they are encountered and mastered.*