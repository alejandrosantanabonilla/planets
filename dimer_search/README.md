# Dimer Search via Basin Hopping with RDKit, ASE, and TBLite

This repository provides Python code for processing molecule dimers, including rotation, joining, and relaxation
using RDKit, ASE, and TBLite.

## Code Structure

This Python code defines a MoleculeProcessor class using the ASE library to:

1. Load a molecule from an input file (e.g., XYZ format).
2. Align the molecule by moving its center of mass to the origin and rotating its principal axes to align with the Cartesian axes.
3. Generate an assembly by creating multiple copies of the aligned molecule, applying specified rotations (yaw, pitch, roll) and translations to each copy.
4. Write the initial, unrelaxed assembly to an output file.
5. Optionally relax the final assembly using an external relax function (imported from utils.py), with support for restarting the relaxation from a previous state.

The main process_molecules method orchestrates these steps.

## Installation

It's highly recommended to use a virtual environment to manage dependencies.  Here's how:

1.  **Create a virtual environment:**

```bash
python3 -m venv dimer_env  # Creates a virtual environment named "dimer_env"
```

2. **Activate the virtual environment:**


```bash
source dimer_env/bin/activate  # On Linux/macOS
```

3. **Clone the repository:**

```bash
git clone https://github.com/alejandrosantanabonilla/dimer_search.git
cd dimer_search
```

4. **Install dimer_search**

```bash
pip install .
```

# Tutorial: Using MoleculeProcessor for Generating and Relaxing Molecular Assemblies

This tutorial explains how to use the `MoleculeProcessor` class, presumably from a library like `dimer_search`, to generate molecular assemblies based on specified orientations and translations. It also covers how to perform structural relaxation (e.g., energy minimization or conformational search) on these assemblies and utilize the restart feature for long calculations. We will walk through the three distinct test cases provided in the example script.

## Prerequisites

* A working Python environment.
* Necessary libraries installed: `numpy` and the library containing `MoleculeProcessor` (referred to as `dimer_search` in the script).
* An input coordinate file named `mol.xyz` located in the same directory as the script. This file should contain the atomic coordinates of the single molecule you wish to assemble into dimers or larger structures.

## Initial Setup

The script begins by importing required libraries and defining the parameters that control how the molecular assemblies are generated and, optionally, relaxed.

```python
import numpy as np
# Assuming MoleculeProcessor is in dimer_search
from dimer_search import *
import glob
import os

print("Running MoleculeProcessor Example...")

# --- Define Orientations ---
# Yaw, Pitch, and Roll angles specify the rotation applied to the second molecule
# relative to the first one. Multiple angles can be provided to generate
# different relative orientations. Angles are converted to radians.
yaw = np.deg2rad([0, 180])
pitch = np.deg2rad([0, 180])
roll = np.deg2rad([0, 0])

# --- Define Translations ---
# These vectors specify the displacement (shift) of the second molecule's
# center relative to the first molecule's center after rotation.
# Each translation vector corresponds to a set of (yaw, pitch, roll) angles.
translations = [
    [10.5, 0.0, 8.0], # Corresponds to yaw[0], pitch[0], roll[0]
    [1.0, 0.0, 9.0]   # Corresponds to yaw[1], pitch[1], roll[1]
]

# --- Define Relaxation Parameters (if relaxation is enabled) ---

# Configuration for the underlying energy calculator (e.g., TBLite)
# These parameters are passed directly to the calculator.
tblite_config = {
    "electronic_temperature": 5500.0, # Example: Fermi-Dirac smearing temperature
    "max_iterations": 300,            # Example: Max SCF cycles
}

# Configuration for the relaxation/search algorithm (e.g., Metropolis-Hastings)
# These parameters control the simulation dynamics.
mh_config = {
    "T0": 1200.0,   # Initial temperature for sampling (e.g., in Kelvin)
    "Ediff0": 0.6,  # Initial energy difference criterion (e.g., in eV)
    "fmax": 0.1     # Force convergence criterion for geometry optimization (e.g., in eV/Angstrom)
}
```

This initial block sets up two distinct relative placements (defined by rotation and translation) for creating dimer structures from the base molecule in mol.xyz. It also prepares
configuration dictionaries (tblite_config, mh_config) that will only be used if the relaxation feature is turned on in the subsequent steps.
