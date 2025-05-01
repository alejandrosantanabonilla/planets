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

```

This initial block sets up two distinct relative placements (defined by rotation and translation) for creating dimer structures from the base molecule in mol.xyz. It also prepares
configuration dictionaries (tblite_config, mh_config) that will only be used if the relaxation feature is turned on in the subsequent steps.

## Test Case 1: Generate Assembly Without Relaxation

### Goal

To generate the initial dimer structures based on the defined orientations and translations without performing any subsequent energy minimization or structural relaxation. This is 
useful for quickly creating and visualizing the starting geometries.

```python
# --- Test Case 1: Process without relaxation ---
print("\n--- Test Case 1: No Relaxation ---")

# Initialize the processor with input and the desired output file for the assemblies
processor1 = MoleculeProcessor(input_file="mol.xyz", output_file="initial_assembly.xyz")

# Process the molecules
result1 = processor1.process_molecules(
    yaw_rad=yaw,
    pitch_rad=pitch,
    roll_rad=roll,
    translation_vector=translations,
    relax_molecule=False # Key parameter: Relaxation is turned OFF
)

# Check if the process completed successfully
if result1:
    print("Test Case 1 completed. Check 'initial_assembly.xyz'.")
```

### Explanation
1. An instance of MoleculeProcessor is created. It's configured to read from mol.xyz and write the generated structures to initial_assembly.xyz.
2. The process_molecules method is called with the orientation angles (yaw_rad, pitch_rad, roll_rad) and translation_vector defined earlier.
3. The crucial parameter here is relax_molecule=False. This instructs the processor to:
    - Read the base molecule from mol.xyz.
    - For each specified placement (combination of rotation and translation):
      - Create a copy of the base molecule.
      - Apply the rotation (yaw, pitch, roll) and translation to the copy.
      - Combine the original molecule and the transformed copy to form a dimer.
4. Write the coordinates of all generated dimer structures sequentially into the specified output_file (initial_assembly.xyz).
   No energy calculations or geometry changes are performed.

### Expected Output
A single file named initial_assembly.xyz. This file will contain the atomic coordinates of the generated dimer structures, concatenated one after another in standard XYZ format. 
The number of structures will match the number of orientation/translation pairs provided (two in this example).

## Test Case 2: Generate Assembly WITH Relaxation (Fresh Start)

### Goal
To generate the initial dimer structures (as in Test Case 1) and then perform a structural relaxation or energy minimization on each generated structure. This test case starts the relaxation process from scratch, ignoring any previous relaxation runs.

```python
# --- Test Case 2: Process WITH relaxation (fresh start) ---
print("\n--- Test Case 2: With Relaxation (Fresh Start) ---")

# Clean up potential relaxation output files from any previous runs
# This ensures that restart_relax=False truly starts fresh.
for f in ['minima.traj', 'qn*.traj', 'hop.log', 'relax_run_*.xyz', 'relax_run_*.png']:
    for ff in glob.glob(f):
        if os.path.exists(ff): os.remove(ff)

# Initialize the processor. Note: output_file still refers to the *initial* assembly.
processor2 = MoleculeProcessor(input_file="mol.xyz", output_file="initial_assembly_for_relax.xyz")

# Process the molecules, enabling relaxation
result2 = processor2.process_molecules(
    yaw_rad=yaw, pitch_rad=pitch, roll_rad=roll, translation_vector=translations,
    relax_molecule=True,          # Key parameter: Relaxation is turned ON
    restart_relax=False,         # Key parameter: Do NOT restart, start fresh
    totalsteps=50,               # Number of steps for the relaxation simulation/optimization
    tblite_params=tblite_config, # Pass the calculator configuration
    mh_params= mh_config,        # Pass the relaxation algorithm configuration
    output_filename_prefix="relax_run" # Prefix for relaxation-specific output files
)

if result2:
    print("Test Case 2 completed. Check 'initial_assembly_for_relax.xyz', 'relax_run_minima.xyz', etc.")
```

### Explanation

1. Cleanup: The code first attempts to delete common output files associated with the relaxation process (like .traj, .log files). This step is crucial to guarantee that setting restart_relax=False results in a completely new relaxation run, unaffected by any previous attempts.
2. An instance of MoleculeProcessor is created. The output_file (initial_assembly_for_relax.xyz) will contain the initial, unrelaxed structures generated before relaxation begins.
3. The process_molecules method is called with several important parameters enabled for relaxation:
   - relax_molecule=True: This activates the relaxation procedure for each dimer generated.
   - restart_relax=False: This explicitly tells the relaxation algorithm not to look for or use any existing state/restart files. Each relaxation starts from the initial generated geometry.
   - totalsteps=50: This defines the duration or extent of each relaxation run (e.g., maximum number of optimization steps, molecular dynamics steps, or search steps). Here, it's set to a small value for a short demonstration run.
   - tblite_params, mh_params: The configuration dictionaries defined earlier are passed to configure the energy calculator and the relaxation algorithm, respectively.
   - output_filename_prefix="relax_run": This string is used as a base name for various output files generated during the relaxation phase. This helps organize outputs when multiple relaxation runs might occur.

### Expected Output


