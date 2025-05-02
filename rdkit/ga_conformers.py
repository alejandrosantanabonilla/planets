#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os # Import os module for path manipulation
from openbabel import pybel, openbabel as ob

# --- The Function (generate_conformers - unchanged from previous answer) ---

def generate_conformers(molecule: pybel.Molecule,
                        num_conformers: int,
                        mutability: int = 5,
                        convergence: int = 5,
                        num_children: int = None) -> pybel.Molecule:
    """
    Generates conformers for a given molecule using Open Babel's
    weighted rotor search.

    Args:
        molecule: A pybel.Molecule object. It MUST have 3D coordinates
                  generated beforehand (e.g., using molecule.make3D()).
        num_conformers: The target number of conformers to generate.
        mutability: Parameter controlling the rotor key search randomization
                    (default: 5). Higher values explore more diverse torsions.
        convergence: Parameter controlling the convergence criteria of the
                     search (default: 5). Higher values mean stricter convergence.
        num_children: Number of children conformers generated in each step
                      (default: num_conformers * 2).

    Returns:
        A *new* pybel.Molecule object containing the generated conformers.
        The original molecule object is not modified. Returns None if
        num_conformers is not positive or if the input molecule is invalid.

    Raises:
        ValueError: If num_conformers is less than 1.
    """
    if num_conformers < 1:
        raise ValueError("Number of conformers (num_conformers) must be at least 1.")

    if not molecule or not molecule.OBMol.Has3D():
         print("Error: Input molecule is invalid or does not have 3D coordinates.", file=sys.stderr)
         print("Hint: Generate 3D coords first, e.g., using molecule.make3D()", file=sys.stderr)
         return None

    # Work on a clone to avoid modifying the original molecule
    mol_clone = pybel.Molecule(molecule.OBMol)

    # Determine the number of children if not specified
    children = num_children if num_children is not None else num_conformers * 2

    # --- Conformer Search Setup ---
    cs = ob.OBConformerSearch()
    setup_successful = cs.Setup(mol_clone.OBMol,
                                num_conformers,  # target number of conformers
                                children,        # numChildren in docs
                                mutability,      # mutability in docs
                                convergence)     # convergence in docs

    if not setup_successful:
        print(f"Error: OBConformerSearch setup failed for molecule: {molecule.title}", file=sys.stderr)
        return None # Or raise an exception

    # --- Perform the Search ---
    print(f"Starting conformer search for {molecule.title or 'molecule'}...")
    cs.Search()
    print("Conformer search finished.")

    # --- Retrieve Conformers ---
    cs.GetConformers(mol_clone.OBMol)

    # Get the actual number found (might be less than requested)
    actual_conformers = mol_clone.OBMol.NumConformers()
    print(f"Found {actual_conformers} conformers.")

    return mol_clone

# --- Example Usage (Modified for separate file saving) ---

if __name__ == "__main__":
    # 1. Create a molecule (e.g., from SMILES)
    smiles = "CCCO"  # Propanol
    mol = pybel.readstring("smi", smiles)
    base_name = "propanol" # Used for filenames
    mol.title = base_name.capitalize()

    # 2. Generate initial 3D coordinates (REQUIRED for conformer search)
    print(f"Generating initial 3D coordinates for {mol.title}...")
    mol.make3D(forcefield="mmff94", steps=50)
    if not mol.OBMol.Has3D():
       print("Failed to generate 3D coordinates. Exiting.", file=sys.stderr)
       sys.exit(1)
    print("Initial 3D coordinates generated.")

    # 3. Generate Conformers using the function
    n_conformers_to_find = 15 # Let's find 15 conformers
    try:
        molecule_with_conformers = generate_conformers(
            molecule=mol,
            num_conformers=n_conformers_to_find
            # You can add mutability, convergence etc. here if needed
        )
    except ValueError as e:
        print(f"Error generating conformers: {e}", file=sys.stderr)
        sys.exit(1)

    # 4. Work with the result and save to SEPARATE files
    if molecule_with_conformers:
        num_found = molecule_with_conformers.OBMol.NumConformers()
        print(f"\nSuccessfully generated {num_found} conformers for {molecule_with_conformers.title}.")

        # Create a directory to store the conformer files (optional but good practice)
        output_dir = f"{base_name}_conformers"
        os.makedirs(output_dir, exist_ok=True) # exist_ok=True prevents error if dir exists
        print(f"Saving individual conformers into directory: '{output_dir}/'")

        # Loop through each found conformer index
        for i in range(num_found):
            # Set the molecule's coordinates to the i-th conformer
            molecule_with_conformers.OBMol.SetConformer(i)

            # Create a *new* molecule object containing ONLY the current conformer's coordinates
            # This is necessary because molecule_with_conformers still holds *all* conformer data internally.
            # Creating a new molecule from the OBMol captures only the *active* coordinate set.
            single_conformer_mol = pybel.Molecule(molecule_with_conformers.OBMol)

            # Optional: Set a unique title for the individual conformer file
            single_conformer_mol.title = f"{molecule_with_conformers.title}_conformer_{i+1}"

            # Define the output filename for this specific conformer
            # Using zfill for padding (e.g., conf_01, conf_02 ... conf_10)
            padding = len(str(num_found)) # Dynamically determine padding based on total number
            output_filename = os.path.join(
                output_dir,
                f"{base_name}_conf_{str(i+1).zfill(padding)}.sdf" # e.g., propanol_conf_01.sdf
            )
            file_format = "sdf" # Or "mol", "pdb", etc.

            # Save the single-conformer molecule to its own file
            try:
                single_conformer_mol.write(file_format, output_filename, overwrite=True)
                # Print progress without newline for cleaner output
                print(f"  Saved: {os.path.basename(output_filename)}", end='\r')
            except Exception as e:
                print(f"\nError saving conformer {i} to file '{output_filename}': {e}", file=sys.stderr)

        print("\nFinished saving all conformers to separate files.") # Print newline after loop

    else:
        print("\nConformer generation failed, no files saved.")
