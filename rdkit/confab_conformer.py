from openbabel import pybel, openbabel 
from typing import Union, Optional, List
import os
import re # For sanitizing filenames

def generate_and_save_conformers_ob(
    molecule_input: Union[str, pybel.Molecule],
    output_dir: str,
    output_format: str = 'sdf',
    base_filename: Optional[str] = None,
    forcefield: str = 'uff',
    nconfs: int = 10000,
    rmsd: float = 0.5,
    energy_gap: float = 50.0, # Changed default, unit is kJ/mol for DiverseConfGen
    initial_opt_steps: int = 50,
    verbose: bool = False
) -> int:
    """
    Generates conformers for a given molecule using Open Babel and saves
    each conformer to a separate file in the specified directory.

    Args:
        molecule_input: The input molecule. Can be a SMILES string or an
                        existing pybel.Molecule object.
        output_dir: The directory where individual conformer files will be saved.
                    It will be created if it doesn't exist.
        output_format: The file format for the output conformers (e.g., 'sdf',
                       'pdb', 'mol2'). Defaults to 'sdf'.
        base_filename: Optional base name for the output files. If None, uses
                       the molecule's title or a sanitized version of the SMILES.
                       Files will be named like:
                       '{output_dir}/{base_filename}_conf_{index}.{output_format}'
        forcefield: The force field to use (e.g., 'uff', 'mmff94', 'ghemical').
                    Defaults to 'uff'.
        nconfs: The target number of conformers to generate in the initial
                search. Defaults to 10000.
        rmsd: The RMSD cutoff (in Angstroms) for considering conformers
              different. Defaults to 0.5.
        energy_gap: The energy window (in kJ/mol) above the minimum energy
                    conformer to keep. Defaults to 50.0 kJ/mol.
        initial_opt_steps: Number of optimization steps for initial 3D structure.
                           Defaults to 50.
        verbose: If True, enables verbose output from the conformer
                 generation process. Defaults to False.

    Returns:
        The number of conformer files successfully saved. Returns 0 if
        critical errors occur before saving.

    Raises:
        ValueError: If input is invalid or molecule creation fails.
        KeyError: If the specified forcefield is not available.
        RuntimeError: If initial 3D structure generation fails.
        OSError: If the output directory cannot be created.
    """
    mol = None
    input_was_smiles = False
    original_smiles = "" # Keep track if input was SMILES

    # --- 1. Input Handling and Molecule Initialization ---
    if isinstance(molecule_input, pybel.Molecule):
        try:
            mol = pybel.Molecule(molecule_input.OBMol) # Clone
            if mol.OBMol is None or mol.OBMol.NumAtoms() == 0:
                raise ValueError("Input pybel.Molecule object is invalid or empty.")
        except Exception as e:
             raise ValueError(f"Failed to clone input pybel.Molecule object: {e}") from e
    elif isinstance(molecule_input, str):
        input_was_smiles = True
        original_smiles = molecule_input
        try:
            mol = pybel.readstring("smi", molecule_input)
            if mol.OBMol is None or mol.OBMol.NumAtoms() == 0:
                raise ValueError(f"Failed to create molecule from SMILES: {molecule_input}")
            if not mol.title: # Set title from SMILES if not present
                mol.title = molecule_input
        except IOError as e:
             raise ValueError(f"Failed to parse SMILES string: '{molecule_input}'") from e
    else:
        raise TypeError("Input must be a SMILES string or a pybel.Molecule object.")

    # --- 2. Prepare Molecule: Add Hydrogens and Generate Initial 3D Coords ---
    mol.addh()
    try:
        mol.make3D(forcefield=forcefield, steps=initial_opt_steps)
    except RuntimeError as e:
        raise RuntimeError(f"Initial 3D structure generation failed using '{forcefield}': {e}") from e
    except Exception as e:
        raise RuntimeError(f"An unexpected error occurred during initial 3D generation: {e}") from e

    # --- 3. Conformer Generation using Force Field ---
    try:
        ff = pybel._forcefields[forcefield]
    except KeyError:
        available_ffs = ', '.join(pybel._forcefields.keys())
        raise KeyError(
            f"Force field '{forcefield}' not found. Available: {available_ffs}"
        )

    if not ff.Setup(mol.OBMol):
        print(
            f"Warning: Force field '{forcefield}' setup failed for molecule "
            f"'{mol.title or 'untitled'}'. "
            f"Cannot generate conformers. Returning 0."
        )
        return 0 # Cannot proceed

    # Perform the diverse conformer generation
    try:
        ff.DiverseConfGen(rmsd, nconfs, energy_gap, verbose)
        ff.GetConformers(mol.OBMol) # Update the OBMol with conformers
    except Exception as e:
        print(f"Warning: An error occurred during conformer generation: {e}")
        print(f"Proceeding with any conformers generated so far (might just be the initial 3D).")
        # Continue, mol.OBMol might still have the initial conformer or some generated ones

    num_conformers = mol.OBMol.NumConformers()
    if num_conformers <= 0:
         print("No conformers generated (or only initial 3D failed). Returning 0.")
         return 0

    print(f"Generated {num_conformers} conformers for molecule '{mol.title or 'untitled'}'.")

    # --- 4. Prepare Output Directory ---
    try:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Saving conformers to directory: {output_dir}")
    except OSError as e:
        raise OSError(f"Failed to create output directory '{output_dir}': {e}") from e

    # --- 5. Determine Base Filename ---
    if base_filename is None:
        if mol.title:
            # Sanitize title for use as filename
            # Remove or replace characters invalid for filenames
            sanitized_title = re.sub(r'[<>:"/\\|?*]', '_', mol.title)
            sanitized_title = re.sub(r'\s+', '_', sanitized_title) # Replace spaces
            base_filename = sanitized_title or "molecule" # Fallback if title is empty after sanitizing
        elif input_was_smiles:
             # Sanitize SMILES for use as filename (can be long/complex)
             sanitized_smiles = re.sub(r'[<>:"/\\|?*()\[\]{}=]', '_', original_smiles)
             sanitized_smiles = re.sub(r'\s+', '_', sanitized_smiles)
             base_filename = sanitized_smiles[:50] or "molecule" # Limit length
        else:
            base_filename = "molecule" # Generic fallback
    else:
        # Sanitize user-provided base_filename as well
        base_filename = re.sub(r'[<>:"/\\|?*]', '_', base_filename)
        base_filename = re.sub(r'\s+', '_', base_filename)

    # --- 6. Iterate and Save Each Conformer ---
    saved_count = 0
    for i in range(num_conformers):
        # Set the molecule object to the i-th conformer
        mol.OBMol.SetConformer(i)

        # Create a *new* molecule object containing only this conformer's coordinates
        # This is crucial because molecule.write() might otherwise write all conformers
        single_conf_mol = pybel.Molecule(mol.OBMol)
        single_conf_mol.title = f"{mol.title or base_filename}_conf_{i+1}" # Add index to title

        # Construct the output filename
        output_filename = f"{base_filename}_conf_{i+1}.{output_format}"
        output_path = os.path.join(output_dir, output_filename)

        # Save the single conformer molecule
        try:
            # Use opt=None with write() to avoid implicit make3D or opt steps
            # overwrite=True is good practice here
            single_conf_mol.write(
                format=output_format,
                filename=output_path,
                overwrite=True,
                opt=None # Prevent re-optimization or coordinate changes
            )
            saved_count += 1
            if verbose and (i + 1) % 50 == 0: # Print progress periodically
                 print(f"  Saved conformer {i + 1}/{num_conformers}")
        except IOError as e:
            print(f"Warning: Failed to write conformer {i+1} to {output_path}: {e}")
        except Exception as e:
            print(f"Warning: An unexpected error occurred writing conformer {i+1}: {e}")

    print(f"Successfully saved {saved_count} out of {num_conformers} conformers.")
    return saved_count

# --- Example Usage ---
if __name__ == "__main__":
    output_base_dir = "conformer_output" # Base directory for examples

    # Example 1: Using SMILES input
    smiles_ethane = "CC"
    ethane_output_dir = os.path.join(output_base_dir, "ethane")
    print(f"\n--- Generating conformers for Ethane ({smiles_ethane}) ---")
    try:
        saved_ethane_count = generate_and_save_conformers_ob(
            smiles_ethane,
            output_dir=ethane_output_dir,
            output_format='sdf',
            forcefield='mmff94',
            rmsd=0.1,
            nconfs=100,
            energy_gap=10.0, # kJ/mol
            verbose=False # Set to True for more output
        )
        print(f"Saved {saved_ethane_count} conformers for Ethane in '{ethane_output_dir}'")

    except (ValueError, KeyError, RuntimeError, TypeError, OSError) as e:
        print(f"Error generating/saving conformers for Ethane: {e}")

    # Example 2: Using a pybel.Molecule object
    smiles_complex = "C1=CC(=C(C=C1Cl)O)C(=O)O" # 3-Chloro-4-hydroxybenzoic acid
    complex_output_dir = os.path.join(output_base_dir, "complex_acid")
    print(f"\n--- Generating conformers for {smiles_complex} ---")
    try:
        mol_obj = pybel.readstring("smi", smiles_complex)
        # Providing a specific base filename
        acid_base_name = "3_Chloro_4_hydroxybenzoic_acid"
        mol_obj.title = acid_base_name # Set title on original object too

        saved_complex_count = generate_and_save_conformers_ob(
            mol_obj,
            output_dir=complex_output_dir,
            output_format='pdb', # Save as PDB this time
            base_filename=acid_base_name, # Explicitly provide base name
            forcefield='uff',
            nconfs=500, # Reduced for faster example
            rmsd=0.5,
            energy_gap=40.0, # kJ/mol
            verbose=False
        )
        print(f"Saved {saved_complex_count} conformers in PDB format in '{complex_output_dir}'")

    except (ValueError, KeyError, RuntimeError, TypeError, OSError) as e:
        print(f"Error generating/saving conformers for complex molecule: {e}")

    # Example 3: Invalid forcefield (should raise KeyError)
    print("\n--- Testing invalid forcefield ---")
    try:
        generate_and_save_conformers_ob("C", output_dir="invalid_ff_test", forcefield='invalid_ff')
    except KeyError as e:
        print(f"Caught expected error: {e}")
    except Exception as e:
        print(f"Caught unexpected error: {e}") # Catch others just in case

    # Example 4: Invalid SMILES (should raise ValueError)
    print("\n--- Testing invalid SMILES ---")
    try:
        generate_and_save_conformers_ob("C(C)(C)(C)(C)C-", output_dir="invalid_smiles_test")
    except ValueError as e:
        print(f"Caught expected error: {e}")
    except Exception as e:
        print(f"Caught unexpected error: {e}")
