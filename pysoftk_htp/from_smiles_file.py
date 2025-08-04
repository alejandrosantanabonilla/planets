import os
from openbabel import pybel # pybel is syntax sugar for openbabel
from openbabel import openbabel as ob
import shutil

from pysoftk.folder_manager.folder_creator import *
from pysoftk.topologies.diblock import *
from pysoftk.format_printers.format_mol import *
from pysoftk.htp_tools.calculator_htp import *


def run_xtb_gfn_calculations(xtb_input_directory, xtb_command="xtb", max_concurrent_jobs=5, cores_per_xtb_job=2, threshold="vtight"):
    """
    Initiates parallel GFN-xTB geometry optimizations on all XYZ files
    within the subdirectories of the specified input directory, using the Htp class.

    Args:
        xtb_input_directory (str): The main directory containing subdirectories,
                                   each with one .xyz file for an xtb job.
        xtb_command (str): The command or path to the xtb executable. Defaults to "xtb".
        max_concurrent_jobs (int): The maximum number of xtb processes to run simultaneously.
                                   This is `max_work` for the Htp class.
        cores_per_xtb_job (int): The number of CPU cores to allocate to each individual
                                   xtb calculation via its `--parallel` flag. This is `num_cores` for Htp.
        threshold (str): The convergence threshold for the xtb optimization (e.g., "crude", "normal", "tight").
    """
    print(f"\n--- Starting XTB GFN-xTB Calculations in '{xtb_input_directory}' ---")
    try:
        htp_manager = Htp(xtb_command=xtb_command)
        htp_manager.htp_xtb_gfn(
            directory=xtb_input_directory,
            max_work=max_concurrent_jobs,
            num_cores=cores_per_xtb_job,
            threshold=threshold
        )
        print("\n--- XTB GFN-xTB Calculations Finished ---")
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("Please ensure xtb is installed and accessible via PATH, or provide its full path.")
    except Exception as e:
        print(f"An unexpected error occurred during XTB calculations: {e}")


def prepare_molecules_for_xtb(smiles_file_path, main_output_directory="xtb_molecules"):
    """
    Reads SMILES strings from a file, canonicalizes them using OpenBabel,
    creates a directory for each molecule, converts the SMILES to an XYZ file
    within that directory, and then initiates XTB GFN-xTB calculations.

    Args:
        smiles_file_path (str): The path to the file containing SMILES strings, one per line.
        main_output_directory (str): The main directory where molecule subdirectories
                                     and XYZ files will be created.
    """
    print(f"\n--- Preparing molecules from '{smiles_file_path}' for XTB calculations ---")

    # Ensure the main output directory exists. If it already exists, new molecules
    # will be added or overwritten if their index matches.
    os.makedirs(main_output_directory, exist_ok=True)
    print(f"Ensured main output directory exists: {main_output_directory}")

    processed_molecules_count = 0
    failed_molecules = []

    try:
        with open(smiles_file_path, 'r') as f:
            for i, line in enumerate(f):
                original_smiles = line.strip()
                if not original_smiles:
                    continue  # Skip empty lines

                print(f"Processing SMILES {i+1}: {original_smiles}")
                try:
                    # Create an OpenBabel molecule object from the SMILES string
                    # This also performs an initial perception of the molecule.
                    mol = pybel.readstring("smi", original_smiles)

                    # Canonicalize the SMILES. This ensures a consistent representation.
                    canonical_smiles = mol.write("can").strip()
                    print(f"  Canonical SMILES: {canonical_smiles}")

                    # Create a subdirectory for the current molecule.
                    # Each molecule gets its own folder for XTB input.
                    mol_dir_name = f"mol_{i}"
                    mol_output_path = os.path.join(main_output_directory, mol_dir_name)
                    os.makedirs(mol_output_path, exist_ok=True)
                    print(f"  Created subdirectory: {mol_output_path}")

                    # Define the XYZ file path within the molecule's subdirectory.
                    # XTB typically takes XYZ files as input.
                    xyz_file_name = f"{mol_dir_name}.xyz"
                    xyz_file_path = os.path.join(mol_output_path, xyz_file_name)

                    # --- Start of new code for 3D coordinate generation and forcefield relaxation ---
                    mol.make3D() # Generate 3D coordinates, handles adding H atoms and a quick forcefield cleanup

                    # Now perform conformer searching and optimization
                    ff = pybel._forcefields["mmff94"]
                    success = False
                    if ff: # Check if MMFF94 is available
                        success = ff.Setup(mol.OBMol)
                        if success:
                            print(f"  Using MMFF94 forcefield for SMILES '{original_smiles}'.")
                        else:
                            print(f"  Warning: MMFF94 setup failed for SMILES '{original_smiles}'. Trying UFF.")

                    if not success: # If MMFF94 failed or wasn't available, try UFF
                        ff = pybel._forcefields["uff"]
                        if ff: # Check if UFF is available
                            success = ff.Setup(mol.OBMol)
                            if success:
                                print(f"  Using UFF forcefield for SMILES '{original_smiles}'.")
                            else:
                                print(f"  Warning: UFF setup also failed for SMILES '{original_smiles}'. Skipping optimization.")
                        else:
                            print(f"  Warning: UFF forcefield not found. Skipping optimization for SMILES '{original_smiles}'.")

                    if success: # Only proceed with optimization if a forcefield setup was successful
                        # Feel free to tweak these to your balance of time / quality
                        ff.ConjugateGradients(100, 1.0e-3)
                        ff.FastRotorSearch(True) # permute central bonds
                        ff.WeightedRotorSearch(100, 25) # 100 cycles, each with 25 forcefield ops
                        # Final optimization
                        ff.ConjugateGradients(250, 1.0e-4)
                        # Update the coordinates in the pybel molecule object
                        ff.GetCoordinates(mol.OBMol)
                        print(f"  Performed forcefield optimization for SMILES '{original_smiles}'.")
                    else:
                        print(f"  No forcefield optimization performed for SMILES '{original_smiles}'. Using initial 3D coordinates.")
                    # --- End of new code ---

                    # Write the molecule to an XYZ file with the (potentially) optimized coordinates.
                    mol.write("xyz", xyz_file_path, overwrite=True)
                    print(f"  Saved XYZ file to: {xyz_file_path}")
                    processed_molecules_count += 1

                except Exception as e:
                    print(f"  Error processing SMILES '{original_smiles}': {e}")
                    failed_molecules.append(original_smiles)
                    continue

    except FileNotFoundError:
        print(f"Error: The SMILES file '{smiles_file_path}' was not found.")
        return
    except Exception as e:
        print(f"An unexpected error occurred while reading the SMILES file: {e}")
        return

    print(f"\n--- Finished preparing molecules. Processed {processed_molecules_count} molecules. ---")
    if failed_molecules:
        print(f"Failed to process {len(failed_molecules)} molecules:")
        for smiles in failed_molecules:
            print(f"  - {smiles}")

    # After preparing all molecules, run the XTB calculations.
    # This assumes that the `run_xtb_gfn_calculations` function will iterate
    # through the subdirectories created in `main_output_directory`.
    if processed_molecules_count > 0:
        xtb_command="/home/alejandro/Documents/xtb-6.6.1-linux-x86_64/xtb-6.6.1/bin/xtb"
        run_xtb_gfn_calculations(main_output_directory, xtb_command)
    else:
        print("No molecules were successfully prepared, skipping XTB calculations.")

# --- Example Usage ---
# 1. Call the function to prepare molecules and initiate XTB calculations.
#    The 'xtb_input_data' directory will be created to store the XYZ files.
prepare_molecules_for_xtb("all_smiles.dat", "xtb_input_data")
