import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

from pysoftk.folder_manager.folder_creator import *
from pysoftk.topologies.diblock import *
from pysoftk.format_printers.format_mol import *
from pysoftk.htp_tools.calculator_htp import *

import os
import itertools
import numpy as np # Added numpy import


def high_throughput_combination_patterning(precursors_folder="precursors", xtb_input_base_dir="xtb_gfn_jobs"):
    """
    Performs high-throughput patterning by creating all possible combinations
    of donor, bridge, and acceptor molecules. Each combination's XYZ file
    is placed into its own unique subdirectory, utilizing the Fld class
    for directory creation, serving as direct input for the Htp class.

    Args:
        precursors_folder (str): The path to the folder containing 'donor', 'bridge', and 'acceptor' subfolders.
        xtb_input_base_dir (str): The base directory where subdirectories for xtb jobs will be created.
    """

    fld_manager = Fld() # Instantiate Fld manager

    # Ensure the base directory for XTB inputs exists using os.makedirs
    # Fld.create does not directly support creating a single directory at an arbitrary absolute path,
    # it primarily makes directories relative to the current working directory.
    if not os.path.exists(xtb_input_base_dir):
        os.makedirs(xtb_input_base_dir)
        print(f"Created base directory for XTB jobs: {xtb_input_base_dir}")

    donor_mols = []
    bridge_mols = []
    acceptor_mols = []

    mol_types = {"donor": donor_mols, "bridge": bridge_mols, "acceptor": acceptor_mols}

    for mol_type, mol_list in mol_types.items():
        folder_path = os.path.join(precursors_folder, mol_type)
        if not os.path.exists(folder_path):
            print(f"Warning: Subfolder '{mol_type}' not found in '{precursors_folder}'. Skipping.")
            continue
        for f_name in os.listdir(folder_path):
            if f_name.endswith(".mol2"):
                mol_path = os.path.join(folder_path, f_name)
                try:
                    mol = Chem.MolFromMol2File(mol_path, removeHs=False, sanitize=False)
                    if mol is None:
                        print(f"Warning: Could not read molecule from {mol_path}. Skipping file: {mol_path}")
                        continue
                    # Store a tuple of (molecule object, molecule file name)
                    mol_list.append((mol, f_name))
                except Exception as e:
                    print(f"Error reading molecule from {mol_path}: {e}")
        if not mol_list:
            print(f"Warning: No valid .mol files found in '{folder_path}'.")

    if not (donor_mols and bridge_mols and acceptor_mols):
        print("Error: Not enough molecules found in all donor, bridge, and acceptor folders to form combinations.")
        return None

    all_combinations = list(itertools.product(donor_mols, bridge_mols, acceptor_mols))
    print(f"Total number of combinations to process: {len(all_combinations)}")

    # Generate all combination subdirectory names using Fld.fxd_name
    # We use a placeholder 'combination_prefix' because fxd_name appends '_i'
    combination_subdir_base_names = fld_manager.fxd_name("combination", len(all_combinations))

    # Store the original working directory
    original_cwd = os.getcwd()

    # Change directory to the base directory where job folders will be created
    # This is necessary because Fld.create makes directories relative to CWD.
    try:
        os.chdir(xtb_input_base_dir)
        print(f"Temporarily changed working directory to: {os.getcwd()}")

        # Use Fld.create to create all the subdirectories for combinations
        # We need to explicitly convert to list then numpy array for fixed_names
        fld_manager.create(fixed_names=np.array(list(combination_subdir_base_names)))

    except Exception as e:
        print(f"Error creating subdirectories with Fld: {e}")
        return None
    finally:
        # Always change back to the original working directory
        os.chdir(original_cwd)
        print(f"Changed working directory back to: {os.getcwd()}")

    # Open log files for successful and failed combinations
    # These files will be created in the directory where the script is run
    successful_log_file = "successful_combinations.log"
    failed_log_file = "failed_combinations.log"

    # Use 'with open' to ensure files are properly closed
    with open(successful_log_file, 'w') as f_success, open(failed_log_file, 'w') as f_failed:
        f_success.write("Successful Combinations Log:\n")
        f_failed.write("Failed Combinations Log:\n")

        string = "DBA"
        combination_counter = 0

        for i, (donor_tuple, bridge_tuple, acceptor_tuple) in enumerate(all_combinations):
            # Unpack the molecule objects and their names
            donor_mol, donor_name = donor_tuple
            bridge_mol, bridge_name = bridge_tuple
            acceptor_mol, acceptor_name = acceptor_tuple

            mols = [donor_mol, bridge_mol, acceptor_mol]
            combination_id = f"Combination_{i+1}"
            # Create a readable label including the molecule names
            combination_label = f"({donor_name}, {bridge_name}, {acceptor_name})"
            
            try:
                patt = Pt(string, mols, "Br").pattern_block_poly(swap_H=True)

                # Construct the full path to the already created subdirectory
                job_subdir_name = combination_subdir_base_names[i]
                job_subdir_path = os.path.join(xtb_input_base_dir, job_subdir_name)

                output_filename = f"initial_structure.xyz"
                output_file_path = os.path.join(job_subdir_path, output_filename)

                Fmt(patt).xyz_print(output_file_path) # Writes directly to the job's subdirectory
                combination_counter += 1
                print(f"Generated {output_filename} in {job_subdir_name}/")
                # Log successful combination with molecule names
                f_success.write(f"SUCCESS: {combination_id} {combination_label} -> {output_file_path}\n")

            # --- Specific Error Handling ---
            except rdkit.Chem.rdchem.AtomValenceException as e:
                error_message = f"RDKit Valence Error processing {combination_id} {combination_label}."
                print(f"ERROR: {error_message}")
                print(f"Details: {e}")
                # Log failed combination with molecule names
                f_failed.write(f"FAILED (ValenceError): {combination_id} {combination_label} | Details: {e}\n")
                continue # Skip to the next combination

            # --- Generic Error Handling ---
            except Exception as e:
                error_message = f"Unexpected Error processing {combination_id} {combination_label}."
                print(f"ERROR: {error_message}")
                print(f"Error type: {type(e).__name__}")
                print(f"Details: {e}")
                # Log failed combination with molecule names
                f_failed.write(f"FAILED (OtherError - {type(e).__name__}): {combination_id} {combination_label} | Details: {e}\n")
                continue # Skip to the next combination

    print(f"\nSuccessfully generated {combination_counter} XYZ files in subdirectories under '{xtb_input_base_dir}'.")
    print(f"Successful combinations logged to: {successful_log_file}")
    print(f"Failed combinations logged to: {failed_log_file}")
    return xtb_input_base_dir

def run_xtb_gfn_calculations(xtb_input_directory, xtb_command="xtb", max_concurrent_jobs=4, cores_per_xtb_job=2, threshold="crude"):
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


if __name__ == "__main__":
    high_throughput_combination_patterning()

    #Define variables needed for the High-throughput
    xtb_identifier = "/home/alejandro/Documents/xtb-6.6.1-linux-x86_64/xtb-6.6.1/bin/xtb"
    precursors_folder = "precursors"
    xtb_job_staging_dir = "xtb_gfn_jobs"

    # --- Step 1: Generate combination structures and prepare them for XTB ---
    print("STEP 1: Generating molecular combinations and preparing for XTB...")
    generated_xtb_input_path = high_throughput_combination_patterning(
        precursors_folder=precursors_folder,
        xtb_input_base_dir=xtb_job_staging_dir
    )

    if generated_xtb_input_path:
        # --- Step 2: Run XTB GFN-xTB calculations on the prepared structures ---
        print("\nSTEP 2: Running GFN-xTB calculations on prepared structures...")
        run_xtb_gfn_calculations(
            xtb_input_directory=generated_xtb_input_path,
            xtb_command=xtb_identifier, # Change this if xtb is not in your system's PATH
            max_concurrent_jobs=3, # Adjust based on your system's capabilities
            cores_per_xtb_job=4, # Adjust based on your system's capabilities
            threshold="tight"
        )
    else:
        print("\nSkipping XTB calculations as no molecular combinations were generated.")
