import ase
from ase.atoms import Atoms
import numpy as np
from ase.io import read, write

from ase.optimize import BFGS
from ase.optimize.minimahopping import MinimaHopping
from tblite.ase import TBLite
from ase.units import Bohr, Rydberg, kJ, kB, fs, Hartree, mol, kcal
from ase.optimize.minimahopping import MHPlot


def relax(atoms, method="GFN2-xTB", Ediff0=1.0, T0=6500.0, totalsteps=165, output_filename_prefix="relaxed"):
    """Relaxes the molecule using TBLite and MinimaHopping with options.

    Args:
        atoms: The ASE Atoms object to relax.
        method (str): The TBLite method to use (default: "GFN2-xTB").
        Ediff0 (float): The energy difference threshold for MinimaHopping (default: 1.0).
        T0 (float): The initial temperature for MinimaHopping (default: 2000.0).
        totalsteps (int): The total number of MinimaHopping steps (default: 10).
        output_filename_prefix (str): Prefix for output filenames.

    Returns:
        ase.Atoms: The relaxed ASE Atoms object.
    """

    #print(f"Relaxing molecule with: method={method}, Ediff0={Ediff0}, T0={T0}, totalsteps={totalsteps}")

    calculator = TBLite(method=method, electronic_temperature=4500.0, max_iterations=1500)
    atoms.set_calculator(calculator)

    hop = MinimaHopping(atoms=atoms, Ediff0=Ediff0, T0=T0)
    hop(totalsteps=totalsteps)

    mhplot = MHPlot()
    mhplot.save_figure(f"{output_filename_prefix}_summary.png")

    traj_filename = str("minima.traj")
    traj = read(traj_filename)
    write(f"{output_filename_prefix}_minima.xyz", traj, format="xyz")

    return atoms


class MoleculeProcessor:
    def __init__(self, input_file=None, atoms=None, output_file=None):
        self.input_file = input_file
        self.output_file = output_file
        self.original_atoms = atoms
        self.aligned_atoms = None

    def load_atoms(self):
        if self.original_atoms is None and self.input_file is not None:
            try:
                self.original_atoms = read(self.input_file)
                return True
            except FileNotFoundError:
                print(f"Error: File '{self.input_file}' not found.")
                return False
            except Exception as e:
                print(f"An error occurred while reading the file: {e}")
                return False
        elif self.original_atoms is not None:
            return True
        else:
            print("Error: No input file or Atoms object provided.")
            return False

    def align_and_translate(self, atoms):
        if len(atoms) <= 1:
            print("Cannot align a single atom or an empty molecule.")
            return atoms.copy()

        new_atoms = atoms.copy()
        center_of_mass = new_atoms.get_center_of_mass()
        new_atoms.translate(-center_of_mass)

        moments_of_inertia = new_atoms.get_moments_of_inertia()
        sorted_indices = np.argsort(moments_of_inertia)[::-1]
        I1, I2, I3 = moments_of_inertia[sorted_indices]

        new_atoms.rotate(I1, 'x')
        new_atoms.rotate(I2, 'y')
        new_atoms.rotate(I3, 'z')

        return new_atoms

    def rotate_and_translate_copies(self, yaw, pitch, roll, translation_vector, num_copies):
        if not (len(yaw) == len(pitch) == len(roll) == len(translation_vector) == num_copies):
            raise ValueError("The number of yaw, pitch, roll, and translation vectors must match the number of copies.")

        all_atoms = Atoms()

        for i in range(num_copies):
            new_atoms = self.aligned_atoms.copy()

            new_atoms.rotate(yaw[i], 'z')
            new_atoms.rotate(pitch[i], 'y')
            new_atoms.rotate(roll[i], 'x')
            new_atoms.translate(translation_vector[i])
            all_atoms.extend(new_atoms)

        return all_atoms

    def process_molecules(self, yaw_rad, pitch_rad, roll_rad, translation_vector, num_copies, relax_molecule=False):
        if not self.load_atoms():
            return None

        self.aligned_atoms = self.align_and_translate(self.original_atoms)
        rotated_atoms = self.rotate_and_translate_copies(yaw_rad, pitch_rad, roll_rad, translation_vector, num_copies)

        if rotated_atoms:
            write(self.output_file, rotated_atoms)
            print(f"Original and translated molecules written to '{self.output_file}'.")

            if relax_molecule:
                relaxed_atoms = relax(rotated_atoms, output_filename_prefix=self.output_file[:-4])
                write(f"{self.output_file[:-4]}_relaxed.xyz", relaxed_atoms)
                print(f"Relaxed molecule written to '{self.output_file[:-4]}_relaxed.xyz'.")
            return True
        else:
            return False

# Example of restarting from a previous (potentially unfinished) calculation
# Assuming cluster_molecules.xyz exists from a previous run
processor_restart = MoleculeProcessor(input_file="mol.xyz", output_file="cluster_molecules_relaxed_restart.xyz")

# Directly relax the loaded atoms
if processor_restart.load_atoms():
    relaxed_atoms_restart = relax(processor_restart.original_atoms, output_filename_prefix="cluster_molecules_relaxed_restart")
    write("cluster_molecules_relaxed_restart.xyz", relaxed_atoms_restart)
    print("Restart relaxation complete.")
else:
    print("Restart failed: could not load atoms from cluster_molecules.xyz")
