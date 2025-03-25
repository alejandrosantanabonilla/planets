import ase
from ase.atoms import Atoms
import numpy as np
from ase.io import read, write

from ase.optimize import BFGS
from ase.optimize.minimahopping import MinimaHopping
from tblite.ase import TBLite

from ase.units import Bohr, Rydberg, kJ, kB, fs, Hartree, mol, kcal
from ase.optimize.minimahopping import MHPlot
from scipy.spatial.transform import Rotation as R

def relax(atoms, method="GFN2-xTB", Ediff0=1.0, T0=6500.0, totalsteps=150, output_filename_prefix="relaxed"):
    """Relaxes the molecule using TBLite and MinimaHopping with options."""

    print(f"Relaxing molecule with: method={method}, Ediff0={Ediff0}, T0={T0}, totalsteps={totalsteps}")

    calculator = TBLite(method=method)
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
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file
        self.original_atoms = None
        self.aligned_atoms = None

    def load_atoms(self):
        try:
            self.original_atoms = read(self.input_file)
            return True
        except FileNotFoundError:
            print(f"Error: File '{self.input_file}' not found.")
            return False
        except Exception as e:
            print(f"An error occurred while reading the file: {e}")
            return False

    def align_and_translate(self, atoms):
        """Aligns the molecule and translates it to the origin."""
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

    def rotate_and_translate_copies(self, yaw, pitch, roll, translation_vector):
        """Rotates and then translates multiple aligned molecules using quaternions."""
        num_copies = len(yaw)
        all_atoms = Atoms()

        for i in range(num_copies):
            new_atoms = self.aligned_atoms.copy()

            rotation = R.from_euler('zyx', [yaw[i], pitch[i], roll[i]])
            quaternion = rotation.as_quat()

            center = new_atoms.get_center_of_mass()
            positions = new_atoms.get_positions() - center
            rotated_positions = rotation.apply(positions) + center

            new_atoms.set_positions(rotated_positions)
            new_atoms.translate(translation_vector[i])
            all_atoms.extend(new_atoms)

        return all_atoms

    def process_molecules(self, yaw_rad, pitch_rad, roll_rad, translation_vector, relax_molecule=False):
        """Processes the molecules: loads, aligns, rotates, translates, and writes."""
        if not self.load_atoms():
            return None  # File error

        self.aligned_atoms = self.align_and_translate(self.original_atoms)
        rotated_atoms = self.rotate_and_translate_copies(yaw_rad, pitch_rad, roll_rad, translation_vector)

        if rotated_atoms:
            write(self.output_file, rotated_atoms)
            print(f"Original and translated molecules written to '{self.output_file}'.")

            if relax_molecule:  # Conditionally relax
                relaxed_atoms = relax(rotated_atoms, output_filename_prefix=self.output_file[:-4])
                write(f"{self.output_file[:-4]}_relaxed.xyz", relaxed_atoms)
                print(f"Relaxed molecule written to '{self.output_file[:-4]}_relaxed.xyz'.")
            return True
        else:
            return False

    def restart_relaxation(self, relaxed_file, yaw_rad):
        """Restarts the relaxation from a given relaxed structure file."""
        try:
            atoms = read(relaxed_file)
            atoms_per_monomer = len(atoms) // len(yaw_rad)

            if len(atoms) % len(yaw_rad) != 0:
                print("Error: Atom count in restart file is not a multiple of the number of monomers. Aborting restart.")
                return False

            print(f"Restarting relaxation from '{relaxed_file}'.")
            relaxed_atoms = relax(atoms, output_filename_prefix=relaxed_file[:-4] + "_restarted")
            write(f"{relaxed_file[:-4]}_restarted_relaxed.xyz", relaxed_atoms)
            print(f"Restarted relaxed molecule written to '{relaxed_file[:-4]}_restarted_relaxed.xyz'.")
            return True
        except FileNotFoundError:
            print(f"Error: Restart file '{relaxed_file}' not found.")
            return False
        except Exception as e:
            print(f"An error occurred during restart: {e}")
            return False


# Example usage:
processor = MoleculeProcessor("mol.xyz", "cluster_molecules.xyz")

yaw_rad = np.deg2rad([0, 180])
pitch_rad = np.deg2rad([0, 0])
roll_rad = np.deg2rad([0, 180])

translation_vectors = [
    [2.5, 3.0, -4.0],
    [-3.0, 4.0, 3.0]
]

success = processor.process_molecules(yaw_rad, pitch_rad, roll_rad, translation_vectors, relax_molecule=True)

#if success is None:
#    print("File error occurred.")
#elif success:
#    print("Molecule processing complete.")
    # Restart the relaxation from the generated relaxed file:
    #restart_success = processor.restart_relaxation("cluster_molecules_relaxed.xyz", yaw_rad)
    #if restart_success:
    #    print("Relaxation restarted successfully.")
    #else:
    #    print("Relaxation restart failed.")

#else:
    #print("Molecule processing failed.")
