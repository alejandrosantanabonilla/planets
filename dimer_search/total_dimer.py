import ase
from ase.atoms import Atoms
import numpy as np
from ase.io import read, write
from scipy.spatial.transform import Rotation as R
from .utils import relax

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

    def process_molecules(self, yaw_rad, pitch_rad, roll_rad, translation_vector, relax_molecule=False, tblite_params=None, mh_params=None, totalsteps=20, output_filename_prefix="relaxed"):
        """Processes the molecules: loads, aligns, rotates, translates, and writes."""
        if not self.load_atoms():
            return None

        self.aligned_atoms = self.align_and_translate(self.original_atoms)
        rotated_atoms = self.rotate_and_translate_copies(yaw_rad, pitch_rad, roll_rad, translation_vector)

        if rotated_atoms:
            write(self.output_file, rotated_atoms)
            print(f"Original and translated molecules written to '{self.output_file}'.")

            if relax_molecule:
                relaxed_atoms = relax(rotated_atoms, output_filename_prefix=output_filename_prefix, tblite_params=tblite_params, mh_params=mh_params, totalsteps=totalsteps)
                write(f"{output_filename_prefix}_relaxed.xyz", relaxed_atoms)
                print(f"Relaxed molecule written to '{output_filename_prefix}_relaxed.xyz'.")
            return True
        else:
            return False

    def restart_relaxation(self, relaxed_file, yaw_rad, tblite_params=None, mh_params=None, totalsteps=20, output_filename_prefix="restarted"):
        """Restarts the relaxation from a given relaxed structure file."""
        try:
            atoms = read(relaxed_file)
            print(f"Atoms read from '{relaxed_file}': {len(atoms)} atoms")  # debug
            print(f"Length of yaw_rad: {len(yaw_rad)}")  # debug

            atoms_per_monomer = len(atoms) // len(yaw_rad)

            if len(atoms) % len(yaw_rad) != 0:
                print("Error: Atom count in restart file is not a multiple of the number of monomers. Aborting restart.")
                return False

            print(f"Restarting relaxation from '{relaxed_file}'.")
            relaxed_atoms = relax(atoms, output_filename_prefix=output_filename_prefix, tblite_params=tblite_params, mh_params=mh_params, totalsteps=totalsteps)
            write(f"{output_filename_prefix}_relaxed.xyz", relaxed_atoms)
            print(f"Restarted relaxed molecule written to '{output_filename_prefix}_relaxed.xyz'.")
            return True
        except FileNotFoundError:
            print(f"Error: Restart file '{relaxed_file}' not found.")
            return False
        except Exception as e:
            print(f"An error occurred during restart: {e}")
            return False
