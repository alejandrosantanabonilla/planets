import numpy as np
from scipy.spatial.transform import Rotation as R
from rdkit import Chem
from rdkit.Chem import AllChem

class CreateDimer:
    """
    Class for creating and manipulating molecule dimers.
    This class focuses on applying transformations to the second molecule (mol2) 
    while keeping the first molecule (mol1) unchanged.
    """
    def __init__(self, mol1, mol2):
        """
        Initializes CreateDimer object with two RDKit molecules.

        Args:
            mol1: The first RDKit molecule object.
            mol2: The second RDKit molecule object.
        """
        self.mol1 = mol1
        self.mol2 = mol2

    def translate_molecule_cartesian(self, mol, translation_vector):
        """
        Translates the given molecule along the Cartesian X, Y, and Z axes.

        Args:
            mol: The RDKit molecule object to translate.
            translation_vector: A numpy array of shape (3,) representing the 
                                    translation vector in Cartesian coordinates (x, y, z).

        Returns:
            An RDKit molecule object with the translated conformation.
        """

        # Get coordinates of the molecule
        conf = mol.GetConformer()
        coords = conf.GetPositions()

        # Apply translation to the molecule's coordinates
        new_coords = coords + translation_vector

        # Set the new coordinates for the molecule
        for i in range(len(coords)):
            conf.SetAtomPosition(i, new_coords[i])

        return mol

    def rotate_molecule_around_plane(self, mol, atom_indices, x_angle=0, y_angle=0, z_angle=0, 
                                                translation_vector=None):
        """
        Rotates the given molecule around a plane defined by the provided atom indices.
        The plane's normal vector is calculated using SVD.
        Additional rotations around X, Y, and Z axes can also be applied using scipy.spatial.transform.Rotation.
        Optionally, the molecule can be translated along the plane's coordinate system.

        Args:
            mol: The RDKit molecule object to rotate.
            atom_indices: A list of atom indices (0-based) defining the plane.
            x_angle: Rotation angle around the X-axis in degrees.
            y_angle: Rotation angle around the Y-axis in degrees.
            z_angle: Rotation angle around the Z-axis in degrees.
            translation_vector: A numpy array of shape (3,) representing the translation vector 
                                in the plane's coordinate system.

        Returns:
            An RDKit molecule object with the rotated and optionally translated conformation.
        """

        # Get coordinates of the molecule
        conf = mol.GetConformer()
        coords = conf.GetPositions()

        # Select atoms defining the plane
        plane_atoms = coords[atom_indices]

        # Compute the best-fit plane using SVD
        centroid = np.mean(plane_atoms, axis=0)
        plane_atoms -= centroid
        _, _, vh = np.linalg.svd(plane_atoms)
        normal_vector = vh[-1, :]  # Normal vector of the plane

        # Create a rotation object using scipy.spatial.transform.Rotation
        rotation = R.from_rotvec(np.pi * normal_vector)  # Rotate by 180 degrees around the normal vector

        # Apply additional rotations around X, Y, and Z axes
        rotation = rotation * R.from_euler('x', x_angle, degrees=True)
        rotation = rotation * R.from_euler('y', y_angle, degrees=True)
        rotation = rotation * R.from_euler('z', z_angle, degrees=True)

        # Apply rotation to the molecule's coordinates
        new_coords = rotation.apply(coords - centroid) + centroid

        # Set the new coordinates for the molecule
        for i in range(len(coords)):
            conf.SetAtomPosition(i, new_coords[i])

        # Apply optional translation along the plane
        if translation_vector is not None:
            mol = self.translate_molecule_cartesian(mol, translation_vector)

        return mol

# Load molecules (replace with your actual molecule loading)
mol1 = Chem.MolFromPDBFile('fragment_1.pdb', removeHs=False)
mol2 = Chem.MolFromPDBFile('fragment_2.pdb', removeHs=False)

# Create a dimer object
dimer = CreateDimer(mol1, mol2)

# Define the translation vector
translation_vector = np.array([0, 0, 6.5])

# Define the atom indices for the plane
atom_indices = [0, 1, 2]  # Example indices

# Perform rotation and translation on mol2
transformed_mol2 = dimer.rotate_molecule_around_plane(
    dimer.mol2, atom_indices, x_angle=45, translation_vector=translation_vector
)

print(Chem.MolToPDBBlock(transformed_mol2))

