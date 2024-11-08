import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.transform import Rotation as R
from openbabel import pybel

class MolFrag:
    """
    Class for handling molecule fragments.
    """
    def __init__(self, mol):
        """
        Initializes MolFrag object with an RDKit molecule.

        Args:
            mol: An RDKit molecule object.
        """
        self.mol = mol

    def write_mol_frags_to_pdb(self, pdb_filename_prefix="mol_"):
        """
        Splits a molecule into fragments, generates 3D coordinates, 
        and returns a list of PDB blocks for each fragment. 
        Also writes each fragment to a separate PDB file with names like 
        "mol_1.pdb", "mol_2.pdb", etc.

        Args:
            pdb_filename_prefix: Prefix for the PDB filenames.

        Returns:
            A list of PDB blocks (strings) for each fragment.
        """

        mol_frags = Chem.GetMolFrags(self.mol, asMols=True)
        pdb_blocks = []

        for i, frag in enumerate(mol_frags):
            # Generate PDB block
            pdb_block = Chem.MolToPDBBlock(frag)
            pdb_blocks.append(pdb_block)

            # Write to PDB file
            pdb_filename = f"{pdb_filename_prefix}{i+1}.pdb"
            Chem.MolToPDBFile(frag, pdb_filename)

        return pdb_blocks


class CreateDimer:
    """
    Class for creating and manipulating molecule dimers.
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

    def translate_molecule_cartesian(self, translation_vector):
        """
        Translates the first molecule along the Cartesian X, Y, and Z axes.

        Args:
            translation_vector: A numpy array of shape (3,) representing the 
                                translation vector in Cartesian coordinates (x, y, z).

        Returns:
            An RDKit molecule object with the translated conformation.
        """

        # Get coordinates of the molecule
        conf = self.mol1.GetConformer()
        coords = conf.GetPositions()

        # Apply translation to the molecule's coordinates
        new_coords = coords + translation_vector

        # Set the new coordinates for the molecule
        for i in range(len(coords)):
            conf.SetAtomPosition(i, new_coords[i])

        return self.mol1

    def rotate_molecule_around_plane(self, atom_indices, x_angle=0, y_angle=0, z_angle=0, 
                                     translation_vector=None):
        """
        Rotates the first molecule around a plane defined by the provided atom indices.
        The plane's normal vector is calculated using SVD.
        Additional rotations around X, Y, and Z axes can also be applied using scipy.spatial.transform.Rotation.
        Optionally, the molecule can be translated along the plane's coordinate system.

        Args:
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
        conf = self.mol1.GetConformer()
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

        # --- Copy coordinates from mol2 ---
        conf2 = self.mol2.GetConformer()
        coords2 = conf2.GetPositions()
        for i in range(len(coords2)):
            conf.SetAtomPosition(i, coords2[i])
        # -----------------------------------

        # Apply optional translation along the plane
        if translation_vector is not None:
            self.mol1 = self.translate_molecule_cartesian(translation_vector)

        return self.mol1

    
# Example usage:
# 1. Load the molecule
mol = Chem.MolFromMolFile("dimer.mol", removeHs=False)

# 2. Use MolFrag to split and write fragments
molfrag = MolFrag(mol)
pdb_blocks = molfrag.write_mol_frags_to_pdb()

# 3. Get the first molecule from the PDB blocks
mol1 = Chem.MolFromPDBBlock(pdb_blocks[0], removeHs=False)
mol2 = Chem.MolFromPDBBlock(pdb_blocks[1], removeHs=False)

# 4. Create a dimer object
dimer = CreateDimer(mol1, mol2)

# 5. Define the indices of the 4 atoms for the plane
atom_indices = [4, 5, 9, 10]  # Example indices, change as needed

# 6. Define the translation vector (6.5 along the normal vector)
translation_vector = np.array([0.0, 0, 5.5]) 

# 7. Rotate and translate the molecule using scipy rotations
rotated_translated_mol = dimer.rotate_molecule_around_plane(
    atom_indices,
    x_angle=0, 
    y_angle=0, 
    z_angle=5, 
    translation_vector=translation_vector
)

# 8. Write the rotated and translated molecule to a PDB file
Chem.MolToPDBFile(rotated_translated_mol, 'rotated_translated_molecule.pdb')


