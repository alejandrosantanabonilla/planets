from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def align_molecule(mol):
    # 1. Calculate Center of Mass
    conf = mol.GetConformer()
    positions = conf.GetPositions()
    masses = []
    for atom in mol.GetAtoms():
        masses.append(atom.GetMass())
    center_mass = np.average(positions, axis=0, weights=masses)

    # 2. Translate Molecule to Origin
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, pos - center_mass)

    # 3. Calculate Inertia Tensor
    inertia = np.zeros((3, 3))
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        mass = masses[i]  # Use the masses list from above
        inertia += mass * (np.dot(pos, pos) * np.eye(3) - np.outer(pos, pos))

    # 4. Calculate Eigenvalues & Eigenvectors
    eigenvals, eigenvecs = np.linalg.eig(inertia)

    # 5. Sort Eigenvalues (Ascending)
    order = np.argsort(eigenvals)
    eigenvals = eigenvals[order]
    eigenvecs = eigenvecs[:, order]

    # 6. Rotation Matrix (Eigenvectors as Columns)
    rotation_matrix = eigenvecs

    # 7. Apply Rotation to All Atoms
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        new_pos = np.dot(rotation_matrix, pos)
        conf.SetAtomPosition(i, new_pos)

    return mol


# Example Usage
mol = Chem.MolFromXYZFile("mol.xyz")  # Replace with your file
aligned_mol = align_molecule(mol)

# Write Aligned Molecule (Optional)
Chem.MolToXYZFile(aligned_mol, "aligned_molecule.xyz")
