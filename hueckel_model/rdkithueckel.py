# Import necessary libraries from RDKit
from rdkit import Chem
from rdkit.Chem import AllChem  # For general molecular operations
from rdkit.Chem import rdEHTTools  # For Extended Hückel Theory calculations
import numpy as np  # For numerical calculations


def set_overlap_populations(m, ehtRes):
    """
    Calculates and sets Mulliken overlap populations for bonds in a molecule.

    This function iterates through each bond in the molecule, excluding bonds involving hydrogen atoms.
    It retrieves the reduced overlap population matrix from the EHT results, and for each bond, 
    it calculates an index based on the atom indices involved. This index is used to fetch the 
    corresponding overlap population from the matrix, which is then stored as a property of the bond.

    Args:
        m: An RDKit Mol object representing the molecule.
        ehtRes: The results object from the Extended Hückel Theory calculation (rdEHTTools.RunMol).
    """
    rop = ehtRes.GetReducedOverlapPopulationMatrix()  # Get overlap matrix
    for bnd in m.GetBonds():
        a1, a2 = bnd.GetBeginAtom(), bnd.GetEndAtom()  # Get atoms involved in the bond

        # Skip bonds involving hydrogen atoms
        if a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1:
            continue

        # Calculate index for the overlap population in the matrix
        i1, i2 = max(a1.GetIdx(), a2.GetIdx()), min(a1.GetIdx(), a2.GetIdx())
        idx = (i1 * (i1 + 1)) // 2 + i2 

        # Set the calculated overlap population as a bond property
        bnd.SetDoubleProp("MullikenOverlapPopulation", rop[idx])


def bonds_and_indexes(m):
    """
    Prints the symbol, hybridization, and bond type for each atom in a molecule.

    Args:
        m: An RDKit Mol object representing the molecule.
    """

    for atom in m.GetAtoms():
        symbol = atom.GetSymbol()              # Get the atom's element symbol
        hybridization = atom.GetHybridization() # Get the atom's hybridization

        # Check if the atom has a bond before accessing bond properties
        bond_kind = ""  
        if atom.GetBonds(): 
            bond_kind = atom.GetBonds()[0].GetBondType()  # Get the type of the first bond

        print(f"{symbol} hybridization: {hybridization}, bond: {bond_kind}")


def bonds_pops(mh):
    """
    Prints the bond indices, atom symbols, bond types, and Mulliken overlap populations for each bond in a molecule.

    This function iterates through each bond in the molecule and checks if it has the "MullikenOverlapPopulation" property.
    If so, it prints information about the bond, including its index, the symbols and indices of the atoms it connects,
    its type, and the calculated overlap population.

    Args:
        mh: An RDKit Mol object representing the molecule with added hydrogens.
    """
    for bnd in list(mh.GetBonds()):
        if not bnd.HasProp("MullikenOverlapPopulation"):  # Check if bond has the property
            continue

        # Print bond information in a formatted way
        print(f'{bnd.GetIdx()} {bnd.GetBeginAtom().GetSymbol()}({bnd.GetBeginAtomIdx()})-'
              f'{bnd.GetEndAtom().GetSymbol()}({bnd.GetEndAtomIdx()}) {bnd.GetBondType()} '
              f'{bnd.GetDoubleProp("MullikenOverlapPopulation"):.3f}')

## Example for using the code

m=Chem.MolFromSmiles('CC')
mh=Chem.AddHs(m)
AllChem.EmbedMultipleConfs(mh, 1)

passed,res=rdEHTTools.RunMol(mh, keepOverlapAndHamiltonianMatrices=True)

mat=res.GetHamiltonian()
ov=res.GetOverlapMatrix()
set_overlap_populations(mh,res)


print (mat.shape)
print (ov.shape)
