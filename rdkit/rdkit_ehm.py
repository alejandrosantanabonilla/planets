from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor

# this is the package including the connection to YAeHMOP
from rdkit.Chem import rdEHTTools
from rdkit.Chem import AllChem
import numpy as np

import rdkit
import time

def set_overlap_populations(m,ehtRes):
    rop = ehtRes.GetReducedOverlapPopulationMatrix()
    for bnd in m.GetBonds():
        a1 = bnd.GetBeginAtom()
        a2 = bnd.GetEndAtom()
        if a1.GetAtomicNum()==1:
            continue
        if a2.GetAtomicNum()==1:
            continue
        # symmetric matrix:
        i1 = max(a1.GetIdx(),a2.GetIdx())
        i2 = min(a1.GetIdx(),a2.GetIdx())
        idx = (i1*(i1+1))//2 + i2
        bnd.SetDoubleProp("MullikenOverlapPopulation",rop[idx])

def bonds_and_indexes(m):
    """
    Prints the symbol, hybridization, and bond type for each atom in a molecule.
    
    Args:
        m: An RDKit Mol object representing the molecule.
    """

    for atom in m.GetAtoms():
        symbol = atom.GetSymbol()
        hybridization = atom.GetHybridization()

        # Check if the atom has a bond before accessing bond properties
        bond_kind = ""  # Default to an empty string if no bond exists
        if atom.GetBonds(): 
            bond_kind = atom.GetBonds()[0].GetBondType()  # Get the first bond

        print(f"{symbol} hybridization: {hybridization}, bond: {bond_kind}")

def bonds_pops(mh):
    for bnd in list(mh.GetBonds()):
      if not bnd.HasProp("MullikenOverlapPopulation"):
          continue
      print(f'{bnd.GetIdx()} {bnd.GetBeginAtom().GetSymbol()}({bnd.GetBeginAtomIdx()})-{bnd.GetEndAtom().GetSymbol()}({bnd.GetEndAtomIdx()}) {bnd.GetBondType()} {bnd.GetDoubleProp("MullikenOverlapPopulation") :.3f}')


## Example for using the code

m=Chem.MolFromSmiles('CC')
mh=Chem.AddHs(m)
AllChem.EmbedMultipleConfs(mh, 1)

passed,res=rdEHTTools.RunMol(mh, keepOverlapAndHamiltonianMatrices=True)

mat=res.GetHamiltonian()
ov=res.GetOverlapMatrix()
set_overlap_populations(mh,res)

print (bonds_and_indexes(mh))
