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
    index = 0 
    while index <= m.GetNumAtoms()-1:
        atomObj = m.GetAtomWithIdx(index)
        bond_kind = m.GetBondWithIdx(index).GetBondType()
        symbol = str(atomObj.GetSymbol())
        Hybridization = str(atomObj.GetHybridization())
        print(symbol + " hybridization: " + Hybridization + "bond: " + str(bond_kind))
        index+=1

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


print (mat.shape)
print (ov.shape)
