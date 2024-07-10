import rdkit
from rdkit import Chem

# create molecule from SMILES
mol = Chem.MolFromSmiles('C1CCCCN1')

# find index of nitrogen
nIdx = mol.GetSubstructMatch(Chem.MolFromSmarts('N'))
# set the nitrogen atom's formal charge to 1
mol.GetAtomWithIdx(*nIdx).SetFormalCharge(1)

# print the protonated SMILES
print(Chem.MolToMolBlock(mol))
