from rdkit import Chem
from rdkit.Chem import AllChem

from pysoftk.linear_polymer.linear_polymer import *
from pysoftk.format_printers.format_mol import *

#5 times the monomer
#a=Chem.MolFromSmiles('c1c([C@H](CBr)C[C@@H](c2ccccc2)C[C@@H](c2ccccc2)C[C@H](c2ccccc2)Br)cccc1')
a=Chem.MolFromSmiles('c1c([C@@H](CBr)Br)cccc1')
#Embedding is needed for being parsed as a pysoftk.object
AllChem.EmbedMolecule(a)

new=Lp(a,"Br", 40, shift=0.95).linear_polymer("MMFF", 20000, 1, no_att=False)
Fmt(new).xyz_print("styrene_pol.xyz")
