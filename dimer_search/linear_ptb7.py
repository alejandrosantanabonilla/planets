from rdkit import Chem
from rdkit.Chem import AllChem

from pysoftk.linear_polymer.linear_polymer import *
from pysoftk.format_printers.format_mol import *

#5 times the monomer
a=Chem.MolFromPDBFile('mol.pdb',removeHs=False)

#Embedding is needed for being parsed as a pysoftk.object
AllChem.EmbedMolecule(a)

new=Lp(a,"Br", 2, shift=1.0).linear_polymer("MMFF", 70000, 45, no_att=True)
Fmt(new).xyz_print("ptb7_pol2.xyz")
