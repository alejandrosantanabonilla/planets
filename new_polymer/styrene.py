from rdkit import Chem
from rdkit.Chem import AllChem

from pysoftk.linear_polymer.linear_polymer import *
from pysoftk.format_printers.format_mol import *
from pysoftk.linear_polymer.mol_conformer import *
from pysoftk.topologies.diblock import *

#5 times the monomer
a=Chem.MolFromSmiles('c1c([C@@H](CBr)Br)cccc1')
b=Chem.MolFromSmiles('[C@@H]1(C(=O)OC(=O)[C@H]1Br)Br')

#Embedding is needed for being parsed as a pysoftk.object
AllChem.EmbedMolecule(a)
AllChem.EmbedMolecule(b)

sty=Lp(a,"Br", 4, shift=1.4).linear_polymer("MMFF", 8500, 150, no_att=False)
mal=Lp(b,"Br", 4, shift=1.4).linear_polymer("MMFF", 6500, 75, no_att=False)

Fmt(sty).mol_print("styrene_pol4.mol")
Fmt(mal).mol_print("mal_pol4.mol")

new_sty=Chem.MolFromMolFile("styrene_pol4.mol",removeHs=False)
new_mal=Chem.MolFromMolFile("mal_pol4.mol",removeHs=False)

AllChem.EmbedMolecule(new_sty)
AllChem.EmbedMolecule(new_mal)

Mcon(new_sty,3).conformer("conformers")

confs=[]
mols=Chem.SDMolSupplier('conformers.sdf')
for i, mol in enumerate(mols):
    confs.append(Chem.MolToMolBlock(mol))

sty_gs=Chem.MolFromMolBlock(confs[0])
AllChem.EmbedMolecule(sty_gs)

mols=[new_mal, sty_gs]
total=Pt('AB', mols, "Br").pattern_block_poly(swap_H=True)
Fmt(total).pdb_print("psam_4pol.pdb")

