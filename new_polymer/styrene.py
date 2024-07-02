from rdkit import Chem
from rdkit.Chem import AllChem

from pysoftk.linear_polymer.linear_polymer import Lp
from pysoftk.format_printers.format_mol import Fmt
from pysoftk.linear_polymer.mol_conformer import Mcon
from pysoftk.topologies.diblock import Pt

from openbabel import pybel as pb 


def pol_psam(mol1_smiles, mol2_smiles, mol1_length, mol2_length):
    """
    Generates a PSAM (Polystyrene-alt-Maleic Anhydride) polymer structure.

    Args:
        mol1_smiles (str): SMILES string of the first monomer (styrene).
        mol2_smiles (str): SMILES string of the second monomer (maleic acid).
        mol1_length (int): Number of repeating units of the first monomer.
        mol2_length (int): Number of repeating units of the second monomer.
    """

    # Create RDKit molecules from SMILES
    mol1 = Chem.MolFromSmiles(mol1_smiles)
    mol2 = Chem.MolFromSmiles(mol2_smiles)

    # Embed molecules for 3D coordinates
    AllChem.EmbedMolecule(mol1)
    AllChem.EmbedMolecule(mol2)

    # Generate linear polymers using pysoftk (force field, steps, temp, no_att=False for minimization)
    sty = Lp(mol1, "Br", mol1_length, shift=1.4).linear_polymer("MMFF", 8500, 175, no_att=False)
    mal = Lp(mol2, "Br", mol2_length, shift=1.4).linear_polymer("MMFF", 6500, 50, no_att=False)

    # Save polymers as .mol files
    Fmt(sty).mol_print(f"styrene_pol{mol1_length}.mol")
    Fmt(mal).mol_print(f"mal_pol{mol2_length}.mol")

    sty.addh()
    sty.make3D()
    new_sty=Chem.MolFromMolBlock(sty.write("mol"))

    mal.addh()
    mal.make3D()
    new_mal=Chem.MolFromMolBlock(mal.write("mol"))

    # Combine the polymers into a diblock copolymer (PSAM)
    total = Pt('AB', [new_mal, new_sty], "Br").pattern_block_poly(swap_H=True)

    total.addh()
    total.make3D()
    #total.write("pdb", f"psam_{mol1_length}_{mol2_length}.pdb")
    # Save the PSAM polymer as a .pdb file
    Fmt(total).pdb_print(f"psam_{mol1_length}_{mol2_length}.pdb")

    

# Example Usage
if __name__ == "__main__":
    mol1_smiles = 'c1c([C@@H](CBr)Br)cccc1'  # Styrene with terminal Br
    mol2_smiles = '[C@@H]1(C(=O)OC(=O)[C@H]1Br)Br'  # Maleic anhydride with terminal Br

    pol_psam(mol1_smiles, mol2_smiles, 40, 40)


