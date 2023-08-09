from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import TorsionFingerprints

import sys
import itertools

def mol_neigh(mol):
    neigh=[(atom.GetIdx(),
            [(nbr.GetIdx())
             for nbr in atom.GetNeighbors()])
           for atom in mol.GetAtoms()]

    return neigh

def mmf94_atomtypes(mol, file_name="MMFF94_atomtypes.dat"):
    myfile = open(str(file_name), 'w')
    properties = AllChem.MMFFGetMoleculeProperties(mol)
    properties.SetMMFFVerbosity(2)

    for atom in range(mol.GetNumAtoms()):
        myfile.write("%s %s %s\n" % (mol.GetAtomWithIdx(atom).GetSymbol(),
                                     int(atom+1), int(properties.GetMMFFAtomType(atom))))

    myfile.close()


def mmmf94_param(mol, filename="MMFF94_param.dat"):

    import sys, os

    if AllChem.MMFFHasAllMoleculeParams(mol) is False:
        raise InputError("RDKit parameters not found for all atom types in molecule.")
    
    #This is done to redirect the output to a file
    sys.stdout.flush()
    file=open(str(filename), "w")
    os.dup2(file.fileno(), 1)
    sys.stdout = file

    properties = AllChem.MMFFGetMoleculeProperties(mol)
    properties.SetMMFFVerbosity(2)

    return AllChem.MMFFGetMoleculeForceField(mol, properties)

if __name__ == "__main__":

    mol=Chem.MolFromPDBFile("test.pdb", removeHs=False)

    #This is for getting a tabulated and customized MMFF94 parametrization for mol. 
    #mmf94_atomtypes(mol)
    #mmmf94_param(mol)

    
    atom_pairs= sum([list(itertools.product([a], b)) for a,b in mol_neigh(mol)], [])

    bond_stretch_param= [(values, AllChem.GetUFFBondStretchParams(mol, values[0], values[1]))
                          for idx, values in enumerate(atom_pairs)]

    vdw_params= [(values, AllChem.GetUFFVdWParams(mol, values[0], values[1]))
                          for idx, values in enumerate(atom_pairs)]

   
    tl, tlr = TorsionFingerprints.CalculateTorsionLists(mol)

    if len(tl) > 1:
        for a,b in tl:
            x,y,z,k=a[0]
