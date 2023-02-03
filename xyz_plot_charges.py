from rdkit import Chem
from rdkit.Chem import AllChem

def show_png(data):
    """Function to plot a RDKit similarity map
    """
    import io
    from PIL import Image
    
    bio = io.BytesIO(data)
    img = Image.open(bio)

    return img.show()

def read_xyz(path_file, charge=0):
    """Function to read an XYZ file. One needs to give the 
       full path.
    """
    from rdkit.Chem import rdDetermineBonds
    
    raw_mol=Chem.MolFromXYZFile(str(path_file))
    conn_mol = Chem.Mol(raw_mol)
    rdDetermineBonds.DetermineBonds(conn_mol, charge=int(charge))

    return conn_mol

def geister_charges(mol):
    """Function to compute the Geister Charges
    """
    from rdkit.Chem import rdPartialCharges
    
    rdPartialCharges.ComputeGasteigerCharges(mol)
    return [x.GetDoubleProp("_GasteigerCharge") for x in mol.GetAtoms()]

def extended_huckel(mol):
    """ Function to compute charges with RDKit Extended Huckel model
    """
    from rdkit.Chem import rdEHTTools
    from rdkit.Chem import rdDistGeom
    
    _,res = rdEHTTools.RunMol(mol)
    static_chgs = res.GetAtomicCharges()[:mol.GetNumAtoms()]

    return list(static_chgs)

def plot_mol(mol, charges, size_image):
    """ Function to plot molecules with attribute charges
    """
    from rdkit.Chem import Draw
    from rdkit.Chem.Draw import SimilarityMaps

    a,b=size_image
    d = Draw.MolDraw2DCairo(int(a), int(b))
    SimilarityMaps.GetSimilarityMapFromWeights(mol, charges, contourLines=10, draw2d=d)
    d.FinishDrawing()
    show_png(d.GetDrawingText())

def random_charges(size):
    """ Function to simulate provide charges
        between -1 and 1
    """
    import random
    
    return [random.uniform(-1, 1) for _ in range(int(size))]

    
if __name__ == "__main__":
    mol = read_xyz("./mol.xyz")
    #charges=geister_charges(mol)
    #charges=extended_huckel(mol)

    charges=random_charges(63)
    
    plot_mol(mol, charges, (1000,1000))
