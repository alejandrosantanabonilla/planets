from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from scipy.spatial.transform import Rotation as R
import parmed as pmd
from io import StringIO

class MolFrag:
    """
    Class for handling molecule fragments.
    """
    def __init__(self, mol, output_format='pdb', write_files=False):
        """
        Initializes MolFrag object with an RDKit molecule.

        Args:
            mol: An RDKit molecule object.
            output_format: The desired output format ('pdb', 'xyz', or 'mol').
            write_files: Whether to write the fragments to files (True/False).
        """
        self.mol = mol
        self.output_format = output_format
        self.write_files = write_files

    def get_mol_frags_blocks(self, pdb_filename_prefix="mol_"):
        """
        Splits a molecule into fragments, generates 3D coordinates,
        and returns a list of PDB blocks for each fragment.

        Args:
            pdb_filename_prefix: Prefix for the PDB filenames.

        Returns:
            A list of PDB blocks (strings) for each fragment.
        """

        mol_frags = Chem.GetMolFrags(self.mol, asMols=True)
        pdb_blocks = []

        # Define a dictionary to map output_format to the corresponding Chem function
        format_to_function = {
            'pdb': Chem.MolToPDBFile,
            'xyz': Chem.MolToXYZFile,
            'mol': Chem.MolToMolFile
        }

        for i, frag in enumerate(mol_frags):
            # Generate PDB block (this remains the same)
            pdb_block = Chem.MolToPDBBlock(frag)  
            pdb_blocks.append(pdb_block)

            if self.write_files:
                # Get the appropriate function from the dictionary
                write_function = format_to_function[self.output_format]  # Directly access using the key

                # Determine file extension
                file_extension = f".{self.output_format}"

                # Call the function to write the file
                write_function(frag, f"{pdb_filename_prefix}{i+1}{file_extension}")

        return pdb_blocks

    def get_mol_frags_parmed(self):
        """
        Splits a molecule into fragments, generates 3D coordinates,
        and returns a list of ParmEd Structure objects for each fragment.

        Returns:
            A list of ParmEd Structure objects for each fragment.
        """

        mol_frags = Chem.GetMolFrags(self.mol, asMols=True)
        parmed_structures = []

        for frag in mol_frags:
            # Create ParmEd Structure from RDKit molecule
            struct = pmd.load_rdkit(frag) 
            parmed_structures.append(struct)

        return parmed_structures

