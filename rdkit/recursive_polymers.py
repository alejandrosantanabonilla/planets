from rdkit import Chem
from openbabel import pybel as pb

def generate_recursive_smiles(smiles, replacements, max_iterations=10, final_replacement="*"):
    """
    Generates a recursive SMILES string with a stopping condition and a final replacement.

    Args:
        smiles: The base SMILES string with a placeholder.
        replacements: A dictionary of replacement strings.
        max_iterations: The maximum number of recursive substitutions.
        final_replacement: The SMILES string or SMART to replace the placeholder after recursion.

    Returns:
        A fully resolved SMILES string.
    """

    for _ in range(max_iterations):
        new_smiles = smiles.format(**replacements)
        if new_smiles == smiles:  # Stop if the SMILES string doesn't change
            break
        smiles = new_smiles

    # Replace the placeholder with the final_replacement
    smiles = smiles.replace("{R}", final_replacement)

    return smiles


# Example to create a chain of 10 benzene rings in a straight line
user_smiles = "c1ccccc1{R}"  # Starting with a benzene ring
user_replacements = {"R": "c2ccc(cc2){R}"}  # Adding another benzene in para position
user_final_replacement = ""  # No final replacement needed

# Generate the recursive SMILES string
smiles_string = generate_recursive_smiles(
    smiles=user_smiles, 
    replacements=user_replacements, 
    max_iterations=100,  # 9 iterations to add 9 more benzene rings
    final_replacement=user_final_replacement
)

# Create the RDKit molecule from the resolved SMILES string
mol = Chem.MolFromSmiles(smiles_string)

# Print the SMILES string of the generated molecule
mol=Chem.MolToSmiles(mol)
mol_pb=pb.readstring('smiles', mol)
mol_pb.make3D() # generate 3D coordinates

# now perform conformer searching
ff = pb._forcefields["mmff94"]
success = ff.Setup(mol_pb.OBMol)
if not success:
    ff = pb._forcefields["uff"]
    success = ff.Setup(mol_pb.OBMol)
    if not success:
        sys.exit("Cannot set up forcefield")

# feel free to tweak these to your balance of time / quality
ff.ConjugateGradients(100, 1.0e-3)
ff.FastRotorSearch(True) # permute central bonds
ff.WeightedRotorSearch(100, 25) # 100 cycles, each with 25 forcefield ops
# final optimization
ff.ConjugateGradients(250, 1.0e-4)
# update the coordinates
ff.GetCoordinates(mol_pb.OBMol) 

mol_pb.write("mol", "optimized.mol", overwrite=True)
