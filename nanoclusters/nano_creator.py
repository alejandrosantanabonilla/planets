from ase.io import read, write
import numpy as np
from ase.visualize import view

def create_nanoparticle_from_cif(cif_file, size=10.0):
    """
    Creates a spherical nanoparticle from a CIF file.

    Args:
        cif_file (str): Path to the CIF file.
        size (float): Radius of the nanoparticle in Angstroms.

    Returns:
        ase.Atoms: The nanoparticle Atoms object.
    """
    try:
        atoms = read(cif_file)
    except FileNotFoundError:
        print(f"Error: CIF file not found at {cif_file}")
        return None

    # Determine a good supercell size.
    lattice_constant = np.max(atoms.get_cell().lengths()) #Aproximate lattice constant.
    supercell_size = int(np.ceil(size * 2 / lattice_constant))

    supercell = atoms * (supercell_size, supercell_size, supercell_size)

    center = supercell.get_center_of_mass()
    positions = supercell.get_positions()
    distances = np.linalg.norm(positions - center, axis=1)
    mask = distances <= size
    nanoparticle = supercell[mask]
    nanoparticle.center()

    return nanoparticle

# Example usage:
cif_file_path = "znHgTe.cif"  # Replace with your CIF file path
nanoparticle = create_nanoparticle_from_cif(cif_file_path, size=75.0)

if nanoparticle is not None:
    write("system_from_cif.xyz", nanoparticle)
    view(nanoparticle)
    print(f"Nanoparticle contains {len(nanoparticle)} atoms.")
