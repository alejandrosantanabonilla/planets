import numpy as np
from ase.io import read
from ase.neighborlist import NeighborList

def generalized_coordination_number(atoms, species="Hg", cutoff=4.0, distance_weighting="inverse"):
    """
    Calculates a generalized coordination number for a specified species.

    Args:
        atoms (ase.Atoms): The Atoms object.
        species (str): The species for which to calculate the coordination number.
        cutoff (float): The cutoff distance for neighbors.
        distance_weighting (str): "inverse" or "gaussian" for distance weighting.

    Returns:
        list: A list of generalized coordination numbers for each atom of the specified species.
    """

    te_indices = [atom.index for atom in atoms if atom.symbol == species]
    neighbor_list = NeighborList([cutoff / 2] * len(atoms), skin=0, bothways=True, self_interaction=False)
    neighbor_list.update(atoms)

    generalized_coordination_numbers = []
    for te_index in te_indices:
        neighbors, offsets = neighbor_list.get_neighbors(te_index)
        te_pos = atoms.positions[te_index]
        generalized_cn = 0.0

        for neighbor_index, offset in zip(neighbors, offsets):
            neighbor_pos = atoms.positions[neighbor_index] + np.dot(offset, atoms.get_cell()) #Apply offset for periodic boundary conditions.
            dist_vec = neighbor_pos - te_pos
            dist = np.linalg.norm(dist_vec)

            if dist <= cutoff:
                if distance_weighting == "inverse":
                    generalized_cn += 1.0 / dist
                elif distance_weighting == "gaussian":
                    sigma = cutoff / 3.0 #Gaussian width
                    generalized_cn += np.exp(-0.5 * (dist / sigma)**2)
                else:
                    raise ValueError("Invalid distance_weighting option.")

        generalized_coordination_numbers.append(generalized_cn)

    return generalized_coordination_numbers

if __name__ == "__main__":
    try:
        atoms = read("last_rel.xyz")
    except FileNotFoundError:
        print("Please provide a valid xyz file named 'your_structure.xyz' in the same directory.")
        exit()

    cutoff = 5.0
    distance_weighting = "inverse" #"gaussian"  or "inverse"

    generalized_cns = generalized_coordination_number(atoms, species="Hg", cutoff=cutoff, distance_weighting=distance_weighting)

    print(f"Generalized Coordination Numbers for Te atoms (cutoff={cutoff} Å, weighting={distance_weighting}):")
    for i, gcn in enumerate(generalized_cns):
        print(f"Te atom {i+1}: {gcn:.3f}")
