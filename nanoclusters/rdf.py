import numpy as np
from ase.io import read
import matplotlib.pyplot as plt

def radial_distribution_function(atoms, r_max, dr, species_i, species_j=None, sigma=0.0):
    """
    Computes the radial distribution function (RDF) for a cluster with optional Gaussian smearing.

    Args:
        atoms (ase.Atoms): The Atoms object.
        r_max (float): Maximum radius for the RDF.
        dr (float): Bin width for the RDF.
        species_i (str): Symbol of the central species.
        species_j (str, optional): Symbol of the neighbor species. If None, it's the same as species_i.
        sigma (float, optional): Standard deviation for Gaussian smearing. Defaults to 0.0 (no smearing).

    Returns:
        tuple: (r_values, g_r) where:
            r_values (numpy.ndarray): Array of radial distances.
            g_r (numpy.ndarray): Array of RDF values.
    """

    if species_j is None:
        species_j = species_i

    indices_i = [atom.index for atom in atoms if atom.symbol == species_i]
    indices_j = [atom.index for atom in atoms if atom.symbol == species_j]

    if not indices_i or not indices_j:
        raise ValueError("Species not found in the Atoms object.")

    n_bins = int(r_max / dr)
    g_r = np.zeros(n_bins)
    r_values = np.linspace(dr / 2, r_max - dr / 2, n_bins)

    for i_index in indices_i:
        for j_index in indices_j:
            if i_index == j_index and species_i == species_j:
                continue
            dist_vec = atoms.positions[j_index] - atoms.positions[i_index]
            dist = np.linalg.norm(dist_vec)

            if dist < r_max:
                if sigma == 0.0:  # No smearing
                    bin_index = int(dist / dr)
                    g_r[bin_index] += 1
                else:  # Gaussian smearing
                    for bin_index, r_bin in enumerate(r_values):
                        g_r[bin_index] += np.exp(-0.5 * ((dist - r_bin) / sigma)**2)

    number_i = len(indices_i)
    number_j = len(indices_j)

    if number_i == 0 or number_j == 0:
        return r_values, g_r

    # Normalization for a cluster.
    for i, r in enumerate(r_values):
        shell_volume = 4 * np.pi * r**2 * dr
        g_r[i] /= number_i * number_j / (r_max**3) * shell_volume # Approximate normalisation for cluster.

    return r_values, g_r

if __name__ == "__main__":
    try:
        atoms = read("no_rel.xyz")
    except FileNotFoundError:
        print("Please provide a valid xyz file named 'last_rel.xyz' in the same directory.")
        exit()

    r_max = 10.0
    dr = 0.1
    species_i = "Te"
    species_j = "Hg"
    sigma = 0.2  # Smearing width

    r_values, g_r = radial_distribution_function(atoms, r_max, dr, species_i, species_j, sigma=sigma)

    plt.plot(r_values, g_r)
    plt.xlabel("r (Å)")
    plt.ylabel("g(r)")
    plt.title(f"Radial Distribution Function ({species_i}-{species_j}) with smearing (sigma={sigma})")
    plt.grid(True)
    plt.show()

    r_values, g_r = radial_distribution_function(atoms, r_max, dr, "Te", sigma=sigma)

    plt.plot(r_values, g_r)
    plt.xlabel("r (Å)")
    plt.ylabel("g(r)")
    plt.title(f"Radial Distribution Function (Te-Te) with smearing (sigma={sigma})")
    plt.grid(True)
    plt.show()

    r_values, g_r = radial_distribution_function(atoms, r_max, dr, "Hg", sigma=sigma)

    plt.plot(r_values, g_r)
    plt.xlabel("r (Å)")
    plt.ylabel("g(r)")
    plt.title(f"Radial Distribution Function (Hg-Hg) with smearing (sigma={sigma})")
    plt.grid(True)
    plt.show()
