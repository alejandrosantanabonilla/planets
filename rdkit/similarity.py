import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.ML.Cluster import Butina
import pathos.multiprocessing as mp

def calculate_fingerprints(smiles):
    """Calculates RDKit fingerprint for a given SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return FingerprintMols.FingerprintMol(mol)

def apply_rdkit_similarity(func_name, fp1, fp2):
    """
    Applies a specified RDKit similarity function to two fingerprints.
    """
    try:
        func = getattr(DataStructs, func_name)
        return func(fp1, fp2)
    except AttributeError:
        raise ValueError(f"Invalid RDKit function name: {func_name}")

def similarity_matrix_parallel(smiles_list, num_cores, similarity_func_name):
    """
    Calculates the similarity matrix in parallel.
    """
    with mp.Pool(num_cores) as pool:
        fingerprints = pool.map(calculate_fingerprints, smiles_list)

    size = len(fingerprints)
    similarity_matrix = np.zeros((size, size))

    for i in range(size):
        for j in range(i, size):
            fp1 = fingerprints[i]
            fp2 = fingerprints[j]
            if fp1 is None or fp2 is None:
                similarity = 0.0
            else:
                similarity = apply_rdkit_similarity(similarity_func_name, fp1, fp2)
                if similarity_func_name == 'BulkTanimotoSimilarity':
                    similarity = similarity[0]

            similarity_matrix[i, j] = similarity
            similarity_matrix[j, i] = similarity

    return similarity_matrix

def cluster_molecules_butina(similarity_matrix, cutoff=0.35):
    """
    Clusters molecules using the Butina algorithm based on a similarity matrix.
    """
    num_molecules = similarity_matrix.shape[0]
    dists = []
    for i in range(1, num_molecules):
        sims = similarity_matrix[i, :i]  # Extract similarities from the matrix
        dists.extend([1 - x for x in sims])

    mol_clusters = Butina.ClusterData(dists, num_molecules, cutoff, isDistData=True)
    cluster_id_list = [0] * num_molecules
    for idx, cluster in enumerate(mol_clusters, 1):
        for member in cluster:
            cluster_id_list[member] = idx
    return cluster_id_list

def create_similar_molecules_list(cluster_labels):
    """
    Creates lists of similar molecules based on cluster labels.
    """
    similar_molecules = {}
    for i, label in enumerate(cluster_labels):
        if label not in similar_molecules:
            similar_molecules[label] = []
        similar_molecules[label].append(i)
    return similar_molecules

if __name__ == "__main__":
    num_cores = 2
    input_file = "total.dat"
    chunksize = 240
    cutoff = 0.35  # Butina cutoff
    similarity_func_name = 'FingerprintSimilarity'

    smiles_list = []
    for chunk in pd.read_csv(input_file, delim_whitespace=True, header=None, 
                             engine='python', chunksize=chunksize, names=['Molecule']):
        smiles_list.extend(chunk['Molecule'].tolist())

    similarity_matrix = similarity_matrix_parallel(smiles_list, num_cores, similarity_func_name)

    cluster_labels = cluster_molecules_butina(similarity_matrix, cutoff=cutoff)
    similar_molecules = create_similar_molecules_list(cluster_labels)

    for label, molecules in similar_molecules.items():
        print(f"Cluster {label}: {molecules}")
