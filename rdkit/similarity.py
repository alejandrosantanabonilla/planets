import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.ML.Cluster import Butina
import pathos.multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor


def read_smiles_from_file(input_file, chunksize=10000):
    """
    Reads SMILES strings from a whitespace-delimited file.

    Args:

      input_file: Path to the input file.

      chunksize: Number of rows to read at a time. Useful for large files.

    Returns:

      A list of SMILES strings.

    """

    smiles_list = []
    for chunk in pd.read_csv(input_file, sep='\s+', header=None,
                             engine='python', chunksize=chunksize, names=['Molecule']):
        smiles_list.extend(chunk['Molecule'].tolist())

    return smiles_list


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


def calculate_fingerprints_threaded(smiles_list, num_cores):
    """Calculates RDKit fingerprints using a thread pool."""
    fingerprints = []
    with ThreadPoolExecutor(max_workers=num_cores) as executor:
        # Use executor.map to apply calculate_fingerprints to each SMILES in parallel
        results = list(executor.map(calculate_fingerprints, smiles_list))
    fingerprints = results  # Assign the results to the fingerprints list
    return fingerprints 

    
def similarity_matrix_parallel(smiles_list, num_cores, similarity_func_name):
    """
    Calculates the similarity matrix in parallel using threads.
    """
    fingerprints = calculate_fingerprints_threaded(smiles_list, num_cores)

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


def create_cluster_dataframe(smiles_list, cluster_labels):
    """
    Creates a Pandas DataFrame with cluster labels and corresponding SMILES strings.
    """

    df = pd.DataFrame({'Cluster': cluster_labels, 'SMILES': smiles_list})

    return df


if __name__ == "__main__":
    num_cores = 1
    input_file = "total.dat"
    chunksize = 240
    cutoff = 0.15  # Butina cutoff
    similarity_func_name = 'FingerprintSimilarity'
    smiles_list = read_smiles_from_file(input_file, chunksize=chunksize)
    similarity_matrix = similarity_matrix_parallel(smiles_list, num_cores, similarity_func_name)
    cluster_labels = cluster_molecules_butina(similarity_matrix, cutoff=cutoff)

    # Create the DataFrame
    df_clusters = create_cluster_dataframe(smiles_list, cluster_labels)

    # Print cluster information
    num_clusters = df_clusters['Cluster'].nunique()

    print(f"Number of clusters found: {num_clusters}")

    cluster_counts = df_clusters['Cluster'].value_counts()
    print("Cluster sizes:")
    print(cluster_counts)
  
    # Save to CSV
    df_clusters.to_csv('clustered_molecules.csv', index=False)
    print("Clustered molecules saved to clustered_molecules.csv")
