import os
import pandas as pd
import numpy as np

from itertools import product

from rdkit import Chem, DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

import pathos.multiprocessing as mp
import parmap
import cProfile


def read_molecules(input_file, chunksize):
    """Reads molecules from a SMILES file in chunks."""
    for chunk in pd.read_csv(input_file, delim_whitespace=True, header=None, 
                              engine='python', chunksize=chunksize, names=['Molecule']):
        yield chunk  # Yield chunks instead of storing them in a list

def calculate_fingerprints(data):
    """Calculates fingerprints for a chunk of molecules."""
    fps = []  # Create a list to store the fingerprints
    for row in data.itertuples():
        mol = Chem.MolFromSmiles(row.Molecule)
        fps.append(FingerprintMols.FingerprintMol(mol))  # Append to the list
    return fps  

def calculate_fingerprints_parallel(molecules, num_proc):
    """Calculates fingerprints in parallel."""
    pool = mp.ProcessingPool(nodes=num_proc)
    results = pool.map(calculate_fingerprints, molecules)
    pool.close()
    pool.join()
    return results

def generate_pairs(total_chunks):
    """Generates all possible pairs of chunk indices."""
    total_indices = range(total_chunks)
    for i in total_indices:
        for j in range(i, total_chunks):
            yield (i, j)

def calculate_similarity(fps_pair):
    """Calculates the similarity between two fingerprints."""
    fp1, fp2 = fps_pair
    return DataStructs.FingerprintSimilarity(fp1, fp2)

def process_chunk_pair(pair, total_fps, num_cores):
    """Calculates similarity for a pair of chunks."""
    i, j = pair
    fps1 = total_fps[i]
    fps2 = total_fps[j]
    results = parmap.map(calculate_similarity, product(fps1, fps2), pm_processes=num_cores)
    return np.array([1 - sim for sim in results])

def calculate_similarity_matrix(input_data, chunksize, num_cores):
    """Calculates the similarity matrix for the input data."""
    molecules = read_molecules(input_data, chunksize)
    total_fps = list(calculate_fingerprints_parallel(molecules, num_cores))
    num_chunks = len(total_fps)
    
    similarity_matrix = np.empty((num_chunks, num_chunks, chunksize, chunksize))

    for pair in generate_pairs(num_chunks):
        i, j = pair
        similarities = process_chunk_pair(pair, total_fps, num_cores)
        similarity_matrix[i, j] = similarities.reshape(chunksize, chunksize)
        similarity_matrix[j, i] = similarities.reshape(chunksize, chunksize).T

    return similarity_matrix


if __name__ == "__main__":
    #similarity_matrix = calculate_similarity_matrix("total.dat", 120, 2)
    cProfile.run('calculate_similarity_matrix("total.dat",1200, 2)', 'profiling_results')
    print ("finished")
    # Access and visualize parts of the matrix as needed
    #import matplotlib.pyplot as plt
    #plt.imshow(similarity_matrix[0, 0], cmap='viridis', interpolation='nearest')
    #plt.colorbar()
    #plt.show()
