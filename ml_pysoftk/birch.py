import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import Birch
from sklearn.manifold import MDS

def mds_birch_clustering(tanimoto_distance_matrix, birch_threshold=0.5):
    """
    Performs MDS followed by BIRCH clustering on a Tanimoto distance matrix.

    Args:
        tanimoto_distance_matrix (numpy.ndarray): The Tanimoto distance matrix (nxn).
        birch_threshold (float): The threshold for BIRCH clustering.

    Returns:
        numpy.ndarray: Cluster labels assigned by BIRCH.
    """

    # MDS projection
    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
    mds_embedding = mds.fit_transform(tanimoto_distance_matrix)

    # BIRCH clustering on MDS embedding
    birch = Birch(n_clusters=None, threshold=birch_threshold, compute_labels=True)
    birch.fit(mds_embedding)
    labels = birch.labels_

    # Visualization (optional)
    plt.figure(figsize=(10, 8))
    plt.scatter(mds_embedding[:, 0], mds_embedding[:, 1], c=labels, cmap='viridis')
    plt.title("MDS + BIRCH Clustering")
    plt.xlabel("MDS Dimension 1")
    plt.ylabel("MDS Dimension 2")
    plt.colorbar(label="Cluster Label")
    plt.show()

    return labels

# Example Usage (replace with your Tanimoto distance matrix)
n_samples = 100
tanimoto_distance_matrix = np.random.rand(n_samples, n_samples)
tanimoto_distance_matrix = (tanimoto_distance_matrix + tanimoto_distance_matrix.T) / 2

cluster_labels = mds_birch_clustering(tanimoto_distance_matrix, birch_threshold=0.3)

# Print cluster labels
print("Cluster Labels:", cluster_labels)
