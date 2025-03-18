import numpy as np
from sklearn.cluster import Birch
from sklearn.metrics import pairwise_distances

def birch_cluster_tanimoto(tanimoto_distance_matrix, threshold=0.5):
    """
    Performs BIRCH clustering based on the indices of a Tanimoto distance matrix.

    Args:
        tanimoto_distance_matrix (numpy.ndarray): The Tanimoto distance matrix.
        threshold (float): The threshold for BIRCH clustering.

    Returns:
        dict: A dictionary where keys are cluster labels and values are lists of indices
              corresponding to the members of each cluster.
    """

    # Assuming the distance matrix is symmetric and represents pairwise distances
    n_samples = tanimoto_distance_matrix.shape[0]

    # Generate a dummy dataset of indices (0 to n_samples - 1)
    data_indices = np.arange(n_samples).reshape(-1, 1) #Reshape to make it 2d for sklearn

    # BIRCH clustering
    birch = Birch(n_clusters=None, threshold=threshold, compute_labels=True)
    birch.fit(data_indices)  # Birch is fit to the indices.
    labels = birch.labels_

    # Create a dictionary to store cluster members
    clusters = {}
    unique_labels = np.unique(labels)

    for cluster_label in unique_labels:
        cluster_members = np.where(labels == cluster_label)[0].tolist()
        clusters[cluster_label] = cluster_members

    return clusters

# Example Usage (replace with your Tanimoto distance matrix)
# Create a dummy tanimoto distance matrix.
data_points = np.random.rand(100, 20) # example data
tanimoto_distance_matrix = pairwise_distances(data_points, metric='jaccard') #jaccard is used as a proxy for tanimoto.

cluster_results = birch_cluster_tanimoto(tanimoto_distance_matrix, threshold=0.3)

# Print the cluster results
for cluster_label, members in cluster_results.items():
    print(f"Cluster {cluster_label}: {members}")
