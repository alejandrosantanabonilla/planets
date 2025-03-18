import numpy as np
import umap
import hdbscan
import matplotlib.pyplot as plt

def umap_hdbscan_clustering_from_distance(distance_matrix, n_neighbors=15, min_dist=0.1, min_cluster_size=10):
    """
    Performs UMAP dimensionality reduction and HDBSCAN clustering using a precomputed distance matrix.

    Args:
        distance_matrix (numpy.ndarray): The precomputed distance matrix (e.g., Tanimoto distance).
        n_neighbors (int): UMAP n_neighbors parameter.
        min_dist (float): UMAP min_dist parameter.
        min_cluster_size (int): HDBSCAN min_cluster_size parameter.

    Returns:
        tuple: (UMAP embedding, HDBSCAN cluster labels)
    """

    # Apply UMAP
    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric="precomputed"  # Crucial for using the distance matrix
    )
    embedding = reducer.fit_transform(distance_matrix)

    # Apply HDBSCAN clustering
    clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size)
    labels = clusterer.fit_predict(embedding)

    return embedding, labels

# Example Usage (replace with your actual Tanimoto distance matrix)
# Generate a dummy tanimoto distance matrix for demonstration.
np.random.seed(42)
num_molecules = 100
dummy_distance_matrix = np.random.rand(num_molecules, num_molecules)
dummy_distance_matrix = (dummy_distance_matrix + dummy_distance_matrix.T) / 2 #ensure symmetry
np.fill_diagonal(dummy_distance_matrix, 0) #ensure diagonal is zero

# Perform clustering
embedding, labels = umap_hdbscan_clustering_from_distance(dummy_distance_matrix)

# Visualize the results
unique_labels = np.unique(labels)
colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))

for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = [0, 0, 0, 1]

    class_member_mask = (labels == k)

    xy = embedding[class_member_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col), markeredgecolor='k', markersize=5)

plt.title('UMAP projection with HDBSCAN clusters')
plt.show()

# Print cluster labels
print("Cluster labels:", labels)
