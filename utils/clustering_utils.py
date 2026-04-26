"""
utils/clustering_utils.py
-------------------------
Hierarchical clustering helpers for the Clustering & Heatmaps lesson.
"""

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, leaves_list, cut_tree


def compute_distance(matrix: np.ndarray, metric: str = "euclidean") -> np.ndarray:
    """
    Compute pairwise distances between rows.

    Parameters
    ----------
    matrix : 2D numpy array (rows = items to cluster)
    metric : 'euclidean', 'pearson', or 'spearman'

    Returns
    -------
    Condensed distance vector (for use with scipy linkage).
    """
    if metric == "euclidean":
        return pdist(matrix, metric="euclidean")

    if metric == "pearson":
        # Convert Pearson r to distance: d = 1 - r
        corr = np.corrcoef(matrix)
        corr = np.clip(corr, -1, 1)
        dist_sq = squareform(1 - corr, checks=False)
        return dist_sq

    if metric == "spearman":
        from scipy.stats import rankdata
        ranked = np.apply_along_axis(rankdata, 1, matrix)
        corr   = np.corrcoef(ranked)
        corr   = np.clip(corr, -1, 1)
        dist_sq = squareform(1 - corr, checks=False)
        return dist_sq

    raise ValueError(f"Unknown metric: {metric}")


def cluster_matrix(matrix: np.ndarray,
                   metric: str = "pearson",
                   method: str = "average") -> np.ndarray:
    """
    Run hierarchical clustering on rows of matrix.

    Returns
    -------
    linkage_matrix : (n-1, 4) array from scipy linkage
    """
    dist_condensed = compute_distance(matrix, metric)
    # Ensure no negative values due to floating point
    dist_condensed = np.clip(dist_condensed, 0, None)
    return linkage(dist_condensed, method=method)


def get_leaf_order(linkage_matrix: np.ndarray) -> np.ndarray:
    """Return the leaf order from a linkage matrix."""
    return leaves_list(linkage_matrix)


def row_zscore(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply row-wise Z-score scaling.
    z = (x - mean) / std
    Rows with zero std are set to 0 to avoid divide-by-zero.
    """
    means = df.mean(axis=1)
    stds  = df.std(axis=1)
    # Avoid division by zero
    stds_safe = stds.copy()
    stds_safe[stds_safe < 1e-10] = 1.0
    scaled = df.subtract(means, axis=0).divide(stds_safe, axis=0)
    # Zero out rows that had no variance
    scaled.loc[stds < 1e-10] = 0.0
    return scaled


def assign_modules(linkage_matrix: np.ndarray,
                   n_clusters: int,
                   labels: list) -> pd.Series:
    """
    Cut dendrogram into k clusters and return module assignments.

    Returns
    -------
    pd.Series indexed by label, values are module numbers (1-indexed).
    """
    assignments = cut_tree(linkage_matrix, n_clusters=n_clusters).flatten()
    return pd.Series(assignments + 1, index=labels, name="module")
