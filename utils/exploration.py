"""
utils/exploration.py
--------------------
Sample-level exploration helpers: distances, correlations, outlier scoring.
"""

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import zscore


def sample_correlation_matrix(log_expr: pd.DataFrame) -> pd.DataFrame:
    """Pearson correlation between all sample pairs."""
    return log_expr.corr(method="pearson")


def sample_distance_matrix(log_expr: pd.DataFrame,
                            metric: str = "euclidean") -> pd.DataFrame:
    """Pairwise distance matrix between samples."""
    dist_arr = squareform(pdist(log_expr.T.values, metric=metric))
    return pd.DataFrame(dist_arr,
                        index=log_expr.columns,
                        columns=log_expr.columns)


def outlier_scores(scores_df: pd.DataFrame,
                   pcs: list = None) -> pd.Series:
    """
    Simple outlier score: z-score of Euclidean distance from centroid in PC space.
    Higher score = more unusual sample.
    """
    if pcs is None:
        pcs = [c for c in scores_df.columns if c.startswith("PC")]
    mat = scores_df[pcs].values
    centroid = mat.mean(axis=0)
    dists = np.linalg.norm(mat - centroid, axis=1)
    z = zscore(dists)
    return pd.Series(np.abs(z), index=scores_df.index, name="outlier_score")
