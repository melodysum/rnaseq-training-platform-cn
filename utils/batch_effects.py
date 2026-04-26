"""
utils/batch_effects.py
----------------------
Batch effect detection and correction helpers.
Educational simplification — not a substitute for limma removeBatchEffect or ComBat.
"""

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression


def to_log_cpm(counts: pd.DataFrame, pseudocount: float = 1.0) -> pd.DataFrame:
    """Convert raw counts to log2(CPM + pseudocount)."""
    lib_sizes = counts.sum(axis=0)
    cpm = counts.divide(lib_sizes, axis=1) * 1e6
    return np.log2(cpm + pseudocount)


def simple_batch_correction(log_expr: pd.DataFrame,
                             metadata: pd.DataFrame,
                             batch_col: str = "batch",
                             group_col: str = "groupA") -> pd.DataFrame:
    """
    Simple educational batch correction via linear regression residualisation.

    For each gene:
      1. Fit: expression ~ batch_dummies + group_dummies
      2. Residuals = expression - batch_contribution
      3. Add back the grand mean so values stay interpretable

    This is conceptually similar to limma::removeBatchEffect.
    It preserves the group (biological) variable.

    Parameters
    ----------
    log_expr  : genes × samples log-CPM matrix
    metadata  : DataFrame with batch_col and group_col
    batch_col : column name for batch
    group_col : column name for biological group

    Returns
    -------
    Corrected genes × samples DataFrame (same shape as input).
    """
    # Align samples
    samples = [s for s in log_expr.columns if s in metadata.index]
    expr    = log_expr[samples]
    meta    = metadata.loc[samples]

    # Build design matrix: dummies for batch + group
    batch_dummies = pd.get_dummies(meta[batch_col], prefix="batch", drop_first=True)
    group_dummies = pd.get_dummies(meta[group_col], prefix="group", drop_first=True)
    X = pd.concat([batch_dummies, group_dummies], axis=1).astype(float).values

    corrected = np.zeros_like(expr.values, dtype=float)
    grand_mean = expr.values.mean(axis=1, keepdims=True)

    for i, gene_expr in enumerate(expr.values):
        if np.std(gene_expr) < 1e-10:
            corrected[i] = gene_expr
            continue
        model = LinearRegression(fit_intercept=True).fit(X, gene_expr)
        # Reconstruct batch-only contribution
        n_batch = batch_dummies.shape[1]
        batch_coefs = model.coef_[:n_batch]
        batch_contrib = batch_dummies.values @ batch_coefs
        corrected[i] = gene_expr - batch_contrib

    # Re-centre to grand mean
    corrected = corrected - corrected.mean(axis=1, keepdims=True) + grand_mean

    return pd.DataFrame(corrected, index=expr.index, columns=expr.columns)


def variance_explained_by(scores: pd.DataFrame,
                           metadata: pd.DataFrame,
                           variable: str) -> float:
    """
    Estimate how much variance in PC1+PC2 is explained by a metadata variable.
    Uses a simple R² from a one-way regression.
    """
    from sklearn.preprocessing import LabelEncoder
    from sklearn.linear_model import LinearRegression

    common = scores.index.intersection(metadata.index)
    if len(common) < 3:
        return np.nan

    y = scores.loc[common, ["PC1", "PC2"]].values.ravel()
    labels = metadata.loc[common, variable]
    x = LabelEncoder().fit_transform(labels).repeat(2).reshape(-1, 1)

    r2 = LinearRegression().fit(x, y).score(x, y)
    return r2


def top_variable_genes(log_expr: pd.DataFrame, n: int = 500) -> pd.DataFrame:
    """Return the top n most variable genes."""
    var = log_expr.var(axis=1)
    return log_expr.loc[var.nlargest(n).index]
