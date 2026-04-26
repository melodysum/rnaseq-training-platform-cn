"""
utils/pca_utils.py
------------------
PCA helper functions for the Batch Correction lesson.
"""

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def run_pca(expr_matrix: pd.DataFrame,
            n_components: int = 10,
            scale: bool = True,
            top_var_genes: int = None) -> tuple:
    """
    Run PCA on an expression matrix (genes × samples).

    Parameters
    ----------
    expr_matrix   : genes × samples DataFrame
    n_components  : number of PCs to compute
    scale         : whether to scale features (samples) before PCA
    top_var_genes : if set, subset to top N most variable genes first

    Returns
    -------
    scores     : pd.DataFrame, shape (samples, n_components)
    explained  : np.array of variance explained ratios
    loadings   : pd.DataFrame, shape (genes, n_components)
    """
    mat = expr_matrix.copy()

    if top_var_genes and top_var_genes < len(mat):
        gene_var = mat.var(axis=1)
        top_idx  = gene_var.nlargest(top_var_genes).index
        mat      = mat.loc[top_idx]

    # PCA is performed on samples, so transpose: rows = samples
    X = mat.T.values
    n_comp = min(n_components, X.shape[0], X.shape[1])

    if scale:
        X = StandardScaler().fit_transform(X)

    pca = PCA(n_components=n_comp)
    scores_arr = pca.fit_transform(X)

    pc_names = [f"PC{i+1}" for i in range(n_comp)]
    scores   = pd.DataFrame(scores_arr, index=mat.columns, columns=pc_names)
    loadings = pd.DataFrame(
        pca.components_.T,
        index=mat.index,
        columns=pc_names,
    )
    return scores, pca.explained_variance_ratio_, loadings


def pca_plot_df(scores: pd.DataFrame,
                metadata: pd.DataFrame,
                pc_x: str = "PC1",
                pc_y: str = "PC2") -> pd.DataFrame:
    """
    Merge PCA scores with metadata for plotting.
    Returns a DataFrame ready for px.scatter.
    """
    df = scores[[pc_x, pc_y]].copy()
    df.index.name = "sample"
    df = df.reset_index()
    if metadata is not None:
        meta = metadata.reset_index().rename(columns={metadata.index.name or "index": "sample"})
        df = df.merge(meta, on="sample", how="left")
    return df


# ══════════════════════════════════════════════════════════════════════════════
# VARIANCE DECOMPOSITION  (PVCA-like analysis)
# ══════════════════════════════════════════════════════════════════════════════

def variance_decomposition(
    scores: pd.DataFrame,
    explained: np.ndarray,
    metadata: pd.DataFrame,
    factors: list = None,
    max_pcs: int = None,
    cumvar_threshold: float = 0.80,
    random_state: int = 42,
) -> pd.DataFrame:
    """
    Quantify the contribution of experimental factors to total expression variance.

    Method (PVCA-like, Bushel et al. 2009):
      1. Select PCs that together explain ≥ cumvar_threshold of variance.
      2. For each selected PC, regress its scores against the provided factors
         using ordinary least squares.
      3. Compute partial R² (partial eta²) for each factor within each PC.
      4. Weight each PC's partial R² by its proportion of variance explained.
      5. Aggregate weighted R² across PCs → per-factor % of total variance.

    Why this matters:
      A simple PCA coloured by batch tells you batch IS a problem if samples
      cluster by batch — but not HOW MUCH variance batch explains relative to
      biology.  This decomposition gives a quantitative answer, e.g.:
        "batch explains 34% of variance; condition explains 41%"
      so you can decide whether batch correction is necessary.

    Parameters
    ----------
    scores            : PC scores (samples × PCs), output of run_pca()
    explained         : array of variance explained ratios, output of run_pca()
    metadata          : sample metadata DataFrame (index = sample IDs)
    factors           : list of metadata column names to test;
                        defaults to all non-numeric columns with < 20 unique values
    max_pcs           : cap on number of PCs to include (ignores cumvar_threshold)
    cumvar_threshold  : include PCs until this cumulative variance is reached
    random_state      : unused; for API consistency

    Returns
    -------
    DataFrame with columns: factor, pct_variance_explained
    Sorted descending by pct_variance_explained.
    Includes a 'residual' row for unexplained variance.
    """
    try:
        import statsmodels.api as sm
    except ImportError:
        raise ImportError(
            "statsmodels is required for variance_decomposition. "
            "It is already in requirements.txt."
        )

    # ── align samples ─────────────────────────────────────────────────────────
    shared = [s for s in scores.index if s in metadata.index]
    if len(shared) < 4:
        raise ValueError(
            f"Only {len(shared)} samples shared between PCA scores and metadata. "
            "Need at least 4."
        )
    scores_  = scores.loc[shared]
    meta_    = metadata.loc[shared]

    # ── select factors ────────────────────────────────────────────────────────
    if factors is None:
        factors = [
            c for c in meta_.columns
            if meta_[c].nunique() < 20 and meta_[c].nunique() > 1
        ]
    factors = [f for f in factors if f in meta_.columns]
    if not factors:
        raise ValueError("No valid factors found in metadata for decomposition.")

    # ── select PCs ────────────────────────────────────────────────────────────
    cumvar = np.cumsum(explained)
    if max_pcs:
        n_pcs = min(max_pcs, len(explained))
    else:
        # PCs until cumvar_threshold is reached
        n_pcs = int(np.searchsorted(cumvar, cumvar_threshold)) + 1
        n_pcs = min(n_pcs, len(explained), scores_.shape[1])

    pc_cols    = scores_.columns[:n_pcs].tolist()
    pc_weights = explained[:n_pcs] / explained[:n_pcs].sum()  # re-normalise

    # ── encode categorical factors as dummy variables ─────────────────────────
    factor_matrices = {}
    for f in factors:
        col = meta_[f]
        if col.dtype.kind in ("O", "b") or col.nunique() <= 10:
            # categorical → one-hot (drop_first to avoid perfect collinearity)
            dummies = pd.get_dummies(col, prefix=f, drop_first=True,
                                     dtype=float)
            factor_matrices[f] = dummies.values
        else:
            # continuous → standardise
            v = col.astype(float).values
            std = v.std()
            factor_matrices[f] = ((v - v.mean()) / std if std > 0 else v).reshape(-1, 1)

    # ── per-PC partial R² ─────────────────────────────────────────────────────
    weighted_r2 = {f: 0.0 for f in factors}

    for pc, weight in zip(pc_cols, pc_weights):
        y = scores_[pc].values.astype(float)

        # Full model R²
        X_all_parts = [factor_matrices[f] for f in factors]
        X_full = np.hstack(X_all_parts)
        X_full = sm.add_constant(X_full, has_constant="add")
        try:
            r2_full = sm.OLS(y, X_full).fit().rsquared
        except Exception:
            r2_full = 0.0

        for f in factors:
            # Reduced model: all factors EXCEPT f
            X_red_parts = [factor_matrices[g] for g in factors if g != f]
            if X_red_parts:
                X_red = np.hstack(X_red_parts)
                X_red = sm.add_constant(X_red, has_constant="add")
                try:
                    r2_red = sm.OLS(y, X_red).fit().rsquared
                except Exception:
                    r2_red = 0.0
            else:
                r2_red = 0.0

            # Partial R²: unique variance explained by f beyond other factors
            partial = max(0.0, r2_full - r2_red)
            weighted_r2[f] += weight * partial

    # ── residual ──────────────────────────────────────────────────────────────
    total_explained = sum(weighted_r2.values())
    residual = max(0.0, 1.0 - total_explained)

    rows = [
        {"factor": f, "pct_variance_explained": round(v * 100, 2)}
        for f, v in weighted_r2.items()
    ]
    rows.append({"factor": "Residual / other", "pct_variance_explained": round(residual * 100, 2)})
    result = pd.DataFrame(rows).sort_values(
        "pct_variance_explained", ascending=False
    ).reset_index(drop=True)
    return result
