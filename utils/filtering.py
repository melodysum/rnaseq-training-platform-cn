"""
utils/filtering.py
------------------
Low-expression filtering logic and related statistics.
"""

import numpy as np
import pandas as pd


def filter_low_expression(counts: pd.DataFrame,
                           min_count: int = 10,
                           min_samples: int = 5) -> pd.DataFrame:
    """Keep genes with >= min_count in at least min_samples samples."""
    mask = (counts >= min_count).sum(axis=1) >= min_samples
    return counts.loc[mask]


def expression_summary(counts: pd.DataFrame) -> pd.DataFrame:
    """Per-gene summary: mean, median, max count, fraction of zero samples."""
    return pd.DataFrame({
        "mean_count":    counts.mean(axis=1),
        "median_count":  counts.median(axis=1),
        "max_count":     counts.max(axis=1),
        "pct_zeros":     (counts == 0).mean(axis=1) * 100,
    })


def check_gene_fate(counts_raw: pd.DataFrame,
                    counts_filtered: pd.DataFrame,
                    gene_list: list) -> pd.DataFrame:
    """
    For each gene in gene_list, report whether it was:
      - retained after filtering
      - removed due to low expression
      - not found in the dataset
    """
    rows = []
    all_genes      = set(counts_raw.index)
    retained_genes = set(counts_filtered.index)

    for gene in gene_list:
        gene = gene.strip()
        if gene not in all_genes:
            status = "❓ Not found in dataset"
            mean_count = None
        elif gene in retained_genes:
            status = "✅ Retained after filtering"
            mean_count = counts_raw.loc[gene].mean()
        else:
            status = "🔴 Removed (low expression)"
            mean_count = counts_raw.loc[gene].mean()

        rows.append({"Gene": gene, "Status": status,
                     "Mean count (all samples)": mean_count})

    return pd.DataFrame(rows)


def normalise_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """Convert raw counts to CPM."""
    lib_sizes = counts.sum(axis=0)
    return counts.divide(lib_sizes, axis=1) * 1e6


def log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """log2(CPM + 1)."""
    return np.log2(normalise_cpm(counts) + 1)


# ══════════════════════════════════════════════════════════════════════════════
# DATA-DRIVEN THRESHOLD SELECTION  (added: filterByExpr-equivalent)
# ══════════════════════════════════════════════════════════════════════════════

def filter_by_expr(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    group_col: str = "groupA",
    min_count: float = 10.0,
    min_total_count: float = 15.0,
) -> pd.DataFrame:
    """
    filterByExpr-equivalent in Python (mirrors edgeR::filterByExpr logic).

    Threshold is derived from library sizes and group sizes — NOT from
    differential expression results.  Using DE signal to choose the filter
    that maximises DE signal is circular reasoning and inflates false positives.

    Logic (edgeR paper, Chen et al. 2016):
      • Compute CPM threshold from median library size:
            cpm_cutoff = min_count / (median_library_size / 1e6)
        This means a gene needs ~min_count reads in the SMALLEST library
        to pass — not a fixed count regardless of sequencing depth.
      • Keep a gene if its CPM exceeds cpm_cutoff in at least
        min_group_size samples, where min_group_size = size of the
        smallest experimental group.
      • Additionally require total counts across all samples ≥ min_total_count.

    Parameters
    ----------
    counts          : raw integer count matrix (genes × samples)
    metadata        : sample metadata DataFrame (index = sample IDs)
    group_col       : column in metadata defining experimental groups
    min_count       : minimum count in the median-sized library (default 10)
    min_total_count : minimum total counts across all samples (default 15)

    Returns
    -------
    Filtered count matrix.
    """
    # Align samples
    shared = [s for s in counts.columns if s in metadata.index]
    counts_  = counts[shared]
    meta_    = metadata.loc[shared]

    # Library sizes and CPM threshold
    lib_sizes   = counts_.sum(axis=0)
    median_lib  = lib_sizes.median()
    cpm_cutoff  = min_count / (median_lib / 1e6)

    # CPM matrix
    cpm_matrix = counts_.divide(lib_sizes, axis=1) * 1e6

    # Smallest group size → minimum samples that must pass
    if group_col in meta_.columns:
        group_sizes  = meta_[group_col].value_counts()
        min_grp_size = int(group_sizes.min())
    else:
        min_grp_size = max(1, counts_.shape[1] // 4)

    # Gene-level filter
    n_passing   = (cpm_matrix > cpm_cutoff).sum(axis=1)
    cpm_pass    = n_passing >= min_grp_size
    total_pass  = counts_.sum(axis=1) >= min_total_count
    keep        = cpm_pass & total_pass

    return counts.loc[keep]


def threshold_sweep_retained(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    group_col: str = "groupA",
    thresholds: list = None,
) -> pd.DataFrame:
    """
    Sweep count thresholds and report genes RETAINED (not DE genes).

    This is the correct diagnostic for threshold sensitivity.
    Plotting retained-gene count vs threshold shows how conservative each
    choice is WITHOUT introducing circularity from DE testing.

    Parameters
    ----------
    counts     : raw count matrix
    metadata   : sample metadata
    group_col  : group column for filterByExpr min_group_size logic
    thresholds : list of min_count values to test; defaults to [1,5,10,20,50]

    Returns
    -------
    DataFrame with columns: threshold, n_retained, pct_retained,
                             n_removed, cpm_cutoff, min_group_size
    """
    if thresholds is None:
        thresholds = [1, 5, 10, 20, 50]

    shared = [s for s in counts.columns if s in metadata.index]
    counts_ = counts[shared]
    meta_   = metadata.loc[shared]

    lib_sizes  = counts_.sum(axis=0)
    median_lib = lib_sizes.median()
    total_genes = len(counts_)

    if group_col in meta_.columns:
        min_grp_size = int(meta_[group_col].value_counts().min())
    else:
        min_grp_size = max(1, counts_.shape[1] // 4)

    cpm_matrix = counts_.divide(lib_sizes, axis=1) * 1e6
    rows = []
    for t in thresholds:
        cpm_cut    = t / (median_lib / 1e6)
        n_pass_cpm = (cpm_matrix > cpm_cut).sum(axis=1)
        cpm_pass   = n_pass_cpm >= min_grp_size
        total_pass = counts_.sum(axis=1) >= max(t + 5, 15)
        n_kept     = int((cpm_pass & total_pass).sum())
        rows.append({
            "threshold":     t,
            "n_retained":    n_kept,
            "pct_retained":  round(n_kept / total_genes * 100, 1),
            "n_removed":     total_genes - n_kept,
            "cpm_cutoff":    round(cpm_cut, 3),
            "min_group_size": min_grp_size,
        })
    return pd.DataFrame(rows)
