"""
utils/fdr_demo.py
-----------------
Simplified FDR / Benjamini-Hochberg teaching utilities.
Also runs a lightweight paired DE analysis on real data for the FDR page.
"""

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


# ── BH simulation (for interactive teaching demo) ────────────────────────────

def simulate_pvalues(n_genes: int = 5000,
                     n_true_de: int = 200,
                     effect_size: float = 2.0,
                     seed: int = 42) -> pd.DataFrame:
    """
    Simulate p-values for n_genes genes.
    n_true_de of them are truly DE (low p-values).
    Returns DataFrame with columns: gene, pvalue, true_de.
    """
    rng = np.random.default_rng(seed)

    # Null genes: p-values uniform under H0
    null_p = rng.uniform(0, 1, n_genes - n_true_de)

    # True DE genes: p-values from a non-central distribution
    true_p = rng.beta(0.5, effect_size, n_true_de)

    pvalues  = np.concatenate([null_p, true_p])
    true_de  = np.array([False] * (n_genes - n_true_de) + [True] * n_true_de)

    genes = [f"Gene_{i+1:05d}" for i in range(n_genes)]
    return pd.DataFrame({"gene": genes, "pvalue": pvalues, "true_de": true_de})


def apply_bh(pvalues: np.ndarray, alpha: float = 0.05):
    """Apply BH correction. Returns (reject, padj)."""
    reject, padj, _, _ = multipletests(pvalues, method="fdr_bh", alpha=alpha)
    return reject, padj


def bh_step_table(pvalues: np.ndarray, alpha: float = 0.05,
                  show_n: int = 20) -> pd.DataFrame:
    """
    Return a small DataFrame showing the BH ranking procedure step-by-step.
    Useful for teaching the BH algorithm.
    """
    n = len(pvalues)
    sorted_idx = np.argsort(pvalues)
    sorted_p   = pvalues[sorted_idx]
    ranks      = np.arange(1, n + 1)
    bh_thresh  = (ranks / n) * alpha
    reject     = sorted_p <= bh_thresh

    # Find the critical rank (largest k where p_k <= (k/n)*alpha)
    critical_rank = 0
    for k in range(n - 1, -1, -1):
        if sorted_p[k] <= bh_thresh[k]:
            critical_rank = k + 1
            break

    df = pd.DataFrame({
        "Rank (k)":         ranks[:show_n],
        "p-value":          sorted_p[:show_n],
        "BH threshold (k/m × α)": bh_thresh[:show_n],
        "Rejected?":        reject[:show_n],
    })
    return df, critical_rank


# ── Real data DE (paired t-test + BH) ────────────────────────────────────────

def run_paired_de(counts_filtered: pd.DataFrame,
                  ctrl_cols: list,
                  treat_cols: list,
                  fdr_cutoff: float = 0.05,
                  lfc_cutoff: float = 1.0) -> pd.DataFrame:
    """
    Paired t-test on log-CPM values + BH correction.
    Educational simplification — not a substitute for DESeq2/edgeR.
    """
    from utils.filtering import log_cpm

    lc = log_cpm(counts_filtered)

    ctrl_mat  = lc[ctrl_cols].values
    treat_mat = lc[treat_cols].values

    diff   = treat_mat - ctrl_mat
    log2fc = diff.mean(axis=1)
    mean_expr = ((ctrl_mat + treat_mat) / 2).mean(axis=1)

    t_stats, p_values = stats.ttest_rel(treat_mat, ctrl_mat, axis=1)
    _, padj, _, _ = multipletests(p_values, method="fdr_bh")

    results = pd.DataFrame({
        "gene":      counts_filtered.index,
        "log2FC":    log2fc,
        "mean_expr": mean_expr,
        "pvalue":    p_values,
        "padj":      padj,
    })
    results["significant"] = (
        (results["padj"] < fdr_cutoff) &
        (results["log2FC"].abs() >= lfc_cutoff)
    )
    results["neg_log10_padj"] = -np.log10(results["padj"].clip(lower=1e-300))
    results["direction"] = "Not significant"
    results.loc[
        results["significant"] & (results["log2FC"] > 0), "direction"
    ] = "Up"
    results.loc[
        results["significant"] & (results["log2FC"] < 0), "direction"
    ] = "Down"

    return results.set_index("gene")
