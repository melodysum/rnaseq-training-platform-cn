"""
utils/stats_utils.py
--------------------
Standardised DE output formatting and statistical validation.

Why this module exists:
  Different pages previously produced DE result tables with inconsistent
  column names.  This module enforces a single schema so every downstream
  step (enrichment, volcano, export) can rely on the same column names.

p-value vs FDR — brief explanation
  p-value: probability of observing data at least as extreme as measured,
           assuming the null hypothesis is true for THAT gene.
  FDR (Benjamini-Hochberg adjusted p-value / padj):
           controls the EXPECTED PROPORTION of false discoveries among all
           genes called significant.  When testing 15,000 genes at α=0.05,
           ~750 are expected to be false positives by chance — FDR controls
           this; raw p-value alone does not.
"""

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


# ── Required columns in every DE result table ─────────────────────────────────
REQUIRED_COLS = ["gene_id", "log2FC", "pvalue", "padj", "significant"]

OPTIONAL_COLS = ["baseMean", "lfcSE", "stat",
                 "mean_ref", "mean_target",
                 "direction", "neg_log10_padj",
                 "ref_group", "target_group", "test_type"]


def format_de_results(
    df: pd.DataFrame,
    gene_col: str = None,
    lfc_col: str = "log2FC",
    pval_col: str = "pvalue",
    fdr_cutoff: float = 0.05,
    lfc_cutoff: float = 1.0,
    recalc_padj: bool = False,
) -> pd.DataFrame:
    """
    Standardise a DE result DataFrame to the canonical schema.

    Parameters
    ----------
    df           : raw DE result table (index or column can be gene IDs)
    gene_col     : column containing gene IDs; if None, uses df.index
    lfc_col      : column name for log2 fold-change
    pval_col     : column name for raw p-value
    fdr_cutoff   : FDR threshold for 'significant' flag
    lfc_cutoff   : |log2FC| threshold for 'significant' flag
    recalc_padj  : if True, recompute BH-FDR from pvalue column
                   (use when the table has no padj or padj is suspect)

    Returns
    -------
    DataFrame with canonical columns; gene_id as index.
    """
    out = df.copy()

    # ── gene_id ───────────────────────────────────────────────────────────────
    if gene_col and gene_col in out.columns:
        out = out.set_index(gene_col)
    out.index.name = "gene_id"

    # ── rename lfc / pvalue if needed ─────────────────────────────────────────
    rename_map = {}
    if lfc_col != "log2FC" and lfc_col in out.columns:
        rename_map[lfc_col] = "log2FC"
    if pval_col != "pvalue" and pval_col in out.columns:
        rename_map[pval_col] = "pvalue"
    if rename_map:
        out = out.rename(columns=rename_map)

    # ── padj ──────────────────────────────────────────────────────────────────
    if recalc_padj or "padj" not in out.columns:
        pvals = out["pvalue"].fillna(1.0).values
        _, padj, _, _ = multipletests(pvals, method="fdr_bh")
        out["padj"] = padj
    # Clamp any floating-point edge cases
    out["padj"] = out["padj"].clip(lower=0.0, upper=1.0)

    # ── significant flag ──────────────────────────────────────────────────────
    out["significant"] = (
        (out["padj"] < fdr_cutoff) &
        (out["log2FC"].abs() >= lfc_cutoff)
    )

    # ── direction ─────────────────────────────────────────────────────────────
    if "direction" not in out.columns:
        out["direction"] = "Not significant"
        out.loc[out["significant"] & (out["log2FC"] > 0), "direction"] = "Up"
        out.loc[out["significant"] & (out["log2FC"] < 0), "direction"] = "Down"

    # ── neg_log10_padj ────────────────────────────────────────────────────────
    out["neg_log10_padj"] = -np.log10(out["padj"].clip(lower=1e-300))

    return out


def validate_de_table(df: pd.DataFrame) -> tuple:
    """
    Check that a DE table satisfies minimum requirements.

    Returns
    -------
    (is_valid: bool, issues: list of str)
    """
    issues = []
    for col in REQUIRED_COLS:
        if col not in df.columns and col != "gene_id":
            issues.append(f"Missing required column: {col}")
    if df.get("padj") is not None:
        n_na = df["padj"].isna().sum()
        if n_na > 0:
            issues.append(f"{n_na} NA values in padj column")
        if (df["padj"] > 1).any():
            issues.append("padj values > 1 detected — check correction method")
    return len(issues) == 0, issues


def summarise_de(df: pd.DataFrame,
                 fdr_cutoff: float = 0.05,
                 lfc_cutoff: float = 1.0) -> dict:
    """
    Return a brief summary dict for a formatted DE result table.

    Useful for displaying quick stats at the top of a results section.
    """
    sig = df[
        (df["padj"] < fdr_cutoff) &
        (df["log2FC"].abs() >= lfc_cutoff)
    ]
    return {
        "total_genes":   len(df),
        "n_significant": len(sig),
        "n_up":          int((sig["log2FC"] > 0).sum()),
        "n_down":        int((sig["log2FC"] < 0).sum()),
        "fdr_cutoff":    fdr_cutoff,
        "lfc_cutoff":    lfc_cutoff,
    }
