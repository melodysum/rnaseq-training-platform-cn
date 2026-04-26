"""
utils/de_analysis.py
--------------------
Educational DE workflow: logCPM + linear model/t-test + BH FDR.
Not a substitute for DESeq2 / edgeR / limma-voom.
"""

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def run_de(counts: pd.DataFrame,
           metadata: pd.DataFrame,
           group_col: str = "groupA",
           ref_group: str = None,
           target_group: str = None,
           batch_col: str = None,
           donor_col: str = None,
           fdr_cutoff: float = 0.05,
           lfc_cutoff: float = 1.0) -> pd.DataFrame:
    """
    Educational DE analysis.

    Steps:
      1. Convert counts to log2(CPM + 1)
      2. If donor exists → paired t-test (treatment - control per donor)
         If no donor    → Welch t-test between the two groups
         If batch given → residualise batch before testing
      3. BH FDR correction

    Returns
    -------
    DataFrame indexed by gene with columns:
        mean_ref, mean_target, log2FC, pvalue, padj,
        significant, direction, neg_log10_padj
    """
    # ── log-CPM ──────────────────────────────────────────────────────────────
    lib_sizes = counts.sum(axis=0)
    cpm       = counts.divide(lib_sizes, axis=1) * 1e6
    log_cpm   = np.log2(cpm + 1)

    # ── subset to the two groups ──────────────────────────────────────────────
    groups = metadata[group_col]
    if ref_group is None:
        ref_group    = sorted(groups.unique())[0]
    if target_group is None:
        target_group = sorted(groups.unique())[1]

    ref_samples    = groups[groups == ref_group].index.tolist()
    target_samples = groups[groups == target_group].index.tolist()

    ref_samples    = [s for s in ref_samples    if s in log_cpm.columns]
    target_samples = [s for s in target_samples if s in log_cpm.columns]

    # ── optional: residualise batch ────────────────────────────────────────────
    if batch_col and batch_col in metadata.columns:
        from utils.batch_effects import simple_batch_correction
        log_cpm = simple_batch_correction(log_cpm, metadata,
                                          batch_col=batch_col,
                                          group_col=group_col)

    # ── DE testing ─────────────────────────────────────────────────────────────
    ref_mat    = log_cpm[ref_samples].values     # genes × ref_samples
    target_mat = log_cpm[target_samples].values  # genes × target_samples

    mean_ref    = ref_mat.mean(axis=1)
    mean_target = target_mat.mean(axis=1)
    log2fc      = mean_target - mean_ref

    # ── Paired alignment by donor ──────────────────────────────────────────────
    # Attempt donor-aligned pairing only if donor_col is available.
    # We find donors with exactly one sample in each group, then align explicitly.
    pairing_status = "unpaired"
    if donor_col is not None and donor_col in metadata.columns:
        meta_sub = metadata.loc[
            metadata[group_col].isin([ref_group, target_group])
        ]
        ref_meta    = meta_sub[meta_sub[group_col] == ref_group]
        target_meta = meta_sub[meta_sub[group_col] == target_group]

        # donors present in both groups exactly once
        ref_donors    = ref_meta[donor_col].value_counts()
        target_donors = target_meta[donor_col].value_counts()
        valid_donors  = sorted(
            d for d in ref_donors.index
            if d in target_donors.index
            and ref_donors[d] == 1
            and target_donors[d] == 1
        )

        if len(valid_donors) >= 2:
            # Re-align ref and target in matched donor order
            ref_aligned    = [ref_meta[ref_meta[donor_col] == d].index[0]
                              for d in valid_donors]
            target_aligned = [target_meta[target_meta[donor_col] == d].index[0]
                              for d in valid_donors]
            # Only keep samples present in log_cpm
            ref_aligned    = [s for s in ref_aligned    if s in log_cpm.columns]
            target_aligned = [s for s in target_aligned if s in log_cpm.columns]

            if len(ref_aligned) >= 2 and len(ref_aligned) == len(target_aligned):
                ref_mat    = log_cpm[ref_aligned].values
                target_mat = log_cpm[target_aligned].values
                mean_ref    = ref_mat.mean(axis=1)
                mean_target = target_mat.mean(axis=1)
                log2fc      = mean_target - mean_ref
                pairing_status = f"paired ({len(ref_aligned)} donors)"
            else:
                pairing_status = "fallback_unpaired (alignment failed)"
        else:
            pairing_status = "fallback_unpaired (insufficient matched donors)"

    use_paired = pairing_status.startswith("paired")

    if use_paired:
        _, pvalues = stats.ttest_rel(target_mat, ref_mat, axis=1)
    else:
        _, pvalues = stats.ttest_ind(target_mat, ref_mat,
                                     axis=1, equal_var=False)

    # BH correction
    _, padj, _, _ = multipletests(pvalues, method="fdr_bh")

    results = pd.DataFrame({
        "gene":        counts.index,
        "mean_ref":    mean_ref,
        "mean_target": mean_target,
        "log2FC":      log2fc,
        "pvalue":      pvalues,
        "padj":        padj,
    }).set_index("gene")

    results["significant"] = (
        (results["padj"] < fdr_cutoff) &
        (results["log2FC"].abs() >= lfc_cutoff)
    )
    results["neg_log10_padj"] = -np.log10(results["padj"].clip(lower=1e-300))

    results["direction"] = "Not significant"
    results.loc[results["significant"] & (results["log2FC"] > 0), "direction"] = "Up"
    results.loc[results["significant"] & (results["log2FC"] < 0), "direction"] = "Down"

    results["ref_group"]      = ref_group
    results["target_group"]   = target_group
    results["test_type"]      = "paired t-test" if use_paired else "Welch t-test"
    results["pairing_status"] = pairing_status

    return results
