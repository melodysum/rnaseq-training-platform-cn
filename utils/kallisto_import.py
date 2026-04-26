"""
utils/kallisto_import.py
------------------------
Bridge between Kallisto pseudo-alignment output and count-based analysis.

Why this module is needed
  Kallisto produces *estimated* transcript-level counts (floating-point).
  DESeq2 and edgeR require *integer* gene-level counts.
  The conversion is NOT simply rounding — it requires length-scaling to
  account for differences in effective transcript length across samples.

  Correct workflow (mirrors tximport in R):
    1. Parse abundance.tsv for each sample
    2. Scale est_counts by mean_eff_length / eff_length  (length correction)
    3. Sum scaled counts over transcripts → gene-level counts
    4. Round to integer

  The resulting count matrix is suitable for DESeq2 / edgeR.
  The TPM matrix is for visualisation only — NEVER use it as input to a
  negative binomial or count-based statistical model.

Reference:
  Soneson, Love, Robinson (2015) — "Differential analyses for RNA-seq:
  transcript-level estimates improve gene-level inferences"
"""

import numpy as np
import pandas as pd
from pathlib import Path


# ── Single-sample parsing ─────────────────────────────────────────────────────

def parse_abundance_tsv(filepath: str) -> pd.DataFrame:
    """
    Parse a single Kallisto abundance.tsv file.

    Kallisto output columns:
      target_id   — transcript ID (e.g. ENST00000456328.2)
      length      — transcript length (bp)
      eff_length  — effective length (accounts for fragment length distribution)
      est_counts  — estimated read counts (float)
      tpm         — transcripts per million

    Parameters
    ----------
    filepath : path to abundance.tsv

    Returns
    -------
    DataFrame with columns: target_id, length, eff_length, est_counts, tpm
    """
    df = pd.read_csv(
        filepath,
        sep="\t",
        dtype={
            "target_id":  str,
            "length":     float,
            "eff_length": float,
            "est_counts": float,
            "tpm":        float,
        },
    )
    required = {"target_id", "length", "eff_length", "est_counts", "tpm"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"abundance.tsv at {filepath} is missing columns: {missing}. "
            "Make sure this is a Kallisto output file."
        )
    # Guard: eff_length of 0 would cause division by zero during scaling
    df["eff_length"] = df["eff_length"].clip(lower=1.0)
    return df


# ── Transcript-to-gene mapping ────────────────────────────────────────────────

def _validate_tx2gene(tx2gene: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure tx2gene has columns 'tx_id' and 'gene_id'.
    Accept common alternative column names.
    """
    col_map = {}
    for col in tx2gene.columns:
        cl = col.lower()
        if cl in ("tx_id", "transcript_id", "txid", "target_id", "tx"):
            col_map[col] = "tx_id"
        elif cl in ("gene_id", "geneid", "gene", "gene_name", "ensembl_gene"):
            col_map[col] = "gene_id"
    tx2gene = tx2gene.rename(columns=col_map)
    if "tx_id" not in tx2gene.columns or "gene_id" not in tx2gene.columns:
        raise ValueError(
            "tx2gene must have columns for transcript ID and gene ID. "
            "Recognised names: tx_id/transcript_id/target_id, "
            "gene_id/gene_name/ensembl_gene."
        )
    return tx2gene[["tx_id", "gene_id"]]


# ── Multi-sample import ────────────────────────────────────────────────────────

def import_kallisto_samples(
    sample_dirs: dict,
    tx2gene: pd.DataFrame = None,
    counts_from_abundance: str = "lengthScaledTPM",
    filename: str = "abundance.tsv",
) -> tuple:
    """
    Import multiple Kallisto samples into a gene-level count matrix.

    Parameters
    ----------
    sample_dirs           : dict {sample_id: directory_path}
                            The directory should contain an abundance.tsv
                            (or abundance.h5 — but this function reads TSV).
    tx2gene               : DataFrame mapping tx_id → gene_id.
                            If None, uses transcript IDs as gene IDs directly
                            (only suitable when Kallisto was run with gene-level
                            decoys or a gene-level reference).
    counts_from_abundance : scaling method; currently only "lengthScaledTPM"
                            is implemented (mirrors tximport default).
                            "lengthScaledTPM":
                              scaled_counts = tpm × mean_eff_length / 1e6
                              ≈ est_counts corrected for length differences
                              across samples — preferred for DESeq2.
    filename              : name of the abundance file inside each directory

    Returns
    -------
    counts_matrix : pd.DataFrame, genes × samples, integer dtype
                    Suitable as input to DESeq2 / edgeR modelling.
    tpm_matrix    : pd.DataFrame, genes × samples, float
                    For visualisation ONLY — not for statistical models.

    Notes
    -----
    The count matrix values are rounded integers.  Small rounding errors
    (~1 count) are expected and do not affect downstream analysis.
    """
    if counts_from_abundance != "lengthScaledTPM":
        raise NotImplementedError(
            "Only 'lengthScaledTPM' is currently implemented. "
            "This mirrors tximport's default and is appropriate for DESeq2."
        )

    raw_counts  = {}   # sample → Series(transcript → scaled_count)
    raw_tpm     = {}   # sample → Series(transcript → tpm)
    eff_lengths = {}   # sample → Series(transcript → eff_length)

    for sample_id, dirpath in sample_dirs.items():
        abund_path = Path(dirpath) / filename
        if not abund_path.exists():
            raise FileNotFoundError(
                f"Abundance file not found: {abund_path}\n"
                "Check that Kallisto was run and the path is correct."
            )
        df = parse_abundance_tsv(str(abund_path))
        df = df.set_index("target_id")
        raw_tpm[sample_id]     = df["tpm"]
        eff_lengths[sample_id] = df["eff_length"]
        raw_counts[sample_id]  = df["est_counts"]

    # ── Align transcripts across all samples ─────────────────────────────────
    all_tx = sorted(
        set.intersection(*[set(s.index) for s in raw_counts.values()])
    )
    if len(all_tx) == 0:
        raise ValueError(
            "No transcripts are shared across all samples. "
            "Check that all samples were quantified against the same reference."
        )

    tpm_df    = pd.DataFrame({s: raw_tpm[s].reindex(all_tx)    for s in sample_dirs})
    eff_df    = pd.DataFrame({s: eff_lengths[s].reindex(all_tx) for s in sample_dirs})
    counts_df = pd.DataFrame({s: raw_counts[s].reindex(all_tx) for s in sample_dirs})

    # ── lengthScaledTPM: scale TPM by mean effective length across samples ────
    #   scaled_count = tpm_s × mean_eff_length_across_samples / 1e6 × lib_size
    #   Simplified: use mean eff_length as transcript length normalisation.
    #   This removes between-sample length differences that inflate variance.
    mean_eff_len = eff_df.mean(axis=1)   # mean eff_length per transcript

    # lib_size per sample (sum of TPM-based counts, proportional to actual depth)
    lib_sizes = {}
    for s in sample_dirs:
        # approximate library size from est_counts
        lib_sizes[s] = counts_df[s].sum()
    mean_lib = np.mean(list(lib_sizes.values()))

    scaled = pd.DataFrame(index=all_tx, columns=list(sample_dirs.keys()), dtype=float)
    for s in sample_dirs:
        # length-scaled count: tpm × (mean_eff_len / 1e6) × (mean_lib)
        scaled[s] = tpm_df[s] * (mean_eff_len / 1e6) * mean_lib

    scaled = scaled.fillna(0.0).clip(lower=0.0)

    # ── Aggregate to gene level ────────────────────────────────────────────────
    if tx2gene is not None:
        tx2gene = _validate_tx2gene(tx2gene)
        tx2g    = tx2gene.set_index("tx_id")["gene_id"]
        # Only keep transcripts in our tx2gene map
        in_map  = scaled.index.isin(tx2g.index)
        if in_map.sum() == 0:
            raise ValueError(
                "No transcripts in your abundance files match the tx2gene table. "
                "Check that transcript IDs (e.g. ENST...) match between Kallisto "
                "index and tx2gene."
            )
        scaled_mapped = scaled.loc[in_map].copy()
        scaled_mapped["gene_id"] = tx2g.reindex(scaled_mapped.index).values
        counts_gene = (
            scaled_mapped
            .groupby("gene_id")[list(sample_dirs.keys())]
            .sum()
        )
        # TPM: aggregate directly for visualisation
        tpm_mapped = tpm_df.loc[in_map].copy()
        tpm_mapped["gene_id"] = tx2g.reindex(tpm_mapped.index).values
        tpm_gene = (
            tpm_mapped
            .groupby("gene_id")[list(sample_dirs.keys())]
            .sum()
        )
    else:
        # No tx2gene: treat transcript IDs as gene IDs
        # (only valid for gene-level references)
        counts_gene = scaled
        tpm_gene    = tpm_df

    # ── Round to integer ──────────────────────────────────────────────────────
    # DESeq2 requires integers.  np.round is more stable than int() for floats.
    counts_int = counts_gene.round(0).astype(int)

    return counts_int, tpm_gene.astype(float)


# ── Utility: summarise import ─────────────────────────────────────────────────

def import_summary(counts_matrix: pd.DataFrame, tpm_matrix: pd.DataFrame) -> str:
    """
    Return a human-readable summary of the imported data.
    """
    n_genes   = len(counts_matrix)
    n_samples = len(counts_matrix.columns)
    lib_sizes = counts_matrix.sum(axis=0)
    return (
        f"Imported {n_genes:,} genes × {n_samples} samples.\n"
        f"Library sizes: {lib_sizes.min():,.0f} – {lib_sizes.max():,.0f} reads "
        f"(median {lib_sizes.median():,.0f}).\n"
        f"Zero-count genes: {(counts_matrix.sum(axis=1) == 0).sum():,}.\n"
        "TPM matrix available for visualisation (do NOT use for DE)."
    )
