"""
utils/enrichment_utils.py
--------------------------
Educational enrichment analysis utilities for the RNA-seq Training Platform.

Provides:
  - TOY_GENE_SETS  : built-in educational pathway library (real human gene symbols)
  - run_ora()      : over-representation analysis (Fisher's exact + BH correction)
  - run_gsea_like(): running-sum enrichment scoring (KS/GSEA-style)
  - demo_ranked_genes(): generate a ranked gene list from built-in demo data

NOT a substitute for clusterProfiler, fgsea, GSEA, or MSigDB-based analysis.
All pathways are educational toy sets for teaching purposes only.
"""

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

# ═══════════════════════════════════════════════════════════════════════════════
# TOY GENE SET LIBRARY
# Real human gene symbols; loosely modelled on GO / MSigDB-style categories.
# Each set contains 15–25 genes chosen to be biologically recognisable.
# ═══════════════════════════════════════════════════════════════════════════════

TOY_GENE_SETS = {

    # ── 1. Interferon response ────────────────────────────────────────────────
    # Anchored to canonical ISG and IRF/STAT signalling genes
    "INTERFERON_RESPONSE": [
        "IFNA1", "IFNB1", "IFNG", "IRF1", "IRF3", "IRF7", "IRF9",
        "STAT1", "STAT2", "MX1", "MX2", "OAS1", "OAS2", "ISG15", "ISG20",
        "IFIT1", "IFIT2", "IFIT3", "RSAD2", "CXCL10", "DDX58", "IFITM1",
        "IFITM3", "BST2", "HERC5",
    ],

    # ── 2. Inflammation / TNF / NF-kB ─────────────────────────────────────────
    # Core inflammatory cytokines and NF-kB pathway members
    "INFLAMMATION_NFKB": [
        "TNF", "IL1A", "IL1B", "IL6", "CXCL8", "CXCL1", "CXCL2",
        "CCL2", "CCL5", "NFKB1", "NFKB2", "RELA", "RELB", "IKBKB",
        "PTGS2", "MMP9", "VCAM1", "ICAM1", "TNFAIP3", "CXCL3",
        "IL18", "NLRP3", "CASP1",
    ],

    # ── 3. Cell cycle / mitosis ────────────────────────────────────────────────
    # Anchored to known mitotic and S-phase regulators
    "CELL_CYCLE_MITOSIS": [
        "MKI67", "PCNA", "CCNA2", "CCNB1", "CCNB2", "CDC20", "CDC25C",
        "CDK1", "AURKA", "AURKB", "BUB1", "BUB1B", "PLK1", "TOP2A",
        "UBE2C", "CENPE", "CENPF", "KIF11", "KIF20A", "NUSAP1",
        "CDKN2A", "E2F1", "MCM2",
    ],

    # ── 4. Mitochondrial / oxidative phosphorylation ──────────────────────────
    # Complex I–V subunits and mitochondrial biogenesis factors
    "OXPHOS_MITOCHONDRIA": [
        "NDUFA1", "NDUFA2", "NDUFS1", "NDUFB1", "SDHA", "SDHB",
        "UQCRC1", "UQCRB", "COX4I1", "COX5A", "ATP5A1", "ATP5B",
        "CYCS", "TFAM", "VDAC1", "TOMM20", "TIMM23", "SLC25A4",
        "IDH2", "FH", "MDH2", "SUCLA2",
    ],

    # ── 5. Ribosome / translation ──────────────────────────────────────────────
    # Ribosomal proteins and translation initiation factors
    "RIBOSOME_TRANSLATION": [
        "RPS3", "RPS4X", "RPS6", "RPS7", "RPS11", "RPS14",
        "RPL3", "RPL4", "RPL6", "RPL10", "RPL13", "RPL18",
        "RPL23", "RPL26", "EIF4E", "EIF4G1", "EIF2S1", "PABPC1",
        "EIF3A", "EIF1AX", "RPLP0",
    ],

    # ── 6. Epithelial markers ──────────────────────────────────────────────────
    # Tight junction, keratin, and epithelial identity genes
    "EPITHELIAL_MARKERS": [
        "EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "KRT7",
        "OCLN", "CLDN3", "CLDN4", "TJP1", "MUC1", "ESRP1",
        "ESRP2", "RAB25", "GRHL2", "OVOL2", "ST14", "DSP",
        "PKP3", "CLDN7",
    ],

    # ── 7. Apoptosis / stress response ────────────────────────────────────────
    # Intrinsic/extrinsic apoptosis and DNA damage response
    "APOPTOSIS_STRESS": [
        "TP53", "BAX", "BAD", "BCL2", "BCL2L1", "CASP3", "CASP7",
        "CASP8", "CASP9", "PARP1", "CYCS", "APAF1", "BID",
        "TNFRSF10A", "TNFRSF10B", "GADD45A", "GADD45B", "CDKN1A",
        "MDM2", "DDIT3", "ATM", "CHEK2",
    ],

    # ── 8. Antigen presentation / immune activation ───────────────────────────
    # MHC class I & II, antigen processing, and immune co-stimulation
    "ANTIGEN_PRESENTATION": [
        "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1",
        "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1", "B2M",
        "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "CD74",
        "CIITA", "NLRC5", "CD80", "CD86",
    ],

    # ── 9. Hypoxia response ────────────────────────────────────────────────────
    # HIF targets and glycolytic shift genes
    "HYPOXIA_RESPONSE": [
        "HIF1A", "EPAS1", "VEGFA", "VEGFC", "SLC2A1", "LDHA",
        "PGAM1", "ENO1", "ALDOA", "PKM", "TPI1", "BNIP3",
        "BNIP3L", "HMOX1", "CA9", "EGLN1", "EGLN3", "PDK1",
        "PFKL", "PGK1",
    ],

    # ── 10. Lipid / metabolic reprogramming ───────────────────────────────────
    # Fatty acid synthesis, cholesterol, and nuclear receptor targets
    "LIPID_METABOLISM": [
        "FASN", "ACACA", "ACACB", "SCD", "ELOVL1", "ELOVL6",
        "DGAT1", "DGAT2", "PPARA", "PPARG", "PPARD", "SREBF1",
        "SREBF2", "INSIG1", "INSIG2", "HMGCR", "LDLR", "ACLY",
        "ACSL4", "FADS1",
    ],
}


# ═══════════════════════════════════════════════════════════════════════════════
# DEMO DATA HELPER
# ═══════════════════════════════════════════════════════════════════════════════

def demo_ranked_genes(counts: pd.DataFrame, metadata: pd.DataFrame,
                      group_col: str = "groupA") -> pd.DataFrame:
    """
    Generate a ranked gene list from built-in demo data using log2FC
    between the first two groups found in metadata[group_col].

    Returns a DataFrame with columns: gene, log2FC, mean_g1, mean_g2
    sorted by log2FC descending (for GSEA-like input).
    """
    groups = metadata[group_col].dropna().unique()
    if len(groups) < 2:
        raise ValueError("Need at least 2 groups in metadata for demo ranking.")

    g1, g2 = sorted(groups)[:2]
    s1 = metadata[metadata[group_col] == g1].index.tolist()
    s2 = metadata[metadata[group_col] == g2].index.tolist()
    s1 = [s for s in s1 if s in counts.columns]
    s2 = [s for s in s2 if s in counts.columns]

    lib = counts.sum(axis=0)
    cpm = counts.divide(lib, axis=1) * 1e6
    lcpm = np.log2(cpm + 1)

    m1 = lcpm[s1].mean(axis=1)
    m2 = lcpm[s2].mean(axis=1)
    lfc = m2 - m1

    df = pd.DataFrame({
        "gene":    counts.index,
        "log2FC":  lfc.values,
        "mean_g1": m1.values,
        "mean_g2": m2.values,
    })
    df["group1"] = g1
    df["group2"] = g2
    return df.sort_values("log2FC", ascending=False).reset_index(drop=True)


# ═══════════════════════════════════════════════════════════════════════════════
# OVER-REPRESENTATION ANALYSIS (ORA)
# ═══════════════════════════════════════════════════════════════════════════════

def run_ora(gene_list: list,
            gene_sets: dict = None,
            universe: list = None) -> pd.DataFrame:
    """
    Over-representation analysis using Fisher's exact test + BH correction.

    Parameters
    ----------
    gene_list : list of gene symbols (query set)
    gene_sets : dict {pathway_name: [gene_symbols]}; defaults to TOY_GENE_SETS
    universe  : background gene list; defaults to all unique genes in gene_sets

    Returns
    -------
    DataFrame with columns:
        pathway, pathway_size, query_size, overlap, gene_ratio,
        pvalue, padj, overlap_genes
    sorted by padj ascending.
    """
    if gene_sets is None:
        gene_sets = TOY_GENE_SETS

    # Build universe from all pathway genes if not supplied
    if universe is None or len(universe) == 0:
        universe = list({g for gs in gene_sets.values() for g in gs})

    universe_set = set(universe)
    query_set    = set(g for g in gene_list if g in universe_set)
    N = len(universe_set)   # total background genes
    K = len(query_set)      # query genes in background

    rows = []
    for pathway, members in gene_sets.items():
        pathway_set   = set(members) & universe_set
        M             = len(pathway_set)           # pathway size in background
        overlap_genes = sorted(query_set & pathway_set)
        k             = len(overlap_genes)         # overlap

        # Fisher's exact test (one-sided: enrichment)
        # Contingency table:
        #             in_pathway | not_in_pathway
        # in_query:      k       |   K - k
        # not_query:   M - k     |   N - M - (K - k)
        table = [
            [k,         K - k],
            [M - k,     N - M - (K - k)],
        ]
        _, pval = stats.fisher_exact(table, alternative="greater")

        rows.append({
            "pathway":      pathway,
            "pathway_size": M,
            "query_size":   K,
            "overlap":      k,
            "gene_ratio":   round(k / M, 4) if M > 0 else 0,
            "pvalue":       pval,
            "overlap_genes": ", ".join(overlap_genes),
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    _, padj, _, _ = multipletests(df["pvalue"].values, method="fdr_bh")
    df["padj"] = padj
    df = df.sort_values("padj").reset_index(drop=True)
    return df


# ═══════════════════════════════════════════════════════════════════════════════
# GSEA-LIKE RUNNING SUM ENRICHMENT
# ═══════════════════════════════════════════════════════════════════════════════

def run_gsea_like(ranked_genes: pd.Series,
                  gene_sets: dict = None,
                  lfc_weights: pd.Series = None) -> pd.DataFrame:
    """
    Running-sum enrichment scoring, conceptually similar to KS/GSEA.

    Parameters
    ----------
    ranked_genes : pd.Series of gene names in ranked order (best to worst).
                   Index should be integer rank position.
    gene_sets    : dict {pathway_name: [gene_symbols]}; defaults to TOY_GENE_SETS
    lfc_weights  : pd.Series {gene_name: abs(log2FC)} for weighted stepping.
                   If None, uses unweighted ±1/N stepping.

    Returns
    -------
    DataFrame with columns:
        pathway, ES, direction, n_hits, pathway_size
    sorted by abs(ES) descending.
    """
    if gene_sets is None:
        gene_sets = TOY_GENE_SETS

    gene_list = list(ranked_genes)
    N         = len(gene_list)
    results   = []

    for pathway, members in gene_sets.items():
        pathway_set = set(members)
        hits        = [g in pathway_set for g in gene_list]
        n_hits      = sum(hits)

        if n_hits == 0:
            results.append({
                "pathway":      pathway,
                "ES":           0.0,
                "direction":    "—",
                "n_hits":       0,
                "pathway_size": len(members),
            })
            continue

        n_miss = N - n_hits

        # Weight for hit positions: abs(log2FC) / sum_abs_lfc_hits
        if lfc_weights is not None:
            lfc_in_path = [abs(lfc_weights.get(g, 0.0)) for g, h in zip(gene_list, hits) if h]
            total_weight = sum(lfc_in_path) or 1.0
            hit_step  = [abs(lfc_weights.get(g, 0.0)) / total_weight if h else 0
                         for g, h in zip(gene_list, hits)]
        else:
            hit_step = [1.0 / n_hits if h else 0 for h in hits]

        miss_step = 1.0 / n_miss if n_miss > 0 else 0.0

        # Walk down the ranked list accumulating the running sum
        running_sum = []
        es = 0.0
        for i, is_hit in enumerate(hits):
            if is_hit:
                es += hit_step[i]
            else:
                es -= miss_step
            running_sum.append(es)

        # Enrichment score = maximum deviation from zero
        max_pos = max(running_sum)
        max_neg = min(running_sum)
        if abs(max_pos) >= abs(max_neg):
            enrichment_score = max_pos
            direction = "Positive (up)"
        else:
            enrichment_score = max_neg
            direction = "Negative (down)"

        results.append({
            "pathway":      pathway,
            "ES":           round(enrichment_score, 4),
            "direction":    direction,
            "n_hits":       n_hits,
            "pathway_size": len(members),
            "_running_sum": running_sum,   # kept for plotting, dropped in return
        })

    df = pd.DataFrame(results)
    if df.empty:
        return df

    df = df.sort_values("ES", key=abs, ascending=False).reset_index(drop=True)
    return df


def get_running_sum(ranked_genes: pd.Series,
                    pathway_genes: list,
                    lfc_weights: pd.Series = None) -> list:
    """
    Return the running-sum vector for a single pathway (for enrichment plots).

    Parameters
    ----------
    ranked_genes   : ordered list/Series of gene names
    pathway_genes  : genes belonging to this pathway
    lfc_weights    : abs(log2FC) per gene for weighted stepping

    Returns
    -------
    list of float — running enrichment score at each rank position
    """
    gene_list   = list(ranked_genes)
    N           = len(gene_list)
    pathway_set = set(pathway_genes)
    hits        = [g in pathway_set for g in gene_list]
    n_hits      = sum(hits)

    if n_hits == 0:
        return [0.0] * N

    n_miss = N - n_hits

    if lfc_weights is not None:
        lfc_in_path  = [abs(lfc_weights.get(g, 0.0)) for g, h in zip(gene_list, hits) if h]
        total_weight = sum(lfc_in_path) or 1.0
        hit_step     = [abs(lfc_weights.get(g, 0.0)) / total_weight if h else 0
                        for g, h in zip(gene_list, hits)]
    else:
        hit_step = [1.0 / n_hits if h else 0 for h in hits]

    miss_step   = 1.0 / n_miss if n_miss > 0 else 0.0
    running_sum = []
    es = 0.0
    for i, is_hit in enumerate(hits):
        es += hit_step[i] if is_hit else -miss_step
        running_sum.append(es)

    return running_sum


# ══════════════════════════════════════════════════════════════════════════════
# MSigDB GMT FILE PARSER
# ══════════════════════════════════════════════════════════════════════════════

# Gene-set source labels (used in the UI)
GENESET_SOURCES = {
    "Toy gene sets — educational only":       "toy",
    "MSigDB Hallmark (50 gene sets)":         "hallmark",
    "MSigDB C2:CP Canonical Pathways (4115 gene sets)": "c2cp",
}

HALLMARK_GMT_PATH = "data/gene_sets/h.all.v2026.1.Hs.symbols.gmt"
C2CP_GMT_PATH     = "data/gene_sets/c2.cp.v2026.1.Hs.symbols.gmt"


def parse_gmt(filepath: str) -> dict:
    """
    Parse an MSigDB GMT file into a gene-set dictionary.

    GMT format (tab-separated):
      pathway_name <TAB> description_or_url <TAB> gene1 <TAB> gene2 ...

    Parameters
    ----------
    filepath : path to .gmt file (relative to app root or absolute)

    Returns
    -------
    dict: {pathway_name: [gene_symbol, ...]}

    Raises
    ------
    FileNotFoundError if the file does not exist.
    """
    import os
    if not os.path.exists(filepath):
        raise FileNotFoundError(
            f"GMT file not found: {filepath}\n"
            "Place the MSigDB GMT file in data/gene_sets/ and restart the app."
        )

    gene_sets = {}
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            name  = parts[0].strip()
            genes = list(dict.fromkeys(           # remove duplicates, preserve order
                g.strip().upper() for g in parts[2:] if g.strip()
            ))
            gene_sets[name] = genes

    return gene_sets


def load_gene_sets(source_key: str) -> tuple:
    """
    Load gene sets based on a source key.

    Parameters
    ----------
    source_key : one of 'toy', 'hallmark', 'c2cp'

    Returns
    -------
    (gene_sets: dict, label: str, n_sets: int)
    """
    if source_key == "toy":
        return TOY_GENE_SETS, "Toy gene sets (educational)", len(TOY_GENE_SETS)

    path_map = {
        "hallmark": (HALLMARK_GMT_PATH, "MSigDB Hallmark"),
        "c2cp":     (C2CP_GMT_PATH,     "MSigDB C2:CP Canonical Pathways"),
    }
    if source_key not in path_map:
        return TOY_GENE_SETS, "Toy gene sets (fallback)", len(TOY_GENE_SETS)

    filepath, label = path_map[source_key]
    gene_sets = parse_gmt(filepath)
    return gene_sets, label, len(gene_sets)


def detect_gene_id_format(gene_index, n_sample: int = 20) -> str:
    """
    Heuristically detect whether gene IDs are HGNC symbols or Ensembl IDs.

    Parameters
    ----------
    gene_index : pandas Index of gene identifiers
    n_sample   : number of genes to sample for detection

    Returns
    -------
    'ensembl'  — most IDs start with ENSG
    'symbol'   — most IDs look like HGNC gene symbols
    'unknown'  — cannot determine
    """
    sample = list(gene_index[:n_sample])
    n_ensg = sum(1 for g in sample if str(g).startswith("ENSG"))
    if n_ensg / len(sample) > 0.6:
        return "ensembl"
    n_sym = sum(1 for g in sample
                if str(g).replace("-", "").replace(".", "").isalnum()
                and not str(g).startswith("ENSG"))
    if n_sym / len(sample) > 0.5:
        return "symbol"
    return "unknown"


# ══════════════════════════════════════════════════════════════════════════════
# GSEA WITH PERMUTATION TESTING (added — proper NES + p-values)
# ══════════════════════════════════════════════════════════════════════════════

def rank_by_statistic(
    de_results: pd.DataFrame,
    stat_col: str = "stat",
    lfc_col: str = "log2FC",
    pval_col: str = "pvalue",
    method: str = "stat",
) -> pd.Series:
    """
    Produce a ranked gene Series for GSEA input from a DE result table.

    Why the choice of ranking metric matters:
      - log2FC alone: ignores variance — noisy genes with small samples
        can rank extremely high/low purely by chance.
      - -log10(p): ignores effect direction and magnitude.
      - Wald / t-statistic: integrates effect size AND precision (SE).
        A gene with large logFC but high variance ranks lower than a gene
        with moderate logFC and low variance.  This is the correct metric.

    Parameters
    ----------
    de_results : DE result DataFrame (index = gene IDs)
    stat_col   : column containing the test statistic (t-stat or Wald stat)
    lfc_col    : column for log2FC (fallback if stat not available)
    pval_col   : column for raw p-value (used for signed-log method)
    method     : one of:
                   "stat"       — use test statistic directly (preferred)
                   "signed_log" — -log10(p) × sign(log2FC)
                   "lfc"        — log2FC (simplest, not recommended)

    Returns
    -------
    pd.Series indexed by gene ID, sorted descending (best → worst).
    """
    df = de_results.copy()
    df.index.name = "gene_id"

    if method == "stat" and stat_col in df.columns:
        metric = df[stat_col].fillna(0.0)
    elif method == "signed_log" and pval_col in df.columns and lfc_col in df.columns:
        sign = np.sign(df[lfc_col].fillna(0.0))
        pval = df[pval_col].clip(lower=1e-300).fillna(1.0)
        metric = -np.log10(pval) * sign
    elif lfc_col in df.columns:
        metric = df[lfc_col].fillna(0.0)
    else:
        raise ValueError(
            f"Could not compute ranking metric. "
            f"Available columns: {list(df.columns)}"
        )

    return metric.sort_values(ascending=False)


def run_gsea_permutation(
    ranked_series: pd.Series,
    gene_sets: dict = None,
    n_permutations: int = 1000,
    random_state: int = 42,
    min_set_size: int = 5,
    max_set_size: int = 500,
) -> pd.DataFrame:
    """
    GSEA with permutation-based normalized enrichment scores (NES) and p-values.

    Key differences from run_gsea_like():
      1. Permutation null: gene labels are shuffled n_permutations times to
         build a null distribution of ES for each pathway size.
      2. NES: ES is normalised by the mean ES of same-size permutations, making
         it comparable across pathways of different sizes.
      3. p-value: proportion of permuted ES ≥ observed ES (one-sided).
      4. FDR: BH correction across all pathways.

    Note on interpretation:
      ORA uses ONLY significant DE genes.
      GSEA uses ALL genes, ranked.
      They answer different questions:
        ORA: "Are my significant hits enriched for a pathway?"
        GSEA: "Is there a systematic tendency for genes in this pathway
               to be coordinately up/downregulated?"
      GSEA finds weaker but consistent signals that ORA misses because ORA
      discards sub-threshold genes entirely.

    Parameters
    ----------
    ranked_series  : pd.Series {gene_id: ranking_metric}, sorted descending
                     Typically output of rank_by_statistic().
    gene_sets      : dict {pathway_name: [gene_ids]}; defaults to TOY_GENE_SETS
    n_permutations : number of gene-label permutations for null distribution
    random_state   : for reproducibility
    min_set_size   : skip pathways with fewer than this many genes in ranked list
    max_set_size   : skip pathways with more than this many genes

    Returns
    -------
    DataFrame with columns:
      pathway, n_genes_in_list, pathway_size, ES, NES,
      pvalue, padj (BH-FDR), direction, leading_edge_genes
    Sorted by NES absolute value descending.
    """
    if gene_sets is None:
        gene_sets = TOY_GENE_SETS

    rng       = np.random.default_rng(random_state)
    gene_list = list(ranked_series.index)
    weights   = np.abs(ranked_series.values)
    N         = len(gene_list)

    def _compute_es(gene_list, pathway_set, weights):
        """Running-sum ES with weighted step (Subramanian 2005 method)."""
        hits     = np.array([g in pathway_set for g in gene_list])
        n_hits   = hits.sum()
        n_miss   = N - n_hits
        if n_hits == 0:
            return 0.0, []

        hit_weights = np.where(hits, weights, 0.0)
        total_w     = hit_weights.sum() or 1.0
        step_hit    = hit_weights / total_w
        step_miss   = 1.0 / n_miss if n_miss > 0 else 0.0

        running = np.cumsum(np.where(hits, step_hit, -step_miss))
        max_dev = running[np.argmax(np.abs(running))]
        # Leading edge: genes before the peak
        peak_idx = int(np.argmax(np.abs(running)))
        leading  = [g for i, (g, h) in enumerate(zip(gene_list, hits))
                    if h and i <= peak_idx]
        return float(max_dev), leading

    results = []
    for pathway, members in gene_sets.items():
        pathway_set = set(members)
        n_in_list   = sum(1 for g in gene_list if g in pathway_set)

        if n_in_list < min_set_size or n_in_list > max_set_size:
            continue

        # Observed ES
        obs_es, leading_edge = _compute_es(gene_list, pathway_set, weights)

        # Permutation null: shuffle gene labels (not sample labels)
        null_es = []
        for _ in range(n_permutations):
            shuffled = rng.permutation(gene_list).tolist()
            es_perm, _ = _compute_es(shuffled, pathway_set, weights)
            null_es.append(es_perm)
        null_es = np.array(null_es)

        # NES: normalise by mean of null ES with same sign
        if obs_es >= 0:
            mean_pos = null_es[null_es >= 0].mean() if (null_es >= 0).any() else 1.0
            nes      = obs_es / mean_pos if mean_pos != 0 else obs_es
            pval     = float((null_es >= obs_es).mean())
        else:
            mean_neg = null_es[null_es <= 0].mean() if (null_es <= 0).any() else -1.0
            nes      = obs_es / abs(mean_neg) if mean_neg != 0 else obs_es
            pval     = float((null_es <= obs_es).mean())

        results.append({
            "pathway":           pathway,
            "n_genes_in_list":   n_in_list,
            "pathway_size":      len(members),
            "ES":                round(obs_es, 4),
            "NES":               round(nes, 4),
            "pvalue":            max(pval, 1 / (n_permutations + 1)),  # floor
            "direction":         "Positive (up)" if obs_es >= 0 else "Negative (down)",
            "leading_edge_genes": ", ".join(leading_edge[:10]),
        })

    if not results:
        return pd.DataFrame()

    df = pd.DataFrame(results)
    _, padj, _, _ = multipletests(df["pvalue"].values, method="fdr_bh")
    df["padj"] = padj
    df = df.sort_values("NES", key=abs, ascending=False).reset_index(drop=True)
    return df
