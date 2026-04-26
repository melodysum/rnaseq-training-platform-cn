"""
utils/simulation.py
-------------------
Simulate RNA-seq count data with batch effects for teaching purposes.
"""

import numpy as np
import pandas as pd


def simulate_batch_data(n_genes: int = 500,
                        n_samples_per_group: int = 5,
                        n_groups: int = 2,
                        n_batches: int = 2,
                        bio_effect: float = 2.0,
                        batch_effect: float = 1.5,
                        confounded: bool = False,
                        seed: int = 42) -> tuple:
    """
    Simulate a log-CPM-like expression matrix with:
      - biological group differences
      - technical batch differences

    Parameters
    ----------
    n_genes              : number of genes to simulate
    n_samples_per_group  : samples per group
    n_groups             : number of biological groups
    n_batches            : number of technical batches
    bio_effect           : standard deviations of biological signal
    batch_effect         : standard deviations of batch shift
    confounded           : if True, each group maps to one batch (danger zone)

    Returns
    -------
    expr     : pd.DataFrame genes × samples (log-CPM scale)
    metadata : pd.DataFrame with sample, group, batch columns
    """
    rng = np.random.default_rng(seed)
    n_samples = n_groups * n_samples_per_group

    # Assign groups
    groups = [f"Group{g+1}" for g in range(n_groups)
              for _ in range(n_samples_per_group)]

    # Assign batches
    if confounded:
        # Each group → one batch (worst case)
        batches = [f"Batch{(g % n_batches) + 1}" for g in range(n_groups)
                   for _ in range(n_samples_per_group)]
    else:
        # Balanced: rotate batches within each group
        batches = []
        for _ in range(n_groups):
            for s in range(n_samples_per_group):
                batches.append(f"Batch{(s % n_batches) + 1}")

    sample_names = [f"S{i+1:02d}" for i in range(n_samples)]
    metadata = pd.DataFrame({
        "sample": sample_names,
        "groupA": groups,
        "batch":  batches,
    }).set_index("sample")

    # Base expression for each gene
    base = rng.normal(5, 1.5, n_genes)

    # Biological effect per group (some genes differ across groups)
    n_de = n_genes // 5
    bio_matrix = np.zeros((n_genes, n_samples))
    group_effects = {g: rng.normal(0, bio_effect, n_de) for g in set(groups)}
    for s_idx, grp in enumerate(groups):
        bio_matrix[:n_de, s_idx] = group_effects[grp]

    # Batch effect per batch (affects all genes)
    batch_offsets = {b: rng.normal(0, batch_effect, n_genes) for b in set(batches)}
    batch_matrix = np.zeros((n_genes, n_samples))
    for s_idx, bat in enumerate(batches):
        batch_matrix[:, s_idx] = batch_offsets[bat]

    # Noise
    noise = rng.normal(0, 0.5, (n_genes, n_samples))

    expr_arr = base[:, None] + bio_matrix + batch_matrix + noise
    gene_names = [f"Gene_{i+1:04d}" for i in range(n_genes)]
    expr = pd.DataFrame(expr_arr, index=gene_names, columns=sample_names)
    return expr, metadata
