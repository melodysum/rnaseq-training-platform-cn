"""
utils/data_loader.py
--------------------
处理 RNA-seq 应用的数据加载、验证和演示数据回退。
"""

import os
import numpy as np
import pandas as pd
import streamlit as st

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")


def load_demo_data():
    """加载内置演示数据集（真实计数 + 元数据）。"""
    counts = pd.read_csv(os.path.join(DATA_DIR, "counts.csv"), index_col=0)
    metadata = pd.read_csv(os.path.join(DATA_DIR, "metadata.csv"), index_col=0)
    return counts, metadata


def validate_and_parse(counts_file, metadata_file):
    """
    解析并验证上传的 CSV 文件，进行全面检查。
    返回 (counts_df, metadata_df, error_message)。
    若全部通过，error_message 为 None。

    输入要求：
    - counts.csv：基因为行，样本为列，原始整数计数
    - metadata.csv：样本为行，必须包含 'groupA' 列
    - 输入必须是原始计数，不能是 TPM / FPKM / CPM
    """
    # ── 1. 读取 CSV ──────────────────────────────────────────────────────────
    try:
        counts = pd.read_csv(counts_file, index_col=0)
    except Exception as e:
        return None, None, f"无法读取 counts.csv：{e}"

    try:
        metadata = pd.read_csv(metadata_file, index_col=0)
    except Exception as e:
        return None, None, f"无法读取 metadata.csv：{e}"

    # ── 2. 空值检查 ───────────────────────────────────────────────────────────
    if counts.empty or counts.isna().all().all():
        return None, None, "counts.csv 为空或仅包含缺失值。"
    if metadata.empty or metadata.isna().all().all():
        return None, None, "metadata.csv 为空或仅包含缺失值。"

    # ── 3. 最小维度检查 ───────────────────────────────────────────────────────
    if counts.shape[1] < 2:
        return None, None, (
            f"counts.csv 必须至少有 2 个样本列。"
            f"当前只有 {counts.shape[1]} 列。请检查基因名是否在第一列（用作行索引），"
            f"其余每列是否为一个样本。"
        )
    if counts.shape[0] < 1:
        return None, None, "counts.csv 必须包含至少 1 行基因。"
    if len(metadata) < 2:
        return None, None, f"metadata.csv 必须有至少 2 个样本。当前只有 {len(metadata)} 个。"

    # ── 4. 必需的元数据列检查 ─────────────────────────────────────────────────
    if "groupA" not in metadata.columns:
        return None, None, (
            "metadata.csv 必须包含 'groupA' 列以标识样本分组。"
            "可选但重要的列：'donor'（配对差异表达）、'batch'（批次校正）、"
            "'sex'、'age'（描述性信息）。"
        )

    # ── 5. 重复项检查 ─────────────────────────────────────────────────────────
    dup_count_cols = counts.columns[counts.columns.duplicated()].tolist()
    if dup_count_cols:
        return None, None, (
            f"counts.csv 列中发现重复样本名：{', '.join(dup_count_cols)}。"
            f"每个样本必须有唯一的名称。"
        )
    dup_meta_idx = metadata.index[metadata.index.duplicated()].tolist()
    if dup_meta_idx:
        return None, None, (
            f"metadata.csv 索引中发现重复样本名：{', '.join(str(d) for d in dup_meta_idx)}。"
        )
    dup_genes = counts.index[counts.index.duplicated()].tolist()
    if dup_genes:
        preview = ', '.join(str(g) for g in dup_genes[:5])
        suffix = '...' if len(dup_genes) > 5 else ''
        return None, None, (
            f"counts.csv 索引中发现重复基因标识符：{preview}{suffix}。"
            f"每个基因只能出现一次。"
        )

    # ── 6. 样本名匹配检查 ────────────────────────────────────────────────────
    count_samples    = set(counts.columns)
    meta_samples     = set(metadata.index)
    missing_in_meta  = count_samples - meta_samples
    missing_in_count = meta_samples  - count_samples

    if missing_in_meta:
        return None, None, (
            f"counts.csv 中以下样本名在 metadata.csv 中缺失：{', '.join(sorted(missing_in_meta))}。"
            "样本名必须完全匹配（区分大小写）。"
        )
    if missing_in_count:
        return None, None, (
            f"metadata.csv 中以下样本在 counts.csv 列中未找到：{', '.join(sorted(missing_in_count))}。"
            "元数据中的所有样本必须存在于计数矩阵中。"
        )

    # ── 7. 数值计数检查 ───────────────────────────────────────────────────────
    try:
        counts = counts.apply(pd.to_numeric, errors="raise")
    except Exception:
        return None, None, (
            "counts.csv 包含非数值。"
            "所有计数值必须是整数或浮点数。"
            "请确保第一列包含基因标识符（行索引），其余所有列包含数值计数。"
        )

    # ── 8. 负值检查 ──────────────────────────────────────────────────────────
    if (counts < 0).any().any():
        return None, None, (
            "counts.csv 包含负值。"
            "原始计数必须为非负整数。"
            "本应用需要原始计数——不能是 TPM、FPKM 或 CPM 值，"
            "这些经 log 变换后可能产生负值。"
        )

    # ── 9. 重新排列元数据以匹配计数列顺序 ────────────────────────────────────
    metadata = metadata.reindex(counts.columns)

    return counts, metadata, None


def get_paired_columns(metadata: pd.DataFrame):
    """返回配对设计的 (donors, ctrl_cols, treat_cols)。"""
    if "donor" not in metadata.columns:
        return None, None, None
    donors = metadata["donor"].unique()
    ctrl_cols  = [f"{d}_control"   for d in donors if f"{d}_control"   in metadata.index]
    treat_cols = [f"{d}_treatment" for d in donors if f"{d}_treatment" in metadata.index]
    return donors, ctrl_cols, treat_cols


def init_session_data():
    """
    若尚未上传数据，用演示数据初始化 session_state。
    在每个页面顶部调用一次。
    """
    if "counts" not in st.session_state:
        counts, metadata = load_demo_data()
        st.session_state["counts"]      = counts
        st.session_state["metadata"]    = metadata
        st.session_state["data_source"] = "demo"
