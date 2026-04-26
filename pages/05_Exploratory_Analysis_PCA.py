"""
pages/05_Exploratory_Analysis_PCA.py — 第5课：探索性分析与 PCA
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from utils.data_loader import init_session_data
from utils.filtering import filter_low_expression, log_cpm
from utils.batch_effects import top_variable_genes
from utils.pca_utils import run_pca, pca_plot_df
from utils.exploration import (
    sample_correlation_matrix,
    sample_distance_matrix,
    outlier_scores,
)
from utils.batch_effects import top_variable_genes as tvg

st.set_page_config(
    page_title="第5课 — 探索性分析",
    page_icon="📊",
    layout="wide",
)

init_session_data()

st.title("📊 第5课 — 探索性分析与 PCA")
st.markdown("""
> **学习目标：** 学会如何检查 RNA-seq 样本结构、检测离群样本和批次效应，
> 并在运行差异表达分析之前正确解读 PCA 结果。
""")

counts_raw = st.session_state["counts"]
metadata   = st.session_state["metadata"]

filter_res = st.session_state.get("filter_results", {})

with st.sidebar:
    st.header("🔧 分析参数")

    use_filtered = st.radio(
        "基因集",
        ["使用过滤后的基因（来自第3课）", "使用所有基因"],
        index=0 if filter_res else 1,
    )
    if use_filtered == "使用过滤后的基因（来自第3课）" and filter_res:
        counts_use = filter_res["counts_filtered"]
        st.caption(f"使用 {len(counts_use):,} 个过滤后的基因。")
    else:
        counts_use = filter_low_expression(counts_raw, 10, 5)
        st.caption(f"使用 {len(counts_use):,} 个基因（默认过滤）。")

    st.divider()
    top_n_genes = st.slider("用于 PCA 的高变异基因数", 100,
                             min(3000, len(counts_use)), 500, 100)
    pc_choices  = [f"PC{i}" for i in range(1, 6)]
    pc_x = st.selectbox("PCA X 轴", pc_choices, index=0)
    pc_y = st.selectbox("PCA Y 轴", pc_choices, index=1)

    meta_cols = [c for c in metadata.columns
                 if metadata[c].nunique() <= 20 or c == "age"]
    color_by = st.selectbox("样本着色依据", meta_cols)
    shape_by = st.selectbox("样本形状依据（可选）",
                             ["无"] + [c for c in meta_cols if c != color_by])
    show_labels = st.checkbox("显示样本标签", value=False)

log_expr = log_cpm(counts_use)
scores, explained, loadings = run_pca(log_expr, n_components=5,
                                       top_var_genes=top_n_genes)
pca_df = pca_plot_df(scores, metadata, pc_x, pc_y)

# ── 第一节 — 为什么探索性分析很重要 ─────────────────────────────────────────
st.subheader("🔭 第一节 — 为什么探索性分析很重要")

c1, c2, c3 = st.columns(3)
for col, icon, title, body in [
    (c1, "🕵️", "检测离群样本",
     "看起来与其他样本差异很大的样本，可能存在质量问题、污染或意外的生物学变异。"),
    (c2, "🔄", "揭示批次效应",
     "如果样本按测序批次或文库制备日期而非生物学分组聚集，则需要对此进行建模处理。"),
    (c3, "✅", "验证研究设计",
     "PCA 确认你的分组是否按预期分离——在运行差异表达分析前给你信心。"),
]:
    col.markdown(f"""
<div style="background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;
padding:1rem;height:150px;">
<div style="font-size:1.5rem">{icon}</div>
<div style="font-weight:600;margin:0.3rem 0">{title}</div>
<div style="font-size:0.85rem;color:#475569">{body}</div>
</div>
""", unsafe_allow_html=True)

st.divider()

# ── 第二节 — 数据集概览 ──────────────────────────────────────────────────────
st.subheader("📋 第二节 — 数据集概览")

col_l, col_r = st.columns(2)
with col_l:
    st.markdown(f"""
| | |
|---|---|
| **使用的基因数** | {len(counts_use):,} |
| **样本数** | {len(counts_use.columns)} |
| **元数据列** | {', '.join(metadata.columns)} |
| **分组** | {', '.join(metadata['groupA'].unique())} |
| **是否应用过滤** | {'是（来自第3课）' if filter_res and use_filtered.startswith('使用过滤') else '默认（最低10计数，5个样本）'} |
""")
with col_r:
    st.info("""
**运行差异表达分析前，始终检查：**
- 同组样本是否聚集在一起？
- 是否存在由技术变量引起的意外分离？
- 是否有明显的离群样本？
    """)

st.divider()

# ── 第三节 — 为什么在 PCA 前要转换计数 ───────────────────────────────────────
st.subheader("🔢 第三节 — 为什么在 PCA 前要转换计数")

col_why, col_effect = st.columns(2)
with col_why:
    st.markdown("""
**原始计数不适合直接用于 PCA，因为：**
- 不同样本的文库大小不同——在 500 万 reads 文库中有 1000 个计数的基因，
  与在 2000 万 reads 文库中有同样计数的基因，表达量并不相同
- 如果直接使用计数，高表达基因会主导 PCA
- 计数数据的方差随均值变化（异方差性）

**log₂(CPM + 1) 通过以下方式解决这一问题：**
- 对文库大小进行归一化（CPM 步骤）
- 压缩动态范围（log 步骤）
- 在不同表达水平上稳定方差

**聚焦于高变异基因：**
- 保留信息量最大的基因
- 减少来自未表达/组成型基因的噪声
- 使 PCA 更快、更清晰
    """)
with col_effect:
    gene_means_raw = counts_use.mean(axis=1)
    gene_vars_raw  = counts_use.var(axis=1)
    gene_means_log = log_expr.mean(axis=1)
    gene_vars_log  = log_expr.var(axis=1)

    sample_idx = np.random.default_rng(42).choice(
        len(counts_use), size=min(1500, len(counts_use)), replace=False)

    fig_mv = go.Figure()
    fig_mv.add_trace(go.Scatter(
        x=np.log2(gene_means_raw.iloc[sample_idx] + 1),
        y=np.log2(gene_vars_raw.iloc[sample_idx] + 1),
        mode="markers", name="原始计数",
        marker=dict(color="#94a3b8", size=3, opacity=0.5),
    ))
    fig_mv.add_trace(go.Scatter(
        x=gene_means_log.iloc[sample_idx],
        y=gene_vars_log.iloc[sample_idx],
        mode="markers", name="log-CPM",
        marker=dict(color="#3b82f6", size=3, opacity=0.5),
    ))
    fig_mv.update_layout(
        xaxis_title="平均表达量（log 尺度）",
        yaxis_title="方差（log 尺度）",
        height=300, margin=dict(t=10),
        legend=dict(x=0.6, y=0.95),
        title="均值-方差关系",
    )
    st.plotly_chart(fig_mv, use_container_width=True)
    st.caption("log-CPM（蓝色）的均值-方差趋势更平坦——更适合 PCA。")

st.divider()

# ── 第四节 — PCA 分析 ────────────────────────────────────────────────────────
st.subheader("🔵 第四节 — PCA 分析")

col_pca, col_scree = st.columns([3, 2])

with col_pca:
    shape_col = None if shape_by == "无" else shape_by
    fig_pca = px.scatter(
        pca_df, x=pc_x, y=pc_y,
        color=color_by,
        symbol=shape_col,
        text="sample" if show_labels else None,
        labels={
            pc_x: f"{pc_x} ({explained[int(pc_x[-1])-1]*100:.1f}%)",
            pc_y: f"{pc_y} ({explained[int(pc_y[-1])-1]*100:.1f}%)",
        },
        title=f"PCA — 按 {color_by} 着色",
        height=430,
    )
    fig_pca.update_traces(marker=dict(size=10), textposition="top center")
    fig_pca.update_layout(margin=dict(t=40))
    st.plotly_chart(fig_pca, use_container_width=True)

with col_scree:
    exp_df = pd.DataFrame({
        "PC":    [f"PC{i+1}" for i in range(len(explained))],
        "方差 %": explained * 100,
        "累积 %": np.cumsum(explained) * 100,
    })
    fig_scree = go.Figure()
    fig_scree.add_trace(go.Bar(
        x=exp_df["PC"], y=exp_df["方差 %"],
        name="每 PC 方差", marker_color="#3b82f6",
    ))
    fig_scree.add_trace(go.Scatter(
        x=exp_df["PC"], y=exp_df["累积 %"],
        name="累积方差", line=dict(color="#ef4444"),
        yaxis="y2", mode="lines+markers",
    ))
    fig_scree.update_layout(
        title="碎石图",
        yaxis=dict(title="解释方差 (%)"),
        yaxis2=dict(title="累积 (%)", overlaying="y",
                    side="right", range=[0, 100]),
        height=430, margin=dict(t=40),
        legend=dict(x=0.5, y=0.95),
    )
    st.plotly_chart(fig_scree, use_container_width=True)

    st.dataframe(
        exp_df.style.format({"方差 %": "{:.2f}", "累积 %": "{:.2f}"}),
        use_container_width=True, hide_index=True, height=180,
    )

with st.expander("💡 如何解读这张 PCA 图"):
    st.markdown(f"""
- **{pc_x}** 解释了 **{explained[int(pc_x[-1])-1]*100:.1f}%** 的方差——
  该数据集中最大的单一表达变异来源。
- **{pc_y}** 解释了 **{explained[int(pc_y[-1])-1]*100:.1f}%**——第二大变异来源。
- 距离近的样本具有**更相似的表达谱**。
- PCA 不会自动将轴标记为"生物学"或"批次"——
  你需要叠加元数据才能解释分离的含义。
- 在侧边栏中切换着色变量，看看不同元数据变量与 PC 结构的对应关系。
    """)

st.divider()

# ── 第五节 — 离群样本检测 ───────────────────────────────────────────────────
st.subheader("🚨 第五节 — 离群样本检测")

outlier_z = outlier_scores(scores, pcs=[pc_x, pc_y])
pca_df["outlier_score"] = pca_df["sample"].map(outlier_z)
threshold = 2.0

col_out1, col_out2 = st.columns(2)

with col_out1:
    pca_df["outlier"] = pca_df["outlier_score"] > threshold
    fig_out = px.scatter(
        pca_df, x=pc_x, y=pc_y,
        color="outlier",
        color_discrete_map={True: "#ef4444", False: "#94a3b8"},
        hover_name="sample",
        hover_data={col: True for col in metadata.columns
                    if col in pca_df.columns},
        labels={
            pc_x: f"{pc_x} ({explained[int(pc_x[-1])-1]*100:.1f}%)",
            pc_y: f"{pc_y} ({explained[int(pc_y[-1])-1]*100:.1f}%)",
        },
        title="潜在离群样本（z 分数 > 2）",
        height=370,
    )
    fig_out.update_traces(marker=dict(size=9))
    fig_out.update_layout(margin=dict(t=40))
    st.plotly_chart(fig_out, use_container_width=True)

with col_out2:
    outlier_df = pca_df[["sample", "outlier_score", "outlier"]].copy()
    if color_by in pca_df.columns:
        outlier_df[color_by] = pca_df[color_by].values
    outlier_df = outlier_df.sort_values("outlier_score", ascending=False)

    st.markdown("**所有样本的离群分数**")
    st.dataframe(
        outlier_df.style
        .format({"outlier_score": "{:.2f}"})
        .apply(lambda col: [
            "background-color: #fee2e2" if v else ""
            for v in col
        ], subset=["outlier"]),
        use_container_width=True, hide_index=True, height=320,
    )

    n_outliers = outlier_df["outlier"].sum()
    if n_outliers > 0:
        st.warning(
            f"**检测到 {n_outliers} 个潜在离群样本。** "
            "在决定是否排除之前，请检查其元数据和原始表达量。"
        )
    else:
        st.success("在 z 分数 > 2 的标准下，未检测到明显的离群样本。")

st.caption("""
⚠️ 此处的离群检测基于 PCA 中到质心的距离。
这是一个探索性标记，不是确定性的排除标准。
在排除任何样本之前，务必检查其生物学背景和 QC 指标。
""")

st.divider()

# ── 第六节 — 样本间距离热图 ─────────────────────────────────────────────────
st.subheader("🗺️ 第六节 — 样本间距离")

dist_mat = sample_distance_matrix(log_expr)
corr_mat = sample_correlation_matrix(log_expr)

tab_dist, tab_corr = st.tabs(["距离热图", "相关性热图"])

with tab_dist:
    fig_dist = px.imshow(
        dist_mat,
        color_continuous_scale="Blues_r",
        title="样本间欧氏距离（log-CPM）",
        height=500,
        labels=dict(color="距离"),
    )
    fig_dist.update_layout(margin=dict(t=40))
    st.plotly_chart(fig_dist, use_container_width=True)
    st.caption(
        "颜色越深 = 越相似。同组样本通常应比不同组样本彼此更相似（颜色更深）。"
    )

with tab_corr:
    fig_corr = px.imshow(
        corr_mat,
        color_continuous_scale="RdBu",
        zmin=0.8, zmax=1.0,
        title="样本间 Pearson 相关性（log-CPM）",
        height=500,
        labels=dict(color="r"),
    )
    fig_corr.update_layout(margin=dict(t=40))
    st.plotly_chart(fig_corr, use_container_width=True)
    st.caption(
        "接近 1.0 的值表示样本高度相似。"
        "相关性较低的样本（异常行/列）可能是离群样本。"
    )

st.divider()

# ── 第七节 — 元数据解读 ──────────────────────────────────────────────────────
st.subheader("🏷️ 第七节 — 在 PCA 上解读元数据")

st.markdown("""
切换以下元数据列，理解是什么驱动了样本分离。
""")

meta_numeric = [c for c in metadata.columns
                if pd.api.types.is_numeric_dtype(metadata[c])]
meta_categ   = [c for c in metadata.columns
                if not pd.api.types.is_numeric_dtype(metadata[c])]

n_cols   = min(len(meta_cols), 4)
tab_list = st.tabs(meta_cols[:n_cols])

for tab, col_name in zip(tab_list, meta_cols[:n_cols]):
    with tab:
        pca_tab = pca_plot_df(scores, metadata, pc_x, pc_y)
        color_seq = px.colors.qualitative.Set2

        if col_name in meta_numeric:
            fig_t = px.scatter(
                pca_tab, x=pc_x, y=pc_y, color=col_name,
                color_continuous_scale="Viridis",
                hover_name="sample", height=360,
                labels={pc_x: f"{pc_x} ({explained[int(pc_x[-1])-1]*100:.1f}%)",
                        pc_y: f"{pc_y} ({explained[int(pc_y[-1])-1]*100:.1f}%)"},
                title=f"按 {col_name} 着色",
            )
        else:
            fig_t = px.scatter(
                pca_tab, x=pc_x, y=pc_y, color=col_name,
                color_discrete_sequence=color_seq,
                hover_name="sample", height=360,
                labels={pc_x: f"{pc_x} ({explained[int(pc_x[-1])-1]*100:.1f}%)",
                        pc_y: f"{pc_y} ({explained[int(pc_y[-1])-1]*100:.1f}%)"},
                title=f"按 {col_name} 着色",
            )
        fig_t.update_traces(marker=dict(size=10))
        fig_t.update_layout(margin=dict(t=40))
        st.plotly_chart(fig_t, use_container_width=True)

        guidance = {
            "groupA": "✅ 如果分组在此处分离，支持存在真实生物学差异。",
            "batch":  "⚠️ 如果批次分离，在差异表达分析前考虑批次校正。",
            "donor":  "ℹ️ 供体效应在配对设计中很常见——使用配对检验。",
            "sex":    "ℹ️ 性别差异可能是真实的生物学信号——检查其是否相关。",
            "age":    "ℹ️ 年龄梯度可以反映生物学——值得在模型中注明。",
        }
        msg = guidance.get(col_name, f"ℹ️ 考虑 {col_name} 是否应该被纳入模型。")
        st.caption(msg)

st.divider()

# ── 第八节 — 常见错误 ────────────────────────────────────────────────────────
st.subheader("⚠️ 第八节 — 常见 PCA 解读错误")

st.error("**假设分离一定意味着生物学差异。** 批次效应、文库大小差异或样本调换都可能导致 PCA 分离。")
st.warning("**忽略元数据。** PCA 轴没有固有的生物学含义，只有叠加元数据后才能解读。")
st.warning("**过度解读弱 PC。** 解释 3% 方差的 PC3 可能只是噪声。")
st.info("**对原始计数运行 PCA。** 始终先进行转换——原始计数因文库大小效应会给出误导性的 PCA 结果。")
st.info("**将 PCA 视为统计检验。** PCA 是探索性的——它产生假说，而非结论。")

st.divider()

# ── 第九节 — 关键总结 ────────────────────────────────────────────────────────
st.subheader("📌 关键总结")

cols = st.columns(4)
msgs = [
    ("🔭", "先做探索性分析", "在运行差异表达分析之前，始终先探索你的数据。"),
    ("🔢", "转换计数", "log-CPM 可稳定 PCA 所需的方差。"),
    ("🏷️", "利用元数据", "PCA 解读需要元数据背景。"),
    ("🚨", "检查离群样本", "在差异表达之前识别异常样本——不要忽视它们。"),
]
for col, (icon, title, body) in zip(cols, msgs):
    col.markdown(f"""
<div style="background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;
padding:0.9rem;text-align:center;">
<div style="font-size:1.5rem">{icon}</div>
<div style="font-weight:600;font-size:0.9rem;margin:0.3rem 0">{title}</div>
<div style="font-size:0.82rem;color:#475569">{body}</div>
</div>
""", unsafe_allow_html=True)

st.divider()

# ══════════════════════════════════════════════════════════════════════════════
# 第9节 — 方差分解（类 PVCA 分析）
# ══════════════════════════════════════════════════════════════════════════════

from utils.pca_utils import variance_decomposition

st.subheader("📊 第9节 — 方差分解")
st.markdown("""
PCA 散点图只能告诉你样本是否按批次或实验条件聚集。
方差分解则能**定量**回答：每个因素解释了多少比例的总表达方差？

**方法（类 PVCA，Bushel et al. 2009）：**
1. 选取能解释 ≥ 80% 方差的主成分（PC）。
2. 对每个 PC，以各实验因素为自变量进行 OLS 回归。
3. 计算每个因素的偏 R²（partial R²）。
4. 按各 PC 的方差贡献加权汇总。

**解读示例：**
- 批次解释 35%，实验条件解释 40% → 需要进行批次校正。
- 批次解释 8%，实验条件解释 55% → 批次影响较小，可直接分析。
""")

with st.expander("⚙️ 配置分解参数", expanded=True):
    vd_factors = []
    meta_cols_vd = [c for c in metadata.columns if metadata[c].nunique() < 20]
    if meta_cols_vd:
        vd_factors = st.multiselect(
            "选择分解因素",
            options=meta_cols_vd,
            default=meta_cols_vd[:3] if len(meta_cols_vd) >= 3 else meta_cols_vd,
            help="选择元数据列，分类变量和连续变量均可。",
        )
    vd_cumvar = st.slider(
        "纳入的累积方差比例（%）",
        50, 95, 80, 5,
        key="vd_cumvar",
    )

if vd_factors:
    with st.spinner("正在运行方差分解…"):
        try:
            vd_result = variance_decomposition(
                scores, explained,
                metadata,
                factors=vd_factors,
                cumvar_threshold=vd_cumvar / 100,
            )

            col_vd1, col_vd2 = st.columns([1, 1])
            with col_vd1:
                st.dataframe(
                    vd_result.style.format({"pct_variance_explained": "{:.1f}%"}),
                    use_container_width=True,
                    hide_index=True,
                )
            with col_vd2:
                fig_vd = px.bar(
                    vd_result,
                    x="factor",
                    y="pct_variance_explained",
                    color="factor",
                    labels={"pct_variance_explained": "方差解释比例（%）",
                            "factor": "因素"},
                    color_discrete_sequence=["#2563eb", "#16a34a", "#dc2626",
                                             "#d97706", "#7c3aed", "#94a3b8"],
                )
                fig_vd.update_layout(
                    height=350,
                    showlegend=False,
                    yaxis=dict(range=[0, 100]),
                )
                st.plotly_chart(fig_vd, use_container_width=True)

            top = vd_result[vd_result["factor"] != "Residual / other"].iloc[0]
            st.info(
                f"**主导因素：** {top['factor']} 解释了 {top['pct_variance_explained']:.1f}% 的总方差。"
                "若批次因素超过 ~30%，建议在差异表达分析前进行批次校正。"
            )
        except Exception as e:
            st.error(f"方差分解失败：{e}")
else:
    st.info("请在上方至少选择一个因素以运行分解分析。")

st.divider()
col_n1, col_n2 = st.columns(2)
with col_n1:
    st.page_link("pages/06_Batch_Correction.py", label="← 第6课：批次校正", icon="🔄")
with col_n2:
    st.page_link("pages/07_Differential_Expression.py", label="第7课：差异表达 →", icon="🧪")
