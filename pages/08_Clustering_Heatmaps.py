"""
pages/08_Clustering_Heatmaps.py — 第8课：聚类与热图
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from utils.data_loader import init_session_data
from utils.filtering import filter_low_expression, log_cpm
from utils.clustering_utils import (
    cluster_matrix,
    get_leaf_order,
    row_zscore,
    assign_modules,
)
from utils.exploration import sample_correlation_matrix

st.set_page_config(
    page_title="第8课 — 聚类与热图",
    page_icon="🗺️",
    layout="wide",
)

init_session_data()

st.title("🗺️ 第8课 — 聚类与热图")
st.markdown("""
> **学习目标：** 理解如何将 RNA-seq 表达模式整理成可解读的基因和样本聚类，
> 以及热图如何概括复杂数据。
""")

counts_raw = st.session_state["counts"]
metadata   = st.session_state["metadata"]
filter_res = st.session_state.get("filter_results", {})
de_results = st.session_state.get("de_results_lesson7")

MAX_GENES_LIVE = 200

# ── 第一节 ────────────────────────────────────────────────────────────────────
st.subheader("🧩 第一节 — 为什么差异表达分析后需要聚类")

c1, c2, c3 = st.columns(3)
for col, icon, title, body in [
    (c1, "📋", "差异表达给出列表",
     "差异表达告诉你*哪些*基因发生变化——聚类展示它们如何在样本间*协同*变化。"),
    (c2, "🔗", "聚类揭示模式",
     "具有相似表达谱的基因可能共享生物学功能或调控通路。"),
    (c3, "👥", "样本结构",
     "样本聚类揭示生物学分组分离、批次效应或意外的离群样本。"),
]:
    col.markdown(f"""
<div style="background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;
padding:1rem;height:130px;">
<div style="font-size:1.4rem">{icon}</div>
<div style="font-weight:600;margin:0.3rem 0">{title}</div>
<div style="font-size:0.85rem;color:#475569">{body}</div>
</div>
""", unsafe_allow_html=True)

st.divider()

sources = []
if de_results is not None and de_results["significant"].sum() > 0:
    sources.append("显著差异表达基因（来自第7课）")
sources.append("高变异基因")
sources.append("内置演示数据基因集")

with st.sidebar:
    st.header("🔧 聚类参数")
    data_source  = st.selectbox("基因来源", sources)
    n_genes      = st.slider("聚类基因数", 10, 500, 50, 10)

    if n_genes > MAX_GENES_LIVE:
        st.warning(f"⚠️ > {MAX_GENES_LIVE} 个基因：请点击 **▶ 运行聚类** 按钮。")

    st.divider()
    st.subheader("数据变换")
    use_scaling = st.checkbox("行 Z 分数标准化", value=True)

    st.divider()
    st.subheader("聚类设置")
    dist_metric  = st.selectbox("距离度量", ["pearson", "spearman", "euclidean"])
    linkage_meth = st.selectbox("连接方法", ["average", "complete", "ward"])
    cluster_genes   = st.checkbox("基因聚类", value=True)
    cluster_samples = st.checkbox("样本聚类", value=True)

    st.divider()
    st.subheader("基因模块")
    n_modules = st.slider("基因模块数（k）", 2, 10, 3)
    run_btn   = st.button("▶ 运行聚类", type="primary")

st.subheader("📥 第二节 — 输入数据与基因选择")

counts_use = filter_res.get("counts_filtered", filter_low_expression(counts_raw, 10, 5))
log_expr   = log_cpm(counts_use)

if data_source.startswith("显著差异表达") and de_results is not None:
    sig_genes  = [g for g in de_results[de_results["significant"]].index if g in log_expr.index]
    gene_pool  = sig_genes
else:
    gene_pool  = log_expr.var(axis=1).nlargest(1000).index.tolist()

selected_genes = gene_pool[:min(n_genes, len(gene_pool))]
expr_subset    = log_expr.loc[selected_genes]

c1, c2, c3, c4 = st.columns(4)
c1.metric("选中基因数", f"{len(selected_genes):,}")
c2.metric("样本数",     f"{expr_subset.shape[1]}")
c3.metric("来源",       data_source.split("（")[0][:20])
c4.metric("差异表达结果", "有" if de_results is not None else "无")

if len(selected_genes) == 0:
    st.error("没有可用基因。请尝试不同的来源或检查数据是否已加载。")
    st.stop()

st.divider()

# ── 第三节 — 行标准化 ────────────────────────────────────────────────────────
st.subheader("📏 第三节 — 数据变换与行标准化")

with st.expander("📖 为什么行标准化很重要", expanded=False):
    col_why, col_ex = st.columns(2)
    with col_why:
        st.markdown("""
**行 Z 分数标准化：** `z = (x − 均值) / 标准差`

- 每个基因以 0 为中心
- ±1 = 比该基因均值高/低一个标准差
- 无论绝对表达量如何，所有基因视觉上可比较
- 使用**发散色标**（蓝—白—红）

| 目标 | 是否标准化 |
|---|---|
| 比较绝对水平 | 否 |
| 比较上调/下调模式 | 是 |
| 大多数论文中的热图 | 是（Z 分数） |
        """)
    with col_ex:
        demo   = pd.DataFrame({"S1":[1000,5,50],"S2":[1200,3,80],"S3":[800,8,30]},
                               index=["高表达","低表达","中表达"])
        demo_z = row_zscore(demo)
        for mat, title, cs, mid in [
            (demo,   "未标准化", "Viridis", None),
            (demo_z, "Z 分数标准化", "RdBu_r",  0),
        ]:
            kw = dict(z=mat.values, x=list(mat.columns), y=list(mat.index),
                      colorscale=cs, showscale=True)
            if mid is not None:
                kw["zmid"] = mid
            fig = go.Figure(go.Heatmap(**kw))
            fig.update_layout(title=title, height=160, margin=dict(t=30,l=80))
            st.plotly_chart(fig, use_container_width=True)

st.divider()

# ── 第四节 — 距离度量 ────────────────────────────────────────────────────────
st.subheader("📐 第四节 — 距离度量与连接方法")

with st.expander("📖 度量和连接选择如何影响聚类结果", expanded=False):
    st.markdown("""
| 度量 | 衡量什么 | 最适用于 |
|---|---|---|
| **Pearson** | 线性表达谱形状 | 共调控基因发现 |
| **Spearman** | 基于排名的表达谱形状 | 稳健、非线性模式 |
| **Euclidean** | 绝对量级差异 | 对量级敏感的分组 |

| 连接方法 | 行为 |
|---|---|
| **Average** | 簇中心之间的距离——均衡 |
| **Complete** | 最大配对距离——紧凑的簇 |
| **Ward** | 最小化簇内方差——通常产生最清晰的模块 |

不同组合可能给出非常不同的结果。尝试多种方法来检验稳健性。
    """)

st.divider()

# ── 第五、六节 — 聚类 + 热图 ─────────────────────────────────────────────────
st.subheader("🔥 第五、六节 — 聚类与交互式热图")

if use_scaling:
    plot_matrix = row_zscore(expr_subset)
    colorscale  = "RdBu_r"
    zmid        = 0
    color_label = "Z 分数"
else:
    plot_matrix = expr_subset.copy()
    colorscale  = "Viridis"
    zmid        = None
    color_label = "log-CPM"

gene_order   = list(plot_matrix.index)
sample_order = list(plot_matrix.columns)
lkg_genes    = None
lkg_samp     = None

should_cluster = run_btn or (n_genes <= MAX_GENES_LIVE)

if should_cluster and (cluster_genes or cluster_samples):
    mat = plot_matrix.values
    if cluster_genes and len(gene_order) > 1:
        try:
            lkg_genes  = cluster_matrix(mat, metric=dist_metric, method=linkage_meth)
            gene_order = [plot_matrix.index[i] for i in get_leaf_order(lkg_genes)]
        except Exception as e:
            st.warning(f"基因聚类失败：{e}")
            lkg_genes = None
    if cluster_samples and len(sample_order) > 1:
        try:
            lkg_samp    = cluster_matrix(mat.T, metric=dist_metric, method=linkage_meth)
            sample_order = [plot_matrix.columns[i] for i in get_leaf_order(lkg_samp)]
        except Exception as e:
            st.warning(f"样本聚类失败：{e}")
            lkg_samp = None
elif n_genes > MAX_GENES_LIVE and not run_btn:
    st.info(f"⚡ 超过 {MAX_GENES_LIVE} 个基因时，聚类已暂停。点击侧边栏中的 **▶ 运行聚类**。")

plot_ordered = plot_matrix.loc[gene_order, sample_order]

hover_text = []
for gene in gene_order:
    row = []
    for sample in sample_order:
        raw_val = expr_subset.loc[gene, sample]
        z_val   = plot_matrix.loc[gene, sample]
        if use_scaling:
            row.append(f"基因：{gene}<br>样本：{sample}<br>log-CPM：{raw_val:.2f}<br>Z 分数：{z_val:.2f}")
        else:
            row.append(f"基因：{gene}<br>样本：{sample}<br>log-CPM：{raw_val:.2f}")
    hover_text.append(row)

hm_kw = dict(
    z=plot_ordered.values, x=sample_order, y=gene_order,
    colorscale=colorscale, text=hover_text, hoverinfo="text",
    colorbar=dict(title=color_label, len=0.8),
)
if zmid is not None:
    hm_kw["zmid"] = zmid

fig_hm = go.Figure(go.Heatmap(**hm_kw))
fig_hm.update_layout(
    height=max(400, min(900, len(gene_order) * 14 + 100)),
    margin=dict(l=130, r=60, t=50, b=110),
    xaxis=dict(tickangle=45, tickfont=dict(size=9 if len(sample_order) > 20 else 11)),
    yaxis=dict(tickfont=dict(size=8 if len(gene_order) > 50 else 10),
               showticklabels=(len(gene_order) <= 60)),
    title=f"{'Z 分数标准化' if use_scaling else '未标准化'}热图 — "
          f"{len(gene_order)} 个基因 × {len(sample_order)} 个样本",
)
st.plotly_chart(fig_hm, use_container_width=True)

if len(gene_order) > 60:
    st.caption("为提高可读性，基因标签已隐藏。将鼠标悬停在格子上可查看基因 ID。")

meta_cols = [c for c in ["groupA", "batch"] if c in metadata.columns]
if meta_cols:
    with st.expander("样本注释"):
        ann_samples = [s for s in sample_order if s in metadata.index]
        st.dataframe(metadata.loc[ann_samples, meta_cols].T, use_container_width=True)

st.divider()

# ── 第七节 ────────────────────────────────────────────────────────────────────
st.subheader("🔄 第七节 — 样本聚类 vs 基因聚类")

st.markdown("""
| 模式 | 揭示什么 |
|---|---|
| **基因聚类 开启** | 共表达基因群——潜在的共调控或共享功能 |
| **样本聚类 开启** | 队列结构——生物学分组、批次效应、离群样本 |
| **两者都开启** | 块状结构：哪些基因群在哪些样本组中上调/下调 |
| **两者都关闭** | 保持原始顺序——可用作参照对比 |

在侧边栏中切换这些选项，探索热图如何变化。
""")

st.divider()

# ── 第八节 — 基因模块 ────────────────────────────────────────────────────────
st.subheader("🧬 第八节 — 识别基因模块")

if lkg_genes is not None and len(gene_order) > 1:
    module_assignments = assign_modules(lkg_genes, n_modules, gene_order)
    mod_summary = module_assignments.value_counts().sort_index().reset_index()
    mod_summary.columns = ["模块", "基因数"]

    col_a, col_b = st.columns([2, 3])
    with col_a:
        fig_mod = px.bar(
            mod_summary, x="模块", y="基因数",
            color="模块", title=f"{n_modules} 个基因模块",
            height=280,
        )
        fig_mod.update_layout(margin=dict(t=40), showlegend=False)
        st.plotly_chart(fig_mod, use_container_width=True)

    with col_b:
        sel_mod   = st.selectbox(
            "查看模块",
            sorted(module_assignments.unique()),
            format_func=lambda x: f"模块 {x}（{(module_assignments==x).sum()} 个基因）",
        )
        mod_genes = module_assignments[module_assignments == sel_mod].index.tolist()
        st.markdown(f"**模块 {sel_mod} — {len(mod_genes)} 个基因：**")
        st.dataframe(pd.DataFrame({"基因": mod_genes}),
                     use_container_width=True, height=200, hide_index=True)
        st.download_button(
            f"⬇️ 模块 {sel_mod} 基因列表（.txt）",
            data="\n".join(mod_genes),
            file_name=f"module_{sel_mod}_genes.txt",
            mime="text/plain",
        )

    st.download_button(
        "⬇️ 所有模块分配（CSV）",
        data=module_assignments.reset_index().rename(columns={"index":"gene"}).to_csv(index=False),
        file_name="gene_module_assignments.csv",
        mime="text/csv",
    )
    st.info("模块基因列表可作为第9课通路富集分析的输入。")
else:
    st.info("启用**基因聚类**并点击 **▶ 运行聚类** 以识别模块。")

st.divider()

# ── 第九节 — 样本相关性 ──────────────────────────────────────────────────────
st.subheader("🗺️ 第九节 — 样本间相关性")

corr_mat = sample_correlation_matrix(log_expr)
corr_ord = corr_mat.loc[sample_order, sample_order]
fig_corr = px.imshow(
    corr_ord, color_continuous_scale="RdBu",
    zmin=0.8, zmax=1.0,
    title="样本 Pearson 相关性（log-CPM，所有过滤后基因）",
    height=460, labels=dict(color="r"),
)
fig_corr.update_layout(margin=dict(t=40))
st.plotly_chart(fig_corr, use_container_width=True)
st.caption("如果样本按批次而非生物学聚集，请回顾第6课——批次校正。")

st.divider()

# ── 第十节 — 常见错误 ────────────────────────────────────────────────────────
st.subheader("⚠️ 第十节 — 常见解读错误")

st.error("**基因太多 → 热图无法阅读。** 坚持使用 50–200 个有意义的基因。")
st.warning("**忘记标准化会改变解读方式。** 标准化和未标准化的热图不可相互比较。")
st.warning("**过度解读颜色块。** 红色 = 该基因相对较高，不一定具有生物学重要性。")
st.info("**假设一个结果就是真相。** 不同的度量、连接方法和 k 值可以给出不同的模式。")
st.info("**对随机基因聚类。** 始终从显著差异表达或高变异基因开始。")

st.divider()

# ── 第十一节 — 关键总结 ──────────────────────────────────────────────────────
st.subheader("📌 关键总结")

cols = st.columns(5)
for col, (icon, title, body) in zip(cols, [
    ("🗺️", "差异表达后的模式", "聚类将基因列表整理成可解读的表达程序。"),
    ("📏", "Z 分数标准化", "将每个基因归中，突出样本间的相对上调/下调。"),
    ("📐", "度量很重要", "Pearson/Spearman 关注形状；Euclidean 关注量级。"),
    ("🧬", "模块 → 富集", "基因模块可直接用于第9课的通路分析。"),
    ("🔍", "仅用于探索", "聚类产生假说——它不是统计检验。"),
]):
    col.markdown(f"""
<div style="background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;
padding:0.9rem;text-align:center;">
<div style="font-size:1.5rem">{icon}</div>
<div style="font-weight:600;font-size:0.9rem;margin:0.3rem 0">{title}</div>
<div style="font-size:0.82rem;color:#475569">{body}</div>
</div>
""", unsafe_allow_html=True)

st.divider()
col_n1, col_n2 = st.columns(2)
with col_n1:
    st.page_link("pages/07_Differential_Expression.py",
                 label="← 第7课：差异表达分析", icon="🧪")
with col_n2:
    st.page_link("pages/09_Functional_Enrichment.py",
                 label="第9课：功能富集分析 →", icon="🧩")
