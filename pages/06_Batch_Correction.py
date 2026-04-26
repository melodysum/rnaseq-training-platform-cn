"""
pages/06_Batch_Correction.py — 第6课：批次效应与校正
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from utils.data_loader import init_session_data
from utils.filtering import filter_low_expression
from utils.batch_effects import (
    to_log_cpm,
    simple_batch_correction,
    variance_explained_by,
    top_variable_genes,
)
from utils.pca_utils import run_pca, pca_plot_df
from utils.simulation import simulate_batch_data

st.set_page_config(
    page_title="第6课 — 批次校正",
    page_icon="🔄",
    layout="wide",
)

init_session_data()

st.title("🔄 第6课 — 批次效应与校正")
st.markdown("""
> **学习目标：** 理解技术性批次效应的来源，如何用 PCA 检测它们，如何进行校正，
> 以及——至关重要地——校正何时会出错。
""")

counts_raw = st.session_state["counts"]
metadata   = st.session_state["metadata"]
has_batch  = "batch" in metadata.columns

filter_res = st.session_state.get("filter_results", {})
if filter_res:
    counts_use = filter_res["counts_filtered"]
else:
    counts_use = filter_low_expression(counts_raw, 10, 5)

log_expr = to_log_cpm(counts_use)

# ── 第一节 — 什么是批次效应？ ────────────────────────────────────────────────
st.subheader("🧩 第一节 — 什么是批次效应？")

col_def, col_diag = st.columns([3, 2])
with col_def:
    st.markdown("""
**批次效应**是样本组之间不需要的技术性差异，与你所关心的生物学无关。

它们来自以下因素：
- 不同的测序批次或测序仪
- 不同的操作员或试剂批号
- 不同的文库制备日期
- 不同的实验室或机构

危险在于批次效应可能让样本**按技术变量而非生物学**分离——
从而干扰你的下游分析。
    """)

with col_diag:
    fig_schematic = go.Figure()
    labels = ["生物学信号", "批次效应", "你观察到的"]
    colors = ["#22c55e", "#f87171", "#3b82f6"]
    x_pos  = [0.15, 0.50, 0.85]
    for x, label, color in zip(x_pos, labels, colors):
        fig_schematic.add_shape(
            type="rect", x0=x-0.12, x1=x+0.12, y0=0.3, y1=0.7,
            fillcolor=color, opacity=0.85, line_width=0,
        )
        fig_schematic.add_annotation(
            x=x, y=0.5, text=label, showarrow=False,
            font=dict(color="white", size=11, family="Inter"),
            align="center",
        )
    fig_schematic.add_annotation(
        x=0.48, y=0.5, ax=0.27, ay=0.5,
        arrowhead=2, arrowcolor="#64748b", arrowwidth=2,
        text="+", font=dict(size=16, color="#64748b"), showarrow=True,
    )
    fig_schematic.add_annotation(
        x=0.73, y=0.5, ax=0.62, ay=0.5,
        arrowhead=2, arrowcolor="#64748b", arrowwidth=2,
        text="=", font=dict(size=16, color="#64748b"), showarrow=True,
    )
    fig_schematic.update_layout(
        height=200, margin=dict(l=0, r=0, t=10, b=10),
        xaxis=dict(visible=False, range=[0, 1]),
        yaxis=dict(visible=False, range=[0, 1]),
        plot_bgcolor="white",
    )
    st.plotly_chart(fig_schematic, use_container_width=True)
    st.caption("观测到的表达量 = 真实生物学 + 技术性批次偏移 + 噪声")

st.divider()

# ── 第二节 — 数据集概览 ──────────────────────────────────────────────────────
st.subheader("📋 第二节 — 数据集概览")

col_meta, col_warn = st.columns([2, 2])
with col_meta:
    n_genes   = len(counts_use)
    n_samples = len(counts_use.columns)
    meta_cols = list(metadata.columns)
    st.markdown(f"""
| | |
|---|---|
| **基因数（过滤后）** | {n_genes:,} |
| **样本数** | {n_samples} |
| **元数据列** | {', '.join(meta_cols)} |
| **是否检测到批次列** | {'✅ 是' if has_batch else '❌ 否'} |
""")

with col_warn:
    if has_batch:
        batch_counts = metadata["batch"].value_counts()
        st.success(
            f"✅ 找到批次列，共 **{len(batch_counts)} 个批次**：" +
            "、".join(f"{b}（{n} 个样本）" for b, n in batch_counts.items())
        )
    else:
        st.warning(
            "⚠️ 元数据中未找到 `batch` 列。"
            "校正演示将使用模拟数据。"
            "上传包含 `batch` 列的元数据以使用你自己的数据。"
        )

st.divider()

# ── 第三节 — PCA：检测批次效应 ──────────────────────────────────────────────
st.subheader("🔭 第三节 — 用 PCA 检测批次效应")

st.markdown("""
PCA 是检查批次效应时首先应使用的工具。
如果样本按批次而非按生物学分组聚集，那就是一个警示信号。
""")

with st.sidebar:
    st.header("🔧 PCA 参数")
    top_genes = st.slider("用于 PCA 的高变异基因数", 100, min(2000, n_genes), 500, 100)
    pc_x = st.selectbox("X 轴", ["PC1", "PC2", "PC3"], index=0)
    pc_y = st.selectbox("Y 轴", ["PC1", "PC2", "PC3"], index=1)
    show_labels = st.checkbox("显示样本标签", value=False)
    color_by = st.selectbox(
        "着色依据",
        [c for c in ["groupA", "batch", "sex", "donor"] if c in metadata.columns],
    )

scores, explained, _ = run_pca(log_expr, n_components=5, top_var_genes=top_genes)
pca_df = pca_plot_df(scores, metadata, pc_x, pc_y)

col_bio, col_batch = st.columns(2)

def make_pca_fig(df, color_col, title, explained, pc_x, pc_y, show_labels):
    x_var = explained[int(pc_x[-1]) - 1] * 100
    y_var = explained[int(pc_y[-1]) - 1] * 100
    fig = px.scatter(
        df, x=pc_x, y=pc_y,
        color=color_col if color_col in df.columns else None,
        text="sample" if show_labels else None,
        title=title,
        labels={
            pc_x: f"{pc_x} ({x_var:.1f}%)",
            pc_y: f"{pc_y} ({y_var:.1f}%)",
        },
        height=380,
    )
    fig.update_traces(marker=dict(size=9), textposition="top center")
    fig.update_layout(margin=dict(t=40))
    return fig

with col_bio:
    fig_bio = make_pca_fig(pca_df, "groupA", "PCA — 按生物学分组着色",
                           explained, pc_x, pc_y, show_labels)
    st.plotly_chart(fig_bio, use_container_width=True)

with col_batch:
    color_col = "batch" if has_batch else "groupA"
    title_str = "PCA — 按批次着色" if has_batch else "PCA（无批次列——按分组着色）"
    fig_bat = make_pca_fig(pca_df, color_col, title_str,
                           explained, pc_x, pc_y, show_labels)
    st.plotly_chart(fig_bat, use_container_width=True)

exp_df = pd.DataFrame({
    "PC": [f"PC{i+1}" for i in range(len(explained))],
    "解释方差 (%)": explained * 100,
})
fig_var = px.bar(
    exp_df, x="PC", y="解释方差 (%)",
    title="每个 PC 解释的方差",
    color_discrete_sequence=["#3b82f6"],
    height=260,
)
fig_var.update_layout(margin=dict(t=40, b=10))
st.plotly_chart(fig_var, use_container_width=True)

if has_batch:
    r2_bio   = variance_explained_by(scores, metadata, "groupA")
    r2_batch = variance_explained_by(scores, metadata, "batch")
    st.info(
        f"PC1+PC2 解释的方差："
        f"**生物学分组 ≈ {r2_bio:.1%}** | "
        f"**批次 ≈ {r2_batch:.1%}**。"
        + ("批次解释更多——校正可能有帮助。" if r2_batch > r2_bio
           else "生物学占主导——批次效应看起来较小。")
    )

st.divider()

# ── 第四节 — 批次效应模拟 ────────────────────────────────────────────────────
st.subheader("🧪 第四节 — 模拟批次效应")

st.markdown("""
使用下方控件模拟具有不同批次效应强度的表达数据。
比较**平衡设计**（每个分组出现在所有批次中）与**混淆设计**（每个分组只对应一个批次）。
""")

sim_col, sim_out = st.columns([1, 2])

with sim_col:
    sim_n_genes   = st.slider("基因数", 100, 1000, 300, 100, key="sim_genes")
    sim_n_samp    = st.slider("每组样本数", 3, 10, 5, 1, key="sim_samp")
    sim_bio       = st.slider("生物学效应强度", 0.0, 5.0, 2.0, 0.5, key="sim_bio")
    sim_batch     = st.slider("批次效应强度",   0.0, 5.0, 2.0, 0.5, key="sim_bat")
    sim_confound  = st.checkbox("混淆设计（危险！）", value=False, key="sim_conf")

sim_expr, sim_meta = simulate_batch_data(
    n_genes=sim_n_genes,
    n_samples_per_group=sim_n_samp,
    bio_effect=sim_bio,
    batch_effect=sim_batch,
    confounded=sim_confound,
)

sim_scores, sim_exp, _ = run_pca(sim_expr, n_components=3, scale=True)
sim_pca_df = pca_plot_df(sim_scores, sim_meta, "PC1", "PC2")

with sim_out:
    tab_bio, tab_bat = st.tabs(["按生物学着色", "按批次着色"])
    with tab_bio:
        fig_s1 = px.scatter(
            sim_pca_df, x="PC1", y="PC2", color="groupA",
            title="模拟 PCA — 生物学分组",
            labels={"PC1": f"PC1 ({sim_exp[0]*100:.1f}%)",
                    "PC2": f"PC2 ({sim_exp[1]*100:.1f}%)"},
            height=350,
        )
        fig_s1.update_traces(marker=dict(size=10))
        st.plotly_chart(fig_s1, use_container_width=True)
    with tab_bat:
        fig_s2 = px.scatter(
            sim_pca_df, x="PC1", y="PC2", color="batch",
            title="模拟 PCA — 批次",
            labels={"PC1": f"PC1 ({sim_exp[0]*100:.1f}%)",
                    "PC2": f"PC2 ({sim_exp[1]*100:.1f}%)"},
            color_discrete_sequence=["#f97316", "#8b5cf6"],
            height=350,
        )
        fig_s2.update_traces(marker=dict(size=10))
        st.plotly_chart(fig_s2, use_container_width=True)

    if sim_confound:
        st.error(
            "⚠️ **混淆设计**：批次与生物学完全对齐。"
            "此处的任何校正也会去除真实的生物学信号。"
            "这就是为什么实验设计比任何计算方法都更重要。"
        )
    else:
        st.success(
            "✅ **平衡设计**：批次分布在各生物学分组中。"
            "校正可以将技术变异与生物学变异分开。"
        )

st.divider()

# ── 第五节 — 校正前后对比 ────────────────────────────────────────────────────
st.subheader("⚖️ 第五节 — 批次校正前后对比")

if has_batch:
    corrected_expr = simple_batch_correction(log_expr, metadata, "batch", "groupA")
    data_label = "真实上传数据"
else:
    st.info("你的元数据中没有批次列——使用模拟数据进行校正演示。")
    sim_demo, meta_demo = simulate_batch_data(
        n_genes=500, n_samples_per_group=5,
        bio_effect=2.0, batch_effect=3.0, confounded=False,
    )
    log_expr_demo   = sim_demo
    corrected_expr  = simple_batch_correction(sim_demo, meta_demo, "batch", "groupA")
    metadata        = meta_demo
    log_expr        = log_expr_demo
    has_batch       = True
    data_label      = "模拟演示数据"

scores_before, exp_before, _ = run_pca(
    top_variable_genes(log_expr, 500), n_components=3)
scores_after,  exp_after,  _ = run_pca(
    top_variable_genes(corrected_expr, 500), n_components=3)

pca_before = pca_plot_df(scores_before, metadata, "PC1", "PC2")
pca_after  = pca_plot_df(scores_after,  metadata, "PC1", "PC2")

col_before, col_after = st.columns(2)

with col_before:
    st.markdown("**校正前**")
    for color_var, tab_label in [("groupA", "分组"), ("batch", "批次")]:
        if color_var not in pca_before.columns:
            continue
        fig_b = px.scatter(
            pca_before, x="PC1", y="PC2", color=color_var,
            labels={"PC1": f"PC1 ({exp_before[0]*100:.1f}%)",
                    "PC2": f"PC2 ({exp_before[1]*100:.1f}%)"},
            height=320, title=f"按 {color_var} 着色",
        )
        fig_b.update_traces(marker=dict(size=9))
        fig_b.update_layout(margin=dict(t=40))
        st.plotly_chart(fig_b, use_container_width=True)

with col_after:
    st.markdown("**校正后**")
    for color_var in ["groupA", "batch"]:
        if color_var not in pca_after.columns:
            continue
        fig_a = px.scatter(
            pca_after, x="PC1", y="PC2", color=color_var,
            labels={"PC1": f"PC1 ({exp_after[0]*100:.1f}%)",
                    "PC2": f"PC2 ({exp_after[1]*100:.1f}%)"},
            height=320, title=f"按 {color_var} 着色",
        )
        fig_a.update_traces(marker=dict(size=9))
        fig_a.update_layout(margin=dict(t=40))
        st.plotly_chart(fig_a, use_container_width=True)

st.caption(
    f"使用{data_label}。"
    "校正应减少批次聚集，同时保留生物学分离。"
)

st.divider()

# ── 第六节 — 批次效应对差异表达的影响 ──────────────────────────────────────
st.subheader("📉 第六节 — 批次效应如何扭曲差异表达结果")

st.markdown("""
批次效应不只是改变样本在 PCA 中的外观——它们还会**夸大或掩盖真实的差异表达基因**。
""")

col_l, col_r = st.columns(2)
with col_l:
    st.markdown("""
**不进行校正/批次建模：**
- 在批次1中表达量高的基因看起来"在A组上调"，即使它并不真正差异表达
- 批次驱动的假阳性膨胀了你的差异表达基因列表

**进行校正或将批次纳入模型：**
- 在进行分组比较之前减去批次贡献
- 真实的生物学差异变得更清晰

**最佳实践：**
在你的差异表达模型中将批次作为协变量（例如在 DESeq2 设计公式中：
`~ batch + groupA`），即使你已经对数据进行了可视化校正。
    """)

with col_r:
    gene_ex = pd.DataFrame({
        "样本": ["A_批次1", "A_批次2", "B_批次1", "B_批次2"],
        "分组":  ["A", "A", "B", "B"],
        "批次":  ["批次1", "批次2", "批次1", "批次2"],
        "未校正": [8.5, 5.5, 6.5, 3.5],
        "校正后":  [7.0, 7.0, 5.0, 5.0],
    })
    fig_ex = go.Figure()
    for grp, color in [("A", "#3b82f6"), ("B", "#ef4444")]:
        sub = gene_ex[gene_ex["分组"] == grp]
        fig_ex.add_trace(go.Bar(
            name=f"分组 {grp}（未校正）",
            x=sub["样本"], y=sub["未校正"],
            marker_color=color, opacity=0.5,
        ))
        fig_ex.add_trace(go.Bar(
            name=f"分组 {grp}（校正后）",
            x=sub["样本"], y=sub["校正后"],
            marker_color=color, opacity=1.0,
        ))
    fig_ex.update_layout(
        barmode="group",
        title="示例基因：表观差异 vs 真实差异",
        yaxis_title="log-CPM",
        height=300,
        margin=dict(t=40),
        legend=dict(font=dict(size=10)),
    )
    st.plotly_chart(fig_ex, use_container_width=True)
    st.caption(
        "未校正数值（浅色）显示 A 组 vs B 组之间存在较大差异。"
        "去除批次偏移后，真实差异更小。"
    )

st.divider()

# ── 第七节 — 重要注意事项 ────────────────────────────────────────────────────
st.subheader("⚠️ 第七节 — 重要注意事项")

st.error("""
**不要盲目进行批次校正。**

如果你的批次与生物学混淆（例如所有对照在批次1，所有处理样本在批次2），
校正会去除真实的生物学信号。
对于设计糟糕的实验，没有任何计算方法可以补救。
""")

st.warning("""
**要避免的常见错误：**

- 在批次与生物学完全混淆时进行校正
- 校正后忘记在差异表达模型中包含批次
- 将校正后的计数作为 DESeq2 的输入（应使用原始计数并在设计中包含批次）
- 过度校正：若批次效应较弱，校正可能增加噪声而非去除它
""")

st.info("""
**最佳实践：**

- 始终在差异表达模型设计公式中将批次作为协变量
- 仅将校正后的数据用于**可视化**（PCA、热图），不作为差异表达工具的原始输入
- 在校正前后都检查 PCA，确认生物学信号得到保留
- 从一开始就设计平衡批次的研究——校正是最后手段，不是计划的一部分
""")

st.divider()

# ── 第八节 — 关键总结 ────────────────────────────────────────────────────────
st.subheader("📌 关键总结")

t1, t2, t3, t4, t5 = st.columns(5)
for col, icon, text in [
    (t1, "🧬", "批次效应是**技术性的**，不是生物学的"),
    (t2, "🔭", "PCA 揭示数据中的批次结构"),
    (t3, "⚖️", "校正可以**改善可视化**，但需谨慎"),
    (t4, "🎯", "实验设计比任何校正方法都重要"),
    (t5, "📊", "始终在**差异表达设计公式**中建模批次"),
]:
    col.markdown(f"""
<div style="background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;
padding:0.9rem;text-align:center;height:120px;display:flex;
align-items:center;justify-content:center;flex-direction:column;">
<div style="font-size:1.6rem">{icon}</div>
<div style="font-size:0.82rem;color:#334155;margin-top:0.3rem">{text}</div>
</div>
""", unsafe_allow_html=True)

st.divider()
col_nav1, col_nav2 = st.columns(2)
with col_nav1:
    st.page_link("pages/04_FDR.py", label="← 第4课：FDR", icon="📐")
with col_nav2:
    st.page_link("Home.py", label="返回课程首页 →", icon="🏠")
