"""
pages/07_Differential_Expression.py — 第7课：差异表达分析
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from utils.data_loader import init_session_data
from utils.filtering import filter_low_expression
from utils.de_analysis import run_de

st.set_page_config(
    page_title="第7课 — 差异表达分析",
    page_icon="🧪",
    layout="wide",
)

init_session_data()

st.title("🧪 第7课 — 差异表达分析")
st.markdown("""
> **学习目标：** 理解如何估计、检验和解读 RNA-seq 数据中
> 两组之间的基因水平差异。
""")

counts_raw = st.session_state["counts"]
metadata   = st.session_state["metadata"]
filter_res = st.session_state.get("filter_results", {})

with st.sidebar:
    st.header("🔧 差异表达参数")

    if filter_res:
        counts_use = filter_res["counts_filtered"]
        st.caption(f"使用 {len(counts_use):,} 个过滤后的基因（来自第3课）。")
    else:
        counts_use = filter_low_expression(counts_raw, 10, 5)
        st.caption(f"使用 {len(counts_use):,} 个基因（默认过滤）。")

    st.divider()

    groups = sorted(metadata["groupA"].unique())
    ref_group    = st.selectbox("参考组", groups, index=0)
    target_group = st.selectbox("目标组", groups, index=min(1, len(groups)-1))

    st.divider()

    has_batch = "batch" in metadata.columns
    has_donor = "donor" in metadata.columns
    use_batch = st.checkbox("校正批次效应", value=False,
                             disabled=not has_batch,
                             help="需要元数据中有 'batch' 列。")

    st.divider()

    fdr_cutoff = st.slider("FDR 阈值",    0.01, 0.20, 0.05, 0.01)
    lfc_cutoff = st.slider("最小 |log₂FC|", 0.0, 3.0,  1.0,  0.1)

    st.caption(
        "⚠️ 教学演示版。使用 log-CPM + t 检验 + BH 校正。"
        "不能替代 DESeq2 / edgeR / limma-voom。"
    )

# ── 第一节 — 什么是差异表达？ ────────────────────────────────────────────────
st.subheader("🧬 第一节 — 什么是差异表达？")

col_a, col_b = st.columns(2)
with col_a:
    st.markdown("""
**差异表达（DE）分析**回答的问题是：

> *"基因在两个（或多个）生物学分组之间是否表现出系统性的表达水平差异？"*

关键点：
- 我们比较的是**分组均值**，而不是单个样本
- 观察到的差异可能只是**噪声**——需要统计检验
- 多个基因同时被检验 → 需要 **FDR 校正**
- 结果取决于所有上游决策：过滤、归一化、批次处理
    """)
with col_b:
    st.info("""
**本课的差异表达流程：**

1. 过滤低表达基因（第3课）
2. 归一化为 log-CPM
3. 可选：校正批次效应
4. 逐基因 t 检验
5. BH/FDR 多重检验校正
6. 按 |log₂FC| 和 FDR 阈值过滤
    """)

st.divider()

# ── 第二节 — 数据集概览 ──────────────────────────────────────────────────────
st.subheader("📋 第二节 — 数据集与比较设置")

c1, c2, c3, c4 = st.columns(4)
c1.metric("被检验基因数", f"{len(counts_use):,}")
c2.metric("样本数",      f"{len(counts_use.columns)}")
c3.metric("参考组",      ref_group)
c4.metric("目标组",      target_group)

if ref_group == target_group:
    st.error("参考组和目标组必须不同。")
    st.stop()

if has_donor:
    st.info(
        "**检测到供体列。** 将使用**配对 t 检验**"
        "（每个供体：处理组 − 对照组），去除供体间噪声。"
    )
elif has_batch and use_batch:
    st.info("**已启用批次校正。** 检验前将去除批次效应的影响。")

st.divider()

@st.cache_data(show_spinner="正在运行差异表达分析…")
def cached_de(counts_key, ref, target, batch, fdr, lfc):
    return run_de(
        counts_use,
        metadata,
        group_col="groupA",
        ref_group=ref,
        target_group=target,
        batch_col="batch" if batch else None,
        donor_col="donor" if has_donor else None,
        fdr_cutoff=fdr,
        lfc_cutoff=lfc,
    )

de_results = cached_de(
    str(counts_use.shape), ref_group, target_group,
    use_batch, fdr_cutoff, lfc_cutoff,
)

# ── 第三节 — 结果摘要 ────────────────────────────────────────────────────────
st.subheader("📊 第三节 — 差异表达结果概览")

n_sig  = de_results["significant"].sum()
n_up   = (de_results["direction"] == "Up").sum()
n_down = (de_results["direction"] == "Down").sum()
test_type = de_results["test_type"].iloc[0]

m1, m2, m3, m4, m5 = st.columns(5)
m1.metric("被检验基因数",             f"{len(de_results):,}")
m2.metric("显著差异表达基因数",        f"{n_sig:,}")
m3.metric(f"⬆ 在 {target_group} 中上调", f"{n_up:,}")
m4.metric(f"⬇ 在 {target_group} 中下调", f"{n_down:,}")
m5.metric("使用的检验方法", test_type)

st.divider()

# ── 第四节 — 火山图 ──────────────────────────────────────────────────────────
st.subheader("🌋 第四节 — 火山图")

col_vol, col_info = st.columns([3, 1])

with col_vol:
    de_plot = de_results.reset_index()
    color_map = {
        "Not significant": "#94a3b8",
        "Up":   "#ef4444",
        "Down": "#3b82f6",
    }
    y_cap = min(de_plot["neg_log10_padj"].quantile(0.999) * 1.1, 50)

    fig_vol = px.scatter(
        de_plot, x="log2FC", y="neg_log10_padj",
        color="direction",
        color_discrete_map=color_map,
        hover_name="gene",
        hover_data={
            "log2FC":         ":.3f",
            "pvalue":         ":.2e",
            "padj":           ":.2e",
            "mean_ref":       ":.2f",
            "mean_target":    ":.2f",
            "direction":      False,
            "neg_log10_padj": False,
        },
        labels={
            "log2FC":         f"log₂FC（{target_group} / {ref_group}）",
            "neg_log10_padj": "−log₁₀（校正 p 值）",
        },
        category_orders={"direction": ["Up", "Down", "Not significant"]},
        height=480,
        opacity=0.65,
    )
    fig_vol.update_traces(marker=dict(size=4))
    fig_vol.add_vline(x= lfc_cutoff, line_dash="dot", line_color="#64748b", opacity=0.7)
    fig_vol.add_vline(x=-lfc_cutoff, line_dash="dot", line_color="#64748b", opacity=0.7)
    fig_vol.add_hline(y=-np.log10(fdr_cutoff), line_dash="dot",
                      line_color="#64748b", opacity=0.7)
    fig_vol.update_layout(
        yaxis_range=[0, y_cap],
        legend_title="方向",
        margin=dict(t=10),
    )
    st.plotly_chart(fig_vol, use_container_width=True)

with col_info:
    st.markdown("""
**如何读懂火山图：**

- **X 轴：** log₂ 倍数变化
  - 右侧 = 在目标组中更高
  - 左侧 = 在参考组中更高

- **Y 轴：** −log₁₀（校正 p 值）
  - 越高 = 越显著

- **虚线：** 你设定的阈值

- **右上 / 左上** = 生物学和统计学上都显著

- **底部** = 无论倍数变化大小，均不显著
    """)

st.divider()

# ── 第五节 — MA 图 ────────────────────────────────────────────────────────────
st.subheader("📈 第五节 — MA 图")

st.markdown("""
**MA 图**展示倍数变化 vs 平均表达量。
它可以揭示倍数变化是否在低表达或高表达水平上存在偏差。
理想情况下，非显著基因的点云应在各表达水平上以 log₂FC = 0 为中心。
""")

de_plot["mean_expr"] = (de_plot["mean_ref"] + de_plot["mean_target"]) / 2

fig_ma = px.scatter(
    de_plot, x="mean_expr", y="log2FC",
    color="direction",
    color_discrete_map=color_map,
    hover_name="gene",
    hover_data={"log2FC": ":.3f", "padj": ":.2e", "direction": False},
    labels={
        "mean_expr": "平均 log-CPM",
        "log2FC":    f"log₂FC（{target_group} / {ref_group}）",
    },
    category_orders={"direction": ["Up", "Down", "Not significant"]},
    height=380, opacity=0.55,
)
fig_ma.update_traces(marker=dict(size=3))
fig_ma.add_hline(y=0,           line_dash="solid", line_color="black", opacity=0.3)
fig_ma.add_hline(y= lfc_cutoff, line_dash="dot",   line_color="#64748b", opacity=0.6)
fig_ma.add_hline(y=-lfc_cutoff, line_dash="dot",   line_color="#64748b", opacity=0.6)
fig_ma.update_layout(margin=dict(t=10))
st.plotly_chart(fig_ma, use_container_width=True)
st.caption(
    "最左侧的基因（极低表达）往往具有不稳定的倍数变化。"
    "这是低表达过滤很重要的原因之一。"
)

st.divider()

# ── 第六节 — 结果表格 ────────────────────────────────────────────────────────
st.subheader("📋 第六节 — 结果表格")

col_tab, col_search = st.columns([3, 1])

with col_search:
    st.markdown("**搜索特定基因**")
    gene_query = st.text_input("基因名", placeholder="例如：HLCS")
    top_n = st.slider("显示前 N 个显著基因", 10, 200, 50, 10)

with col_tab:
    if gene_query.strip():
        query = gene_query.strip()
        if query in de_results.index:
            row = de_results.loc[[query]]
            st.success(f"找到：**{query}**")
            st.dataframe(
                row[["mean_ref", "mean_target", "log2FC", "pvalue", "padj",
                     "significant", "direction"]]
                .rename(columns={
                    "mean_ref":    f"均值（{ref_group}）",
                    "mean_target": f"均值（{target_group}）",
                    "log2FC":      "log₂FC",
                    "pvalue":      "p 值",
                    "padj":        "校正 p 值",
                })
                .style.format({
                    f"均值（{ref_group}）":    "{:.2f}",
                    f"均值（{target_group}）": "{:.2f}",
                    "log₂FC":   "{:.3f}",
                    "p 值":     "{:.2e}",
                    "校正 p 值": "{:.2e}",
                }),
                use_container_width=True,
            )
        else:
            st.warning(f"在结果中未找到基因 '{query}'。")
    else:
        sig = de_results[de_results["significant"]].sort_values("padj")
        if len(sig) == 0:
            st.warning("当前阈值下没有显著差异表达基因。请尝试调整参数。")
        else:
            display = sig.head(top_n)[
                ["mean_ref", "mean_target", "log2FC", "pvalue", "padj", "direction"]
            ].rename(columns={
                "mean_ref":    f"均值（{ref_group}）",
                "mean_target": f"均值（{target_group}）",
                "log2FC":      "log₂FC",
                "pvalue":      "p 值",
                "padj":        "校正 p 值",
                "direction":   "方向",
            })
            st.dataframe(
                display.style.format({
                    f"均值（{ref_group}）":    "{:.2f}",
                    f"均值（{target_group}）": "{:.2f}",
                    "log₂FC":      "{:.3f}",
                    "p 值":        "{:.2e}",
                    "校正 p 值":   "{:.2e}",
                }),
                use_container_width=True,
                height=380,
            )

st.download_button(
    "⬇️ 下载完整差异表达结果（CSV）",
    data=de_results.reset_index().to_csv(index=False),
    file_name="de_results_lesson7.csv",
    mime="text/csv",
)

st.divider()

# ── 第七节 — 正确解读差异表达 ───────────────────────────────────────────────
st.subheader("🎓 第七节 — 正确解读差异表达结果")

col_l, col_r = st.columns(2)
with col_l:
    st.markdown("""
**统计显著性 ≠ 生物学重要性**

| 情景 | 解读 |
|---|---|
| 小 FDR + 大 log₂FC | 强有力的候选基因 |
| 小 FDR + 微小 log₂FC | 统计显著但生物学意义可能可忽略 |
| 大 log₂FC + 大 FDR | 有趣但不可靠——常见于低表达基因 |
| FDR 和 FC 都显著 | 最具说服力 |

**FDR 和 log₂FC 回答不同的问题：**
- **FDR** → "我们有多确信这个差异是真实的？"
- **log₂FC** → "差异有多大？"

你需要两者才能做出良好的生物学解读。
    """)

with col_r:
    st.markdown("""
**上游决策影响差异表达结果**

| 决策 | 对差异表达的影响 |
|---|---|
| 过滤阈值 | 改变被检验基因数 → 改变 FDR 校正 |
| 归一化方法 | 改变检验的输入 |
| 批次校正 | 可以减少或夸大表观差异 |
| 配对 vs 非配对检验 | 配对减少供体噪声 → 更高功效 |

这就是为什么本课程的每一课都很重要——它们都影响最终差异表达结果的质量。
    """)

st.divider()

# ── 第八节 — 常见错误 ────────────────────────────────────────────────────────
st.subheader("⚠️ 第八节 — 常见错误")

st.error("**使用未校正的原始 p 值。** 检验数千个基因时，原始 p < 0.05 会产生数百个假阳性。")
st.warning("**忽视低表达问题。** 极低计数时的倍数变化不稳定。先进行过滤。")
st.warning("**忽略批次效应。** 未建模的批次效应会大幅膨胀差异表达基因列表。")
st.info("**将校正后的计数输入 DESeq2。** 批次校正后的数值仅用于可视化——应使用原始计数并在设计公式中包含批次。")
st.info("**把火山图当作最终答案。** 它是摘要，不是结论。对排名靠前的基因应进行功能背景的后续分析。")

st.divider()

# ── 第九节 — 关键总结 ────────────────────────────────────────────────────────
st.subheader("📌 关键总结")

cols = st.columns(5)
msgs = [
    ("🧬", "比较分组而非样本", "差异表达比较分组均值，而非单个测量值。"),
    ("🔢", "两个指标都重要", "FDR 和 log₂FC 结合才能使结果可解读。"),
    ("📊", "火山图是摘要", "图表帮助探索——它们不是统计证明。"),
    ("🔗", "上游决策很重要", "过滤、归一化和批次处理都塑造了差异表达结果。"),
    ("🎯", "设计是关键", "配对设计、良好元数据和平衡批次提高统计功效。"),
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
st.page_link("pages/05_Exploratory_Analysis_PCA.py",
             label="← 第5课：探索性分析", icon="📊")
st.page_link("Home.py", label="返回课程首页 →", icon="🏠")
