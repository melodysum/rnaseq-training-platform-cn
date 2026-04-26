"""
pages/04_FDR.py — 第4课：多重检验与 FDR
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from utils.data_loader import init_session_data, get_paired_columns
from utils.filtering import filter_low_expression
from utils.fdr_demo import (
    simulate_pvalues,
    apply_bh,
    bh_step_table,
    run_paired_de,
)

st.set_page_config(
    page_title="第4课 — FDR",
    page_icon="📐",
    layout="wide",
)

init_session_data()

st.title("📐 第4课 — 多重检验与 FDR")
st.markdown("""
> **学习目标：** 理解 Benjamini-Hochberg 校正的工作原理，
> 以及——至关重要地——你在第3课中的*过滤决策*如何改变 FDR 结果。
""")

counts_raw = st.session_state["counts"]
metadata   = st.session_state["metadata"]

with st.sidebar:
    st.header("🔧 参数设置")

    st.subheader("过滤参数（与第3课联动）")
    prev = st.session_state.get("filter_results", {})
    min_count = st.slider(
        "每样本最低计数",
        0, 100, prev.get("min_count", 10), 1,
    )
    min_samples = st.slider(
        "达到阈值的最少样本数",
        1, len(counts_raw.columns),
        prev.get("min_samples", max(1, len(counts_raw.columns) // 4)), 1,
    )

    st.subheader("差异表达显著性")
    fdr_cutoff = st.slider("FDR 阈值", 0.01, 0.20, 0.05, 0.01)
    lfc_cutoff = st.slider("最小 |log₂FC|", 0.0, 3.0, 1.0, 0.1)

    st.divider()
    st.subheader("BH 模拟面板")
    sim_n_genes  = st.slider("模拟总基因数", 500, 20000, 5000, 500)
    sim_n_de     = st.slider("模拟中的真实差异表达基因数", 0, 500, 100, 10)
    sim_fdr      = st.slider("FDR 阈值（模拟）", 0.01, 0.20, 0.05, 0.01)

counts_filtered = filter_low_expression(counts_raw, min_count, min_samples)
n_before = len(counts_raw)
n_after  = len(counts_filtered)
n_removed_filter = n_before - n_after

donors, ctrl_cols, treat_cols = get_paired_columns(metadata)

de_results = None
if ctrl_cols and treat_cols and n_after > 0:
    de_results = run_paired_de(
        counts_filtered, ctrl_cols, treat_cols, fdr_cutoff, lfc_cutoff
    )

# ── 第一节：多重检验问题 ──────────────────────────────────────────────────────
st.subheader("🎯 多重检验问题")

col_a, col_b = st.columns([3, 2])
with col_a:
    st.markdown(f"""
当你对 **{n_after:,} 个基因**进行差异表达检验时，
你实际上运行了 **{n_after:,} 次统计检验**——每个基因一次。

即使没有任何基因真正存在差异，p 值阈值 0.05 也会*随机*产生大约：

> **{int(n_after * 0.05):,} 个假阳性**（5% × {n_after:,} 次检验）

这就是多重检验问题。检验越多，仅靠偶然性就能预期产生的假阳性就越多。

**通过过滤减少被检验的基因数量，可以直接缓解这一问题。**
    """)

with col_b:
    thresholds = list(range(0, 51, 5))
    n_tested   = [
        len(filter_low_expression(counts_raw, t, min_samples))
        for t in thresholds
    ]
    fig_burden = px.line(
        x=thresholds, y=n_tested,
        labels={"x": "最低计数阈值", "y": "被检验基因数"},
        title="检验负担 vs 过滤阈值",
        markers=True,
    )
    fig_burden.add_vline(
        x=min_count, line_dash="dot", line_color="#ef4444",
        annotation_text=f"当前：{min_count}", annotation_position="top right",
    )
    fig_burden.update_layout(height=300, margin=dict(t=40, b=10))
    st.plotly_chart(fig_burden, use_container_width=True)

st.divider()

# ── 第二节：BH 校正原理 ──────────────────────────────────────────────────────
st.subheader("📊 Benjamini-Hochberg 校正的工作原理")

with st.expander("BH 逐步说明", expanded=True):
    col_exp, col_table = st.columns([2, 3])
    with col_exp:
        st.markdown("""
**Benjamini-Hochberg（BH）** 方法控制*假发现率（FDR）*——
即显著结果中假阳性所占的预期比例。

**算法步骤：**
1. 将所有 p 值从小到大排序（第1名到第 m 名）。
2. 对于每个排名 k，计算 BH 阈值：**(k / m) × α**
3. 找到 p 值 ≤ BH 阈值的**最大排名**。
4. 拒绝该排名及以前的所有假设。

**核心洞见：** 随着 m 的增大，每次单独检验的 BH 阈值变得*更严格*。
检验 10,000 个基因时，排名第 50 的 p 值 0.001 需与
(50/10000) × 0.05 = 0.00025 比较——比只检验 1,000 个基因时严格得多。

这就是**过滤为何重要**：基因越少 → 校正越宽松 → 统计功效越高。
        """)
    with col_table:
        sim_df = simulate_pvalues(sim_n_genes, sim_n_de)
        step_df, critical_rank = bh_step_table(sim_df["pvalue"].values, sim_fdr, show_n=15)
        st.markdown(f"**BH 步骤（前15行，m = {sim_n_genes:,}，α = {sim_fdr}）**")
        st.dataframe(
            step_df.style
            .format({"p-value": "{:.4f}", "BH threshold (k/m × α)": "{:.4f}"})
            .apply(lambda col: [
                "background-color: #dcfce7" if v else "background-color: #fee2e2"
                for v in col
            ], subset=["Rejected?"])
            .set_properties(**{"font-size": "0.82rem"}),
            use_container_width=True,
            hide_index=True,
        )
        st.caption(f"临界排名：**{critical_rank}** — 该排名及以前的所有基因均被拒绝（即显著）。")

st.divider()

# ── 第三节：交互式 BH 模拟 ───────────────────────────────────────────────────
st.subheader("🧪 交互式 FDR 模拟")

sim_df         = simulate_pvalues(sim_n_genes, sim_n_de)
reject, padj   = apply_bh(sim_df["pvalue"].values, sim_fdr)

sim_df["padj"]     = padj
sim_df["rejected"] = reject
sim_df["neg_log10_p"]    = -np.log10(sim_df["pvalue"].clip(lower=1e-300))
sim_df["neg_log10_padj"] = -np.log10(sim_df["padj"].clip(lower=1e-300))

n_sig_sim    = reject.sum()
n_tp_sim     = (reject & sim_df["true_de"]).sum()
n_fp_sim     = (reject & ~sim_df["true_de"]).sum()
fdr_actual   = n_fp_sim / max(n_sig_sim, 1)

m1, m2, m3, m4 = st.columns(4)
m1.metric("被检验基因数",    f"{sim_n_genes:,}")
m2.metric("真实差异表达基因数",   f"{sim_n_de:,}")
m3.metric("BH 显著基因数", f"{n_sig_sim:,}")
m4.metric("实际 FDR",      f"{fdr_actual:.1%}")

col_p, col_padj = st.columns(2)

with col_p:
    fig_p = px.histogram(
        sim_df, x="pvalue", nbins=50,
        color="true_de",
        color_discrete_map={True: "#ef4444", False: "#94a3b8"},
        labels={"pvalue": "原始 p 值", "true_de": "真实差异表达基因"},
        title="原始 p 值分布",
    )
    fig_p.add_vline(x=sim_fdr, line_dash="dot", line_color="black",
                    annotation_text=f"α = {sim_fdr}")
    fig_p.update_layout(height=320, margin=dict(t=40))
    st.plotly_chart(fig_p, use_container_width=True)

with col_padj:
    fig_padj = px.histogram(
        sim_df, x="padj", nbins=50,
        color="true_de",
        color_discrete_map={True: "#ef4444", False: "#94a3b8"},
        labels={"padj": "BH 校正 p 值", "true_de": "真实差异表达基因"},
        title="校正后 p 值分布",
    )
    fig_padj.add_vline(x=sim_fdr, line_dash="dot", line_color="black",
                       annotation_text=f"FDR = {sim_fdr}")
    fig_padj.update_layout(height=320, margin=dict(t=40))
    st.plotly_chart(fig_padj, use_container_width=True)

st.caption(
    "红色 = 真实差异表达基因（ground truth），灰色 = 零假设基因。"
    "注意校正后许多小的原始 p 值如何被推到阈值以上。"
)

st.divider()

# ── 第四节：基因命运——两种丢失方式 ─────────────────────────────────────────
st.subheader("🗺️ 基因消失的两种方式")

fate_col, fate_viz = st.columns([2, 3])

with fate_col:
    st.markdown(f"""
经过所有分析步骤后，一个基因最终落入以下三种结果之一：

| 命运 | 原因 | 阶段 |
|---|---|---|
| ✅ **显著差异表达** | 通过过滤 + FDR | 两步均通过 |
| 🔴 **被过滤去除** | 低表达 | 第3课 |
| 🟡 **检验但不显著** | FDR 未通过 | 第4课 |

当前设置下：
- **{n_removed_filter:,}** 个基因被低表达过滤去除
- **{n_after:,}** 个基因进入统计检验
""")
    if de_results is not None:
        n_sig  = de_results["significant"].sum()
        n_fail = n_after - n_sig
        st.markdown(f"""
- **{n_sig:,}** 个基因显著差异表达
- **{n_fail:,}** 个基因已检验但不显著
        """)

with fate_viz:
    if de_results is not None:
        n_sig  = int(de_results["significant"].sum())
        n_fail = n_after - n_sig

        fig_fate = px.funnel(
            x=[n_before, n_after, n_sig],
            y=["数据集中所有基因",
               "低表达过滤后",
               "显著差异表达基因"],
            color_discrete_sequence=["#64748b", "#3b82f6", "#22c55e"],
        )
        fig_fate.update_layout(height=320, margin=dict(t=10))
        st.plotly_chart(fig_fate, use_container_width=True)

        fig_sankey = go.Figure(go.Sankey(
            node=dict(
                label=["所有基因", "过滤去除",
                       "进入检验", "显著", "不显著"],
                color=["#94a3b8", "#f87171", "#60a5fa", "#4ade80", "#fbbf24"],
            ),
            link=dict(
                source=[0, 0, 2, 2],
                target=[1, 2, 3, 4],
                value=[n_removed_filter, n_after, n_sig, n_fail],
                color=["#fecaca", "#bfdbfe", "#bbf7d0", "#fde68a"],
            ),
        ))
        fig_sankey.update_layout(height=280, margin=dict(t=10, b=10))
        st.plotly_chart(fig_sankey, use_container_width=True)
        st.caption("追踪基因从完整数据集 → 过滤 → 统计检验 → 显著性的流向。")

st.divider()

# ── 第五节：过滤 → FDR 的联系 ───────────────────────────────────────────────
st.subheader("🔗 过滤如何改变 FDR 结果")

st.markdown("""
下图展示了改变过滤阈值如何影响被检验的基因数量，
以及由此带来的显著结果数量变化。
""")

thresholds_scan = list(range(0, 51, 5))
scan_rows = []
for t in thresholds_scan:
    cf = filter_low_expression(counts_raw, t, min_samples)
    n_t = len(cf)
    n_sig_t = 0
    if ctrl_cols and treat_cols and n_t > 0:
        try:
            de_t = run_paired_de(cf, ctrl_cols, treat_cols, fdr_cutoff, lfc_cutoff)
            n_sig_t = int(de_t["significant"].sum())
        except Exception:
            pass
    scan_rows.append({"最低计数": t, "被检验基因数": n_t, "显著差异表达基因数": n_sig_t})

scan_df = pd.DataFrame(scan_rows)

fig_scan = go.Figure()
fig_scan.add_trace(go.Scatter(
    x=scan_df["最低计数"], y=scan_df["被检验基因数"],
    name="被检验基因数", line=dict(color="#3b82f6"), mode="lines+markers",
))
fig_scan.add_trace(go.Scatter(
    x=scan_df["最低计数"], y=scan_df["显著差异表达基因数"],
    name="显著差异表达基因数", line=dict(color="#22c55e"), mode="lines+markers",
    yaxis="y2",
))
fig_scan.add_vline(x=min_count, line_dash="dot", line_color="#ef4444",
                   annotation_text=f"当前：{min_count}")
fig_scan.update_layout(
    xaxis_title="最低计数阈值",
    yaxis=dict(title="被检验基因数", color="#3b82f6"),
    yaxis2=dict(title="显著差异表达基因数", color="#22c55e",
                overlaying="y", side="right"),
    height=380,
    margin=dict(t=20),
    legend=dict(x=0.6, y=0.95),
)
st.plotly_chart(fig_scan, use_container_width=True)
st.caption(
    "随着过滤阈值升高，被检验基因数减少（蓝色），"
    "但 BH 校正变得不那么严格——通常在一定范围内会产生*更多*显著基因（绿色）。"
)

st.divider()

if de_results is not None:
    st.subheader("🌋 火山图 — 真实数据差异表达结果")

    de_plot = de_results.reset_index()
    color_map = {
        "Not significant": "#94a3b8",
        "Up":              "#ef4444",
        "Down":            "#3b82f6",
    }
    fig_vol = px.scatter(
        de_plot, x="log2FC", y="neg_log10_padj",
        color="direction",
        color_discrete_map=color_map,
        hover_name="gene",
        hover_data={"log2FC": ":.3f", "padj": ":.2e", "direction": False,
                    "neg_log10_padj": False},
        labels={"log2FC": "log₂FC（处理组 / 对照组）",
                "neg_log10_padj": "−log₁₀（校正 p 值）"},
        category_orders={"direction": ["Up", "Down", "Not significant"]},
        height=480, opacity=0.65,
    )
    fig_vol.update_traces(marker=dict(size=4))
    fig_vol.add_vline(x= lfc_cutoff, line_dash="dot", line_color="#64748b", opacity=0.7)
    fig_vol.add_vline(x=-lfc_cutoff, line_dash="dot", line_color="#64748b", opacity=0.7)
    fig_vol.add_hline(y=-np.log10(fdr_cutoff), line_dash="dot",
                      line_color="#64748b", opacity=0.7)
    fig_vol.update_layout(margin=dict(t=10))
    st.plotly_chart(fig_vol, use_container_width=True)

    st.subheader("📋 显著差异表达基因前30名")
    sig = de_results[de_results["significant"]].sort_values("padj")
    if len(sig) > 0:
        st.dataframe(
            sig.head(30)[["log2FC", "mean_expr", "pvalue", "padj", "direction"]]
            .rename(columns={"log2FC": "log₂FC", "mean_expr": "平均 log-CPM",
                             "pvalue": "p 值", "padj": "校正 p 值",
                             "direction": "方向"})
            .style.format({"log₂FC": "{:.3f}", "平均 log-CPM": "{:.2f}",
                            "p 值": "{:.2e}", "校正 p 值": "{:.2e}"}),
            use_container_width=True, height=380,
        )
        st.download_button(
            "⬇️ 下载完整差异表达结果（CSV）",
            data=de_results.reset_index().to_csv(index=False),
            file_name="de_results.csv", mime="text/csv",
        )
    else:
        st.warning("当前阈值下没有显著差异表达基因。请尝试调整侧边栏参数。")

st.divider()
st.info("""
**本课关键收获：**

过滤和 FDR 校正不是独立的步骤——它们直接相连。
被检验基因越少 → BH 校正越宽松 → 检测真实差异的统计功效越高。
两种基因丢失方式（过滤去除 vs FDR 不通过）在生物学和统计学上截然不同，
你应该始终在报告中说明各自的数量。
""")
st.page_link("Home.py", label="← 返回课程首页", icon="🏠")
