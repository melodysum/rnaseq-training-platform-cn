"""
pages/03_Filtering.py — 第3课：低表达基因过滤
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from utils.data_loader import init_session_data
from utils.filtering import (
    filter_low_expression,
    expression_summary,
    check_gene_fate,
    log_cpm,
)

st.set_page_config(
    page_title="第3课 — 过滤",
    page_icon="🔍",
    layout="wide",
)

init_session_data()

st.title("🔍 第3课 — 低表达基因过滤")
st.markdown("""
> **学习目标：** 理解*为什么*在分析前需要去除低计数基因，
> 以及阈值的选择如何影响哪些基因被保留——哪些被丢弃。
""")

counts_raw = st.session_state["counts"]
metadata   = st.session_state["metadata"]
N_GENES    = len(counts_raw)
N_SAMPLES  = len(counts_raw.columns)

with st.sidebar:
    st.header("🔧 过滤参数")
    min_count = st.slider(
        "每个样本的最低计数",
        0, 100, 10, 1,
        help="一个基因在足够多的样本中至少需要这么多 reads 才能被保留。",
    )
    min_samples = st.slider(
        "达到阈值的最少样本数",
        1, N_SAMPLES, max(1, N_SAMPLES // 4), 1,
        help="有多少个样本必须超过计数阈值。",
    )
    st.divider()
    st.caption(
        "💡 试着将最低计数设为 0（保留所有基因），"
        "然后逐渐提高——观察有多少基因消失。"
    )

counts_filtered = filter_low_expression(counts_raw, min_count, min_samples)
n_kept    = len(counts_filtered)
n_removed = N_GENES - n_kept
pct_kept  = n_kept / N_GENES * 100

st.session_state["filter_results"] = {
    "counts_filtered": counts_filtered,
    "min_count":       min_count,
    "min_samples":     min_samples,
    "n_kept":          n_kept,
    "n_removed":       n_removed,
}

c1, c2, c3, c4 = st.columns(4)
c1.metric("总基因数",      f"{N_GENES:,}")
c2.metric("保留基因数",   f"{n_kept:,}",    delta=f"{pct_kept:.1f}% 保留")
c3.metric("去除基因数",    f"{n_removed:,}", delta=f"-{100-pct_kept:.1f}%", delta_color="inverse")
c4.metric("样本数",          f"{N_SAMPLES}")

st.divider()

with st.expander("📖 为什么要过滤低表达基因？", expanded=True):
    col_a, col_b = st.columns(2)
    with col_a:
        st.markdown("""
**极低计数的问题**

一个在各样本中计数为 0、0、1、0、2 的基因，
给不出可靠的生物学信号——它基本上是测序过程中的噪声。

这类基因会造成两个问题：
1. **膨胀方差**——倍数变化变得巨大且毫无意义（例如从 0 到 1 在技术上是无限大）。
2. **增加检验负担**——每个被检验的基因都是一次产生假阳性的机会。
        """)
    with col_b:
        st.markdown("""
**过滤的作用**

过滤去除那些表达量过低、无法可靠测量的基因。

- 这**不是任意删除**——而是一个有据可依的预处理步骤。
- 阈值由你决定：过严会丢失真实生物学信号，过宽则保留了噪声。
- 常用规则：在至少**与最小组大小相同的样本数**中，保留计数 **≥ 10** 的基因。

过滤后，每个剩余基因都有足够的数据可以被公平检验。
        """)

st.divider()

st.subheader("📊 表达分布：过滤前 vs 过滤后")

tab1, tab2 = st.tabs(["直方图", "每样本箱线图"])

with tab1:
    rng = np.random.default_rng(42)
    n_plot = min(3000, N_GENES)
    idx_b  = rng.choice(N_GENES,  size=n_plot, replace=False)
    idx_a  = rng.choice(n_kept,   size=min(n_plot, n_kept), replace=False)

    vals_before = np.log2(counts_raw.iloc[idx_b].values.flatten() + 1)
    vals_after  = np.log2(counts_filtered.iloc[idx_a].values.flatten() + 1)

    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=vals_before, name="过滤前",
        opacity=0.55, nbinsx=80, marker_color="#64748b",
    ))
    fig.add_trace(go.Histogram(
        x=vals_after, name="过滤后",
        opacity=0.65, nbinsx=80, marker_color="#3b82f6",
    ))
    fig.update_layout(
        barmode="overlay",
        xaxis_title="log₂(计数 + 1)",
        yaxis_title="观测数",
        height=360,
        margin=dict(t=10),
        legend=dict(x=0.7, y=0.95),
    )
    st.plotly_chart(fig, use_container_width=True)
    st.caption(
        "左侧接近零的峰代表低/零计数基因。"
        "过滤去除了这个峰，留下更干净的分布。"
    )

with tab2:
    sample_subset = list(counts_raw.columns[:10])
    lc_before = log_cpm(counts_raw)[sample_subset]
    lc_after  = log_cpm(counts_filtered)[sample_subset] if n_kept > 0 else lc_before

    fig2 = go.Figure()
    for col in sample_subset:
        grp = metadata.loc[col, "groupA"] if col in metadata.index else "unknown"
        color = "#3b82f6" if grp == "control" else "#ef4444"
        fig2.add_trace(go.Box(
            y=lc_before[col], name=col, marker_color=color,
            opacity=0.4, showlegend=False,
        ))
        fig2.add_trace(go.Box(
            y=lc_after[col], name=col + "（过滤后）",
            marker_color=color, opacity=0.9, showlegend=False,
        ))
    fig2.update_layout(
        yaxis_title="log₂(CPM + 1)",
        height=380, margin=dict(t=10),
    )
    st.plotly_chart(fig2, use_container_width=True)
    st.caption("颜色较深 = 过滤后。去除噪声低计数基因后，分布变得更紧致。")

st.divider()

st.subheader("🗺️ 被去除的基因去哪了？")

col_l, col_r = st.columns(2)

with col_l:
    st.markdown("""
基因可能在**两个不同阶段**离开你的分析：

| 去除阶段 | 原因 | 时机 |
|---|---|---|
| 🔴 **低表达过滤** | 在样本间计数过低 | 现在，本步骤 |
| 🟡 **FDR 不显著** | 已检验但未通过校正 | 后续，在差异表达步骤 |

理解这一区别很重要：
- 因过滤而丢失的基因**从未被统计检验过**。
- 因 FDR 而丢失的基因**被检验过**，但没有显示出显著差异。
- 这在生物学和统计学上是完全不同的结果。
    """)

with col_r:
    if n_kept > 0:
        fig_pie = px.pie(
            values=[n_kept, n_removed],
            names=["保留用于检验", "去除（低表达）"],
            color_discrete_sequence=["#3b82f6", "#f87171"],
            hole=0.45,
        )
        fig_pie.update_layout(height=280, margin=dict(t=10, b=10))
        st.plotly_chart(fig_pie, use_container_width=True)

st.divider()

st.subheader("🔬 追踪特定基因")

col_track, col_gt = st.columns(2)

with col_track:
    st.markdown("**搜索任意基因**")
    gene_input = st.text_area(
        "输入基因名（每行一个或用逗号分隔）",
        placeholder="例如：\nHLCS\nNKX2-1\nSNU13",
        height=120,
        key="gene_search",
    )
    if gene_input.strip():
        query_genes = [
            g.strip()
            for g in gene_input.replace(",", "\n").split("\n")
            if g.strip()
        ]
        fate_df = check_gene_fate(counts_raw, counts_filtered, query_genes)
        st.dataframe(
            fate_df.style.apply(
                lambda col: [
                    "background-color: #dcfce7" if "Retained" in v
                    else "background-color: #fee2e2" if "Removed" in v
                    else "background-color: #fef9c3"
                    for v in col
                ],
                subset=["Status"],
            ),
            use_container_width=True,
            hide_index=True,
        )

with col_gt:
    st.markdown("**已知生物学相关基因**")
    st.caption(
        "粘贴一份你预期具有生物学意义的基因列表。"
        "这样你可以检查过滤阈值是否意外地去除了已知的重要基因。"
    )
    gt_input = st.text_area(
        "粘贴已知基因列表",
        placeholder="例如：\nHLCS\nTBC1D13\nH2AZ2\nCLNS1A",
        height=120,
        key="gt_genes",
    )
    if gt_input.strip():
        gt_genes = [
            g.strip()
            for g in gt_input.replace(",", "\n").split("\n")
            if g.strip()
        ]
        gt_fate = check_gene_fate(counts_raw, counts_filtered, gt_genes)

        retained = (gt_fate["Status"].str.contains("Retained")).sum()
        removed  = (gt_fate["Status"].str.contains("Removed")).sum()
        notfound = (gt_fate["Status"].str.contains("Not found")).sum()

        if removed > 0:
            st.warning(
                f"⚠️ 你的已知基因中有 {removed} 个正在被当前阈值**去除**。"
                f"考虑降低过滤条件。"
            )
        if retained > 0:
            st.success(f"✅ {retained} 个已知基因已保留。")

        st.dataframe(
            gt_fate.style.apply(
                lambda col: [
                    "background-color: #dcfce7" if "Retained" in v
                    else "background-color: #fee2e2" if "Removed" in v
                    else "background-color: #fef9c3"
                    for v in col
                ],
                subset=["Status"],
            ),
            use_container_width=True,
            hide_index=True,
        )

st.divider()

# ══════════════════════════════════════════════════════════════════════════════
# 第5节 — 数据驱动的阈值选择（filterByExpr 方法）
# ══════════════════════════════════════════════════════════════════════════════

from utils.filtering import filter_by_expr, threshold_sweep_retained

st.subheader("📐 第5节 — 数据驱动的阈值选择（filterByExpr 方法）")
st.markdown("""
上面的手动阈值适合学习演示，但在真实分析中，阈值应由数据本身决定——而非凭经验随意指定。

**为什么任意阈值会有问题：**  
无论测序深度如何，一律使用 `min_count = 10` 是不一致的。  
在 100 万 reads 的文库中检测到 10 条 reads 代表良好表达；  
而在 1000 万 reads 的文库中，同样的 10 条 reads 则代表极低表达。

**filterByExpr 逻辑（edgeR，Chen et al. 2016）：**  
阈值随文库大小动态调整：

```
CPM 截止值 = min_count / (中位文库大小 / 1,000,000)
```

只有在至少 `min_group_size` 个样本中 CPM 超过该截止值的基因才会被保留。
""")

with st.expander("⚠️ 为什么不能用 DE 结果来优化阈值？", expanded=False):
    st.markdown("""
一个常见错误是：测试多个阈值 → 在每个阈值下跑 DE → 选择 DE 显著基因最多的那个阈值。

**这是循环论证：**  
你在用答案（DE 结果）来决定输入（哪些基因参与检验）。  
这会系统性地膨胀假阳性率，因为你实际上是在针对数据的噪声结构调参。

**正确做法：** 仅根据文库大小和分组大小确定阈值——与 DE 结果完全无关。
""")


with st.expander("filterByExpr 背后的统计原理", expanded=False):
    st.markdown(
        "filterByExpr 只保留表达量达到可进行统计检验水平的基因。"
        "它从两个维度评判每个基因：(1) 基于中位文库大小动态缩放的 CPM 截止值；"
        "(2) 样本支持数：基因必须在至少 min_group_size 个样本中超过 CPM 截止值。"
        "简单计数过滤不考虑文库大小和实验设计。"
        "低计数基因会在负二项模型中导致方差膨胀，使 logFC 被放大、p 值不可靠。"
        "过滤还能减少多重检验负担：检验的基因越少，FDR 校正惩罚越低，剩余基因的统计功效越高。"
    )

col_fe1, col_fe2 = st.columns(2)
with col_fe1:
    fe_min_count = st.slider(
        "filterByExpr min_count",
        1, 50, 10, 1,
        help="中位文库大小下的最低计数。edgeR 默认推荐值为 10。",
        key="fe_min_count",
    )
    fe_min_total = st.slider(
        "filterByExpr min_total_count",
        5, 50, 15, 1,
        help="所有样本合计的最低计数。",
        key="fe_min_total",
    )

counts_fe = filter_by_expr(
    counts_raw, metadata,
    group_col="groupA",
    min_count=fe_min_count,
    min_total_count=fe_min_total,
)

lib_sizes_fe = counts_raw.sum(axis=0)
median_lib_fe = lib_sizes_fe.median()
cpm_cut_fe = fe_min_count / (median_lib_fe / 1e6)

with col_fe2:
    st.metric("filterByExpr 保留基因数", f"{len(counts_fe):,}")
    st.metric("已移除基因数", f"{N_GENES - len(counts_fe):,}")
    st.caption(
        f"推导的 CPM 截止值：**{cpm_cut_fe:.3f}** "
        f"（来自中位文库大小 {median_lib_fe:,.0f}）"
    )

st.markdown("#### 阈值敏感性分析：不同阈值下保留基因数")
st.caption(
    "该图展示各阈值下存活的基因数量——**不涉及 DE 分析**。"
    "用于判断基因保留数量对阈值选择的敏感程度，不用于最大化 DE 信号。"
)

sweep_df = threshold_sweep_retained(counts_raw, metadata, group_col="groupA")

fig_sweep = go.Figure()
fig_sweep.add_trace(go.Scatter(
    x=sweep_df["threshold"],
    y=sweep_df["n_retained"],
    mode="lines+markers",
    marker=dict(size=8, color="#2563eb"),
    line=dict(color="#2563eb", width=2),
    text=[f"CPM 截止: {r:.3f}<br>保留: {n:,} ({p}%)"
          for r, n, p in zip(sweep_df["cpm_cutoff"],
                             sweep_df["n_retained"],
                             sweep_df["pct_retained"])],
    hovertemplate="%{text}<extra></extra>",
))
rec_row = sweep_df.iloc[(sweep_df["threshold"] - fe_min_count).abs().argsort()[:1]]
fig_sweep.add_trace(go.Scatter(
    x=rec_row["threshold"],
    y=rec_row["n_retained"],
    mode="markers",
    marker=dict(size=14, color="red", symbol="diamond"),
    name=f"当前阈值 ({fe_min_count})",
))
fig_sweep.update_layout(
    xaxis_title="min_count 阈值",
    yaxis_title="保留基因数",
    height=350,
    legend=dict(orientation="h", y=-0.25),
)
st.plotly_chart(fig_sweep, use_container_width=True)
st.dataframe(sweep_df, use_container_width=True, hide_index=True)

st.divider()

st.info("""
**本课关键收获：**

过滤是一个预处理决策，而非分析步骤。它决定了哪些基因将*进入*你的统计检验。
你在这里选择的阈值将直接影响被检验的基因数量——进而改变多重检验的负担，
最终影响你的 FDR 结果。

👉 **前往第4课 — FDR**，看看这一切是如何展开的。
""")

st.page_link("pages/04_FDR.py", label="继续第4课 → FDR 与多重检验", icon="📐")
