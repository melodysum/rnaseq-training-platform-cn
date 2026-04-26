"""
pages/01_RNAseq_Foundations.py — 第1课：RNA-seq 基础与实验设计
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import streamlit as st

st.set_page_config(
    page_title="第1课 — RNA-seq 基础",
    page_icon="🔬",
    layout="wide",
)

st.title("🔬 第1课 — RNA-seq 基础与实验设计")
st.markdown("""
> **学习目标：** 理解 RNA-seq 数据的来源、核心测序概念的含义，以及为什么良好的实验设计必须在任何分析开始*之前*就确定下来。
""")

# ── 第一节 — 为什么基础很重要 ────────────────────────────────────────────────
st.subheader("🏗️ 第一节 — 为什么基础很重要")

c1, c2, c3 = st.columns(3)
for col, icon, title, body in [
    (c1, "📐", "设计在先",
     "RNA-seq 分析不从统计开始——它从样本的设计、采集和测序方式开始。"),
    (c2, "🔗", "环环相扣",
     "对 FASTQ 质量、重复数或测序深度的误解，会在后续的过滤、PCA 和差异表达分析中引发连锁问题。"),
    (c3, "⚠️", "亡羊补牢有限",
     "糟糕的实验设计决策往往无法通过计算手段完全弥补，无论下游方法多么复杂。"),
]:
    col.markdown(f"""
<div style="background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;
padding:1rem;height:140px;">
<div style="font-size:1.4rem">{icon}</div>
<div style="font-weight:600;margin:0.3rem 0">{title}</div>
<div style="font-size:0.85rem;color:#475569">{body}</div>
</div>
""", unsafe_allow_html=True)

st.divider()

# ── 第二节 — 什么是 RNA-seq？ ─────────────────────────────────────────────────
st.subheader("🧬 第二节 — 什么是 RNA-seq？")

st.markdown("""
**RNA-seq**（RNA 测序）用于测量生物样本中 RNA 转录本的丰度，
提供全基因组范围内基因活跃程度的快照。

- **Bulk RNA-seq** 测量数千乃至数百万细胞的*平均*表达水平
- reads 比对到参考基因组后进行计数，生成**表达矩阵**
- 表达矩阵（基因 × 样本）是所有下游分析的起点
""")

st.markdown("**RNA-seq 工作流程：**")
stages = [
    ("🧫", "生物\n样本", "组织、细胞或类器官"),
    ("🔬", "建库\n制备", "RNA 提取、片段化、cDNA 合成、接头连接"),
    ("💻", "测序", "生成短 reads（FASTQ 文件）"),
    ("📁", "FASTQ\n文件", "带质量分数的原始测序数据"),
    ("🗺️", "定量", "reads 比对并按基因/转录本计数"),
    ("📊", "计数\n矩阵", "基因 × 样本的原始计数表"),
    ("🔍", "下游\n分析", "过滤、归一化、差异表达、富集分析"),
]

fig_flow = go.Figure()
n = len(stages)
for i, (icon, label, tooltip) in enumerate(stages):
    x = i / (n - 1)
    fig_flow.add_trace(go.Scatter(
        x=[x], y=[0.5],
        mode="markers+text",
        marker=dict(size=40, color="#3b82f6", symbol="circle"),
        text=[icon],
        textposition="middle center",
        textfont=dict(size=18),
        hovertext=tooltip,
        hoverinfo="text",
        showlegend=False,
    ))
    fig_flow.add_annotation(
        x=x, y=0.18, text=label.replace("\n", "<br>"),
        showarrow=False, font=dict(size=10, color="#334155"),
        align="center",
    )
    if i < n - 1:
        fig_flow.add_annotation(
            x=(i + 0.5) / (n - 1), y=0.5,
            text="→", showarrow=False,
            font=dict(size=20, color="#94a3b8"),
        )

fig_flow.update_layout(
    height=180, margin=dict(l=10, r=10, t=10, b=40),
    xaxis=dict(visible=False, range=[-0.05, 1.05]),
    yaxis=dict(visible=False, range=[0, 1]),
    plot_bgcolor="white",
)
st.plotly_chart(fig_flow, use_container_width=True)
st.caption("将鼠标悬停在每个阶段上可查看更多详情。")

st.divider()

# ── 第三节 — FASTQ 文件结构 ────────────────────────────────────────────────────
st.subheader("📄 第三节 — FASTQ 文件结构")

col_fastq, col_qual = st.columns([3, 2])

with col_fastq:
    st.markdown("""
每条测序 read 以 **4 行 FASTQ 格式** 存储：
    """)
    st.code("""@SRR123456.1 read_1 length=150          ← Read 标识符
ATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATG  ← 核苷酸序列
+                                                  ← 分隔符（始终为 +）
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII!!!!        ← 质量字符串（Phred）
""", language="text")

    with st.expander("🔍 质量字符代表什么？"):
        st.markdown("""
每个字符编码一个 **Phred 质量分数**（Q）：

| Q 分数 | 错误概率 | ASCII 字符（示例） |
|---|---|---|
| Q10 | 1/10 | + |
| Q20 | 1/100 | 5 |
| Q30 | 1/1,000 | ? |
| Q40 | 1/10,000 | I |

- 质量编码方式：`Q = ASCII值 - 33`
- `!` = Q0（最差），`I` = Q40（优秀）
- 质量分数通常**在 reads 末端下降**，原因是测序化学反应的累积误差
        """)

with col_qual:
    st.markdown("**质量分数交互探索**")
    q_score = st.slider("Phred Q 分数", min_value=0, max_value=40, value=30)
    error_prob = 10 ** (-q_score / 10)
    accuracy   = (1 - error_prob) * 100

    st.metric("错误概率", f"1/{int(1/error_prob):,}" if error_prob > 0 else "< 1/10,000")
    st.metric("碱基识别准确率", f"{accuracy:.4f}%")

    q_vals = list(range(0, 41))
    probs  = [10 ** (-q / 10) for q in q_vals]
    fig_q = px.line(
        x=q_vals, y=probs,
        labels={"x": "Q 分数", "y": "错误概率"},
        title="Q 分数 → 错误概率",
        log_y=True,
    )
    fig_q.add_vline(x=q_score, line_dash="dot", line_color="#ef4444",
                    annotation_text=f"Q{q_score}", annotation_position="top right")
    fig_q.update_layout(height=260, margin=dict(t=40, b=10))
    st.plotly_chart(fig_q, use_container_width=True)

st.divider()

# ── 第四节 — 核心测序概念 ────────────────────────────────────────────────────
st.subheader("📚 第四节 — 核心测序概念")

concepts = {
    "读长（Read length）": {
        "definition": "每条 read 测序的碱基对数量（如 75 bp、150 bp）。",
        "why": "更长的 reads 可提高比对置信度和剪接位点检测能力。",
        "misconception": "更长并不总是更好——对于基因水平的差异表达分析，150 bp 通常已足够。",
    },
    "测序深度（Sequencing depth）": {
        "definition": "每个样本生成的 reads 总数（如 2000 万、5000 万 reads）。",
        "why": "更大的深度可改善低表达基因的检测，并降低抽样噪声。",
        "misconception": "更高的深度无法替代生物学重复——对 2 个样本进行更深的测序，你仍然只有 n=2。",
    },
    "文库（Library）": {
        "definition": "从一个样本制备的 DNA 片段集合，准备好用于测序。",
        "why": "文库质量（RNA 完整性、接头连接）直接影响数据质量。",
        "misconception": "技术重复（从同一 RNA 制备的两个文库）不等同于生物学重复。",
    },
    "单端 vs 双端测序": {
        "definition": "单端测序仅读取每个片段的一端；双端测序读取两端。",
        "why": "双端测序可提高比对精度，更适合转录本亚型分析。",
        "misconception": "双端测序不会使生物学重复数翻倍——它只是为每个片段提供了两条 reads。",
    },
    "生物学重复（Biological replicates）": {
        "definition": "来自相同条件的独立生物样本（如同一组的 4 只独立小鼠）。",
        "why": "对于估计组内变异和任何有效统计检验来说都是必不可少的。",
        "misconception": "每组 n=1 或 n=2 不足以进行可靠的差异表达分析——无法估计变异性。",
    },
    "技术重复（Technical replicates）": {
        "definition": "对同一生物样本的重复测量（如对同一 RNA 测序两次）。",
        "why": "用于评估平台可重复性，但不应替代生物学重复。",
        "misconception": "技术重复无法反映个体之间的生物学变异。",
    },
}

tabs = st.tabs(list(concepts.keys()))
for tab, (concept, info) in zip(tabs, concepts.items()):
    with tab:
        col_def, col_misc = tab.columns(2)
        with col_def:
            st.markdown(f"**定义：** {info['definition']}")
            st.markdown(f"**为什么重要：** {info['why']}")
        with col_misc:
            st.warning(f"⚠️ **常见误区：**\n\n{info['misconception']}")

st.divider()

# ── 第五节 — 实验设计模拟器 ──────────────────────────────────────────────────
st.subheader("🧪 第五节 — 实验设计模拟器")

st.markdown("探索设计选择如何影响研究的统计功效和成本。")

col_ctrl, col_output = st.columns([1, 2])

with col_ctrl:
    n_groups     = st.slider("生物学分组数",  2, 6, 2)
    n_reps       = st.slider("每组生物学重复数", 2, 10, 4)
    reads_per_m  = st.slider("每样本测序深度（百万 reads）", 5, 100, 20)
    cost_per_m   = st.number_input("每百万 reads 估算成本（£）", 1.0, 20.0, 4.0, 0.5)
    has_batch    = st.checkbox("添加批次结构（2 个批次）")

n_samples     = n_groups * n_reps
total_reads_m = n_samples * reads_per_m
total_cost    = total_reads_m * cost_per_m

if has_batch:
    samples_b1 = n_samples // 2
    samples_b2 = n_samples - samples_b1
    batch_balanced = (samples_b1 == samples_b2) and (n_reps % 2 == 0)

with col_output:
    m1, m2, m3 = st.columns(3)
    m1.metric("总样本数",    f"{n_samples}")
    m2.metric("总 reads（百万）",  f"{total_reads_m:,}")
    m3.metric("估算成本",        f"£{total_cost:,.0f}")

    issues = []
    if n_reps < 3:
        issues.append("⚠️ **重复数过少。** 每组 n < 3 会严重限制统计功效和方差估计。")
    if reads_per_m < 10:
        issues.append("⚠️ **深度不足。** < 1000 万 reads 可能漏检低表达基因。")
    if reads_per_m > 60 and n_reps < 4:
        issues.append("⚠️ **高深度但低重复。** 考虑将测序预算重新分配给更多重复。")
    if has_batch and not batch_balanced:
        issues.append("⚠️ **批次结构不平衡。** 部分分组可能与批次存在混淆。")
    if not issues:
        st.success("✅ **设计良好。** 重复充足、深度合理、结构平衡。")
    for issue in issues:
        st.warning(issue)

    rep_range   = list(range(2, 11))
    depth_range = [max(5, int(total_reads_m / (n_groups * r))) for r in rep_range]
    fig_td = go.Figure()
    fig_td.add_trace(go.Scatter(
        x=rep_range, y=depth_range,
        mode="lines+markers", name="每样本深度",
        line=dict(color="#3b82f6"),
        marker=dict(size=8),
    ))
    fig_td.add_vline(x=n_reps, line_dash="dot", line_color="#ef4444",
                     annotation_text=f"当前：{n_reps} 重复", annotation_position="top right")
    fig_td.add_hline(y=20, line_dash="dot", line_color="#22c55e",
                     annotation_text="推荐 2000 万 reads", annotation_position="right")
    fig_td.update_layout(
        xaxis_title="每组重复数",
        yaxis_title="每样本 reads 数（百万）",
        title="重复数与测序深度的权衡（预算固定）",
        height=300, margin=dict(t=40),
    )
    st.plotly_chart(fig_td, use_container_width=True)

st.divider()

# ── 第六节 — 重复 vs 深度 ────────────────────────────────────────────────────
st.subheader("⚖️ 第六节 — 重复数与测序深度")

col_a, col_b = st.columns(2)
with col_a:
    st.markdown("""
**方案 A — 高深度，少重复**
- 每组 2 个重复
- 每样本 5000 万 reads
- 总成本：与方案 B 相近

❌ 问题：
- 无法可靠估计组内变异
- 统计检验功效极低
- 一个离群样本就会破坏整个比较
- 无法判断差异是真实的还是噪声
    """)
    st.error("**不推荐用于差异表达分析**")

with col_b:
    st.markdown("""
**方案 B — 均衡重复**
- 每组 4–5 个重复
- 每样本 1500–2000 万 reads
- 总成本相同或更低

✅ 优势：
- 每个基因的方差估计稳健
- 差异表达统计功效更高
- 对离群样本更具抵抗力
- 适用于大多数标准差异表达流程
    """)
    st.success("**大多数 RNA-seq 研究的推荐方案**")

st.info("""
**核心洞见：** 多增加一个生物学重复，通常比将测序深度翻倍更能改善差异表达分析结果。
深度对于检测稀有转录本有意义，但重复数才是统计推断的根本。
""")

st.divider()

# ── 第七节 — 常见设计错误 ────────────────────────────────────────────────────
st.subheader("⚠️ 第七节 — 常见设计错误")

failures = [
    ("🚫", "无生物学重复",
     "每组只有一个样本，在统计上根本无法估计方差或进行有效检验。",
     "error"),
    ("🚫", "没有合适的对照组",
     "没有适当的未处理或基线对照，倍数变化就没有参考基准。",
     "error"),
    ("⚠️", "批次与处理条件混淆",
     "所有处理样本在批次1，所有对照在批次2——计算校正无法解开这种混淆。",
     "warning"),
    ("⚠️", "元数据缺失或不完整",
     "测序后无法分辨哪个样本是哪个，数据将无从解读。",
     "warning"),
    ("ℹ️", "仅用 RNA-seq 验证 RNA-seq",
     "显著的差异表达结果最好通过正交方法验证（如 qPCR、蛋白质检测）。",
     "info"),
    ("ℹ️", "将技术重复当生物学重复",
     "这会人为夸大表观 n 值，产生虚假的统计把握度——它们并非独立观测。",
     "info"),
]

for icon, title, body, level in failures:
    getattr(st, level)(f"**{icon} {title}**\n\n{body}")

st.divider()

# ── 第八节 — 关键总结 ─────────────────────────────────────────────────────────
st.subheader("📌 关键总结")

cols = st.columns(5)
msgs = [
    ("🧬", "RNA-seq 测量转录本", "提供基因活跃程度的全基因组快照。"),
    ("📄", "FASTQ = 序列 + 质量", "质量分数指导过滤和修剪决策。"),
    ("👥", "重复是必须的", "没有生物学重复，统计推断就无效。"),
    ("⚖️", "重复 > 深度", "适度深度下增加重复，通常优于少重复下加大深度。"),
    ("🏗️", "设计无法事后弥补", "良好的实验设计比任何计算校正都更有力。"),
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
col_n1, col_n2 = st.columns(2)
with col_n1:
    st.page_link("Home.py", label="← 返回课程首页", icon="🏠")
with col_n2:
    st.page_link("pages/02_Quantification_Import_Annotation.py",
                 label="第2课：定量与注释 →", icon="📦")
