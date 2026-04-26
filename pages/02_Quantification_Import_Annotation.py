"""
pages/02_Quantification_Import_Annotation.py — 第2课：定量、导入与注释
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from utils.data_loader import init_session_data

st.set_page_config(
    page_title="第2课 — 定量与注释",
    page_icon="📦",
    layout="wide",
)

init_session_data()

st.title("📦 第2课 — 定量、导入与注释")
st.markdown("""
> **学习目标：** 理解测序 reads 是如何变成表达矩阵的，
> 以及为什么元数据和注释在任何下游分析之前都是不可或缺的。
""")

# ── 第一节 — 从 FASTQ 到计数矩阵 ────────────────────────────────────────────
st.subheader("🗺️ 第一节 — 从 FASTQ 到表达矩阵")

stages = [
    ("📁", "FASTQ", "来自测序仪的原始 reads 及质量分数"),
    ("📚", "参考\n转录组", "该物种的已知转录本序列"),
    ("🔢", "定量", "将 reads 比对到转录本；估计丰度"),
    ("🧮", "转录本\n丰度", "每个转录本的计数、TPM 及有效长度"),
    ("🔗", "tx2gene\n汇总", "将转录本聚合为基因水平数值"),
    ("📊", "基因水平\n矩阵", "基因 × 样本计数矩阵——可供分析"),
]

fig = go.Figure()
n = len(stages)
for i, (icon, label, tip) in enumerate(stages):
    x = i / (n - 1)
    fig.add_trace(go.Scatter(
        x=[x], y=[0.55], mode="markers+text",
        marker=dict(size=44, color="#6366f1", symbol="circle"),
        text=[icon], textposition="middle center",
        textfont=dict(size=18),
        hovertext=tip, hoverinfo="text", showlegend=False,
    ))
    fig.add_annotation(
        x=x, y=0.15, text=label.replace("\n", "<br>"),
        showarrow=False, font=dict(size=10, color="#334155"), align="center",
    )
    if i < n - 1:
        fig.add_annotation(
            x=(i + 0.5) / (n - 1), y=0.55, text="→",
            showarrow=False, font=dict(size=20, color="#94a3b8"),
        )
fig.update_layout(
    height=170, margin=dict(l=10, r=10, t=10, b=50),
    xaxis=dict(visible=False, range=[-0.05, 1.05]),
    yaxis=dict(visible=False, range=[0, 1]),
    plot_bgcolor="white",
)
st.plotly_chart(fig, use_container_width=True)
st.caption("将鼠标悬停在每个阶段上可查看详情。")

st.divider()

# ── 第二节 — 传统比对 vs 伪比对 ──────────────────────────────────────────────
st.subheader("⚡ 第二节 — 传统比对 vs 伪比对")

col_al, col_ps = st.columns(2)

with col_al:
    st.markdown("""
**传统比对** *（如 STAR、HISAT2）*

将每条 read 逐碱基比对到参考基因组。

✅ 精确的位置信息  
✅ 可检测新型剪接位点  
✅ 基因组水平输出（适用于变异检测）  
❌ 速度慢，内存占用高  
❌ 需要完整的基因组参考序列  
    """)

with col_ps:
    st.markdown("""
**伪比对** *（如 kallisto、salmon）*

判断每条 read 与哪些转录本*兼容*——无需碱基级别比对。

✅ 非常快（快 10–100 倍）  
✅ 内存需求低  
✅ 内置转录本水平估计  
❌ 无精确基因组位置  
❌ 需要转录组参考序列  
    """)

comparison_df = pd.DataFrame({
    "特性":             ["速度", "内存", "输出水平", "新型剪接检测", "最适用于"],
    "STAR/HISAT2":     ["慢", "高（>30 GB）", "基因组", "是", "全面表征"],
    "kallisto/salmon": ["极快", "低（< 4 GB）", "转录本", "否", "以定量为核心的差异表达"],
})
st.dataframe(comparison_df, use_container_width=True, hide_index=True)

st.divider()

# ── 第三节 — 计数、TPM 与有效长度 ───────────────────────────────────────────
st.subheader("🔢 第三节 — 计数、TPM 与有效长度")

col_def, col_interactive = st.columns([2, 3])

with col_def:
    st.markdown("""
**计数（Counts）** — 分配到某基因/转录本的 reads 数量。
- 整数值，反映原始 reads 支持度
- DESeq2/edgeR 直接使用（这些工具对计数分布建模）
- *不对基因长度或文库大小进行归一化*

**TPM**（每百万转录本）
- 同时对基因长度和文库大小进行归一化
- 适合在样本内或样本间比较表达水平
- *不推荐作为基于计数的差异表达工具的输入*

**有效长度（Effective length）**
- 针对片段大小分布进行调整
- 较短的有效长度意味着片段可能起源的位置更少
- 影响丰度估计：`估计计数 = TPM × 有效长度 / 1e6 × 文库大小`
    """)

with col_interactive:
    st.markdown("**交互演示：长度和深度如何影响 TPM**")

    col_s1, col_s2 = st.columns(2)
    with col_s1:
        tx_len_A = st.slider("TxA 长度（bp）", 200, 3000, 1000, 100)
    with col_s2:
        tx_len_B = st.slider("TxB 长度（bp）", 200, 3000, 500, 100)

    frag_len = st.slider("平均片段长度（bp）", 50, 400, 180, 10)

    lengths = [max(1, tx_len_A - frag_len + 1),
               max(1, tx_len_B - frag_len + 1),
               max(1, 250 - frag_len + 1)]
    counts  = [500, 200, 50]
    rpk     = [c / (l / 1000) for c, l in zip(counts, lengths)]
    tpm     = [r / sum(rpk) * 1e6 for r in rpk]

    result_df = pd.DataFrame({
        "转录本":        ["TxA", "TxB", "TxC"],
        "原始计数":      counts,
        "有效长度（bp）": [int(l) for l in lengths],
        "TPM":           [f"{t:.1f}" for t in tpm],
    })
    st.dataframe(result_df, use_container_width=True, hide_index=True)
    st.caption(
        "调整转录本长度和片段长度，观察即使原始计数不变，TPM 如何随之变化。"
    )

st.divider()

# ── 第四节 — 导入核查清单 ────────────────────────────────────────────────────
st.subheader("✅ 第四节 — 导入定量数据")

col_check, col_example = st.columns([2, 3])

with col_check:
    st.markdown("""
分析之前，数据导入应满足以下条件：
    """)

    checklist = [
        ("✅", "所有样本的定量文件均已就位"),
        ("✅", "文件中的样本名与元数据匹配"),
        ("✅", "元数据文件包含 groupA 和 batch 列"),
        ("✅", "已记录参考转录组/注释版本"),
        ("✅", "文件路径结构清晰、可重现"),
        ("⚠️", "若样本名不一致，分析将静默失败或产生错误结果"),
    ]
    for icon, item in checklist:
        st.markdown(f"{icon} {item}")

with col_example:
    st.markdown("**示例：样本名不一致导致分析失败**")

    ok_df = pd.DataFrame({
        "counts.csv 列名":    ["D01_control", "D01_treatment", "D02_control"],
        "metadata sample_name": ["D01_control", "D01_treatment", "D02_control"],
        "匹配？": ["✅", "✅", "✅"],
    })
    bad_df = pd.DataFrame({
        "counts.csv 列名":    ["sample1", "sample2", "sample3"],
        "metadata sample_name": ["D01_control", "D01_treatment", "D02_control"],
        "匹配？": ["❌", "❌", "❌"],
    })

    st.markdown("**✅ 正确——名称匹配：**")
    st.dataframe(ok_df, use_container_width=True, hide_index=True)
    st.markdown("**❌ 错误——名称不匹配：**")
    st.dataframe(bad_df, use_container_width=True, hide_index=True)
    st.error("名称不匹配会导致元数据与样本脱钩——分析结果将失去生物学意义。")

st.divider()

# ── 第五节 — 转录本到基因的汇总 ─────────────────────────────────────────────
st.subheader("🔗 第五节 — 转录本到基因的汇总")

st.markdown("""
一个基因通常有**多个转录本**（剪接变体）。
对于标准的差异表达分析，我们通常在**基因水平**工作——需要进行聚合。
""")

col_tx, col_gene = st.columns(2)

with col_tx:
    st.markdown("**转录本水平（来自 kallisto/salmon）：**")
    tx_df = pd.DataFrame({
        "转录本 ID":   ["ENST001.1", "ENST001.2", "ENST001.3"],
        "基因 ID":     ["ENSG001"] * 3,
        "基因名":      ["BRCA1"] * 3,
        "TPM":         [120.3, 45.2, 8.1],
        "估计计数":    [540, 210, 38],
    })
    st.dataframe(tx_df, use_container_width=True, hide_index=True)

with col_gene:
    st.markdown("**经 tximport / 基因水平聚合后：**")
    gene_df = pd.DataFrame({
        "基因名":   ["BRCA1"],
        "基因 ID":  ["ENSG001"],
        "总计数":   [788],
        "总 TPM":   [173.6],
    })
    st.dataframe(gene_df, use_container_width=True, hide_index=True)
    st.info(
        "tximport 将转录本水平数值聚合到基因水平，"
        "并在汇总过程中考虑有效长度。"
    )

st.divider()

# ── 第六节 — 注释 ────────────────────────────────────────────────────────────
st.subheader("🏷️ 第六节 — 基因注释与 ID")

col_ann, col_map = st.columns([2, 3])

with col_ann:
    st.markdown("""
**基因身份的层次：**

| ID 类型 | 示例 | 说明 |
|---|---|---|
| 转录本 ID | ENST00000367770.8 | 特定亚型 + 版本号 |
| 基因 ID | ENSG00000139618.19 | 稳定的 Ensembl 基因 ID + 版本号 |
| 基因名 | BRCA2 | 人类可读，但不总是稳定 |

**为什么注释很重要：**
- 版本不一致会导致不同文件之间的 ID 无法匹配
- 基因名可能对应多个 ID，或随时间变化
- 非模式生物的注释可能不完整甚至缺失
    """)

with col_map:
    st.markdown("**注释映射示例：**")
    map_df = pd.DataFrame({
        "转录本 ID":    ["ENST00000367770.8", "ENST00000544455.6"],
        "→ 基因 ID":    ["ENSG00000139618.19", "ENSG00000139618.19"],
        "→ 基因名":     ["BRCA2", "BRCA2"],
        "→ 染色体":     ["13", "13"],
        "注释来源":     ["Ensembl v109", "Ensembl v109"],
    })
    st.dataframe(map_df, use_container_width=True, hide_index=True)

    st.warning("""
**注释常见陷阱：**
- 基因名并非总是唯一——某些基因名对应多个基因
- 务必记录所使用的 Ensembl 版本
- 在跨数据库匹配 ID 时，需去掉版本后缀（如 `.8`、`.19`）
    """)

st.divider()

# ── 第七节 — 元数据 ──────────────────────────────────────────────────────────
st.subheader("📋 第七节 — 元数据：分析的骨架")

metadata = st.session_state.get("metadata")

if metadata is not None and len(metadata) > 0:
    st.success("✅ 你上传的元数据如下所示。")
    st.dataframe(metadata.head(10), use_container_width=True)
else:
    st.info("未上传数据——显示示例元数据结构。")
    example_meta = pd.DataFrame({
        "sample_name":   ["D01_control", "D01_treatment", "D02_control", "D02_treatment"],
        "groupA":        ["control", "treatment", "control", "treatment"],
        "donor":         ["D01", "D01", "D02", "D02"],
        "batch":         ["batch1", "batch1", "batch2", "batch2"],
        "sex":           ["F", "F", "M", "M"],
        "age":           [37, 37, 40, 40],
    })
    st.dataframe(example_meta, use_container_width=True, hide_index=True)

st.markdown("""
| 列名 | 用于 | 相关课程 |
|---|---|---|
| `groupA` | 分组比较 / 差异表达 | 第3、4、7课 |
| `batch` | 批次建模 / 校正 | 第4、6课 |
| `donor` | 配对设计 | 第4、7课 |
| `sex`、`age` | 协变量 / 解读 | 第5、7课 |
""")

st.divider()

# ── 第八节 — 常见陷阱 ────────────────────────────────────────────────────────
st.subheader("⚠️ 第八节 — 常见陷阱")

st.error("**样本名不匹配。** 若名称不一致，计数矩阵和元数据会静默脱钩。")
st.error("**将 TPM 用作差异表达输入。** DESeq2 和 edgeR 需要原始计数——TPM 移除了这些工具所依赖的计数分布属性。")
st.warning("**忽略转录本聚合。** 在未聚合的情况下对转录本水平数据进行差异表达分析，会不恰当地混合亚型。")
st.warning("**忽视注释版本。** 来自不同 Ensembl 版本的 ID 可能无法匹配。")
st.info("**假设基因名唯一。** 部分基因名对应旁系同源基因或假基因。")

st.divider()

# ── 第九节 — 关键总结 ────────────────────────────────────────────────────────
st.subheader("📌 关键总结")

cols = st.columns(5)
msgs = [
    ("🗺️", "定量连接 reads 与分析", "FASTQ → 计数矩阵需要经过几个关键步骤。"),
    ("🔢", "计数 ≠ TPM", "它们回答不同的问题，适用场景也不同。"),
    ("🔗", "将转录本聚合到基因", "基因水平分析是差异表达流程的标准做法。"),
    ("🏷️", "注释版本很重要", "务必记录你使用的参考序列和版本。"),
    ("📋", "元数据就是一切", "良好的元数据使每个下游解读成为可能。"),
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
    st.page_link("pages/01_RNAseq_Foundations.py",
                 label="← 第1课：RNA-seq 基础", icon="🔬")
with col_n2:
    st.page_link("pages/03_Filtering.py",
                 label="第3课：过滤 →", icon="🔍")
