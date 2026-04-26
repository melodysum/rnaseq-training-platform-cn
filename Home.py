"""
Home.py — RNA-seq 互动训练平台首页
运行方式: streamlit run Home.py
"""

import streamlit as st
import pandas as pd
from utils.data_loader import validate_and_parse, init_session_data, load_demo_data

st.set_page_config(
    page_title="RNA-seq 训练平台",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap');
    html, body, [class*="css"] { font-family: 'Inter', sans-serif; }
    .hero-title { font-size: 2.6rem; font-weight: 700; color: #0f172a; line-height: 1.2; margin-bottom: 0.4rem; }
    .hero-sub { font-size: 1.1rem; color: #475569; max-width: 680px; line-height: 1.6; }
    .lesson-card { background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 12px; padding: 1.2rem 1.4rem; margin-bottom: 0.8rem; transition: box-shadow 0.2s; }
    .lesson-card:hover { box-shadow: 0 4px 16px rgba(0,0,0,0.07); }
    .lesson-title { font-size: 1.05rem; font-weight: 600; color: #1e293b; }
    .lesson-desc  { font-size: 0.88rem; color: #64748b; margin-top: 0.2rem; }
    .badge-available { background: #dcfce7; color: #166534; border-radius: 999px; padding: 2px 10px; font-size: 0.75rem; font-weight: 600; }
    .badge-coming { background: #f1f5f9; color: #94a3b8; border-radius: 999px; padding: 2px 10px; font-size: 0.75rem; font-weight: 600; }
    .path-step { display: flex; align-items: flex-start; margin-bottom: 0.7rem; }
    .step-num { background: #3b82f6; color: white; border-radius: 50%; width: 28px; height: 28px; display: flex; align-items: center; justify-content: center; font-size: 0.8rem; font-weight: 700; flex-shrink: 0; margin-right: 0.7rem; margin-top: 2px; }
    .step-text { font-size: 0.92rem; color: #334155; }
    .divider { border: none; border-top: 1px solid #e2e8f0; margin: 1.5rem 0; }
    .section-label { font-size: 0.75rem; font-weight: 700; text-transform: uppercase; letter-spacing: 0.08em; color: #94a3b8; margin-bottom: 0.6rem; }
</style>
""", unsafe_allow_html=True)

init_session_data()

st.markdown("""
<div class="hero-title">🧬 RNA-seq 互动训练平台</div>
<div class="hero-sub">
    一个循序渐进的互动学习环境，带你理解转录组数据分析的完整流程——从原始计数矩阵到生物学解读。
</div>
""", unsafe_allow_html=True)

st.markdown('<hr class="divider">', unsafe_allow_html=True)

col_left, col_right = st.columns([3, 2], gap="large")

with col_left:
    st.markdown('<div class="section-label">课程内容</div>', unsafe_allow_html=True)

    lessons = [
        {"icon": "🔬", "title": "第1课 — RNA-seq 基础与实验设计",
         "desc": "了解 RNA-seq 数据是如何产生的、关键测序概念的含义，以及为什么良好的实验设计对下游分析至关重要。",
         "page": "pages/01_RNAseq_Foundations.py", "available": True},
        {"icon": "📦", "title": "第2课 — 定量、导入与注释",
         "desc": "学习 RNA-seq reads 如何转换为转录本和基因水平的表达数据，以及元数据和注释如何使下游分析成为可能。",
         "page": "pages/02_Quantification_Import_Annotation.py", "available": True},
        {"icon": "🔍", "title": "第3课 — 低表达基因过滤",
         "desc": "理解为什么在分析前需要去除低计数基因，以及过滤阈值的选择如何影响保留的基因数量。",
         "page": "pages/03_Filtering.py", "available": True},
        {"icon": "📐", "title": "第4课 — 多重检验与 FDR 校正",
         "desc": "学习 Benjamini-Hochberg 校正的工作原理，理解过滤如何改变 FDR 结果，以及如何区分统计意义与生物学意义上的基因丢失。",
         "page": "pages/04_FDR.py", "available": True},
        {"icon": "📊", "title": "第5课 — 探索性分析与 PCA",
         "desc": "可视化样本间关系、识别离群样本，并解读 RNA-seq 数据中的主成分分析结果。",
         "page": "pages/05_Exploratory_Analysis_PCA.py", "available": True},
        {"icon": "🔄", "title": "第6课 — 批次效应校正",
         "desc": "探索技术性批次效应的来源，学习如何用 PCA 检测批次效应，以及如何在差异表达分析前进行校正。",
         "page": "pages/06_Batch_Correction.py", "available": True},
        {"icon": "🧪", "title": "第7课 — 差异表达分析",
         "desc": "进行完整的差异表达分析，理解统计模型，并正确解读分析结果。",
         "page": "pages/07_Differential_Expression.py", "available": True},
        {"icon": "🗺️", "title": "第8课 — 聚类与热图",
         "desc": "学习如何对基因和样本进行聚类、构建表达热图，并从差异表达结果中识别基因模块。",
         "page": "pages/08_Clustering_Heatmaps.py", "available": True},
        {"icon": "🧩", "title": "第9课 — 功能富集分析（GO / GSEA）",
         "desc": "使用内置的教学通路数据集进行过度代表性分析（ORA）和 GSEA 类排序富集分析。理解 ORA 与 GSEA 各自衡量什么，以及如何解读富集图。",
         "page": "pages/09_Functional_Enrichment.py", "available": True},
        {"icon": "🌐", "title": "第10课 — 公共数据与可重复性分析",
         "desc": "探索真实的 GEO 公共数据集，学习如何选取适合重分析的文件，完成可重复性核查清单，并理解网页端、本地端和 HPC 分析的区别。",
         "page": "pages/10_Public_Data.py", "available": True},
        {"icon": "🔬", "title": "第11课 — 单细胞 RNA-seq（入门模块）",
         "desc": "单细胞 RNA-seq 概念入门：它与 bulk RNA-seq 的区别、聚类和 UMAP 能揭示什么，以及真正的单细胞分析需要哪些工具。包含模拟 UMAP 可视化示例。",
         "page": "pages/11_Single_Cell_RNAseq.py", "available": True},
    ]

    for lesson in lessons:
        badge = (
            '<span class="badge-available">已开放</span>'
            if lesson["available"]
            else '<span class="badge-coming">即将上线</span>'
        )
        st.markdown(f"""
        <div class="lesson-card">
            <div class="lesson-title">{lesson["icon"]} {lesson["title"]} &nbsp; {badge}</div>
            <div class="lesson-desc">{lesson["desc"]}</div>
        </div>
        """, unsafe_allow_html=True)
        if lesson["page"]:
            label = f"进入 {lesson['title']} →" if lesson["available"] else "预览（即将上线）"
            st.page_link(lesson["page"], label=label)

    st.markdown('<hr class="divider">', unsafe_allow_html=True)

    st.markdown('<div class="section-label">推荐学习路径</div>', unsafe_allow_html=True)
    steps = [
        ("从这里开始", "阅读 <b>第1课 — 基础</b>，了解 RNA-seq 测量的是什么，以及实验设计如何影响所有下游分析。"),
        ("数据是怎么来的", "学习 <b>第2课 — 定量</b>，理解 FASTQ reads 是如何变成计数矩阵的。"),
        ("清洗数据", "打开 <b>第3课 — 过滤</b>，探索阈值选择如何影响保留的基因数量。"),
        ("理解 FDR", "进入 <b>第4课 — FDR</b>，看清过滤与多重检验之间的联系。"),
        ("探索数据结构", "使用 <b>第5课 — 探索性分析</b>，检查 PCA 结果，发现批次效应或离群样本。"),
        ("校正批次效应", "如果存在技术变异，打开 <b>第6课 — 批次校正</b>。"),
        ("运行差异表达", "完成 <b>第7课 — 差异表达</b>，比较分组并解读结果。"),
    ]
    for i, (label, text) in enumerate(steps, 1):
        st.markdown(f"""
        <div class="path-step">
            <div class="step-num">{i}</div>
            <div class="step-text"><b>{label}：</b>{text}</div>
        </div>
        """, unsafe_allow_html=True)

with col_right:
    st.markdown('<div class="section-label">你的数据</div>', unsafe_allow_html=True)

    source = st.session_state.get("data_source", "demo")
    if source == "demo":
        st.info(
            "**当前使用内置演示数据。**\n\n"
            "上传你自己的文件以分析你的数据集。",
            icon="📂",
        )
    else:
        st.success("✅ 正在使用你上传的数据。", icon="✅")
        if st.button("切换回演示数据"):
            counts, metadata = load_demo_data()
            st.session_state["counts"]   = counts
            st.session_state["metadata"] = metadata
            st.session_state["data_source"] = "demo"
            st.session_state.pop("filter_results", None)
            st.rerun()

    with st.expander("📁 上传你自己的数据", expanded=(source == "demo")):
        st.markdown("""
**counts.csv** — 原始基因水平计数矩阵
- 基因为行，样本为列
- 数值必须是**原始整数计数**——不能是 TPM、FPKM 或 CPM
- 第一列 = 基因标识符（用作行索引）

```
gene_symbol, D01_control, D01_treatment, ...
IFNA6,       25,          40,            ...
CAD,          0,           9,            ...
```

**metadata.csv** — 每行对应一个样本
- 必须包含的列：**groupA**（定义比较分组）
- 可选但重要：**donor**（配对差异表达）、**batch**（批次校正）
- 可选描述性列：**sex**、**age** 等注释信息

```
sample_name,   groupA,    donor, batch,  sex, age
D01_control,   control,   D01,   batch1, F,   37
D01_treatment, treatment, D01,   batch1, F,   37
```

⚠️ 两个文件中的样本名必须完全一致（区分大小写）。  
⚠️ 基因 ID 重复、元数据列缺失或计数值非数字均会被拒绝并显示明确的错误信息。  
⚠️ 文件过大（>5,000 基因 × 200+ 样本）在网页端可能运行缓慢——大型数据集建议在本地分析。
        """)
        st.info("输入数据必须是**原始计数值**，不能是 TPM/FPKM/CPM。归一化数值会导致验证错误或分析结果不正确。", icon="⚠️")

        counts_file   = st.file_uploader("上传 counts.csv",   type="csv", key="upload_counts")
        metadata_file = st.file_uploader("上传 metadata.csv", type="csv", key="upload_meta")

        if counts_file and metadata_file:
            counts, metadata, err = validate_and_parse(counts_file, metadata_file)
            if err:
                st.error(f"❌ {err}")
            else:
                st.session_state["counts"]   = counts
                st.session_state["metadata"] = metadata
                st.session_state["data_source"] = "uploaded"
                st.session_state.pop("filter_results", None)
                st.success(
                    f"✅ 已加载 **{len(counts):,} 个基因** × **{len(counts.columns)} 个样本**。"
                )

    st.markdown('<hr class="divider">', unsafe_allow_html=True)

    st.markdown('<div class="section-label">当前数据集</div>', unsafe_allow_html=True)
    counts   = st.session_state["counts"]
    metadata = st.session_state["metadata"]

    groups = metadata["groupA"].value_counts()
    st.markdown(f"""
| | |
|---|---|
| **基因数** | {len(counts):,} |
| **样本数** | {len(counts.columns)} |
| **分组** | {', '.join(f'{g} ({n})' for g, n in groups.items())} |
| **数据来源** | {'演示数据' if source == 'demo' else '已上传'} |
""")

    if "donor" in metadata.columns:
        st.caption(f"检测到配对设计：共 {metadata['donor'].nunique()} 个供体。")
    if "batch" in metadata.columns:
        st.caption(f"批次信息可用：共 {metadata['batch'].nunique()} 个批次。")

    st.markdown('<hr class="divider">', unsafe_allow_html=True)

    st.markdown('<div class="section-label">跳转到课程</div>', unsafe_allow_html=True)
    st.page_link("pages/01_RNAseq_Foundations.py",               label="🔬 第1课 — 基础",           icon="▶")
    st.page_link("pages/02_Quantification_Import_Annotation.py", label="📦 第2课 — 定量",            icon="▶")
    st.page_link("pages/03_Filtering.py",                         label="🔍 第3课 — 过滤",            icon="▶")
    st.page_link("pages/04_FDR.py",                               label="📐 第4课 — FDR",             icon="▶")
    st.page_link("pages/05_Exploratory_Analysis_PCA.py",          label="📊 第5课 — 探索性分析",      icon="▶")
    st.page_link("pages/06_Batch_Correction.py",                  label="🔄 第6课 — 批次校正",        icon="▶")
    st.page_link("pages/07_Differential_Expression.py",           label="🧪 第7课 — 差异表达",        icon="▶")
    st.page_link("pages/08_Clustering_Heatmaps.py",               label="🗺️ 第8课 — 聚类与热图",      icon="▶")
    st.page_link("pages/09_Functional_Enrichment.py",             label="🧩 第9课 — 功能富集分析",    icon="▶")
    st.page_link("pages/10_Public_Data.py",                       label="🌐 第10课 — 公共数据",       icon="▶")
    st.page_link("pages/11_Single_Cell_RNAseq.py",                label="🔬 第11课 — 单细胞 RNA-seq", icon="▶")
