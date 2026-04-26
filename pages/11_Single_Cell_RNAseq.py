"""
pages/11_Single_Cell_RNAseq.py
第11课 — 单细胞 RNA-seq（入门模块）

仅为概念性介绍。本应用未实现真正的单细胞分析流程。
所有可视化均使用模拟数据作为示意。
"""

import numpy as np
import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go

st.set_page_config(
    page_title="第11课 — 单细胞 RNA-seq",
    page_icon="🔬",
    layout="wide",
)

st.title("🔬 第11课 — 单细胞 RNA-seq（入门模块）")
st.info(
    "**这仅为概念性介绍。** "
    "本应用未实现真正的单细胞分析流程。"
    "所有可视化均使用模拟数据作为示意。",
    icon="ℹ️",
)

# ── 第一节：Bulk vs 单细胞 ───────────────────────────────────────────────────
st.subheader("📊 Bulk RNA-seq vs 单细胞 RNA-seq")

comparison_data = {
    "属性": [
        "测量单元",
        "生物学分辨率",
        "典型输出",
        "稀疏性",
        "常用可视化",
        "样本规模",
        "计算复杂度",
        "典型应用场景",
    ],
    "Bulk RNA-seq": [
        "数千至数百万细胞的平均值（每样本）",
        "样本水平；掩盖细胞间异质性",
        "基因 × 样本计数矩阵（密集）",
        "低——每个样本大多数基因均可检测",
        "PCA、热图、火山图、MA 图",
        "5–100+ 个样本（个体、时间点、条件）",
        "中等——标准 R/Python 流程",
        "组织水平差异表达、人群研究、生物标志物发现",
    ],
    "单细胞 RNA-seq": [
        "单个细胞的转录组",
        "细胞水平；揭示亚群和罕见细胞类型",
        "基因 × 细胞计数矩阵（极度稀疏）",
        "高——每个细胞大多数基因计数为零（dropout）",
        "UMAP、t-SNE、小提琴图、点图、轨迹图",
        "数千至数百万个细胞；样本数较少",
        "高——大型矩阵、基于图的聚类、整合分析",
        "细胞类型发现、轨迹分析、细胞通讯",
    ],
}

st.dataframe(
    pd.DataFrame(comparison_data).set_index("属性"),
    use_container_width=True,
)

st.divider()

# ── 第二节：关键概念 ─────────────────────────────────────────────────────────
st.subheader("💡 关键单细胞概念")

col1, col2 = st.columns(2)
with col1:
    st.markdown("""
### 细胞 vs 样本
在 bulk RNA-seq 中，计数矩阵的每一**列**代表一个生物学**样本**
（如来自一个患者的组织活检）。测量到的信号是该组织所有细胞表达量的平均值。

在 scRNA-seq 中，每一**列**代表单个**细胞**。典型的 10x Chromium 实验
捕获 2,000–10,000 个细胞，产生一个庞大的矩阵——
通常是 20,000 个基因 × 5,000 个细胞——其中大多数值为零。

### 稀疏性与 Dropout
每个细胞的测序深度远低于 bulk RNA-seq。许多转录本即使存在于细胞中也无法被捕获。
这些 **dropout** 导致矩阵高度稀疏。

这种稀疏性并非纯粹的噪声——它反映了 mRNA 捕获的随机性。
专门设计的归一化方法（如 scran、SCTransform）用于处理这种情况，
而不是简单地应用 CPM 或 log 归一化。
    """)
with col2:
    st.markdown("""
### 聚类与细胞类型
由于组织内细胞的异质性，scRNA-seq 的首要目标通常是
**按转录相似性对细胞进行分组**。基于共享最近邻图的图聚类
（Leiden、Louvain）是标准方法。

每个簇（理想情况下）对应一个细胞群——但簇本身没有标签。
你通过检查其**标记基因**来识别每个簇的身份：
即在该簇中相对于其他所有簇高度且特异性表达的基因。

### UMAP
UMAP（均匀流形近似与投影）将高维转录组空间（数千个基因）
压缩为 2D 以供可视化。UMAP 空间中靠近的点在转录上相似。
UMAP 用于可视化簇结构、轨迹和标记基因表达。

**重要：** UMAP 坐标不能解读为距离——只有邻域关系才有意义。

### 批次效应与整合
当来自多个实验、供体或测序批次的细胞合并时，
技术变异（批次效应）可能主导生物学信号。
整合方法（Harmony、Seurat CCA、scVI）对数据集进行对齐，
使来自不同批次的相同细胞类型聚集在一起，而非按批次分组。
    """)

st.divider()

# ── 第三节：模拟 UMAP 可视化 ─────────────────────────────────────────────────
st.subheader("🗺️ 概念性 UMAP 可视化")
st.caption("**仅使用模拟数据作示意。** 非真实 scRNA-seq 数据。")

rng = np.random.default_rng(42)

CLUSTER_DEFS = [
    ("T 细胞",      (-4.0, 3.5),   350,  "#2563eb"),
    ("单核细胞",    (4.5, 2.0),    280,  "#dc2626"),
    ("NK 细胞",     (-3.5, -3.0),  180,  "#16a34a"),
    ("B 细胞",      (1.5, -4.5),   220,  "#9333ea"),
    ("上皮细胞",    (6.0, -1.5),   160,  "#ea580c"),
    ("巨噬细胞",    (2.5, 4.5),    200,  "#0891b2"),
    ("成纤维细胞",  (-1.0, 0.5),    90,  "#854d0e"),
]

rows = []
for name, (cx, cy), n, color in CLUSTER_DEFS:
    xs = rng.normal(cx, 0.9, n)
    ys = rng.normal(cy, 0.9, n)
    for x, y in zip(xs, ys):
        rows.append({"UMAP_1": x, "UMAP_2": y, "细胞簇": name})

umap_df = pd.DataFrame(rows)

fig_umap = px.scatter(
    umap_df, x="UMAP_1", y="UMAP_2", color="细胞簇",
    title="模拟 UMAP — 7 个示意细胞簇<br><sup>模拟数据，仅作示意</sup>",
    color_discrete_map={name: color for name, _, _, color in CLUSTER_DEFS},
    opacity=0.7,
    labels={"UMAP_1": "UMAP 1", "UMAP_2": "UMAP 2"},
    height=500,
)
fig_umap.update_traces(marker=dict(size=4))
fig_umap.update_layout(legend=dict(title="细胞簇", itemsizing="constant"))
st.plotly_chart(fig_umap, use_container_width=True)

st.markdown("#### 每个簇的示例标记基因")
st.caption("固定的示意示例——真实标记基因来自已发表文献。")

marker_table = pd.DataFrame([
    {"细胞簇": "T 细胞",     "关键标记基因": "CD3D, CD3E, IL7R, LTB, TRAC",
     "说明": "泛 T 细胞标记；IL7R 标记初始/记忆 T 细胞"},
    {"细胞簇": "单核细胞",   "关键标记基因": "LST1, S100A8, S100A9, CTSS, FCN1",
     "说明": "经典单核细胞标记；S100A8/A9 是危险信号"},
    {"细胞簇": "NK 细胞",    "关键标记基因": "GNLY, NKG7, GZMB, FCER1G, KLRD1",
     "说明": "细胞毒性 NK 细胞标记；GNLY/NKG7 标记细胞毒性颗粒"},
    {"细胞簇": "B 细胞",     "关键标记基因": "CD79A, MS4A1, BANK1, CD22, HLA-DQA1",
     "说明": "B 细胞身份标记；MS4A1 = CD20"},
    {"细胞簇": "上皮细胞",   "关键标记基因": "EPCAM, KRT8, KRT18, KRT19, CLDN4",
     "说明": "上皮细胞身份；角蛋白定义上皮起源"},
    {"细胞簇": "巨噬细胞",   "关键标记基因": "CD68, MRC1, MARCO, APOE, C1QA",
     "说明": "组织驻留巨噬细胞标记；C1Q 标记成熟巨噬细胞"},
    {"细胞簇": "成纤维细胞", "关键标记基因": "COL1A1, COL1A2, DCN, LUM, THY1",
     "说明": "间质成纤维细胞标记；胶原蛋白定义产生 ECM 的细胞"},
])
st.dataframe(marker_table.set_index("细胞簇"), use_container_width=True)

st.divider()

# ── 第四节：典型 scRNA-seq 分析流程 ─────────────────────────────────────────
st.subheader("🔧 典型单细胞分析流程")
st.caption("概念概述——本应用中未实现。")

workflow_steps = [
    ("1. 原始数据与比对",
     "FASTQ → Cell Ranger（10x）或 STARsolo → filtered_feature_bc_matrix。"
     "输出：细胞 × 基因 UMI 计数矩阵。"),
    ("2. 质量控制",
     "按以下标准过滤细胞：检测到的最少基因数、最大线粒体 reads 比例（标记垂死细胞）、"
     "最少 UMI 计数。使用 DoubletFinder 或 Scrublet 等工具去除双细胞。"),
    ("3. 归一化",
     "按文库大小归一化（scran 合并法，或简单 CPM）。Log 变换。"
     "Seurat 的 SCTransform 使用负二项回归建模技术变异。"),
    ("4. 特征选择",
     "选择高变异基因（HVG）——通常 2,000–5,000 个基因——"
     "以减少噪声并聚焦于信息丰富的变异。"),
    ("5. 降维",
     "对 HVG 进行 PCA（通常保留 10–50 个 PC）。然后对 PCA 嵌入"
     "进行 UMAP 或 t-SNE，得到 2D 可视化结果。"),
    ("6. 聚类",
     "在 PCA 坐标上构建 k 最近邻图。"
     "应用 Leiden 或 Louvain 社区检测。分辨率参数控制聚类粒度。"),
    ("7. 标记基因识别",
     "找到每个簇相对于其他所有簇差异表达的基因"
     "（Wilcoxon 秩和检验）。排名靠前的标记基因定义簇的身份。"),
    ("8. 细胞类型注释",
     "将标记基因与文献或参考数据集中的已知细胞类型特征匹配"
     "（如 HCA、CellTypist、SingleR）。"),
    ("9. 下游分析",
     "差异丰度分析、轨迹分析（伪时间）、细胞通讯、"
     "多样本整合，或与空间数据整合。"),
]

for step, desc in workflow_steps:
    with st.expander(step):
        st.markdown(desc)

st.divider()

# ── 第五节：为什么单细胞分析未完整实现 ─────────────────────────────────────
st.subheader("🔬 第五节 — 为什么单细胞分析未在本应用中完整实现")

st.info(
    "**第11课仅为入门概念模块。** "
    "以下单细胞功能在本网页应用中不可用——下表解释了每项功能未实现的具体原因。",
    icon="ℹ️",
)


# 直接用 st.dataframe 渲染表格
not_implemented_table = pd.DataFrame({
    "功能": [
        "细胞水平 QC（线粒体 reads 过滤、双细胞去除）",
        "单细胞归一化（scran、SCTransform）",
        "高变异基因（HVG）选择",
        "真实 scRNA-seq 数据的 PCA / SNN 图",
        "真实单细胞矩阵的 UMAP",
        "Leiden / Louvain 图聚类",
        "簇水平差异表达分析",
        "基于标记基因的细胞类型注释",
        "多样本整合（Harmony、scVI、Seurat CCA）",
        "轨迹 / 伪时间分析（Monocle、scVelo）",
        "AnnData（.h5ad）/ Seurat（.rds）文件支持",
        "10x Cell Ranger 输出解析",
    ],
    "未实现原因": [
        "需要 10x Cell Ranger 格式的逐细胞 barcode 级别数据，而非 CSV 文件",
        "专为稀疏 UMI 计数分布设计的方法，不能与 bulk 归一化互换使用",
        "需要跨数千个细胞进行逐细胞方差建模",
        "20,000 基因 × 50,000 细胞的矩阵超出浏览器内存限制（Streamlit Cloud 约 1 GB）",
        "计算量极大；对大型矩阵运行真实 UMAP 需数分钟，在网页应用中会超时",
        "依赖 `leidenalg` 和 `igraph`，含复杂 C 扩展的依赖项与简单 Streamlit 部署不兼容",
        "依赖已完成的聚类结果；没有真实的细胞类型分配就没有意义",
        "需要参考图谱或精选标记基因数据库，未随本应用打包",
        "需要 GPU 或大内存服务器；scVI 依赖 PyTorch",
        "需要 RNA velocity 或超出计数矩阵的伪时间特定数据结构",
        "二进制 HDF5 格式；文件通常达数 GB，无法通过浏览器上传",
        "需要读取稀疏矩阵市场格式（barcodes / features / matrix 三元组），而非 CSV",
    ],
})

st.dataframe(
    not_implemented_table,
    use_container_width=True,
    hide_index=True,
    height=430,
)

st.warning(
    "**一句话总结：** scRNA-seq 分析需要专用文件格式、大内存（8–64 GB RAM）、"
    "复杂的库依赖，有时还需要 GPU 资源。这些约束使其不适合轻量级的浏览器端教学应用。",
    icon="⚠️",
)

st.markdown("""
**进行真正单细胞分析，请使用以下工具：**

| 工具 | 语言 | 说明 |
|------|------|------|
| [Scanpy](https://scanpy.readthedocs.io/) | Python | 标准单细胞分析工具包 |
| [Seurat](https://satijalab.org/seurat/) | R | 领域内广泛使用 |
| [Bioconductor scran / scater](https://bioconductor.org/books/release/OSCA/) | R | 严格的统计框架 |

在配备 **≥16 GB RAM 的本地机器**上运行，或使用 **HPC 集群**。
""")

st.divider()

col_n1, col_n2 = st.columns(2)
with col_n1:
    st.page_link("pages/10_Public_Data.py",
                 label="← 第10课：公共数据与可重复性", icon="🌐")
with col_n2:
    st.markdown("*更多课程即将上线。*")
