"""
pages/10_Public_Data.py
第10课 — 公共数据与可重复性分析
"""

import streamlit as st
import pandas as pd
from utils.data_loader import init_session_data

st.set_page_config(
    page_title="第10课 — 公共数据与可重复性",
    page_icon="🌐",
    layout="wide",
)
init_session_data()

st.title("🌐 第10课 — 公共数据与可重复性分析")
st.caption(
    "公共 RNA-seq 数据集的结构、如何选择正确的文件，以及如何使你的分析可重现。"
)

# ── 第一节：概念说明 ──────────────────────────────────────────────────────────
with st.expander("📖 什么是 GEO 和 ArrayExpress？", expanded=True):
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("""
### NCBI GEO（基因表达综合数据库）
GEO 是由 NCBI 维护的主要高通量基因表达数据公共存储库。
它接受来自微阵列、RNA-seq、ChIP-seq 及许多其他基因组学平台的数据。

每次提交获得一个 **GSE 登录号**（如 `GSE148036`）代表整个研究，
**GSM 登录号**代表单个样本。**GPL 登录号**描述测序平台（如 Illumina HiSeq 2500）。

当论文说"数据在 GEO 中，登录号为 GSE…"时，你就到这里找。
可在 [ncbi.nlm.nih.gov/geo](https://www.ncbi.nlm.nih.gov/geo/) 搜索。
        """)
    with col2:
        st.markdown("""
### EBI ArrayExpress / BioStudies
ArrayExpress（现已整合进 EBI BioStudies）是 GEO 的欧洲对应库。
许多研究同时在两个存储库中提交数据。登录号格式为 `E-MTAB-XXXX`。

**序列读取档案（SRA）：** 提交到 GEO 的原始 FASTQ reads 存储在 SRA 中，
可通过 SRA Toolkit 的 `prefetch` / `fastq-dump` 下载。
欧洲数据镜像在 EBI 的 ENA（欧洲核苷酸档案）中。

两个存储库均免费开放，可按物种、组织、疾病或关键词搜索。
        """)

    st.markdown("---")
    st.markdown("""
### 通常可获取哪些文件？

| 文件类型 | 内容 | 使用场景 |
|----------|------|----------|
| **原始 FASTQ** | 测序 reads；通过 SRA/ENA 获取 | 从头重新比对；最严格 |
| **计数矩阵** | 基因 × 样本 reads 计数；在补充文件中 | 从计数开始重新分析；更快 |
| **处理/归一化值** | TPM、FPKM、RPKM 表；在补充文件中 | 可视化、探索（不用于差异表达） |
| **Series matrix** | GEO 格式的元数据 + 表达摘要 | 微阵列；也包含样本元数据 |
| **SOFT / MINiML** | 结构化格式的 GEO 元数据 | 自动化元数据解析 |

对于 bulk RNA-seq **重新分析**，补充文件中的基因水平计数矩阵通常是最佳起点——
它避免了重新运行完整比对流程，同时保留了 DESeq2/edgeR 所需的整数计数值。
    """)

st.divider()

# ── 第二节：案例研究 ──────────────────────────────────────────────────────────
st.subheader("🔬 真实公共数据集案例研究")
st.caption(
    "两个真实 bulk RNA-seq GEO 研究，为教学目的整理摘要。"
    "详情存储在本地——无需联网。"
)

CASE_STUDIES = [
    {
        "accession": "GSE148036",
        "title": "结核分枝杆菌感染期间人肺泡巨噬细胞的转录组分析",
        "organism": "智人（Homo sapiens）",
        "study_question": (
            "人肺泡巨噬细胞对结核分枝杆菌感染的转录应答是什么，"
            "在活动性结核与潜伏性结核期间哪些宿主通路发生了失调？"
        ),
        "experimental_design": (
            "从健康供体（n=10）和活动性肺结核患者（n=10）中分离的原代人肺泡巨噬细胞进行 bulk RNA-seq。"
            "每个供体两个条件：未感染对照和体外结核分枝杆菌感染。"
            "配对设计；跨采集点存在批次效应的 2×2 因子结构。"
        ),
        "sample_groups": (
            "对照（未感染巨噬细胞）vs 结核分枝杆菌感染巨噬细胞，"
            "来自健康供体和活动性结核患者。"
            "关键比较：每个供体组内感染 vs 未感染；"
            "健康 vs 结核患者的背景效应。"
        ),
        "input_file_types": (
            "原始 FASTQ reads 存入 SRA（双端，Illumina HiSeq）。"
            "补充文件提供制表符分隔的计数矩阵文本文件（基因水平，STAR 比对至 GRCh38）。"
            "样本元数据可从 GEO Series matrix 和补充 Excel 表中获取。"
        ),
        "recommended_start": (
            "下载补充计数矩阵（GSE148036_raw_counts.txt.gz）。"
            "这直接给你 STAR 生成的基因水平计数，可直接用于 DESeq2/edgeR。"
            "将样本名与 GEO 元数据 SOFT 文件匹配以恢复分组/批次注释。"
            "不要用 FPKM 表进行差异表达分析。"
        ),
        "downstream_workflow": (
            "1. 将计数矩阵和元数据加载到 R。"
            "2. 过滤低计数基因（≥50% 样本中 ≥10 计数）。"
            "3. 以 ~ donor + infection_status 为设计运行 DESeq2 以考虑配对结构。"
            "4. 使用 apeglm 进行 LFC 收缩。"
            "5. 使用 clusterProfiler 进行通路富集（GO、KEGG）。"
            "6. 使用 MSigDB Hallmark 基因集进行 GSEA 以识别宏观通路主题。"
        ),
        "reproducibility_notes": (
            "记录使用的 STAR 和 featureCounts 的确切版本。"
            "参考基因组和 GTF 注释版本（Ensembl release 100 vs 104 可能改变基因计数）。"
            "供体匿名化意味着 GEO 中的样本 ID 可能与论文中的不同——"
            "务必通过 Methods 部分交叉核对。"
            "该研究与结核免疫学相关，是配对巨噬细胞感染设计的良好参考模型。"
        ),
    },
    {
        "accession": "GSE167232",
        "title": "单细胞和空间转录组学揭示结核肉芽肿的架构和细胞状态",
        "organism": "智人 / 恒河猴（Homo sapiens / Macaca mulatta）",
        "study_question": (
            "结核肉芽肿内部的空间组织和转录细胞状态是什么，"
            "稳定型和进展型肉芽肿的免疫和间质细胞龛有何不同？"
        ),
        "experimental_design": (
            "结合 bulk RNA-seq（肉芽肿组织）、scRNA-seq（10x Chromium）"
            "和空间转录组学（Visium），对来自非人灵长类动物（恒河猴）结核模型"
            "和人体尸检组织的结核肉芽肿样本进行分析。"
            "每只动物多种肉芽肿类型：稳定型（无菌）vs 进展型（培养阳性）。"
            "无配对设计；按动物 ID 进行批次校正。"
        ),
        "sample_groups": (
            "稳定型肉芽肿 vs 进展型肉芽肿。在 scRNA-seq 中："
            "巨噬细胞亚群（泡沫型、炎症型、过渡型），T 细胞亚群（CD4、CD8、调节性），"
            "中性粒细胞、成纤维细胞、上皮样细胞。"
            "空间数据通过去卷积推断每个点位的细胞类型组成。"
        ),
        "input_file_types": (
            "Bulk RNA-seq：补充 CSV 中的原始计数。"
            "scRNA-seq：每个样本的 Cell Ranger 输出（barcodes、features、matrix）——"
            "可直接加载到 Seurat 或 Scanpy。"
            "空间 Visium：包含组织图像和点位水平计数矩阵的 10x Space Ranger 输出文件夹，"
            "可在 Seurat 或 Squidpy 中加载。"
        ),
        "recommended_start": (
            "bulk 重新分析：下载补充计数矩阵，应用标准 DESeq2 流程，"
            "以 ~ granuloma_type 作为主效应。"
            "scRNA-seq：如提供则下载处理好的 Seurat RDS 对象，"
            "或重新加载原始 Cell Ranger 输出并用 Seurat/Scanpy 重新处理。"
            "空间数据需要 Space Ranger 输出文件夹结构——不适合基于浏览器的重新分析。"
        ),
        "downstream_workflow": (
            "Bulk：DESeq2 → 通路富集 → 识别炎症 vs 消退特征。"
            "scRNA-seq：QC → 归一化 → HVG 选择 → PCA → UMAP → 聚类 → "
            "标记基因识别 → 细胞类型注释 → 差异丰度分析。"
            "空间：RCTD 或 cell2location 去卷积 → 空间相关性 → "
            "CellChat 或 NicheNet 配体-受体相互作用分析。"
        ),
        "reproducibility_notes": (
            "多模态数据集需要仔细追踪版本：Seurat v4 vs v5 更改了默认归一化流程。"
            "Cell Ranger 版本影响 barcode 识别和比对。"
            "空间去卷积工具对所用的参考单细胞图谱非常敏感。"
            "该研究展示了现代复杂数据集如何组合多种数据模态，"
            "需要不同的计算基础设施——不适合在没有 HPC 资源的情况下进行简单网页重新分析。"
        ),
    },
]

study_options = {f"{cs['accession']} — {cs['title'][:50]}…": i
                 for i, cs in enumerate(CASE_STUDIES)}
selected_label = st.selectbox("选择一项研究：", list(study_options.keys()))
cs = CASE_STUDIES[study_options[selected_label]]

st.markdown(f"### {cs['accession']}：{cs['title']}")
st.markdown(f"**物种：** {cs['organism']}")

tabs = st.tabs([
    "🔬 研究设计", "📁 文件与起点",
    "🔧 分析流程", "♻️ 可重复性说明"
])

with tabs[0]:
    st.markdown(f"**研究问题：**\n\n{cs['study_question']}")
    st.markdown(f"**实验设计：**\n\n{cs['experimental_design']}")
    st.markdown(f"**样本分组：**\n\n{cs['sample_groups']}")

with tabs[1]:
    st.markdown(f"**通常可获取的文件类型：**\n\n{cs['input_file_types']}")
    st.info(f"💡 **推荐起始文件：**\n\n{cs['recommended_start']}", icon="📥")

with tabs[2]:
    st.markdown(f"**典型下游分析流程：**\n\n{cs['downstream_workflow']}")
    st.warning(
        "公共数据在使用前通常需要清洗："
        "GEO 元数据中的样本名可能与论文不同，"
        "分组标签可能使用缩写，批次/供体注释可能需要从补充表格中手动整理。",
        icon="⚠️",
    )

with tabs[3]:
    st.markdown(f"{cs['reproducibility_notes']}")

st.divider()

# ── 第三节：计数 vs TPM/FPKM 说明 ────────────────────────────────────────────
st.subheader("📏 原始计数 vs TPM / FPKM / CPM")

col_a, col_b = st.columns(2)
with col_a:
    st.markdown("""
**原始计数**是比对到每个基因的测序 reads 数量。
它们是整数，直接反映测序深度和基因表达量。

✅ 使用原始计数用于：
- 差异表达（DESeq2、edgeR、limma-voom）
- 任何显式对计数分布建模的方法

❌ 不要单独用原始计数进行：
- 样本间归一化比较
  （不同样本的文库大小不同）
    """)
with col_b:
    st.markdown("""
**TPM（每百万转录本）**、**FPKM** 和 **CPM** 是归一化值，
考虑了文库大小以及（对于 FPKM/TPM）基因长度。

✅ 使用归一化值用于：
- 可视化（热图、PCA、箱线图）
- 不进行正式检验的跨样本表达比较
- 单基因探索性图表

❌ 不要将 TPM/FPKM 用于：
- DESeq2、edgeR 的输入——这些模型需要原始计数
- 当原始计数也可获得时的重新分析
    """)

st.info(
    "**本应用期望输入原始计数。** 如果你上传 TPM 或 FPKM 值，"
    "统计模型（log-CPM 归一化、t 检验）仍会运行，"
    "但结果的意义会降低，而且 log 转换后可能出现负值并触发验证警告。",
    icon="⚠️",
)

st.divider()

# ── 第四节：元数据规范 ───────────────────────────────────────────────────────
st.subheader("🧹 元数据规范")
st.markdown("""
良好的元数据与良好的数据同等重要。重复使用公共数据集时的常见问题：

- **样本名不匹配** — GEO 样本 ID（GSM…）很少与计数矩阵中的样本名匹配；
  务必检查补充元数据文件。
- **分组标签不一致** — "control"、"ctrl"、"Control"、"CTRL" 是四个不同的字符串。
  分析前进行标准化。
- **批次或供体信息缺失** — GEO 元数据通常不完整；检查论文本身的 Methods 和补充表格。
- **时间点或条件模糊** — 特别是在纵向或多因子设计中；在对分组进行编码前明确你想要的比较。
- **重复或交换的样本** — 特别是在大型多中心研究中；PCA 和样本相关性检查可以揭示明显异常。

**本应用所需的最少元数据：**
- `groupA` — 必须；定义比较分组
- `donor` — 可选；启用配对差异表达分析
- `batch` — 可选；启用批次校正流程
""")

st.divider()

# ── 第五节：可重复性核查清单 ─────────────────────────────────────────────────
st.subheader("✅ 可重复性核查清单")
st.caption(
    "逐项勾选以追踪你的分析进度。"
    "这是教学工具——状态不会在会话之间保存。"
)

checklist_items = [
    ("📋 研究背景", [
        ("已记录登录号（如 GSE…、E-MTAB-…）", False),
        ("已阅读并引用对应论文", False),
        ("已记录研究设计和分组定义", False),
        ("已确认物种和组织类型", False),
    ]),
    ("📁 数据来源", [
        ("已区分原始文件和处理后文件", False),
        ("已验证下载文件的校验和（若提供）", False),
        ("已记录基因注释版本（如 Ensembl release、GTF 日期）", False),
        ("已记录参考基因组版本（如 GRCh38、hg19）", False),
    ]),
    ("🧬 元数据", [
        ("已检查并清洗样本元数据", False),
        ("计数矩阵和元数据中的样本名已匹配", False),
        ("分组标签已标准化且无歧义", False),
        ("已识别批次和供体列（若存在）", False),
    ]),
    ("💻 计算环境", [
        ("已记录软件版本（R、Python、关键包）", False),
        ("已记录环境（conda env、renv、requirements.txt）", False),
        ("已在适用处设置随机种子", False),
        ("分析脚本已进行版本控制（git）", False),
    ]),
    ("📊 分析输出", [
        ("分析步骤已用说明性注释记录", False),
        ("输出文件命名一致且具描述性", False),
        ("图形可追溯到生成它们的代码", False),
        ("过滤和参数选择已在笔记中说明", False),
        ("已记录重新分析的局限性", False),
    ]),
]

all_checked = 0
all_total   = 0
for section, items in checklist_items:
    st.markdown(f"**{section}**")
    for label, _ in items:
        checked = st.checkbox(label, key=f"chk_{label}")
        if checked:
            all_checked += 1
        all_total += 1

progress = all_checked / all_total if all_total > 0 else 0
st.progress(progress, text=f"核查清单进度：{all_checked}/{all_total} 项")
if all_checked == all_total:
    st.success("🎉 所有核查项目已完成——你的分析记录完善！", icon="✅")

st.divider()

# ── 第六节：网页端 vs 本地 vs HPC 对比 ──────────────────────────────────────
st.subheader("🖥️ 网页端 vs 本地分析 vs HPC")

comparison_data = {
    "特性": [
        "推荐用途",
        "优势",
        "局限",
        "典型数据规模",
        "最适合示例",
    ],
    "本网页应用": [
        "教学、探索、小型数据集",
        "无需配置；交互式；非常适合学习核心概念",
        "计算资源有限；无文件持久化；浏览器内存限制；无法自定义安装包",
        "< 5,000 基因 × < 100 样本可获得流畅体验",
        "学习 RNA-seq、小型试验数据集、教学演示",
    ],
    "本地分析（R/Python）": [
        "研究级分析；完整包访问",
        "完整的 DESeq2/edgeR/limma；完整包生态系统；可重现脚本",
        "需要 R/Python 配置；大型数据集在笔记本上可能较慢",
        "数千基因 × 数百样本；受 RAM 限制（通常 8–64 GB）",
        "标准 bulk RNA-seq；中等规模队列的 scRNA-seq；复现已发表流程",
    ],
    "HPC / 服务器": [
        "大规模基因组学；比对流程；多样本队列",
        "并行计算；大内存；FASTQ 比对；工作流管理器（Snakemake、Nextflow）",
        "需要 HPC 访问权限、SLURM 知识、环境管理",
        "数百到数千个样本；完整比对 + 定量流程",
        "FASTQ → BAM → 计数流程；大型 scRNA-seq 队列；空间转录组学",
    ],
}

st.dataframe(
    pd.DataFrame(comparison_data).set_index("特性"),
    use_container_width=True,
)

st.info(
    "本应用专为教学和中等规模数据集设计。"
    "对于真实队列的研究级分析，请在本地或 HPC 集群上使用 DESeq2/edgeR/limma-voom，"
    "配合适当的工作流管理。",
    icon="ℹ️",
)

st.divider()

col_n1, col_n2 = st.columns(2)
with col_n1:
    st.page_link("pages/09_Functional_Enrichment.py",
                 label="← 第9课：功能富集分析", icon="🧩")
with col_n2:
    st.page_link("pages/11_Single_Cell_RNAseq.py",
                 label="第11课：单细胞 RNA-seq →", icon="🔬")
