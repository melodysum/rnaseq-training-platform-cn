# 🧬 RNA-seq 互动训练平台（中文版）

一个交互式、循序渐进的 Streamlit 网页应用，用于学习 RNA-seq 数据分析——从原始计数到生物学解读。专为初学转录组分析的学生、研究人员和生物信息学工作者设计。

> **教学声明：** 本应用不能替代 DESeq2、edgeR、limma-voom、clusterProfiler、fgsea 或完整的 GSEA 流程。所有分析均以透明的教学模型实现，旨在教授核心概念。对于研究级分析，请使用经过验证的 R/Python 生物信息学工具。

---

## 📋 更新日志

### v2.1 — 2026-04-30

**第9课 — 真实通路数据库支持**

*目的：* 此前实现使用的是教学用玩具通路集，概念演示有效，但缺乏真实生物学意义。本次更新引入 MSigDB 标准通路数据库，使第9课的富集分析结果具备真正的生物学解读价值，而不仅仅是方法论示范。

- 第9课新增通路数据库选择器（单选：教学通路 / Hallmark / C2:CP）
- utils/enrichment_utils.py 新增 parse_gmt() 和 load_gene_sets() — 解析本地 MSigDB GMT 文件为基因集字典
- 新增 detect_gene_id_format() — 自动检测输入基因 ID 格式（HGNC 符号或 Ensembl ID），若不匹配 MSigDB 则发出警告
- data/gene_sets/ 内置两个 MSigDB GMT 文件：

  h.all.v2026.1.Hs.symbols.gmt — MSigDB Hallmark 基因集
  - 是什么：包含 50 个基因集的精选库，每个基因集代表一种特定、明确的生物学状态或过程（例如：干扰素应答、炎症反应、氧化磷酸化、缺氧）。Hallmark 基因集由计算方法从 MSigDB 其他数千个基因集中提取重叠基因衍生而来，仅保留协调表达的基因，具有低冗余、高内聚的特点。
  - 为什么选它：Hallmark 是 RNA-seq 富集分析的标准起点。50 条通路数量小、结果易于解读，适合教学演示和数据集的初步生物学探索，在结核病和免疫学研究中被广泛使用。
  - 来源：Broad Institute 分子特征数据库（MSigDB），人类集合 H。msigdb.org
  - 版本：v2026.1.Hs（人类，HGNC 基因符号）
  - 下载时间：2026-04-30

  c2.cp.v2026.1.Hs.symbols.gmt — MSigDB C2:CP 标准通路
  - 是什么：包含 4115 个基因集，来自多个权威通路数据库的标准生物学通路集合。C2:CP 子集包括：Reactome（1839 条，详细机制通路）、WikiPathways（925 条，社区整理）、KEGG MEDICUS（658 条，代谢与信号通路）、BioCarta（292 条）、PID（196 条）。与 Hallmark 不同，C2:CP 通路直接来源于专家整理的通路数据库，而非从表达数据中计算衍生。
  - 为什么选它：C2:CP 提供了 Hallmark 所不具备的机制深度。Hallmark 能告知干扰素应答被激活，C2:CP 能进一步区分具体是 REACTOME_INTERFERON_GAMMA_SIGNALING 还是 REACTOME_INTERFERON_ALPHA_BETA_SIGNALING。这种分辨率对结核病研究尤为重要，前缀过滤功能允许聚焦于单一数据库（如仅用 REACTOME_）以避免冗余。
  - 来源：Broad Institute 分子特征数据库（MSigDB），人类集合 C2，子集 CP。msigdb.org
  - 版本：v2026.1.Hs（人类，HGNC 基因符号）
  - 下载时间：2026-04-30

- C2:CP 模式支持前缀过滤（如 REACTOME_、KEGG_、WP_），无需修改代码即可聚焦特定子集
- 所有 ORA、GSEA-like 和置换检验 GSEA 均自动使用当前选中数据库
- 教学通路集保留并明确标注为仅供教学使用
- 新增诊断提示：当使用真实数据库但 demo 数据无显著富集结果时，说明这是预期行为而非代码错误
- 选择 C2:CP 时新增运行时间提示（4,115 条基因集，预计需 30-60 秒）
- 首页上传区域新增清除缓存说明（Mac：Cmd+Shift+R，Windows：Ctrl+Shift+R）

### v2.0 — 2026-04-27

**新增工具模块**

- `utils/stats_utils.py` — 统一 DE 输出格式，在代码注释中详细说明 p 值与 FDR 的区别，以及为何需要多重检验校正
- `utils/kallisto_import.py` — Kallisto 伪比对输出 → 基因级整数计数矩阵（长度缩放 TPM 方法，对应 R 中的 tximport）；解决了伪比对工作流与计数模型（DESeq2 / edgeR）之间的衔接问题

**扩展现有模块** *（所有原有函数均保留）*

- `utils/filtering.py` — 新增 `filter_by_expr()`（edgeR filterByExpr 等价实现：阈值由文库大小和分组大小推导，不依赖 DE 结果）和 `threshold_sweep_retained()`（诊断性阈值扫描，展示不同阈值下保留基因数，无循环论证）
- `utils/pca_utils.py` — 新增 `variance_decomposition()`（类 PVCA 方法：通过加权偏 R² 量化批次、实验条件等因素对总方差的贡献百分比）
- `utils/enrichment_utils.py` — 新增 `rank_by_statistic()`（用检验统计量对基因排序，用于 GSEA 输入）和 `run_gsea_permutation()`（完整 GSEA：基因标签置换、NES 归一化、置换 p 值、BH-FDR 校正）

**更新页面** *（所有原有内容均保留）*

- **第3课 — 低表达基因过滤**：新增"数据驱动阈值（filterByExpr 方法）"模块——解释为何用 DE 结果优化阈值是循环论证；展示阈值 vs 保留基因数的诊断图
- **第5课 — 探索性分析与 PCA**：新增"方差分解"模块——各因素方差贡献柱状图与表格；辅助判断是否需要批次校正
- **第9课 — 功能富集分析**：新增"基于置换检验的 GSEA（NES + p 值）"模块——置换 NES、p 值、FDR；ORA vs GSEA 对比表；NES 条形图

---

## 🚀 快速开始

### 本地运行

```bash
git clone https://github.com/melodysum/rnaseq-training-platform-cn.git
cd rnaseq-training-platform-cn
pip install -r requirements.txt
streamlit run Home.py
```

### 在线访问
[rnaseq-training-platform-cn.streamlit.app](https://rnaseq-training-platform-cn.streamlit.app/)

---

## 📚 课程内容

| # | 课程 | 状态 |
|---|------|------|
| 1 | RNA-seq 基础与实验设计 | ✅ 已开放 |
| 2 | 定量、导入与注释 | ✅ 已开放 |
| 3 | 低表达基因过滤 | ✅ 已开放 |
| 4 | 多重检验与 FDR | ✅ 已开放 |
| 5 | 探索性分析与 PCA | ✅ 已开放 |
| 6 | 批次效应校正 | ✅ 已开放 |
| 7 | 差异表达分析 | ✅ 已开放 |
| 8 | 聚类与热图 | ✅ 已开放 |
| 9 | 功能富集分析（GO / GSEA） | ✅ 已开放 |
| 10 | 公共数据与可重复性分析 | ✅ 已开放 |
| 11 | 单细胞 RNA-seq（入门模块） | ✅ 已开放 |

---

## 📁 输入文件格式

本应用接受两个 CSV 文件：

### counts.csv — 原始基因水平计数矩阵

- **基因为行**，样本为列
- 第一列 = 基因标识符（用作行索引）
- 数值必须是**原始整数计数**——不能是 TPM、FPKM 或 CPM
- 至少需要 2 个样本列和 1 行基因

```
gene_symbol, D01_control, D01_treatment, D02_control, D02_treatment
IFNA6,       25,          40,            15,          48
CAD,          0,           9,            12,          10
```

### metadata.csv — 样本注释表

- **样本为行**，注释为列
- 样本名必须与 counts.csv 的列名完全匹配（区分大小写）
- 至少需要 2 个样本

```
sample_name,   groupA,    donor, batch,  sex, age
D01_control,   control,   D01,   batch1, F,   37
D01_treatment, treatment, D01,   batch1, F,   37
D02_control,   control,   D02,   batch2, F,   40
D02_treatment, treatment, D02,   batch2, F,   40
```

---

## 🏷️ 元数据列要求

| 列名 | 是否必需 | 用途 |
|------|----------|------|
| `groupA` | **必需** | 定义所有分析的比较分组 |
| `donor` | 可选 | 启用配对差异表达分析（第7课） |
| `batch` | 可选 | 启用批次校正流程（第6课） |
| `sex`、`age` 等 | 可选 | 描述性信息；在数据集摘要中显示 |

**验证规则：**
- 两个文件中的重复样本名将被拒绝
- 重复的基因标识符将被拒绝
- 非数值计数将被拒绝
- 负值将被拒绝（计数必须 ≥ 0）
- 两个文件中的样本名必须完全匹配

---

## 🔬 功能范围

### 第1–8课（Bulk RNA-seq 核心）
- 原始计数加载与演示数据回退
- 低表达基因过滤
- 多重检验 / FDR 校正演示
- 基于 PCA 的探索性分析
- 简单线性模型批次校正
- 教学版差异表达：log-CPM + 配对/非配对 t 检验 + BH FDR
  - 配对差异表达在检验前按供体对齐样本
  - 若供体匹配不足，回退到 Welch t 检验
- 层次聚类和表达热图

### 第9课 — 功能富集分析
- 内置教学玩具通路库（10 条通路，真实人类基因名）
- 过度代表性分析（ORA）：Fisher's 精确检验 + BH 校正
- GSEA 类累积求和富集评分
- 富集图：累积求和 vs 基因排名，通路命中标记
- 使用会话中的差异表达结果或演示数据运行（若无前序差异表达结果）

### 第10课 — 公共数据与可重复性
- 真实 GEO 案例研究（GSE148036、GSE167232）本地存储
- 交互式可重复性核查清单
- 计数 vs TPM/FPKM/CPM 说明
- 网页端 vs 本地 vs HPC 对比表
- 元数据规范指导

### 第11课 — 单细胞 RNA-seq（入门）
- 概念介绍：细胞 vs 样本、稀疏性、UMAP、聚类、标记基因
- 模拟 UMAP 散点图（7 个彩色教学簇）
- 每簇示例标记基因表
- Bulk vs scRNA-seq 结构化对比表
- 明确的"本应用未实现"说明

---

## ⚠️ 教学局限性

本应用专为**教学和中等规模数据集**设计：

- 统计模型已简化（对 log-CPM 进行 t 检验，而非负二项 GLM）
- 不能替代 DESeq2、edgeR 或 limma-voom 进行差异表达分析
- 功能富集使用教学用玩具通路，不是真正的 GO/KEGG/MSigDB 数据库
- 第11课使用模拟 UMAP 数据；未实现真正的 scRNA-seq 流程
- 会话之间无文件持久化
- 大型数据集（>5,000 基因 × 200+ 样本）可能运行缓慢或超出浏览器内存
- 交互式聚类和热图应限制在合理的基因数量内

对于研究级分析，请使用：
- **差异表达：** DESeq2、edgeR、limma-voom
- **富集分析：** clusterProfiler、fgsea、GSEA desktop
- **单细胞：** Scanpy（Python）、Seurat（R）
- **比对：** STAR、HISAT2、Salmon + tximeta
- **工作流：** HPC 上的 Snakemake、Nextflow

---

## 📂 项目结构

```
rnaseq-training-platform-cn/
├── Home.py                          # 首页 + 数据上传
├── requirements.txt
├── data/
│   ├── counts.csv                   # 内置演示计数矩阵（10,000 基因 × 40 样本）
│   └── metadata.csv                 # 内置演示元数据（配对设计，20 个供体）
├── pages/
│   ├── 01_RNAseq_Foundations.py
│   ├── 02_Quantification_Import_Annotation.py
│   ├── 03_Filtering.py
│   ├── 04_FDR.py
│   ├── 05_Exploratory_Analysis_PCA.py
│   ├── 06_Batch_Correction.py
│   ├── 07_Differential_Expression.py
│   ├── 08_Clustering_Heatmaps.py
│   ├── 09_Functional_Enrichment.py
│   ├── 10_Public_Data.py
│   └── 11_Single_Cell_RNAseq.py
└── utils/
    ├── data_loader.py
    ├── de_analysis.py
    ├── enrichment_utils.py
    ├── batch_effects.py
    ├── clustering_utils.py
    ├── exploration.py
    ├── fdr_demo.py
    ├── filtering.py
    ├── pca_utils.py
    └── simulation.py
```

---

## 📄 许可证

MIT License。详见 LICENSE 文件。

---

## 🔗 相关链接

- **英文原版：** [rnaseq-training-platform.streamlit.app](https://rnaseq-training-platform.streamlit.app/)
- **英文 GitHub：** [melodysum/rnaseq-training-platform](https://github.com/melodysum/rnaseq-training-platform)
