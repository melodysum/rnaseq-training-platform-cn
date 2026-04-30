"""
pages/09_Functional_Enrichment.py
第9课 — 功能富集分析（GO / GSEA）

使用内置教学通路数据集的教育性富集分析。
不能替代 clusterProfiler、fgsea 或基于 MSigDB 的流程。
"""

import numpy as np
import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go

from utils.data_loader import init_session_data, load_demo_data
from utils.enrichment_utils import (
    TOY_GENE_SETS, run_ora, run_gsea_like, get_running_sum, demo_ranked_genes,
    GENESET_SOURCES, load_gene_sets, detect_gene_id_format,
    run_gsea_permutation, rank_by_statistic,
)

st.set_page_config(
    page_title="第9课 — 功能富集分析",
    page_icon="🧩",
    layout="wide",
)
init_session_data()

st.title("🧩 第9课 — 功能富集分析（GO / GSEA）")
st.caption(
    "使用内置教学通路数据集的教育性实现。"
    "不能替代 clusterProfiler、fgsea 或基于 MSigDB 的分析。"
)

# ── 第一节：概念说明 ──────────────────────────────────────────────────────────
with st.expander("📖 什么是功能富集分析？", expanded=True):
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("""
### 过度代表性分析（ORA）
ORA 回答的问题是：*我的差异表达基因列表中，来自某个生物学通路的基因
是否比随机情况下出现得更频繁？*

**工作原理：**
1. 获取你的显著差异表达基因列表。
2. 对数据库中的每条通路（如 GO、KEGG），统计有多少通路基因出现在你的列表中。
3. 使用 Fisher's 精确检验，判断该重叠是否比背景基因组所预期的更大。
4. 进行多重检验校正（Benjamini-Hochberg FDR）。

**优势：** 简单、快速、可解读。  
**局限：** 对所有显著基因同等对待——忽略倍数变化的大小。
背景基因组的选择会显著影响结果。
        """)
    with col2:
        st.markdown("""
### 基因集富集分析（GSEA）
GSEA 回答的问题是：*即使单个基因未通过显著性截断，
某条通路的基因是否倾向于出现在排名基因列表的顶部（或底部）？*

**工作原理：**
1. 按某统计量（如 log2FC、带符号的 p 值）对**所有**基因进行排名。
2. 沿排名列表向下走；当遇到通路中的基因时向上走，否则向下走。
3. 富集分数（ES）是这个累积求和的最大偏差。
4. 基因集中在顶部的通路为正向富集；集中在底部的为负向富集。

**优势：** 使用完整基因列表——比 ORA 更灵敏。  
**局限：** 解读更复杂；基于置换检验的显著性计算量大（教学版使用更简单的评分方式）。
        """)
    st.info("""
**为什么基因组（universe）很重要？**  
在 ORA 中，背景（universe）定义了"随机情况下的期望"是什么。
只使用你实际测量到的基因（而非所有 20,000 个人类基因）可以得到更真实的富集统计量。
过窄的 universe 会夸大，过宽的 universe 会低估富集 p 值。
    """)

st.divider()

# ── 第二节：数据来源 ──────────────────────────────────────────────────────────
st.subheader("📂 数据来源")

counts   = st.session_state["counts"]
metadata = st.session_state["metadata"]
source   = st.session_state.get("data_source", "demo")
de_res   = st.session_state.get("de_results", None)

have_de = (
    de_res is not None
    and isinstance(de_res, pd.DataFrame)
    and "log2FC" in de_res.columns
    and len(de_res) > 0
)

if have_de:
    st.success(
        f"✅ 已在会话中找到差异表达结果（{len(de_res):,} 个基因）。"
        "可在下方将其用作基因列表来源。",
        icon="🧪",
    )
else:
    st.info(
        "会话状态中未找到差异表达结果。"
        f"将使用{'内置演示' if source == 'demo' else '已上传'}数据中的 log2FC "
        "生成排名基因列表。",
        icon="📊",
    )
    try:
        de_res = demo_ranked_genes(counts, metadata)
        have_de = True
        st.caption(
            f"演示排名：{de_res['group2'].iloc[0]} vs {de_res['group1'].iloc[0]} "
            f"（{len(de_res):,} 个基因按 log2FC 排名）。"
        )
    except Exception as e:
        st.error(f"无法生成演示排名：{e}")
        have_de = False

# ── 基因 ID 格式检测 ─────────────────────────────────────────────────────────
gene_fmt = detect_gene_id_format(counts.index)
if gene_fmt == "ensembl":
    st.warning(
        "⚠️ 你的计数矩阵似乎使用的是 **Ensembl ID**（ENSG...）。"
        "MSigDB GMT 文件使用 HGNC 基因符号——对真实通路数据库的富集分析将几乎无法匹配。"
        "请先进行基因 ID 转换。"
    )
elif gene_fmt == "symbol":
    st.success(
        "✅ 基因 ID 看起来是 **HGNC 基因符号**——与 MSigDB GMT 文件兼容。",
        icon="🧬",
    )

# ── 通路数据库选择器 ───────────────────────────────────────────────────────────
st.markdown("#### 🗂️ 选择通路数据库")
st.caption(
    "**教学通路集**仅用于演示 ORA/GSEA 原理，不具备生物学验证意义。"
    "**MSigDB Hallmark** 包含 50 个简洁、内聚的生物学过程特征基因集。"
    "**MSigDB C2:CP** 包含 4115 条来自 Reactome、WikiPathways、KEGG 等的标准通路。"
)

gs_label = st.radio(
    "通路数据库",
    options=list(GENESET_SOURCES.keys()),
    index=0,
    horizontal=False,
    label_visibility="collapsed",
)
gs_key = GENESET_SOURCES[gs_label]

try:
    active_gene_sets, gs_display_name, n_gs = load_gene_sets(gs_key)
    if gs_key == "toy":
        st.caption(
            f"📚 **{gs_display_name}** — {n_gs} 条通路，仅供教学。"
        )
    else:
        st.caption(
            f"🧬 **{gs_display_name}** — 已从本地 GMT 文件加载 {n_gs:,} 条基因集。"
        )
        if gs_key == "c2cp":
            filter_prefix = st.text_input(
                "按前缀过滤通路（可选）",
                value="",
                placeholder="例：REACTOME_ 或 KEGG_ 或 WP_",
                help="留空则使用全部 C2:CP 通路。输入前缀可聚焦于某一子集。"
            )
            if filter_prefix.strip():
                prefix = filter_prefix.strip().upper()
                active_gene_sets = {k: v for k, v in active_gene_sets.items()
                                    if k.upper().startswith(prefix)}
                st.caption(f"已过滤至 **{len(active_gene_sets):,}** 条匹配 '{prefix}' 的通路。")
except FileNotFoundError as e:
    st.error(str(e))
    active_gene_sets = TOY_GENE_SETS
    gs_display_name  = "教学通路集（后备）"
    n_gs = len(TOY_GENE_SETS)

st.divider()

if not have_de:
    st.error("没有基因列表无法运行富集分析。请先完成第7课——差异表达分析。")
    st.stop()

# ── 第三节：基因列表来源选择 ─────────────────────────────────────────────────
st.subheader("🎯 选择基因列表来源")

gene_list_options = []
if have_de and "log2FC" in de_res.columns:
    gene_list_options += [
        "显著差异表达基因（padj < 0.05，|log2FC| ≥ 1）",
        "前100个上调基因",
        "前100个下调基因",
    ]
gene_list_options.append("手动输入基因列表（在下方输入）")

source_choice = st.radio("基因列表来源：", gene_list_options, horizontal=True)

if source_choice == "显著差异表达基因（padj < 0.05，|log2FC| ≥ 1）":
    if "padj" in de_res.columns and "significant" in de_res.columns:
        sig = de_res[de_res["significant"]].copy()
    elif "padj" in de_res.columns:
        sig = de_res[(de_res["padj"] < 0.05) & (de_res["log2FC"].abs() >= 1)].copy()
    else:
        sig = de_res[de_res["log2FC"].abs() >= 1].head(200).copy()
    gene_list_for_ora = sig.index.tolist() if sig.index.name == "gene" else sig["gene"].tolist() if "gene" in sig.columns else sig.index.tolist()
    st.caption(f"已选择 {len(gene_list_for_ora)} 个显著差异表达基因。")

elif source_choice == "前100个上调基因":
    top = de_res.sort_values("log2FC", ascending=False).head(100)
    gene_list_for_ora = top.index.tolist() if top.index.name == "gene" else top["gene"].tolist() if "gene" in top.columns else top.index.tolist()
    st.caption(f"按 log2FC 排名前 {len(gene_list_for_ora)} 个上调基因。")

elif source_choice == "前100个下调基因":
    top = de_res.sort_values("log2FC", ascending=True).head(100)
    gene_list_for_ora = top.index.tolist() if top.index.name == "gene" else top["gene"].tolist() if "gene" in top.columns else top.index.tolist()
    st.caption(f"按 log2FC 排名前 {len(gene_list_for_ora)} 个下调基因。")

else:
    manual_input = st.text_area(
        "输入基因名（每行一个或用逗号分隔）：",
        placeholder="STAT1\nMX1\nISG15\nIFIT1",
        height=120,
    )
    gene_list_for_ora = [
        g.strip().upper()
        for g in manual_input.replace(",", "\n").splitlines()
        if g.strip()
    ]
    st.caption(f"已输入 {len(gene_list_for_ora)} 个基因。")

if len(gene_list_for_ora) < 3:
    st.warning("请至少选择或输入 3 个基因才能运行富集分析。")
    st.stop()

if "gene" in de_res.columns:
    ranked_series = de_res.set_index("gene")["log2FC"].sort_values(ascending=False)
elif de_res.index.name == "gene" or de_res.index.dtype == object:
    ranked_series = de_res["log2FC"].sort_values(ascending=False)
else:
    ranked_series = de_res["log2FC"].sort_values(ascending=False)

lfc_weights = ranked_series.abs()

st.divider()

# ── 第四节：ORA ──────────────────────────────────────────────────────────────
st.subheader("📊 过度代表性分析（ORA）")
st.caption(
    "Fisher's 精确检验 + Benjamini-Hochberg FDR 校正。"
    "背景基因组 = 当前数据集中测量到的所有基因。"
)

universe = counts.index.tolist()

with st.spinner("正在运行 ORA…"):
    ora_res = run_ora(gene_list_for_ora, active_gene_sets, universe=universe)

if ora_res.empty:
    st.warning("未找到通路重叠。请尝试更大的基因列表或不同的来源。")
else:
    display_cols = ["pathway", "overlap", "pathway_size", "query_size",
                    "gene_ratio", "pvalue", "padj", "overlap_genes"]
    display_cols = [c for c in display_cols if c in ora_res.columns]
    st.dataframe(
        ora_res[display_cols].style.format({
            "gene_ratio": "{:.3f}",
            "pvalue":     "{:.4f}",
            "padj":       "{:.4f}",
        }),
        use_container_width=True,
        height=280,
    )

    plot_df = ora_res.copy()
    plot_df["neg_log10_padj"] = -np.log10(plot_df["padj"].clip(lower=1e-10))
    plot_df["label"] = plot_df["pathway"].str.replace("_", " ").str.title()
    plot_df = plot_df.sort_values("neg_log10_padj", ascending=True).tail(10)

    fig_ora = px.bar(
        plot_df,
        x="neg_log10_padj",
        y="label",
        orientation="h",
        color="gene_ratio",
        color_continuous_scale="Blues",
        labels={
            "neg_log10_padj": "-log₁₀（校正 p 值）",
            "label":          "通路",
            "gene_ratio":     "基因比例",
        },
        title="ORA 富集最高通路 — 教学通路数据集",
    )
    fig_ora.add_vline(x=-np.log10(0.05), line_dash="dash", line_color="red",
                      annotation_text="padj = 0.05")
    fig_ora.update_layout(height=400)
    st.plotly_chart(fig_ora, use_container_width=True)
    if gs_key == "toy":
        st.caption(
            "⚠️ 这些是教学用教学通路，不是真正的 GO 或 KEGG 数据库。"
            "结果仅作说明之用。"
        )
    else:
        st.info(
            f"ℹ️ 当前使用 **{gs_display_name}**。"
            "若没有通路通过 FDR < 0.05，这很可能是因为当前 demo 数据为模拟数据，"
            "不包含能够富集这些真实通路的生物学信号。"
            "如需有生物学意义的结果，请上传真实实验数据（例如结核病数据集 GSE167232）。"
            "这**不代表**代码存在错误。"
        )

st.divider()

# ── 第五节：GSEA 类排序富集 ──────────────────────────────────────────────────
st.subheader("📈 排序富集分析（GSEA 类）")
st.caption(
    "基因按 log2FC 排名。累积求和富集分数沿排名列表向下走。"
    "标记为**教学性 GSEA 类实现**——非 GSEA 置换检验。"
)

with st.spinner("正在运行 GSEA 类评分…"):
    gsea_res = run_gsea_like(ranked_series.index, active_gene_sets, lfc_weights)

if gsea_res.empty:
    st.warning("未产生富集结果。")
else:
    display_cols_g = ["pathway", "ES", "direction", "n_hits", "pathway_size"]
    display_cols_g = [c for c in display_cols_g if c in gsea_res.columns]

    gsea_display = gsea_res[display_cols_g].copy()
    gsea_display["pathway"] = gsea_display["pathway"].str.replace("_", " ").str.title()

    st.dataframe(
        gsea_display.style.format({"ES": "{:.4f}"}),
        use_container_width=True,
        height=300,
    )

    st.markdown("#### 富集图 — 选择一条通路")
    st.markdown(
        "富集图展示了沿排名基因列表向下走时的累积求和。"
        "在顶部出现峰值表明通路基因集中在上调程度最高的基因中。"
    )

    pathway_options = gsea_res["pathway"].tolist()
    pathway_labels  = [p.replace("_", " ").title() for p in pathway_options]
    selected_label  = st.selectbox(
        "选择要绘图的通路：", pathway_labels, index=0
    )
    selected_pathway = pathway_options[pathway_labels.index(selected_label)]
    pathway_genes    = active_gene_sets[selected_pathway]

    rs = get_running_sum(ranked_series.index, pathway_genes, lfc_weights)
    ranks = list(range(1, len(rs) + 1))

    fig_es = go.Figure()
    fig_es.add_trace(go.Scatter(
        x=ranks, y=rs,
        mode="lines",
        name="累积富集分数",
        line=dict(color="#2563eb", width=2),
    ))
    fig_es.add_hline(y=0, line_dash="dash", line_color="gray", line_width=1)

    hit_positions = [i + 1 for i, g in enumerate(ranked_series.index)
                     if g in set(pathway_genes)]
    for pos in hit_positions:
        fig_es.add_vline(x=pos, line_color="rgba(220,38,38,0.25)", line_width=0.8)

    rs_array = np.array(rs)
    if abs(rs_array.max()) >= abs(rs_array.min()):
        peak_idx = int(np.argmax(rs_array))
    else:
        peak_idx = int(np.argmin(rs_array))
    fig_es.add_trace(go.Scatter(
        x=[ranks[peak_idx]], y=[rs_array[peak_idx]],
        mode="markers",
        name=f"ES = {rs_array[peak_idx]:.3f}",
        marker=dict(color="red", size=10, symbol="diamond"),
    ))

    fig_es.update_layout(
        title=(
            f"富集图：{selected_label}<br>"
            "<sup>教学性 GSEA 类实现 — 演示数据</sup>"
        ),
        xaxis_title="基因排名（log2FC 降序）",
        yaxis_title="累积富集分数",
        height=420,
        legend=dict(orientation="h", y=-0.2),
    )
    st.plotly_chart(fig_es, use_container_width=True)

    es_val = gsea_res[gsea_res["pathway"] == selected_pathway]["ES"].values[0]
    direction = gsea_res[gsea_res["pathway"] == selected_pathway]["direction"].values[0]
    n_hits_val = gsea_res[gsea_res["pathway"] == selected_pathway]["n_hits"].values[0]

    st.info(
        f"**{selected_label}** — ES = {es_val:.4f} | 方向：{direction} | "
        f"在排名列表中找到的通路基因：{n_hits_val}/{len(pathway_genes)}\n\n"
        "竖向红线显示通路基因在排名列表中的位置。"
        "集中在顶部或底部的簇会产生清晰的 ES 峰值。"
    )

st.divider()

# ── 第六节：解读说明 ──────────────────────────────────────────────────────────
st.subheader("💡 解读富集分析结果")
st.markdown("""
**ORA 结果：**
- 较低的校正 p 值表明该通路在你的基因列表中过度代表。
- 基因比例（重叠数 / 通路大小）表示有多少比例的通路成员在你的差异表达基因列表中——
  即使重叠数较少，小通路的高比例仍可能有意义。
- 始终检查重叠基因——有时少数高连接的枢纽基因会驱动多条通路的富集。

**GSEA 类结果：**
- 正 ES 表示通路基因集中在排名列表顶部（在目标组中与上调相关）。
- 负 ES 表示通路基因集中在底部（下调）。
- 富集图的形状很重要：接近第1名的尖锐峰值比列表中部的平缓峰值更有力。

**一般注意事项：**
- 这些教学通路仅用于教育目的。真实分析需要 GO、KEGG、Reactome 或 MSigDB 数据库。
- 富集分析对基因组和排名统计量高度敏感。
- 生物学合理性应该指导解读，而不只是 p 值。
""")

st.divider()

# ══════════════════════════════════════════════════════════════════════════════
# 进阶：基于置换检验的 GSEA（NES + p 值）
# ══════════════════════════════════════════════════════════════════════════════

from utils.enrichment_utils import run_gsea_permutation, rank_by_statistic
from utils.de_analysis import run_de

st.subheader("🧬 进阶：基于置换检验的 GSEA（NES + p 值）")

with st.expander("📖 与上方 GSEA-like 分析的区别", expanded=True):
    st.markdown("""
| | GSEA-like（上方） | 置换检验 GSEA（本节） |
|---|---|---|
| **基因排序指标** | log2FC | 检验统计量（t 统计量）|
| **零假设分布** | 无 | 基因标签置换 |
| **得分** | 仅 ES | ES + **NES**（归一化富集分数）|
| **p 值** | 无 | 置换检验 |
| **FDR 校正** | 无 | Benjamini-Hochberg |
| **跨通路可比性** | ❌ | ✅ |

**为什么用检验统计量排序而不是 log2FC？**  
t 统计量同时整合了效应量和方差。  
3 个样本 logFC = 3 的基因，其 t 统计量低于 20 个样本 logFC = 2 的基因——后者更可靠。  
仅用 log2FC 排序会把噪声大、样本少的基因错误地推至顶端。

**为什么 ORA 和 GSEA 结果不同？**  
ORA 丢弃了所有未达显著性阈值的基因。若某通路有 12 个基因都接近 FDR 0.05 但未超过，  
ORA 找不到这条通路；而 GSEA 能看到这 12 个基因都集中在排序列表顶部，从而正确识别该通路。
""")

with st.spinner("正在运行差异表达分析以生成 GSEA 排序…"):
    de_for_gsea = run_de(counts, metadata, fdr_cutoff=0.05, lfc_cutoff=0.5)

if "stat" not in de_for_gsea.columns:
    de_for_gsea["stat"] = de_for_gsea["log2FC"] / (
        de_for_gsea["log2FC"].abs().mean() /
        (-np.log10(de_for_gsea["pvalue"].clip(1e-300))).clip(lower=0.01)
    ).clip(lower=0.01)

ranked_for_perm = rank_by_statistic(
    de_for_gsea,
    stat_col="stat",
    lfc_col="log2FC",
    pval_col="pvalue",
    method="signed_log",
)

col_pg1, col_pg2 = st.columns([2, 1])
with col_pg1:
    n_perm = st.select_slider(
        "置换次数",
        options=[100, 500, 1000],
        value=500,
        help="置换次数越多，p 值越稳定。标准为 1000 次。",
    )
with col_pg2:
    min_set = st.number_input("通路最小基因数", 3, 20, 5, key="min_set_perm")

run_perm = st.button("▶ 运行置换检验 GSEA", type="primary")

if run_perm or "perm_gsea_result" in st.session_state:
    if run_perm:
        with st.spinner(f"正在运行 {n_perm} 次置换…"):
            perm_result = run_gsea_permutation(
                ranked_for_perm,
                gene_sets=active_gene_sets,
                n_permutations=n_perm,
                random_state=42,
                min_set_size=min_set,
            )
            st.session_state["perm_gsea_result"] = perm_result
    else:
        perm_result = st.session_state["perm_gsea_result"]

    if perm_result.empty:
        st.warning("没有通路通过最小基因数筛选。")
    else:
        display_perm = perm_result[[
            "pathway", "n_genes_in_list", "ES", "NES",
            "pvalue", "padj", "direction", "leading_edge_genes"
        ]].copy()
        display_perm["pathway"] = display_perm["pathway"].str.replace("_", " ").str.title()

        def _style_sig(row):
            if row["padj"] < 0.05:
                return ["background-color: #dcfce7"] * len(row)
            return [""] * len(row)

        st.dataframe(
            display_perm.style
                .apply(_style_sig, axis=1)
                .format({"ES": "{:.4f}", "NES": "{:.4f}",
                         "pvalue": "{:.4f}", "padj": "{:.4f}"}),
            use_container_width=True,
            height=350,
            hide_index=True,
        )

        sig_perm = (perm_result["padj"] < 0.05).sum()
        st.info(
            f"置换检验后共有 **{sig_perm}** 条通路 FDR < 0.05（绿色高亮）。"
            "NES > 0 = 通路在目标组中上调；NES < 0 = 下调。"
        )

        st.markdown("#### NES 条形图")
        fig_nes = go.Figure()
        colors = ["#16a34a" if v > 0 else "#dc2626" for v in perm_result["NES"]]
        pway_labels = [p.replace("_", " ").title() for p in perm_result["pathway"]]
        fig_nes.add_trace(go.Bar(
            y=pway_labels,
            x=perm_result["NES"],
            orientation="h",
            marker_color=colors,
            text=[f"padj={p:.3f}" for p in perm_result["padj"]],
            textposition="outside",
        ))
        fig_nes.add_vline(x=0, line_color="black", line_width=1)
        fig_nes.update_layout(
            xaxis_title="NES（归一化富集分数）",
            yaxis=dict(autorange="reversed"),
            height=max(300, len(perm_result) * 35),
        )
        st.plotly_chart(fig_nes, use_container_width=True)
else:
    st.info("点击 **▶ 运行置换检验 GSEA** 以计算 NES 和置换 p 值。")

st.divider()

col_n1, col_n2 = st.columns(2)
with col_n1:
    st.page_link("pages/08_Clustering_Heatmaps.py",
                 label="← 第8课：聚类与热图", icon="🗺️")
with col_n2:
    st.page_link("pages/10_Public_Data.py",
                 label="第10课：公共数据与可重复性 →", icon="🌐")
