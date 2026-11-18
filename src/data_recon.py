# /// script
# requires-python = ">=3.12"
# dependencies = ["pandas==3.1.1", "openpyxl>=3.1.5", "matplotlib>=3.9"]
# ///
"""一次性侦察。仓库根目录执行 uv run --script src/data_recon.py。"""
import hashlib
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[1]
docs = root / "docs"
book = root / "data/raw/41588_2021_861_MOESM3_ESM.xlsx"
source = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-021-00861-8/MediaObjects/" + book.name
paper = "https://pmc.ncbi.nlm.nih.gov/articles/PMC8628368/"
code = "https://github.com/YichaoOU/ABE_NonCoding_functional_score/blob/cef8dbbd74e5f9359feb6cf709fccaac1127cd80/per_A_base_score"
guide = pd.read_excel(book, sheet_name="Supplementary_table _4", header=2)
base = pd.read_excel(book, sheet_name="Supplementary_table _6", header=2)
score = base["BPRSHbF"]
sites = base.Adenine.str.extract(r"^(?P<chrom>chr\d+):(?P<start>\d+)-(?P<end>\d+)(?P<strand>[+-])$")
assert sites.notna().all().all()
sites = sites.astype({"start": int, "end": int})
assert (sites.end - sites.start).eq(1).all() and not sites.duplicated(["chrom", "start"]).any()
loci = {"chr2": "BCL11A", "chr11": "HBG / β-globin cluster", "chr19": "KLF1", "chr6": "MYB–HBS1L"}
assert set(sites.chrom) == set(loci)
report = []
add = report.append


def table(headers, rows):
    add("| " + " | ".join(headers) + " |")
    add("| " + " | ".join(["---"] * len(headers)) + " |")
    for row in rows:
        add("| " + " | ".join(map(str, row)) + " |")
    add("")


add("# HbF 筛选数据侦察\n")
add(f"Cheng 等（2021）的出版社原版合订工作簿含表 4 和表 6。两页第 1 行是标题，第 2 行空白，第 3 行是字段名。下文只统计数据行。"
    f"[论文与方法]({paper}) · [下载原始 XLSX]({source})\n")
add(f"本地文件 `data/raw/{book.name}`，SHA-256 为 `{hashlib.sha256(book.read_bytes()).hexdigest()}`。"
    "重新下载到该路径后，运行 `uv run --script src/data_recon.py` 可重建报告及图；脚本联网读取 UCSC hg19 参考区段。\n")
add("## 工作表与缺失值\n")
meaning = {
    "sgRNA": "guide 的序列标识，每行一条",
    "Element": "筛选靶调控元件的区间标签，或对照类别；不是单碱基坐标",
    "HbFlow (normalized read counts)": "HbF 低组两个重复的归一化 guide 读数，以 / 分隔",
    "HbFhigh (normalized read counts)": "HbF 高组两个重复的归一化 guide 读数，以 / 分隔",
    "LFC(Log fold change)": "guide 在高组相对低组的对数倍数变化；正值表示高组富集",
    "FDR": "guide 层面差异丰度检验的多重检验校正值",
    "Adenine": "可编辑 A 的 hg19 单碱基区间及其所在链",
    "BPRSHbF": "该 A 的筛选统计分数；数值越大，支持差异富集的证据越强",
}
add(f"dtype 按 pandas {pd.__version__} 读取结果。缺失按空单元格或读入 NA 统计，数值 0 计作观测值。\n")
for i, df in [(4, guide), (6, base)]:
    add(f"### Supplementary_table _{i}（{len(df):,} 行）\n")
    table(["实际列名", "dtype", "字段解释", "缺失数", "缺失比例"], [
        [f"`{col}`", df[col].dtype, meaning[col], df[col].isna().sum(), f"{df[col].isna().mean():.2%}"]
        for col in df])
controls = guide.Element[~guide.Element.str.startswith("chr")].value_counts()
add(f"表 4 有 {guide.sgRNA.nunique():,} 条不同 guide，其中 {len(guide)-controls.sum():,} 行带染色体区间，"
    f"其余 {controls.sum()} 行为对照（" + "、".join(f"{k} × {v}" for k,v in controls.items()) + "）。"
    f"区间标签中 {guide.Element.str.fullmatch(r'chr\d+\.\d+\.\d+').sum()} 行采用点分格式。"
    "两个读数列应分别保留重复；0/0 表示两个重复均为零。\n")
add("## 位点、组装与等位基因\n")
add("表 6 没有 locus 列。按论文四个筛选区域与染色体一一对应分组，统计整个靶区，含远端元件。"
    "区间全部长 1 bp，按染色体与起点检查无重复。下表范围给出最外侧被测碱基的 hg19 一基坐标，中间并非连续覆盖。\n")
rows = []
for chrom, locus in loci.items():
    g = sites[sites.chrom == chrom]
    start, end = g.start.min(), g.end.max()
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom={chrom};start={start};end={end}"
    response = subprocess.check_output(["curl", "-fsSL", "--retry", "2", "--max-time", "60", url])
    ref = json.loads(response)
    assert (ref["genome"], ref["chrom"], ref["start"], ref["end"]) == ("hg19", chrom, start, end)
    dna = ref["dna"].upper()
    matched = sum(dna[r.start-start] == ("A" if r.strand == "+" else "T") for r in g.itertuples())
    rows.append([locus, chrom, len(g), f"{start+1:,}–{end:,}", f"[{matched}/{len(g)}]({url})"])
table(["筛选区域", "染色体", "独立位置数", "最小–最大位置", "正链参考核对"], rows)
add(f"判定为 **hg19（GRCh37）**。作者生成可编辑区间时调用 hg19.fa，并用 BED 起点和 gRNA 内偏移生成单碱基区间。"
    f"上表进一步逐行比较 UCSC hg19 序列与链标记。[组装来源]({code}/get_editable_A_positioin/main.py#L60) · "
    f"[坐标算法]({code}/get_editable_A_positioin/get_editable_bp_per_gRNA.py#L20-L39)\n")
add("给定的 GRCh38 BCL11A TSS 60,553,654 和 +58 kb DHS 60,494,920–60,496,208 符合负链基因下游坐标递减的关系。"
    "这两个锚点属于另一套组装，不能直接拿来筛表 6。BCL11A 的转录方向也不能代替每行 A 的链标记。\n")
add("两表均无独立 REF/ALT 列。表 6 的 `chrom:start-end±` 是零基半开区间。若表示为 hg19 正链、一基位置，"
    "`pos = start + 1 = end`；行末 + 对应 A→G，− 对应 T→C。要构造目标 GRCh38 的四元组，还需做组装转换，"
    "检查唯一的单碱基映射、映射方向（翻转时取互补）及 GRCh38 REF，并按模型接口的位置约定编码。这里仅核对了源组装。\n")
add("## BPRS 分布、符号与高分尾部\n")
table(["分位点", "BPRSHbF（线性插值）"], [[f"{p:.0%}", f"{score.quantile(p):.6f}"] for p in (0,.01,.05,.1,.25,.5,.75,.9,.95,.99,1)])
fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), layout="constrained")
for ax, bins, title in [(axes[0], np.arange(0,np.ceil(score.max()/5)*5+1,5), "Full range · bin width 5"),
                        (axes[1], np.linspace(0,10,41), "Detail [0, 10] · bin width 0.25")]:
    ax.hist(score, bins=bins, color="#306b85")
    ax.set(title=title, xlabel="BPRS (HbF)", ylabel="Number of sites", xlim=(bins[0],bins[-1]))
fig.savefig(docs / "bprs_hist.png", dpi=160)
plt.close(fig)
add("![BPRS 全范围直方图与低分区放大图](bprs_hist.png)\n")
add(f"{score.eq(0).sum():,} 个零值（{score.eq(0).mean():.2%}），{score.lt(0).sum()} 个负值。分布在零和低分区集中，右尾很长。"
    "在图示箱宽下看不到两个清晰分离的主体峰；这是图形描述，未做模态检验。\n")
add(f"作者从 guide FDR 出发，除以对应位置的经验编辑频率，上限截为 1；多个 guide 覆盖一个 A 时用经验 Brown 法合并，"
    f"最后取 −log10。计算不带入 LFC 符号，正 BPRS 既不能解释为 HbF 升高，也不能解释为降低。"
    f"截至 1 可产生零分；代码给原始 FDR 加 10⁻³⁰⁰，极端高分受数值处理影响。[分数实现]({code}/define_functional_score/get_functional_score.py#L46-L79)\n")
add("表 4 的 LFC 保留 guide 层面的方向，高组相对低组富集为正，耗减为负。举例取该列最大、最小值两行查看原始重复。\n")
table(["表 4 Excel 行", "HbF 低组", "HbF 高组", "LFC"], [
    [idx+4, guide.iloc[idx,2], guide.iloc[idx,3], guide.iloc[idx,4]]
    for idx in [guide.iloc[:,4].idxmax(),guide.iloc[:,4].idxmin()]])
add(f"你尚未指定高效应阈值。下表仅复核论文预测任务使用的高分界限 t = 30（含端点），不设本项目的分类阈值。"
    f"全表 ≥30 有 {score.ge(30).sum()} 个，≤10 有 {score.le(10).sum():,} 个，与论文方法中的两组规模一致。"
    f"分数绝对值等于原值，高分计数不能解释为上调或下调的变异数。[论文 Methods]({paper})\n")
table(["区域", "总数", "BPRS ≥ +30", "BPRS ≤ −30", "双尾合计"], [
    [loci[c], len(s), s.ge(30).sum(), s.le(-30).sum(), s.abs().ge(30).sum()]
    for c in ["chr2","chr11"] for s in [score[sites.chrom == c]]])
add("## 后续解释的边界\n")
add("HbF 高低分选是在细胞表型层做的测量；表 6 给每个碱基分配的分数，是从 guide 丰度和编辑窗口推断的证据强度。"
    "一个 guide 能覆盖多个 A，一个 A 也能被多个 guide 覆盖，相邻行可能共享筛选信号。"
    "直接把 guide 的 LFC 赋给每个 A 会引入未经验证的方向。"
    "因此，模型与 BPRS 排序一致时，所支持的是分子预测与 HbF 筛选信号的关联；单碱基 HbF 改变量、方向和中间调控步骤仍需各自的证据。\n")
(docs / "data_recon.md").write_text("\n".join(report), encoding="utf-8")
