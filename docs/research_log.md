# Research Log — Phase 1: AlphaGenome 2 BCL11A 增强子变异效应复现

本文档记录 Phase 1（复现 Google DeepMind 2026 AlphaGenome 2 技术报告 Figure 4j–k 中 `BCL11A` +58 kb 红系增强子变异 `chr2:60495266:G>A` 预测效应）的全部关键科研与工程决策。每项决策均按统一结构说明：**选择了什么、为什么选择、依据来源、替代方案、独立验证方法**，并明确记录原论文、官方教程与已安装 SDK（`alphagenome==2.0.0`）之间的差异及不确定性。

---

## 决策总览与论文-教程差异对照表

在编写分析代码前，我们核对了四份一手材料：
1. AlphaGenome 2 技术报告（`alphagenome-2.pdf`，重点 p. 10 Figure 4h–k、p. 11 正文、pp. 19–35 模型与评分方法、pp. 47–48 *Locus analysis at BCL11A*、pp. 52–53 数据与代码可用性）；
2. AlphaGenome v2 官方 Getting Started 文档；
3. 官方 `enhancer_to_gene_bcl11a.ipynb` 教程；
4. 本地虚拟环境中已安装的 `alphagenome==2.0.0` 源码（`.venv/Lib/site-packages/alphagenome/`）。

核查发现原论文 Figure 4j–k 与官方 `enhancer_to_gene_bcl11a` 教程存在 3 处实质差异，均已**优先按原论文设定复现**并完整记录：

| 维度 | AlphaGenome 2 论文（Fig. 4j–k & pp. 47–48） | 官方 v2 `enhancer_to_gene_bcl11a` 教程 | 本复现采用的设定 |
| :--- | :--- | :--- | :--- |
| 分析内容 | 包含 E2G（Fig. 4h–i）、`chr2:60495266:G>A` 饱和 ISM（Fig. 4j）与 REF/ALT `scATAC` 预测对比（Fig. 4k） | 仅演示 E2G 输入梯度关联（复现 Fig. 4h–i），未包含变异推理、评分或 ISM | 严格按论文 pp. 47–48 与官方 `predict_variant` / `score_variant` / `score_ism_variants` 接口实现 Fig. 4j–k |
| 对照 B 细胞条件 | **外周血初始 B 细胞** `blood naive B cell`（`UBERON:0000178_CL:0000788`） | **骨髓初始 B 细胞** `bone marrow naive B cell`（`UBERON:0002371_CL:0000788`） | 采用论文指定的 `UBERON:0000178_CL:0000788`，并在元数据中同时记录教程 ID |
| +58 kb 增强子区间 | **490-bp 核心增强子** `chr2:60495049-60495539`（Bauer et al. 与 Javierre et al. 交集，Fig. 4k 绘图窗口） | **1,288-bp DHS 宽峰** `chr2:60494920-60496208` | Fig. 4k 对照与核心聚焦使用论文 `490-bp` 核心区；同时保留 `1201-bp` 评分窗口与教程 DHS 边界作多尺度标注 |

---

## Decision 1：变异编码与参考基因组正负链规范

- **选择了什么**：使用 GRCh38 / hg38 参考基因组正链（`+` strand）1-based 坐标构建 `genome.Variant(chromosome='chr2', position=60495266, reference_bases='G', alternate_bases='A')`，对应 0-based 半开区间 `[60495265, 60495266)`；同时在基因结构注释中按 GENCODE v46 真实负链（`-` strand）标注 `BCL11A` 转录方向（TSS 位于 `chr2:60553654`，转录由右向左 `5' -> 3'`）。
- **为什么选择**：基因组 VCF 与 AlphaGenome API 均以参考基因组正链定义 `REF`/`ALT`，与受调控靶基因位于正链还是负链无关。`chr2:60495266` 在 GRCh38 正链上的碱基为 `G`，处于正链红系转录因子 GATA1 结合基序 `GATAA`（`chr2:60495266-60495270`）的第 1 位；若误按 `BCL11A` 负链互补为 `C>T`，将导致输入序列与参考基因组不匹配并产生错误推理。
- **依据来源**：AlphaGenome 2 论文 p. 10 Fig. 4j–k 图注与 p. 48 Methods（明确写出 `chr2:60495266:G>A`）；`alphagenome.data.genome.Variant` 文档（明确 `position` 为 1-based，而 `genome.Interval` 为 0-based 半开区间）。
- **替代方案**：按基因转录链（负链）表示碱基（`C>T`）或混用 0-based 坐标传入 `Variant(position=60495265)`——两者均违反 API 规范，会在参考碱基校验阶段失败。
- **如何独立验证**：代码在发起任何模型推理前，通过 UCSC REST API 拉取 `hg38` `chr2:60495262-60495285` 真实序列（`TGTGATAAAAGCAACTGTTAGCT`），断言 `position_1based=60495266`（0-based `60495265`）处碱基严格等于 `G`，且 `[60495265, 60495270)` 为 `GATAA`；同时由 GENCODE v46 核实 `BCL11A` MANE Select 转录本 `ENST00000642384.2` 位于负链，TSS `60553654 - 60495266 = +58,388 bp`（即经典的 `+58 kb` 内含子增强子）。见 `src/bcl11a_repro/validation.py` 与单元测试 `tests/test_bcl11a_repro.py`。

---

## Decision 2：模型版本与推理检查点（`ModelVersion.ALL_FOLDS`）

- **选择了什么**：显式通过 `alphagenome.models.v2.dna_client.create(..., model_version=dna_client.ModelVersion.ALL_FOLDS)` 调用 AlphaGenome 2 全折蒸馏模型（服务端 gRPC 包路径为 `google.gdm.gdmscience.alphagenome.v2main`）。
- **为什么选择**：AlphaGenome 2 论文指出下游变异效应预测、单细胞轨迹展示与案例分析默认使用由多折教师模型蒸馏得到的单学生模型（all-folds distilled student），以兼顾全基因组泛化能力与跨基因座一致性。显式指定 `v2.dna_client` 与 `ModelVersion.ALL_FOLDS` 可杜绝误调 AlphaGenome 1 或单折交叉验证权重。
- **依据来源**：论文 p. 19 *Model architecture and training*、p. 47 *Locus analysis at BCL11A*，以及 `alphagenome==2.0.0` 中 `v2/dna_client.py` 的默认参数定义。
- **替代方案**：使用交叉验证单折模型（`FOLD_0`–`FOLD_3`）。单折模型适合评估未见序列区间的泛化基准（out-of-fold benchmark），但并非论文 Figure 4 案例展示所使用的主模型。
- **如何独立验证**：检查生成的 `data/processed/bcl11a_ag2/run_metadata.json`，确认 `alphagenome_package_version == "2.0.0"`、`api_client_class == "alphagenome.models.v2.dna_client.DnaClient"`、`grpc_service_package == "google.gdm.gdmscience.alphagenome.v2main"` 且 `model_version == "ALL_FOLDS"`。

---

## Decision 3：主要研究细胞与对照细胞条件的精确 Identifier

- **选择了什么**：
  - **Condition A（主要条件）**：`UBERON:0002371_CL:0000765`（`bone marrow erythroblast`，`biosample_type = primary_cell`）；
  - **Condition B（对照条件）**：`UBERON:0000178_CL:0000788`（`blood naive B cell`，`biosample_type = primary_cell`）。
- **为什么选择**：`BCL11A` 在红系前体细胞与 B 淋巴细胞中均高表达且对 B 细胞发育必需，但其 `+58 kb` 增强子具有严格的红系特异性（正是 CRISPR 疗法 `exa-cel` 选择编辑该增强子而非敲除 `BCL11A` 编码区以避免免疫缺陷的原因）。官方教程使用的是骨髓初始 B 细胞（`UBERON:0002371_CL:0000788`），但原论文 Figure 4k 与 p. 47–48 Methods 明确使用的是外周血初始 B 细胞（`UBERON:0000178_CL:0000788`）。按任务要求，存在差异时优先遵循原论文设定。
- **依据来源**：论文 p. 10 Fig. 4k、pp. 47–48 Methods，以及通过 `dna_client.supported_conditions(Organism.HOMO_SAPIENS)` 实时检索得到的 2,799 个人类条件元数据表。
- **替代方案**：使用官方教程中的 `UBERON:0002371_CL:0000788`（`bone marrow naive B cell`）。两者同属初始 B 细胞谱系，但组织来源本体（UBERON）不同；我们未作名称猜测，而是通过官方 `supported_conditions()` 精确核验并将两者的元数据一并留档。
- **如何独立验证**：在 `data/raw/bcl11a_ag2/supported_conditions_selected.csv` 与 `run_metadata.json` 中核对 `condition`、`biosample_name` 与 `sc_atac == "observed"` 字段。

---

## Decision 4：预测模态（`SC_ATAC`）与可比主轨道（`observed` 10x nuclear）选择

- **选择了什么**：调用 `requested_outputs=[dna_client.OutputType.SC_ATAC]` 获取单细胞染色质开放性预测；在返回的每个细胞条件 5 条单细胞核小体/转座酶轨道中，以**唯一具有真实实验观测支撑（`imputation_mode == "observed"`）的 `scATAC-seq 10X-cellranger-nuclear` 轨道**作为主展示与跨细胞可比轨道，同时完整保留其余 4 条 `imputed_supported` 轨道（`EasySci-ATAC`、`SHARE-seq`、`scATAC-seq 10X-bowtie2`、`sciATAC-seq`）用于稳健性检验。
- **为什么选择**：
  1. **为何选 `SC_ATAC` 而非 bulk `ATAC` 或 `RNA_SEQ`**：原发性骨髓成红细胞（erythroblast）深嵌于异质性骨髓组织中，传统 bulk 测序难以在不引入体外培养伪影的前提下分离纯净原代亚群；AlphaGenome 2 的核心突破之一正是引入单细胞多组学图谱（`scATAC-seq` 与 `scRNA-seq`），而论文 Figure 4j–k 直接评估的物理量正是 `scATAC` 局部染色质开放性（enhancer accessibility），这是 GATA1 结合丧失最直接的顺式（*cis*）表观遗传读出。
  2. **为何聚焦 `observed` 10x 轨道**：API 查询证实，在 `UBERON:0002371_CL:0000765` 与 `UBERON:0000178_CL:0000788` 两个条件下，均仅有 `scATAC-seq 10X-cellranger-nuclear` 为 `observed`（真实训练数据存在），其余 4 条均为跨平台填补轨道（`imputed_supported`），无 `imputed_exploratory` 轨道。在两种细胞间锁定完全相同的测序技术栈（10x Genomics + Cell Ranger + nuclear suspension）与 `observed` 状态，消除了实验平台归一化尺度差异带来的混淆。
- **依据来源**：论文 p. 10 Fig. 4j–k、p. 21–24 单细胞轨道训练与 `imputation_mode` 定义、p. 48 Methods，以及 API 返回的 `TrackData.metadata`（共 10 行轨道元数据）。
- **替代方案**：对每个条件的 5 条轨道取算术平均，或同时请求 `SC_GENE_EXPRESSION` / `RNA_SEQ` 等全部模态。对不同实验协议（如 `EasySci` 与 `10X`）直接求平均会混入跨协议尺度差异；而不必要地请求全部模态既消耗配额，也偏离 Phase 1 核心问题。我们在 Figure C panel d 中将 4 条 `imputed_supported` 轨道作为独立散点叠加在 `observed` 柱状图上，兼顾了主指标严谨性与跨协议一致性展示。
- **如何独立验证**：查看 `data/processed/bcl11a_ag2/variant_effect_scores_comparison.csv` 中全部 10 条轨道的 `Assay title`、`data_pipeline`、`imputation_mode` 与评分值（成红细胞 5 条轨道 $\Delta_{\text{ATAC}}$ 介于 `-2.0383` 至 `-2.2192`，B 细胞介于 `-0.0024` 至 `-0.0120`，结论在所有轨道上完全一致）。

---

## Decision 5：1-Mb 模型输入窗口与多尺度增强子区间选择

- **选择了什么**：
  1. **模型推理输入窗口（`1 Mb = 1,048,576 bp`）**：通过 `variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)` 构建以变异为中心的 0-based 半开区间 `chr2:59970978-61019554`；
  2. **局部轨迹展示窗口**：同时展示 `1201-bp` 官方评分窗口（`chr2:60494665-60495866`）与论文 Figure 4k 使用的 **490-bp `+58 kb` 核心增强子窗口**（`chr2:60495049-60495539`），并在基因座总览图（Figure A）中额外标出教程使用的 `1,288-bp` DHS 窗口（`chr2:60494920-60496208`）以及相邻 `+62 kb`（`60490400-60491800`）和 `+55 kb`（`60497800-60499200`）增强子。
- **为什么选择**：AlphaGenome 2 的标准全分辨率上下文长度为 $2^{20} = 1,048,576\text{ bp}$（`SEQUENCE_LENGTH_1MB`），足以完整覆盖 `BCL11A` 全基因体（`chr2:60450519-60554467`）及两侧超远端调控环境；而在结果可视化上，既展示完整的 `1201-bp` 评分积分域（防止局部截断造成视觉夸大），又按论文 p. 47 明确给出的 `490-bp` 核心增强子坐标复现 Figure 4k。
- **依据来源**：论文 p. 19（1-Mb 模型输入）、p. 47 *Locus analysis at BCL11A*（明确给出 `+58 kb` 490-bp 区间 `chr2:60495049-60495539` 为 Bauer et al. 与 Javierre et al. 交集），以及 `dna_client.SEQUENCE_LENGTH_1MB` 常量定义。
- **替代方案**：使用较短输入窗口（如 `16 kb` / `100 kb`）进行推理，或仅绘制教程中的 `1,288-bp` DHS 窗口。较短输入窗口会损失远端染色质三维上下文并偏离论文默认推理设置。
- **如何独立验证**：在 `run_metadata.json` 的 `windows_0based_half_open` 中核查 `61019554 - 59970978 == 1048576`，并验证变异 0-based 坐标 `60495265` 位于所有上述增强子与评分窗口内部。

---

## Decision 6：1201-bp `CenterMaskScorer` 定量评分与独立数组交叉验证（含 1-bp 居中细节）

- **选择了什么**：使用官方推荐的 `variant_scorers.RECOMMENDED_VARIANT_SCORERS['SC_ATAC']`，其实例为 `CenterMaskScorer(requested_output=SC_ATAC, width=1201, aggregation_type=DIFF_LOG2_SUM)`；同时从 `predict_variant()` 返回的 1-bp 分辨率 `REF`/`ALT` 原始预测数组中，严格截取以变异 0-based 坐标 `60495265` 为中心的对称窗口 `[60494665, 60495866)`（`pos_0 - 600` 至 `pos_0 + 601`，共 `1201 bp`），按以下公式独立复算：
  $$\Delta_{\text{ATAC}} = \log_2\!\left(\sum_{i \in W_{1201}} \hat{y}_{\text{ALT},i} + 1\right) - \log_2\!\left(\sum_{i \in W_{1201}} \hat{y}_{\text{REF},i} + 1\right)$$
- **为什么选择**：
  1. **为何采用 `1201-bp` 窗口**：单个转录因子基序（如 5-bp `GATAA`）的破坏不仅改变突变位点单个碱基的 Tn5 切割概率，而是导致局部约 `0.5–1.0 kb` 核小体排斥区（nucleosome-depleted region, NDR）整体闭合。论文 Table 5（p. 35）将 `SC_ATAC` 推荐评分器统一设定为变异两侧各 `600 bp`（总宽 `1201 bp`）的局部积分对数比值（`DIFF_LOG2_SUM`，伪计数 `+1`），既捕获完整的开放染色质峰塌缩，又抑制单碱基切割偏好的高频噪声。
  2. **1-bp 居中边界细节**：由于 `1,048,576` 为偶数，`variant.reference_interval.resize(2**20)` 得到的窗口 `[59970978, 61019554)` 几何中点在 `60495266`；若直接对 1-Mb `TrackData` 调用 `.resize(1201)`，会得到 `[60494666, 60495867)`（相对变异位点右偏 `1 bp`，与服务端评分存在约 $6.85 \times 10^{-5}$ 的微小偏差）。而服务端 `CenterMaskScorer` 实际以变异本身 `variant.reference_interval.resize(1201)` 即 `[60494665, 60495866)` 为掩码中心。采用精确对称切片后，全部 10 条轨道的手工复算值与服务端 `score_variant()` 返回值最大绝对误差降至 $5.77 \times 10^{-7}$（完全处于 IEEE 754 `float32` 舍入精度内）。
- **依据来源**：论文 p. 35 Table 5（Variant scoring）与 p. 48 Methods（明确写出 *"a 1,201 bp window centered at the variant using the `CenterMaskScorer` (`DIFF_LOG2_SUM`)"*），以及 `alphagenome.models.variant_scorers` 源码。
- **替代方案**：使用 `501-bp` 窗口、`DIFF_MEAN` 线性差值，或直接调用 `track_data.resize(1201)`。线性差值对实验测序深度敏感，而 `track_data.resize(1201)` 会引入上述 1-bp 偶数区间舍入偏移。
- **如何独立验证**：运行单元测试 `TestCenterMaskManualCrossValidation`（包含偶数长度区间 1-bp 偏移的专门验证）并检查 `data/processed/bcl11a_ag2/variant_effect_scores_comparison.csv` 中的 `api_raw_score`、`manual_log2_score` 与 `abs_diff` 列。

---

## Decision 7：Figure 4j 饱和突变（ISM）窗口与 SDK `v2.0.0` 元数据兼容处理

- **选择了什么**：对论文 p. 48 指定的 23-bp 区间 `chr2:60495262-60495285`（0-based 半开，对应 1-based `60,495,263–60,495,285`，完整覆盖治疗性 `sg1617` gRNA 靶点及 `GATAA` 基序）分别在 `REF (G)` 与 `ALT (A)` 序列背景下执行饱和单碱基突变评分（共 $2 \times 23 \times 3 = 138$ 个单碱基替换变异），并使用 `alphagenome.interpretation.ism.ism_matrix(..., multiply_by_sequence=True)` 计算输入序列归因矩阵。针对 `alphagenome==2.0.0` 中 `dna_client.score_ism_variants()` 返回的 `AnnData` 缺少 `adata.uns['variant_scorer']` 键导致 `variant_scorers.tidy_scores()` 抛出 `ValueError` 的已知问题，编写了显式解析函数 `ism_anndata_list_to_tidy_df()` 直接提取 `adata.uns['variant']`、`adata.var` 与 `adata.X`。
- **为什么选择**：Figure 4j 是论文解释“为何 `chr2:60495266:G>A` 会导致红系染色质关闭”的核心机制证据（展示 `REF` 上完整的 `GATAA` 正向归因在 `ALT` 背景下整体消融）。在 ALT 背景计算 ISM 时，我们将 1-Mb 参考序列在 `60495266` 处替换为 `A`，通过 `score_ism_variants` 的自定义 `sequence` 参数直接由真实 API 推理，未作任何近似。
- **依据来源**：论文 p. 10 Fig. 4j 与 p. 48 *In silico mutagenesis*（明确给出 `chr2:60495262-60495285` 区间与 `CenterMaskScorer` 设置），以及本地 `alphagenome/interpretation/ism.py` 与 `alphagenome/models/v2/dna_client.py` 源码核查。
- **替代方案**：跳过 Figure 4j 仅复现 Figure 4k，或手工修补第三方包源码。修改 `.venv` 第三方包会破坏可复现性；在项目自有模块 `src/bcl11a_repro/scoring.py` 中安全解析并辅以离线回归测试更为稳健。
- **如何独立验证**：运行回归测试 `test_ism_anndata_list_to_tidy_df_handles_missing_uns_scorer`，并核查 `data/processed/bcl11a_ag2/ism_attribution_chr2_60495262_60495285.csv`（46 行逐位归因数值）。

---

## Decision 8：原论文精确数值不可得时的比较边界

- **选择了什么**：对原论文 Figure 4j–k 进行严格的**定性一致性与量级核对**（坐标区间、模型版本、细胞条件、`GATAA` 归因正负号与约 `+2` 量级、成红细胞 `>20` 峰值塌缩、B 细胞近零平坦基线），明确将“与论文精确浮点数值逐位比对”标记为**因论文未公开原始数值表而无法执行**，绝不通过目测或像素数字化量取论文 PDF 曲线来声称“数值精确复现”。
- **为什么选择**：核查 `alphagenome-2.pdf` 全文及公开数据下载说明（pp. 52–53）确认，论文仅公开了 GWAS 表达评分汇总表，未附带 Figure 4j–k 单基因座轨迹或单变异分数的精确浮点表格。此外，论文成稿时的内部推理检查点与当前公开 gRPC 服务 `v2.0.0` 在浮点实现细节上可能存在细微版本演进，将像素估读值当作金标准有违科学诚信。
- **依据来源**：论文 pp. 10, 47–48, 52–53，以及任务说明第五、七、九节的严格科研规范。
- **替代方案**：从 PDF 矢量图提取像素坐标估算数值误差——已被任务规范明确禁止且不具科学严谨性。
- **如何独立验证**：直接对比论文 p. 10 Figure 4j–k 原图与本项目生成的 `docs/figures/bcl11a_ag2/fig_d_paper_comparison.png` 及 `data/processed/bcl11a_ag2/` 下的真实数据文件，所有我方数值均可由缓存或实时 API 100% 确定性重算。
