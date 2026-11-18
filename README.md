# AlphaSickleGenome

基于 Google DeepMind AlphaGenome 2 的镰刀型细胞贫血（SCD）与胎儿血红蛋白（HbF）调控机制研究项目。

仓库：https://github.com/KenXiao1/AlphaSickleGenome

## 1. 环境配置

- 使用 `uv` 管理 Python 3.12 与依赖（由 `uv.lock` 严格锁定，`alphagenome==2.0.0`）。
- 同步虚拟环境：
  ```powershell
  uv sync --locked
  ```

## 2. API Key 配置

首次克隆后，将 `.env.example` 复制为 `.env`，从 https://alphagenome.google/ 获取 API Key 并填入：

```dotenv
ALPHAGENOME_API_KEY=your_api_key_here
```

- `.env` 已被 `.gitignore` 忽略，**切勿将真实 API Key 提交至 Git 或打印到日志**。
- 运行分析脚本时，可通过 `uv run --env-file .env ...` 显式注入环境变量，也可直接调用脚本（入口程序会通过 `python-dotenv` 安全加载本地 `.env`，但从不输出密钥内容）。

## 3. Phase 1：AlphaGenome 2 `BCL11A` +58 kb 增强子变异效应复现

复现 AlphaGenome 2 技术报告 Figure 4j–k 中 `GRCh38 chr2:60495266:G>A` 在原代骨髓成红细胞（`UBERON:0002371_CL:0000765`）与外周血初始 B 细胞（`UBERON:0000178_CL:0000788`）中的单细胞染色质开放性（`scATAC-seq`）预测、`1201-bp` `CenterMaskScorer` 评分交叉验证及局部饱和突变（ISM）归因。

### 运行命令

1. **运行完整分析流水线（支持本地缓存，命中缓存时不重复消耗 API 配额）**：
   ```powershell
   uv run --env-file .env python src/reproduce_bcl11a.py
   ```
2. **纯离线校验与重生成（不读取 API Key、不发起任何网络请求，从本地缓存复算所有表格与图表）**：
   ```powershell
   uv run python src/reproduce_bcl11a.py --offline
   ```
3. **仅从已保存的处理后数据重新渲染全部图表（Figures A–D，PNG + PDF）**：
   ```powershell
   uv run python src/reproduce_bcl11a.py --plots-only
   ```
4. **强制忽略本地缓存重新请求 AlphaGenome 2 API（会消耗配额，仅在需要刷新服务端输出时使用）**：
   ```powershell
   uv run --env-file .env python src/reproduce_bcl11a.py --refresh-cache
   ```

### 离线单元测试（不依赖 API Key）

单元测试使用显式标记为 `SYNTHETIC_*` 的模拟数据验证基因组坐标校验、链方向检查、`1201-bp` 评分复算（含偶数长度窗口 `1-bp` 边界校验）、SDK `AnnData` 解析及图表渲染：

```powershell
uv run python -m unittest discover -s tests -v
```

## 4. 项目目录结构

```text
AlphaSickleGenome/
├── src/
│   ├── reproduce_bcl11a.py          # Phase 1 CLI 主入口
│   ├── bcl11a_repro/                # Phase 1 模块化实现（常量、校验、评分、绘图、缓存流水线）
│   └── data_recon.py                # 既有 Cheng et al. (2021) 数据侦察脚本
├── tests/
│   └── test_bcl11a_repro.py         # 不依赖 API Key 的离线单元测试（模拟数据与真实数据严格分离）
├── data/
│   ├── cache/bcl11a_ag2/            # 完整 1-Mb 原始二进制缓存与 GENCODE feather（已被 Git 忽略）
│   ├── raw/
│   │   ├── bcl11a_ag2/              # 真实 API 紧凑原始导出（10-kb 原始 npz、条件表、官方评分、ISM 矩阵）
│   │   └── 41588_2021_861_MOESM3_ESM.xlsx
│   └── processed/bcl11a_ag2/        # 处理后的轨迹表、交叉验证评分表、ISM 逐位表与运行元数据（含 SHA-256）
└── docs/
    ├── bcl11a_reproduction_report.md # Phase 1 科学复现报告
    ├── research_log.md               # Phase 1 关键科研决策日志
    ├── figures/bcl11a_ag2/           # Figure A–D（300 DPI PNG 与矢量 PDF）
    └── data_recon.md                 # 既有 Cheng et al. (2021) 数据侦察文档
```
