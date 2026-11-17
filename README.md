# AlphaSickleGenome

使用 AlphaGenome 的研究项目。目前仅配置环境，不包含业务代码。

仓库：https://github.com/KenXiao1/AlphaSickleGenome

## 环境

- 使用 uv 管理 Python 3.12 和依赖。
- 运行 `uv sync --locked`，根据 `uv.lock` 安装依赖到 `.venv`。
- 已安装 `alphagenome` 和 `python-dotenv`。

## API Key

首次克隆后，将 `.env.example` 复制为 `.env`；当前工作目录已创建 `.env`。
从 https://alphagenome.google/ 获取 API Key，填入 `.env` 的
`ALPHAGENOME_API_KEY`。

`.env` 已被 Git 忽略。后续运行程序时，可通过 `uv run --env-file .env …`
加载环境变量；仅安装 `python-dotenv` 不会自动加载 `.env`。
