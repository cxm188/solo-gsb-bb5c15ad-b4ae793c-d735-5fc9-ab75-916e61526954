# MSPC — 质谱同位素分布预计算 API

基于 **Python + FastAPI + SQLite** 的质谱预计算服务。客户端提交含元素计数的分子式（可带同位素标记），服务按各元素天然同位素丰度将分布表示为"质量-强度多项式"，逐元素卷积得到完整同位素包络，合并质量差小于 0.0001 Da 的峰，按基峰归一化后返回前 N 个峰及分子单同位素质量。

## 快速启动（单条命令）

```bash
docker compose up --build
```

启动过程会自动完成依赖安装（镜像构建时）与 SQLite 数据库初始化（服务启动时建表），无需任何手工迁移。

## 访问地址

| 地址 | 说明 |
|---|---|
| http://localhost:8000 | 服务根地址（端口 8000） |
| http://localhost:8000/docs | Swagger 交互式文档 |
| http://localhost:8000/health | 健康检查 |

## API 用法

### `POST /api/v1/envelope`

请求体：

```json
{
  "formula": "C6H12O6",
  "top_n": 5
}
```

- `formula`（必填）：分子式，如 `C6H12O6`；同位素标记用方括号表示，如 `C6H6[13C]2`（2 个纯 ¹³C 原子）、`C2H5[15N]O[18O]`。
- `top_n`（可选，默认 20，最大 200）：返回强度最高的前 N 个峰。

响应：

```json
{
  "formula": "C6H12O6",
  "composition": {"C": 6, "H": 12, "O": 6},
  "labels": {},
  "monoisotopic_mass": 180.063388,
  "base_peak_intensity": 100.0,
  "peak_count": 5,
  "peaks": [
    {"mass": 180.063388, "intensity": 100.0},
    {"mass": 181.066743, "intensity": 6.489437}
  ],
  "cached": false
}
```

- `monoisotopic_mass`：单同位素质量（各元素取天然丰度最高的同位素，标记原子取标记同位素质量）。
- `peaks`：基峰强度归一化为 100，按质量升序返回前 N 个峰。
- `cached`：结果是否来自 SQLite 缓存（相同 `formula` + `top_n` 的重复请求直接命中缓存）。

### `POST /api/v1/enrichment-fit`

估计分子中单个同位素标记的丰度（富集度）：在 0%～100% 范围内按给定步长网格搜索候选丰度，对每个候选值复用卷积生成归一化包络，与观测峰配对比较，返回均方误差最小的丰度。

请求体：

```json
{
  "formula": "C6",
  "label": "13C",
  "peaks": [
    {"mass": 72.0, "intensity": 5.0},
    {"mass": 73.003355, "intensity": 30.0},
    {"mass": 74.00671, "intensity": 75.0},
    {"mass": 75.010065, "intensity": 100.0},
    {"mass": 76.013419, "intensity": 75.0},
    {"mass": 77.016774, "intensity": 30.0},
    {"mass": 78.020129, "intensity": 5.0}
  ],
  "step_percent": 0.5,
  "mass_tolerance": 0.05
}
```

- `formula`（必填）：分子式，规则同 `/api/v1/envelope`（可含 `[13C]` 等纯同位素标记，标记贡献固定质量位移）。
- `label`（必填）：待估计的单个同位素标签，格式为"质量数+元素符号"，如 `13C`、`15N`、`2H`。
- `peaks`（必填，1～200 个）：观测峰，质量必须**严格递增**，强度非负且至少一个为正。
- `step_percent`（可选，默认 1.0，范围 0.5～10）：丰度搜索步长（百分点）。候选点为 0、step、2·step、…、100（100% 端点始终包含）。
- `mass_tolerance`（必填，>0）：观测峰与预测峰配对的质量容差（Da）。

拟合模型：候选丰度为 p 时，被标记元素的同位素分布改为"标记同位素占 p，其余天然同位素按比例分享 1−p"（p 等于天然丰度时即还原为天然分布），再与其他元素的天然分布卷积得到预测包络；预测包络与观测峰均按各自基峰归一化为 100 后比较。

响应：

```json
{
  "formula": "C6",
  "label": "13C",
  "best_abundance_percent": 50.0,
  "mse": 0.0,
  "search_points": 201,
  "pair_count": 7,
  "pairs": [
    {"observed_mass": 72.0, "observed_intensity": 5.0, "predicted_mass": 72.0, "predicted_intensity": 5.0},
    {"observed_mass": 73.003355, "observed_intensity": 30.0, "predicted_mass": 73.003355, "predicted_intensity": 30.0},
    {"observed_mass": 74.00671, "observed_intensity": 75.0, "predicted_mass": 74.00671, "predicted_intensity": 75.0},
    {"observed_mass": 75.010065, "observed_intensity": 100.0, "predicted_mass": 75.010065, "predicted_intensity": 100.0}
  ]
}
```

- `best_abundance_percent`：最佳丰度（百分比）。MSE 相同时取较低丰度。
- `mse`：最佳丰度下全部配对的强度均方误差。
- `search_points`：搜索点数（0%～100% 按步长划分的候选数）。
- `pairs`：最佳丰度下的逐峰配对，按质量升序。每个观测峰与容差内最近的未占用预测峰配对；未匹配的观测峰 `predicted_*` 为 `null`（预测强度按 0 计），未匹配的预测峰 `observed_*` 为 `null`（观测强度按 0 计）。两侧强度均为基峰=100 的归一化值。

配对与误差细节：所有容差内的（观测, 预测）候选对按质量差从小到大贪心占用，每个峰至多参与一对；MSE 在全部配对（含未匹配峰）上计算。

### 错误定位

非法分子式或标记返回 HTTP 422，并指出出错字符的 **0 基位置**：

```json
{
  "detail": {
    "error": "invalid_formula",
    "message": "Element 'C' has no known isotope with mass number 20; check the label mass",
    "position": 2
  }
}
```

可检测的错误包括：未知元素符号、元素计数为 0、未闭合的 `[`、标记格式错误、未知标记质量数、孤立数字、非法字符等。

`/api/v1/enrichment-fit` 的 422 错误同样定位明确：

```json
{
  "detail": {
    "error": "invalid_label",
    "message": "Element 'C' has no known isotope with mass number 20; check the label mass",
    "position": 0
  }
}
```

- `invalid_formula`：分子式非法（`position` 为分子式中的 0 基字符位置）。
- `invalid_label`：标签非法——格式错误（应为 `13C` 形式）、未知元素符号、未知质量数（`position` 为标签字符串中的 0 基字符位置）。
- `invalid_peaks`：观测峰列表为空、质量非递增或重复（`index` 为出错峰的 0 基索引，空列表时为 `null`）。
- 参数越界（`step_percent` 超出 0.5～10、`mass_tolerance` 非正、峰数超过 200、强度为负等）由请求校验返回标准 422 并指明字段位置。

### 其他端点

- `GET /api/v1/elements`：列出全部支持的元素及其天然同位素质量/丰度。
- `GET /health`：健康检查。

## 配置（环境变量）

均可在 `docker-compose.yml` 的 `environment` 段覆盖：

| 变量 | 默认值 | 说明 |
|---|---|---|
| `MSPC_DB_PATH` | 容器内 `/data/mspc.db`，本地 `./mspc.db` | SQLite 数据库文件路径 |
| `MSPC_MERGE_TOLERANCE` | `0.0001` | 峰合并质量容差（Da） |
| `MSPC_PRUNE_THRESHOLD` | `1e-12` | 卷积过程中按相对强度剪枝的阈值 |
| `MSPC_TOP_N_MAX` | `200` | `top_n` 上限 |
| `MSPC_TOP_N_DEFAULT` | `20` | `top_n` 默认值 |

## 本地开发（不用 Docker）

```bash
pip install -r requirements-dev.txt
uvicorn app.main:app --reload          # 服务: http://127.0.0.1:8000
pytest                                 # 运行测试
```

## 实现说明

- **多项式卷积**：每种元素的同位素模式表示为多项式 `Σ aᵢ·x^mᵢ`（丰度 × 质量的幂），分子包络为各元素多项式（按原子数取幂，使用快速幂）之积；多项式乘法即 (质量, 强度) 峰列表的卷积。纯同位素标记贡献固定质量位移。
- **峰合并**：排序后相邻质量差 < 0.0001 Da 的峰融合为一个，质量取强度加权平均，强度取和。
- **单同位素质量**：按 IUPAC 惯例取各元素天然丰度最高的同位素质量之和。
- **富集拟合**：对每个候选丰度 p（0%～100% 按步长网格搜索），将被标记元素的同位素分布替换为"标记同位素占 p、其余天然同位素按比例分享 1−p"，复用同一套卷积生成归一化包络（基峰=100，观测峰同样归一化）；观测峰与容差内最近的预测峰贪心配对，未匹配峰强度按 0 计，以全部配对的强度均方误差最小者为最佳丰度，同分取较低丰度。除被拟合元素外其余元素的包络与候选丰度无关，只计算一次。
- **同位素数据**：`app/isotopes.py` 内置约 50 种常见元素的 AME/CIAAW 同位素质量与天然丰度，可按需扩充。
- **持久化**：SQLite 存计算结果缓存（`envelope_cache`）与请求日志（`request_log`），schema 随服务启动自动创建。
