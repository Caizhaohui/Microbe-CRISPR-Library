# Cas9_knockout_designer_v9 使用说明

本说明文档对应脚本：

- `Cas9_knockout_designer_v9.py`

该版本在 `v8` 基础上增加了前置策略路由参数 `--dele_model`，用于在同一套命令行入口下切换两类删除策略：

- `Mt`：针对 *M. thermophila* 的定向删除策略（按 PAM 方向定义 cut window）
- `normal`：通用删除策略（V7 风格，走 `legacy_length`）

---

## 1. 运行环境

建议环境：

- Python 3.8 及以上
- 依赖库：
  - `pandas`
  - `biopython`
  - `gffutils`

安装示例：

```bash
pip install pandas biopython gffutils
```

---

## 2. 核心新增参数：`--dele_model`

### 2.1 参数定义

```text
--dele_model {normal, Mt}
```

- 默认值：`normal`
- 含义：
  - `Mt`：强制使用 Mt 删除策略（`cut_window`）
  - `normal`：强制使用通用删除策略（`legacy_length`）

### 2.2 启动审计日志

脚本启动后会明确打印当前策略，便于批量任务审计：

- `[V9审计] 当前采用 Mt 删除策略`
- `[V9审计] 当前采用 normal 删除策略`

---

## 3. 两种删除策略的差异

## 3.1 `--dele_model Mt`

- 内部自动注入：
  - `--deletion_mode cut_window`
- 删除窗口按 PAM 方向定义：
  - 上游 `cut_window_upstream`（默认 20 bp）
  - 下游 `cut_window_downstream`（默认 100 bp）
- 删除长度优先选择窗口内最大可用且非 3 倍数（优先 frameshift）

适用场景：

- 明确采用 Mt 经验窗口（切点上游 20 / 下游 100）
- 需要与 V9 Mt 优化实验策略对齐

## 3.2 `--dele_model normal`

- 内部自动注入：
  - `--deletion_mode legacy_length`
- 删除长度遵循 V7/V8 legacy 逻辑：
  - 通过 `--del_length_per` / `--del_length_bp` 控制长度约束

适用场景：

- 非 Mt 物种或希望保持跨物种通用参数化长度约束

---

## 4. 输入模式

脚本继承 v8 的输入方式，支持：

- FASTA + GFF3：
  - `--input_fna`
  - `--input_gff`
- GBFF：
  - `--input_gbff`

二选一模式，不可混用。

---

## 5. 常用参数速查

- 基础参数：
  - `--output`：输出 CSV 路径
  - `--species`：物种名称（默认 `M_thermophila`）
  - `--synthesis_template`：合成模板文件（必填）
  - `--sgRNA_num`：每基因目标设计数（默认 2）
  - `--barcode_len`：条形码长度
  - `--restriction_site`：需规避的酶切位点列表
  - `--max_oligo_length`：目标 oligo 长度（固定长度设计）
  - `--num_workers`：并行线程数

- 与删除策略相关：
  - `--dele_model {normal, Mt}`
  - `--cut_window`：仅 Mt/cut_window 模式使用，例如 `20:100`
  - `--del_length_per`：仅 normal/legacy 模式常用，例如 `10%:80%`
  - `--del_length_bp`：仅 normal/legacy 模式常用，例如 `300:1000`

---

## 6. 命令示例

## 6.1 Mt 策略（推荐用于 Mt 数据）

```bash
python Cas9_knockout_designer_v9.py \
  --dele_model Mt \
  --input_gbff Mt_genomic.gbff \
  --output Mt_V9_Mt_KO.csv \
  --synthesis_template Mt_knockout_library_oligo_template.txt \
  --species M_thermophila \
  --barcode_len 11 \
  --max_oligo_length 300 \
  --restriction_site GGTCTC GAAGAC
```

## 6.2 normal 策略（通用 V7 风格）

```bash
python Cas9_knockout_designer_v9.py \
  --dele_model normal \
  --input_gbff Mt_genomic.gbff \
  --output Mt_V9_normal_KO.csv \
  --synthesis_template Mt_knockout_library_oligo_template.txt \
  --species M_thermophila \
  --barcode_len 11 \
  --max_oligo_length 300 \
  --restriction_site GGTCTC GAAGAC \
  --del_length_per 10%:80% \
  --del_length_bp 300:1000
```

---

## 7. 输出文件说明

运行成功后会输出两类文件：

- 主结果文件：`<output>.csv`
  - 成功设计条目
- 失败汇总文件：`<output>_failed.csv`
  - `Failed`：0 条设计
  - `Partial`：设计条数小于 `--sgRNA_num`

主结果常见字段：

- `Gene Id`
- `Status`
- `Sgrna Seq`
- `Sgrna Pam`
- `Sgrna Cut Site`
- `Deletion Start`
- `Deletion End`
- `Deletion Length`
- `Final Oligo For Synthesis`

---

## 8. 参数组合建议

- 若选择 `--dele_model Mt`：
  - 优先与 Mt 实验经验窗口一起使用
  - 通常不再需要 `--del_length_per/--del_length_bp`

- 若选择 `--dele_model normal`：
  - 建议总是显式提供 `--del_length_per` 或 `--del_length_bp`（或同时提供交集约束）
  - 对短基因场景可适当放宽 `--del_length_bp` 下限

---

## 9. 常见问题

### Q1：为什么 `normal` 模式报缺少删除长度参数？

因为 `normal` 会强制进入 `legacy_length` 逻辑，该模式需要长度约束参数驱动，请补充：

- `--del_length_per`
- 或 `--del_length_bp`

### Q2：如何确认实际运行的是哪种策略？

看启动日志中的审计行：

- `当前采用 Mt 删除策略`
- `当前采用 normal 删除策略`

---

## 10. 版本说明

`v9` 相对 `v8` 的核心增强：

- 新增前置参数路由 `--dele_model`
- `Mt` 模式下使用 PAM 方向感知的删除窗口
- 增加策略审计日志，便于批量任务追踪
