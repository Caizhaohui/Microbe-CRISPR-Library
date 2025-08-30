# Bact-CRISPR-Library: 细菌CRISPR文库自动化设计工具

`Bact-CRISPR-Library` 是一个功能强大的、基于命令行的Python工具，专门用于在细菌基因组中进行大规模寡核苷酸（Oligo）文库的自动化设计。它支持多种CRISPR基因编辑应用，并内置了多种智能优化策略，以确保设计出的文库具有高质量和高成功率。

## ✨ 主要功能

- **多种设计模式**: 支持三种核心的基因编辑文库设计：
  - `knockout`: 基因敲除文库
  - `promoter_replace`: 启动子替换文库
  - `C_fusion`: 蛋白C端融合标签文库
- **智能sgRNA筛选**: 集成了基于模型的sgRNA活性评分，并能优先选择最理想的切割位点。
- **灵活的同源臂设计**: 支持自定义同源臂长度范围，并为启动子替换模式提供了独特的安全检查，以保护上游基因。
- **支持环状与线性基因组**: 通过`--genome_type`参数，可正确处理环状（如多数细菌）和线性基因组，避免因基因组边界问题导致的设计失败。
- **先进的酶切位点规避**: 可在设计的可变区（sgRNA, 同源臂, barcode）中规避指定的限制性酶切位点，同时豁免固定的克隆位点。
- **自动化Barcode生成**: 为每个设计生成符合GC含量和序列复杂性要求的独特条形码（Barcode）。
- **模板化Oligo合成**: 最终的Oligo序列是基于用户提供的模板文件生成的，具有极高的灵活性和可定制性。

## ⚙️ 系统要求与安装

**1. Python版本**:
   - Python 3.7 或更高版本

**2. 依赖库**:
   - `pandas`
   - `biopython`
   - `gffutils`

**3. 安装**:
   使用pip安装所有依赖库：
   ```bash
   pip install pandas biopython gffutils
   ```

## 🚀 使用方法

脚本通过命令行运行，基本结构如下：

```bash
python Bact-CRISPR-Library <mode> [common_options] [mode_specific_options]
```

### 通用参数

这些参数适用于所有设计模式：

| 参数 | 描述 | 示例 |
| --- | --- | --- |
| `--input_fna` | **[必需]** 输入的基因组FASTA文件路径。 | `--input_fna genome.fna` |
| `--input_gff` | **[必需]** 输入的GFF3基因注释文件路径。 | `--input_gff annotation.gff` |
| `--output` | **[必需]** 输出结果的CSV文件路径。 | `--output knockout_library.csv` |
| `--synthesis_template` | **[必需]** 定义最终Oligo结构的文本模板文件。| `--synthesis_template oligo_template.txt`|
| `--genome_type` | 指定基因组类型，`circle`为环状，`linear`为线性。| `--genome_type circle` |
| `--pam` | PAM序列，N代表任意碱基。 | `--pam NGG` (默认) |
| `--HR_len` | 同源臂长度。格式为'首选:最小'或单个数字。 | `--HR_len 50:40` (默认) |
| `--sgRNA_num` | 为每个基因生成的设计方案数量。 | `--sgRNA_num 2` |
| `--restriction_site` | 需要在可变区中避免的一个或多个酶切位点。 | `--restriction_site gaagac ggtctc` |
| `--cloning_site` | 一个特殊的酶切位点，通常用于后续克隆，它在最终检查中被豁免。| `--cloning_site AAGGTCTCA` |

### 设计模式详解与示例

#### 1. 基因敲除 (`knockout`)

此模式用于设计在每个基因编码区内部产生小片段删除的文库。

**特有参数**:
- `--del_length`: 删除片段的长度范围，格式为 '最小:最大'。 (默认: `50:100`)
- `--promoter_region`: 定义启动子保护区的大小，避免删除。 (默认: `50:150`)
- `--ko_search_range`: 在CDS中搜索sgRNA的百分比范围。 (默认: `5:80`)

**示例命令**:
```bash
python Bact-CRISPR-Library knockout \
  --input_fna ecoli.fna \
  --input_gff ecoli.gff \
  --output ecoli_knockout_lib.csv \
  --genome_type circle \
  --pam NGG \
  --HR_len 46:40 \
  --del_length 20:100 \
  --sgRNA_num 2 \
  --restriction_site gaagac ggtctc \
  --cloning_site AAGGTGAGACCAAGGTCTCTCGAC \
  --synthesis_template knockout_template.txt
```

#### 2. 启动子替换 (`promoter_replace`)

此模式用于在每个基因的起始密码子上游替换天然启动子。最终合成的Oligo将包含一个克隆位点，用于后续插入具体的启动子序列。

**特有参数**:
- `--promoter_search_size`: 在起始密码子（ATG）上游搜索sgRNA的区域大小(bp)。 (默认: `150`)

**示例命令**:
```bash
python Bact-CRISPR-Library promoter_replace \
  --input_fna ecoli.fna \
  --input_gff ecoli.gff \
  --output ecoli_promoter_lib.csv \
  --genome_type circle \
  --HR_len 50:45 \
  --restriction_site gaagac \
  --cloning_site GGTCTCAAGCTT \
  --synthesis_template promoter_template.txt
```

#### 3. C端融合 (`C_fusion`)

此模式用于在每个蛋白的C末端（终止密码子处）插入一个标签序列（如荧光蛋白、降解决定子等）。

**特有参数**:
- `--insert_sequence`: **[必需]** 要插入的序列（例如，Linker-mCherry），必须包含终止密码子。
- `--insert_name`: **[必需]** 插入序列的名称，将记录在输出文件中。
- `--sgrna_search_range`: 在终止密码子周围搜索sgRNA的范围(bp)。 (默认: `50`)

**示例命令**:
```bash
python Bact-CRISPR-Library C_fusion \
  --input_fna ecoli.fna \
  --input_gff ecoli.gff \
  --output ecoli_Cfusion_lib.csv \
  --genome_type circle \
  --insert_sequence "GGAGGCGGATCGTAGTAGTAA" \
  --insert_name "Linker-3xFLAG" \
  --sgrna_search_range 40 \
  --restriction_site gaagac \
  --cloning_site GGTCTC \
  --synthesis_template C_fusion_template.txt
```

## 📄 输入文件格式

#### 1. 基因组序列 (`--input_fna`)
标准的FASTA格式文件。对于环状基因组，文件内容仍是线性的，程序会通过`--genome_type circle`参数来正确处理。

#### 2. 基因注释 (`--input_gff`)
标准的GFF3格式。文件中必须包含`gene`和`CDS`类型的特征（feature），并且`gene`特征需要有`locus_tag`属性作为基因的唯一标识符。

#### 3. Oligo合成模板 (`--synthesis_template`)
一个纯文本文件，定义了最终输出的Oligo的结构。文件中可以使用以下占位符，程序会自动替换它们：

- `{sgRNA_fwd}`: sgRNA序列 (20nt)
- `{sgRNA_rc}`: sgRNA的反向互补序列
- `{upstream_arm}`: 上游同源臂
- `{downstream_arm}`: 下游同源臂
- `{barcode}`: 随机生成的条形码
- `{exempt_restriction_site}`: 用于后续克隆的酶切位点（由`--cloning_site`参数指定）
- `{insert}`: 要插入的序列 (仅用于`C_fusion`模式)

**模板示例 (`promoter_template.txt`)**:
```
GACGACTTGCTATTTCTAGCTCTAAAAC{sgRNA_fwd}GTTAAACCCTATAGTGAGTCGTATTAC{barcode}AGGTGAGACC{exempt_restriction_site}GAC{upstream_arm}{downstream_arm}ACTCGCACTGCTGGTCACTTGGTCTGGA
```

## 📤 输出文件格式

程序会生成两个CSV文件：

1.  **`<output_name>.csv`**: 包含所有成功设计的结果。关键列包括：
    - `Gene_ID`: 目标基因的locus_tag。
    - `Status`: 设计状态 ("Success")。
    - `sgRNA_Sequence`: 设计的sgRNA序列。
    - `Barcode`: 生成的条形码。
    - `Final_Oligo_for_Synthesis`: 最终用于合成的完整Oligo序列。
    - `Arm_Length`: 同源臂长度。
    - `Upstream_Arm`, `Downstream_Arm`: 同源臂序列。
    - `Design_Scenario`: （仅promoter_replace模式）显示采用的设计策略和参数。
    - ... 以及其他详细信息。

2.  **`<output_name>_failed.csv`**: 包含所有设计失败的基因ID列表。

## 📜 许可

本项目采用 [MIT License](LICENSE)。

## 📧 联系方式

- **作者**: [Zhaohui Cai]
- **邮箱**: [cai_zhaohui@163.com]
