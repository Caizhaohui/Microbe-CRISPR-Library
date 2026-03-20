# Microbe-CRISPR-Library

Microbe-CRISPR-Library is a Python toolkit for CRISPR library design across microbial genomes, including bacterial and fungal workflows.
It contains multiple design scripts for Cas9 and CASTs applications, with versioned pipelines that preserve reproducibility while enabling iterative optimization.

---

## Overview

This repository is designed for **batch library generation** rather than one-by-one guide picking.
Typical outputs are:

- a success CSV containing final oligo-ready designs
- a failure/partial CSV summarizing genes that did not meet target design count

The toolkit supports:

- knockout library design
- knockdown library design
- promoter replacement workflows
- C-terminal fusion workflows
- CASTs insertion workflows

---

## Repository Structure and Design Logic

The codebase uses a modular, script-per-mode architecture.

### Core design principle

- Keep each biological task in a dedicated script
- Keep versioned files to preserve old behavior
- Add new constraints as forward-compatible layers

### Important scripts

- `Bact-CRISPR-Library.py`
  Main dispatcher for multi-mode usage
- `Cas9_knockout_designer_v11.py`
  Current Cas9 knockout entry point (V11, dynamic spacing adaptive version)
- `Cas9_knockdown_designer_v1.py`
- `Cas9_PromoterChange_designer_v2.py`
- `Cas9_Cfusion_designer_v1.py`
- `CASTs_designer_v3.py`
- `CRISPR_knockin_v6_standalone.py`
  Self-contained dual-mode knockin designer (`N_start` and `C_stop`) with inlined logic and no dependency on external `knockin_J23119RBS_V*.py` files

For fungal knockout library generation, use:

- `Cas9_knockout_designer_v11.py`

---

## CRISPR_knockin_v6_standalone.py (Standalone Dual-Mode Knockin)

`CRISPR_knockin_v6_standalone.py` is a fully self-contained knockin designer that supports two insertion models in one script:
- `--model N_start`
  - insert payload before the gene start codon
  - V46-equivalent start-codon targeting workflow
- `--model C_stop`
  - insert payload before the gene stop codon
  - intended for C-terminal fusion design (tag fusion use case)

The script accepts `--payload` as either a literal sequence or a FASTA file, and writes oligo-ready CSV output directly from a single standalone file.

### Core design logic

Both modes share the same high-level pipeline:

1. parse CDS features from GBFF
2. build strand-aware junction-centered sequence context
3. scan PAMs around `junction ± search_window`
4. generate candidates by strategy priority
5. apply arm sanitization and RE filtering
6. perform mutation-aware HA balancing and oligo assembly
7. rank and select top designs per gene

#### N_start mode

- Junction: gene start codon boundary
- LHA: upstream region ending at the start-codon boundary
- RHA: coding-side region starting from start codon
- Strategy family:
  - Priority1 (deletion)
  - Priority2 (bridge)
  - Priority3 (RHA mutation)

#### C_stop mode

- Junction: stop codon boundary
- LHA: coding-tail region ending before stop codon
- RHA: downstream region after stop codon
- Strategy family:
  - `CStop_P1_Del_Downstream`
  - `CStop_P2_Bridge`
  - `CStop_P2_Bridge_Mut`
  - `CStop_P3_Mut_LHA`
- Additional post-CDS validation:
  - stop codon must be one of `TAA/TAG/TGA`
- Overlap-aware handling keeps the fusion open while protecting neighboring CDS sequence when adjacent genes share bases

### Mutation strategy

Mutation logic is applied to both LHA and RHA in both modes:

- Level 1: silent/synonymous mutation first
- Level 2: conservative amino-acid-group substitution fallback

This keeps PAM-breaking behavior consistent while minimizing coding impact.

### Standalone V6 updates

- all knockin logic is inlined into one file for easier reuse and repository distribution
- genome-wide sgRNA specificity filtering is available through `--max_offtargets`
- `--rank2_sim_max` controls diversity-aware Rank2 selection in `C_stop` mode
- `--barcode_seed` controls deterministic barcode generation
- same inputs + same seed produce reproducible outputs
- balanced HA trimming, flexible HA lengths, and restriction-site exclusion are retained in the standalone workflow

### Usage examples

#### 1) N_start (start-codon insertion)

```bash
python CRISPR_knockin_v6_standalone.py \
  --model N_start \
  --gbff MG1655_genomic.gbff \
  --payload J23119_RBS \
  --template Knockin_J23100RBS_library_oligo_template.fasta \
  --output CRISPR_Nstart_v6.csv \
  --num_designs 2 \
  --lha_len 70 --rha_len 70 \
  --barcode_seed 42 \
  --max_offtargets 0 \
  --restriction_site GGTCTC --restriction_site GAAGAC
```

#### 2) C_stop (C-terminal fusion insertion before stop codon)

```bash
python CRISPR_knockin_v6_standalone.py \
  --model C_stop \
  --gbff MG1655_genomic.gbff \
  --payload J23119_RBS \
  --template Cfusion_library_oligo_template.fasta \
  --output CRISPR_Cstop_v6.csv \
  --num_designs 2 \
  --lha_len 70 --rha_len 70 \
  --barcode_seed 42 \
  --rank2_sim_max 50 \
  --max_offtargets 0 \
  --restriction_site GGTCTC --restriction_site GAAGAC
```

#### 3) Debug a single gene

```bash
python CRISPR_knockin_v6_standalone.py \
  --model C_stop \
  --target_gene b0002 \
  --gbff MG1655_genomic.gbff \
  --payload J23119_RBS \
  --template Cfusion_library_oligo_template.fasta \
  --output debug_b0002_cstop.csv
```

---

## Cas9 Knockout Pipeline (Conceptual Flow)

`Cas9_knockout_designer_v11.py` wraps and extends the evolutionary chain from v9/v10.

1. **Parse input genome and annotations**
   - FASTA+GFF3 or GBFF mode
2. **Build gene/CDS coordinate model**
   - strand-aware 5' information
3. **Enumerate sgRNA candidates**
   - PAM scanning
   - strand-aware reverse complement handling
   - restriction-site filtering
4. **Compute cut site per sgRNA**
5. **Generate deletion candidates**
   - strategy depends on `--dele_model`
6. **Build homology arms**
   - constrained by sequence availability and oligo budget
7. **Generate barcodes**
   - uniqueness + GC + repeat + restriction constraints
8. **Assemble final synthesis oligo**
   - via template placeholders
9. **Rank candidates and select per-gene top designs**
10. **Write success and failure CSV outputs**

---

## What V9 Adds

V9 introduces a strategy router:

```text
--dele_model {normal, Mt}
```

- default: `normal`
- `Mt`: force Mt deletion strategy (PAM-direction cut-window logic)
- `normal`: force legacy/V7-style length-constrained deletion logic

Startup audit lines are printed for traceability:

- `[V9审计] 当前采用 Mt 删除策略`
- `[V9审计] 当前采用 normal 删除策略`

---

## What V11 Adds

V11 builds on the V10 global-candidate architecture with additional optimizations.

### 1. Dynamic design spacing

V11 introduces CDS-length-aware spacing between multiple designs for the same gene:

- CDS length `>= 2 × min_design_spacing`: keep the standard spacing constraint
- CDS length `< 2 × min_design_spacing`: automatically switch spacing to `0 bp`

With the default `--min_design_spacing 100`, this means:

- CDS `>= 200 bp`: target `100 bp` spacing between designs
- CDS `< 200 bp`: allow overlapping designs to maximize two-design coverage

### 2. Gradient fallback

For genes that still cannot satisfy the preferred spacing, V11 keeps the staged fallback logic:

- first pass: target dynamic spacing
- second pass: target half spacing
- final pass: `0 bp`

### 3. Better behavior for dense candidate pools

When spacing falls back to `0`, V11 still tries to maximize dispersion across already selected cut sites rather than simply taking adjacent top-ranked candidates.

---

## Installation

Recommended:

- Python 3.8+
- Dependencies:
  - pandas
  - biopython
  - gffutils

```bash
pip install pandas biopython gffutils
```

---

## Input Modes

Choose one mode:

- FASTA + GFF3
  - `--input_fna`
  - `--input_gff`
- GBFF
  - `--input_gbff`

Do not mix FASTA/GFF with GBFF in the same run.

---

## Output Files

For `--output X.csv`, the pipeline writes:

- `X.csv`
  successful designs
- `X_failed.csv`
  failed/partial genes

Status definitions:

- `Failed`: no valid design found
- `Partial`: fewer designs than `--sgRNA_num`

---

## Key Parameters (V11 Knockout)

- `--output`
- `--species`
- `--synthesis_template`
- `--sgRNA_num`
- `--barcode_len`
- `--restriction_site`
- `--max_oligo_length`
- `--num_workers`
- `--dele_model {normal,Mt}`
- `--deletion_mode {auto,cut_window,legacy_length}`
- `--min_design_spacing` (V11: short CDS genes auto-relax to `0 bp`)

Deletion-specific:

- Mt / cut_window route:
  - `--cut_window` (e.g. `20:100`)
- normal / legacy_length route:
  - `--del_length_per` (e.g. `10%:80%`)
  - `--del_length_bp` (e.g. `300:1000`)

---

## Usage Examples

### 1) Mt strategy (PAM-direction cut-window)

```bash
python Cas9_knockout_designer_v11.py \
  --dele_model Mt \
  --input_gbff Mt_genomic.gbff \
  --output Mt_V11_Mt_KO.csv \
  --synthesis_template Mt_knockout_library_oligo_template.txt \
  --species M_thermophila \
  --barcode_len 11 \
  --max_oligo_length 300 \
  --restriction_site GGTCTC GAAGAC
```

### 2) normal strategy (legacy length constraints)

```bash
python Cas9_knockout_designer_v11.py \
  --dele_model normal \
  --input_gbff Mt_genomic.gbff \
  --output Mt_V11_normal_KO.csv \
  --synthesis_template Mt_knockout_library_oligo_template.txt \
  --species M_thermophila \
  --barcode_len 11 \
  --max_oligo_length 300 \
  --restriction_site GGTCTC GAAGAC \
  --del_length_per 10%:80% \
  --del_length_bp 300:1000
```

### 3) Override V11 base spacing

```bash
python Cas9_knockout_designer_v11.py \
  --input_gbff Mt_genomic.gbff \
  --output Mt_V11_spacing80.csv \
  --synthesis_template Mt_knockout_library_oligo_template.txt \
  --species M_thermophila \
  --dele_model Mt \
  --min_design_spacing 80
```

With `--min_design_spacing 80`, genes with CDS shorter than `160 bp` will automatically run with `0 bp` spacing.
