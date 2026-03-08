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
- Keep versioned files (v2/v3/.../v9) to preserve old behavior
- Add new constraints as forward-compatible layers

### Important scripts

- `Bact-CRISPR-Library.py`  
  Main dispatcher for multi-mode usage
- `Cas9_knockout_designer_v2.py` ... `Cas9_knockout_designer_v9.py`  
  Cas9 knockout evolution line
- `Cas9_knockdown_designer_v1.py`
- `Cas9_PromoterChange_designer_v2.py`
- `Cas9_Cfusion_designer_v1.py`
- `CASTs_designer_v3.py`

For fungal knockout library generation, use:

- `Cas9_knockout_designer_v9.py`

---

## Cas9 Knockout Pipeline (Conceptual Flow)

`Cas9_knockout_designer_v9.py` wraps and extends the v8 engine.

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

## Key Parameters (V9 Knockout)

- `--output`
- `--species`
- `--synthesis_template`
- `--sgRNA_num`
- `--barcode_len`
- `--restriction_site`
- `--max_oligo_length`
- `--num_workers`
- `--dele_model {normal,Mt}`

Deletion-specific:

- Mt route:
  - `--cut_window` (e.g. `20:100`)
- normal route:
  - `--del_length_per` (e.g. `10%:80%`)
  - `--del_length_bp` (e.g. `300:1000`)

---

## Usage Examples

### 1) Mt strategy (PAM-direction cut-window)

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

### 2) normal strategy (legacy length constraints)

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

### 3) FASTA + GFF mode

```bash
python Cas9_knockout_designer_v9.py \
  --dele_model normal \
  --input_fna genome.fna \
  --input_gff genome.gff \
  --output KO_from_fna_gff.csv \
  --synthesis_template Mt_knockout_library_oligo_template.txt \
  --species E_coli \
  --del_length_per 20%:80% \
  --barcode_len 10
```

---

## Practical Guidance

- Use `--dele_model Mt` when you want the Mt window strategy and PAM-direction-aware deletion behavior.
- Use `--dele_model normal` for cross-species workflows with explicit length constraints.
- In `normal` mode, always pass at least one of:
  - `--del_length_per`
  - `--del_length_bp`
- If many short genes are filtered out, reduce the minimum deletion pressure (especially in `--del_length_bp`).

---

## Troubleshooting

### Why does `normal` mode fail when I omit deletion bounds?

Because `normal` routes to legacy-length logic, which requires deletion constraints.

### How can I verify active strategy quickly?

Check startup audit output:

- `当前采用 Mt 删除策略`
- `当前采用 normal 删除策略`

### Why do I see many `Partial` genes?

Typical reasons:

- strict restriction-site filtering
- too narrow deletion constraints
- oligo-length budget leaves insufficient homology-arm space

---

## Notes

This toolkit generates computational design candidates.  
All designs should be experimentally validated in the target strain and editing protocol.
