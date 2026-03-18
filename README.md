# Microbe-CRISPR-Library

Microbe-CRISPR-Library is a Python toolkit for CRISPR library design across microbial genomes. The current Cas9 knockout workflow in this repository is centered on V11, with older Cas9 knockout version files removed so the maintained entry point is unambiguous.

---

## Current Cas9 Knockout Entry Point

Use `Cas9_knockout_designer_v11.py` for knockout library design.

V11 combines the V10 global-candidate architecture with an additional dynamic spacing policy for short CDS genes, while preserving:

- dual input modes: FASTA+GFF3 or GBFF
- Mt-aware deletion strategy selection via `--dele_model`
- cut-window and legacy-length deletion modes
- thread-safe global barcode allocation
- final oligo assembly from a synthesis template
- success and failure/partial CSV outputs

---

## What V11 Adds

### 1. Dynamic design spacing

V11 introduces CDS-length-aware spacing between multiple designs for the same gene:

- CDS length `>= 2 x min_design_spacing`: keep the standard spacing constraint
- CDS length `< 2 x min_design_spacing`: automatically switch spacing to `0 bp`

With the default `--min_design_spacing 100`, this means:

- CDS `>= 200 bp`: target `100 bp` spacing
- CDS `< 200 bp`: allow overlapping designs to maximize two-design coverage

### 2. Gradient fallback remains in place

For genes that still cannot satisfy the preferred spacing, V11 keeps the staged fallback logic:

- first pass: target dynamic spacing
- second pass: target half spacing
- final pass: `0 bp`

### 3. Better behavior for dense candidate pools

When spacing falls back to `0`, V11 still tries to maximize dispersion across already selected cut sites rather than simply taking adjacent top-ranked candidates.

---

## Other Maintained Designer

The repository also contains `CRISPR_knockin.py`, a dual-mode knockin designer for start-codon (`N_start`) and stop-codon (`C_stop`) insertion workflows.

- `N_start`: insert payload before the gene start codon
- `C_stop`: insert payload before the gene stop codon for C-terminal fusion style designs
- supports mutation-aware homology-arm balancing, restriction-site filtering, and deterministic barcode generation via `--barcode_seed`

Example:

```bash
python CRISPR_knockin.py \
  --model C_stop \
  --gbff MG1655_genomic.gbff \
  --template Cfusion_library_oligo_template.fasta \
  --output CRISPR_Cstop_v1.csv \
  --lha_len 70 --rha_len 70 \
  --barcode_seed 42 \
  --restriction_site GGTCTC --restriction_site GAAGAC
```

---

## Installation

Recommended:

- Python 3.8+
- pandas
- biopython
- gffutils

```bash
pip install pandas biopython gffutils
```

---

## Input Modes

Choose one of the following:

- FASTA + GFF3
  - `--input_fna`
  - `--input_gff`
- GBFF
  - `--input_gbff`

Do not mix FASTA/GFF input with GBFF input in the same run.

---

## Output Files

For `--output X.csv`, V11 writes:

- `X.csv`: successful designs
- `X_failed.csv`: failed or partial genes

Successful output includes V11 multi-design tracking fields such as:

- `Design Index`
- `Cut Site Spacing`
- `Deletion Warning`

Status definitions:

- `Failed`: no valid design found
- `Partial`: fewer designs than `--sgRNA_num`

---

## Key Parameters

- `--output`: output CSV path
- `--species`: currently supports `M_thermophila` and `E_coli`
- `--synthesis_template`: oligo template file
- `--sgRNA_num`: designs per gene, default `2`
- `--barcode_len`: barcode length
- `--restriction_site`: restriction sites to avoid
- `--max_oligo_length`: total oligo length cap
- `--num_workers`: thread count
- `--dele_model {normal,Mt}`: deletion strategy
- `--deletion_mode {auto,cut_window,legacy_length}`: deletion candidate generation mode
- `--min_design_spacing`: V11 base spacing; short CDS genes can auto-relax to `0 bp`

Deletion-specific parameters:

- `cut_window` mode:
  - `--cut_window`, e.g. `20:100`
- `legacy_length` mode:
  - `--del_length_per`, e.g. `10%:80%`
  - `--del_length_bp`, e.g. `300:1000`

---

## Usage Examples

### 1. M. thermophila GBFF mode with Mt strategy

```bash
python Cas9_knockout_designer_v11.py \
  --input_gbff Mt_genomic.gbff \
  --output Mt_V11_KO.csv \
  --synthesis_template Mt_knockout_library_oligo_template.txt \
  --species M_thermophila \
  --dele_model Mt \
  --deletion_mode cut_window \
  --cut_window 20:100 \
  --barcode_len 11 \
  --max_oligo_length 300 \
  --restriction_site GGTCTC GAAGAC \
  --num_workers 8
```

### 2. FASTA + GFF mode with explicit legacy deletion bounds

```bash
python Cas9_knockout_designer_v11.py \
  --input_fna genome.fna \
  --input_gff genome.gff \
  --output KO_from_fna_gff.csv \
  --synthesis_template Mt_knockout_library_oligo_template.txt \
  --species E_coli \
  --dele_model normal \
  --deletion_mode legacy_length \
  --del_length_per 20%:80% \
  --del_length_bp 300:1000 \
  --barcode_len 10
```

### 3. Override V11 base spacing

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

---

## Practical Guidance

- Use `--dele_model Mt` for the M. thermophila strand-aware deletion strategy.
- Use `--deletion_mode auto` if you want the species-default mode selected automatically.
- In `legacy_length` mode, you must provide at least one of `--del_length_per` or `--del_length_bp`.
- If many genes remain `Partial`, first inspect restriction-site filters, oligo-length budget, and deletion bounds.
- If you want stricter multi-design diversity on long genes, increase `--min_design_spacing`.

---

## Conceptual Pipeline

1. Parse genome and annotations from FASTA+GFF3 or GBFF.
2. Build CDS-aware gene models with strand-specific 5' coordinates.
3. Enumerate sgRNA candidates and score them.
4. Compute cut sites using the selected deletion strategy.
5. Generate deletion candidates in cut-window or legacy-length mode.
6. Build a global candidate pool per gene.
7. Select multiple designs using V11 dynamic spacing and fallback rules.
8. Construct homology arms and synthesis oligos.
9. Generate globally unique barcodes across threads.
10. Write success and failure CSV outputs.

---

## Notes

This repository generates computational design candidates. Experimental validation is still required for each target strain, genome annotation set, and editing workflow.
