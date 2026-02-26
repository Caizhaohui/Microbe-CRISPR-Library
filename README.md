# CRISPR Knockout Library Designer for *Myceliophthora thermophila*

## Overview

This project implements a CRISPR-Cas9 knockout library designer for *Myceliophthora thermophila* (SpCas9 with NGG PAM). The designer generates multiple homology arm-based knockout constructs for each target gene, with synthesis quality constraints.

## Features

### V2 (Production-Ready)
- **Multi-contig genome support**: Handles genomes split across multiple contigs
- **GBFF input mode**: Load genome features directly from GenBank format files
- **Multi-guide design**: Generate 2 independent gRNAs per gene with automated fallback
- **Performance optimized**: PAM search via string matching, reverse complement via string translation
- **Flexible parameters**: 
  - `--HR_len` (preferred:minimum homology arm lengths)
  - `--del_length` (deletion size range)
  - `--barcode_len` (barcode length, default 8)
  - `--synthesis_template` (custom oligo template)

### V3 (Oligo Length Constraint)
- **All V2 features**, plus:
- **Exact oligo length targeting**: All synthesized oligos are exactly `--max_oligo_length` bp (default 300)
- **Variable homology arms**: Arm lengths automatically adjusted and maximized within oligo budget
- **Commercial synthesis optimization**: Designed for fixed-length oligo library synthesis

## Usage

### V2 Example
```bash
python Cas9_knockout_designer_v2.py \
  --input_gbff Mt_genomic.gbff \
  --output knockout_library.csv \
  --synthesis_template Mt_knockout_library_oligo_template.txt \
  --species M_thermophila \
  --HR_len 75:65 \
  --del_length 100:220 \
  --barcode_len 10
```

### V3 Example (with oligo length constraint)
```bash
python Cas9_knockout_designer_v3.py \
  --input_gbff Mt_genomic.gbff \
  --output knockout_library_300bp.csv \
  --synthesis_template Mt_knockout_library_oligo_template.txt \
  --species M_thermophila \
  --HR_len 75:65 \
  --del_length 100:220 \
  --barcode_len 10 \
  --max_oligo_length 300
```

## Input Files

- **`Mt_genomic.gbff`**: GenBank format genome file with gene features
- **`Mt_knockout_library_oligo_template.txt`**: Oligo synthesis template with placeholders:
  - `{sgRNA_rc}`: Reverse complement of 20bp guide
  - `{upstream_arm}`: Variable-length homology arm
  - `{barcode}`: Barcode sequence
  - `{downstream_arm}`: Variable-length homology arm

## Output

CSV file with columns:
- `Gene Id`: Target gene identifier
- `Gene Name`: Target gene name
- `sgRNA sequence`: 20bp guide RNA sequence
- `Upstream Arm Length`, `Downstream Arm Length`: Final arm lengths
- `Final Oligo For Synthesis`: Complete synthesis oligo
- Other metadata (barcode, deletion boundaries, etc.)

## Parameters

| Parameter | V2 Default | V3 Default | Description |
|-----------|-----------|-----------|-------------|
| `--HR_len` | `75:65` | `75:65` | Preferred and minimum homology arm lengths (preferred:minimum format) |
| `--del_length` | `200:300` | `200:300` | Deletion size range (min:max) |
| `--barcode_len` | `8` | `8` | Barcode sequence length |
| `--max_oligo_length` | N/A | `300` | Target oligo length (V3 only) |
| `--synthesis_template` | required | required | Path to oligo template file |
| `--species` | `M_thermophila` | `M_thermophila` | Species identifier |
| `--input_gbff` | - | - | GenBank format genome file |

## Performance

- **Speed**: ~3 minutes for full *M. thermophila* genome (9,097 genes)
- **Coverage**: ~9,100 genes, ~18,200 total designs (2 per gene)
- **Success rate**: >99% design success

## Version History

- **V1**: Initial single-contig implementation
- **V2**: Multi-contig support, GBFF input, performance optimization
- **V3**: Oligo length constraint, variable homology arms

## Requirements

- Python 3.6+
- BioPython
- gffutils
- pandas

## License

[To be determined]

## Author

Caizhaohui
