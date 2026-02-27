# Bact-CRISPR-Library Toolkit Overview

The Bact-CRISPR-Library is a versatile Python-based toolkit tailored for designing CRISPR-based gene editing libraries in bacteria, such as *Escherichia coli* and *Streptomyces*. It supports both Cas9-mediated precise editing and CASTs (CRISPR-associated transposase) systems for transposon-based insertions. The central script serves as a dispatcher, invoking specialized sub-scripts for various modes while preserving modularity.

This toolkit streamlines the process of generating sgRNA libraries for applications like gene knockdown, promoter replacement, C-terminal fusions, and knockouts. Designs are customizable, with outputs in CSV format for easy integration into downstream workflows. While effective for bacterial genomes, all designs should undergo experimental validation to account for off-target effects or strain-specific variations.

### Key Features
- **Modular Architecture**: A single entry-point script dispatches tasks to independent sub-scripts.
- **Dual Systems Support**: Cas9 for homology-directed repair and CASTs for efficient insertions.
- **Input Handling**: Processes FASTA genomes (linear/circular) and GFF3 annotations.
- **Output Formats**: CSV files detailing successful designs and failures, including sgRNAs, barcodes, and synthesized oligos.
- **Customization Options**: Adjustable parameters for sgRNA count, barcode length, restriction site avoidance, and mode-specific constraints.

### Installation Guide
Clone the repository and install dependencies using pip:
```
pip install pandas biopython gffutils
```
Ensure Python 3.8+ is installed. No additional tools are required beyond these libraries.

### Quick Usage Examples
Run the toolkit with a mode and required inputs:
```
python Bact-CRISPR-Library.py Knockdown_Cas9 --input_fna ecoli.fna --input_gff ecoli.gff --output knockdown_results.csv --synthesis_template template.txt --sgRNA_num 2
```
For CASTs-based knockout:
```
python Bact-CRISPR-Library.py Knockout_CASTs --input_fna strepto.fna --input_gff strepto.gff --output casts_ko.csv --synthesis_template template.txt --target_cds_range 10:70
```
View general help: `python Bact-CRISPR-Library.py -h`. Mode-specific help: `python Bact-CRISPR-Library.py PromoterChange_Cas9 -h`.

### Common Parameters
These apply across all modes:
- `--input_fna`: Genome FASTA file (required).
- `--input_gff`: GFF3 annotation file (required).
- `--output`: Path for CSV output (required).
- `--genome_type`: 'linear' or 'circle' (default: linear).
- `--sgRNA_num`: Number of designs per gene (default: 1).
- `--barcode_len`: Barcode sequence length (default: 10).
- `--restriction_site`: Space-separated sites to avoid.
- `--synthesis_template`: Oligo structure template file (required).
- `--cloning_site`: Exempt restriction site.
- `--relax`: Bypass checks for adjacent gene disruptions (where supported).

---

## Fungal Knockout Designer (V2 & V3) - Specialized Implementation

This repository provides optimized V2 and V3 implementations of the Cas9-based knockout designer with enhanced performance and synthesis constraints, specifically tailored for fungal genome applications.

### V2 (Production-Ready)
**File**: `Cas9_knockout_designer_v2.py`

Features:
- **Multi-contig genome support**: Handles genomes split across multiple contigs
- **GBFF input mode**: Load genome features directly from GenBank format files
- **Multi-guide design**: Generate 2 independent gRNAs per gene with automated fallback
- **Performance optimized**: PAM search via string matching, reverse complement via string translation
- **Flexible parameters**: 
  - `--HR_len` (preferred:minimum homology arm lengths, e.g., `75:65`)
  - `--del_length` (deletion size range, e.g., `100:220`)
  - `--barcode_len` (barcode length, default 8)
  - `--synthesis_template` (custom oligo template)

**Usage**:
```bash
python Cas9_knockout_designer_v2.py \
  --input_gbff genome.gbff \
  --output knockout_library.csv \
  --synthesis_template knockout_library_oligo_template.txt \
  --species fungal_species \
  --HR_len 75:65 \
  --del_length 100:220 \
  --barcode_len 10
```

### V3 (Oligo Length Constraint)
**File**: `Cas9_knockout_designer_v3.py`

All V2 features, plus:
- **Exact oligo length targeting**: All synthesized oligos are exactly `--max_oligo_length` bp (default 300bp)
- **Variable homology arms**: Arm lengths automatically adjusted and maximized within oligo budget
- **Commercial synthesis optimization**: Designed for fixed-length oligo library synthesis

**Usage**:
```bash
python Cas9_knockout_designer_v3.py \
  --input_gbff genome.gbff \
  --output knockout_library_300bp.csv \
  --synthesis_template knockout_library_oligo_template.txt \
  --species fungal_species \
  --HR_len 75:65 \
  --del_length 100:220 \
  --barcode_len 10 \
  --max_oligo_length 300
```

### Output Format

CSV file with columns:
- `Gene Id`: Target gene identifier
- `Gene Name`: Target gene name
- `sgRNA sequence`: 20bp guide RNA sequence
- `Upstream Arm Length`, `Downstream Arm Length`: Final arm lengths
- `Final Oligo For Synthesis`: Complete synthesis oligo
- `Barcode`: Unique barcode sequence
- Other metadata (deletion boundaries, etc.)

### V2 vs V3 Comparison

| Aspect | V2 | V3 |
|--------|----|----|
| **Oligo Length** | Variable (≤template) | Exactly `max_oligo_length` bp |
| **Arm Lengths** | Fixed preferred/min search | Dynamic, maximized within budget |
| **Use Case** | General knockout library | Commercial synthesis (fixed-length pools) |
| **Speed** | Scales with genome size | Scales with genome size |
| **Coverage** | 2 designs per gene | 2 designs per gene |

---

## Fungal Knockout Designer (V4–V6) - Current Workflow

V6 is the maintained version for fungal knockout library design in this workspace.

**File**: `Cas9_knockout_designer_v6.py`

Key improvements since V3:
- **Deletion length by CDS percent and/or bp**: `--del_length_per` and `--del_length_bp` (supports single value as max).
- **Start-codon-aware scoring**: prioritizes deletions that remove the start codon and cause frameshift.
- **Barcode GC range**: widened to 30%-80% by default.
- **Exact oligo length**: set `--max_oligo_length 300` to keep all oligos fixed-length.

**Example (V6)**:
```bash
python Cas9_knockout_designer_v6.py \
  --input_gbff Mt_genomic.gbff \
  --output Mt_KO_library_v6_delPct10_80_delBp300_1000_bc11.csv \
  --synthesis_template Mt_knockout_library_oligo_template.txt \
  --species M_thermophila \
  --barcode_len 11 \
  --del_length_per 10%:80% \
  --del_length_bp 300:1000 \
  --max_oligo_length 300 \
  --restriction_site GGTCTC GAAGAC
```

---

The Bact-CRISPR-Library toolkit represents a significant advancement in computational support for bacterial CRISPR engineering, bridging traditional Cas9 methodologies with emerging CASTs technologies. Developed to address the need for scalable, customizable library designs, it integrates genome processing, sgRNA optimization, and oligo synthesis into a unified framework. This comprehensive guide expands on the toolkit's architecture, implementation details, mode-specific workflows, potential extensions, and practical considerations for researchers in synthetic biology and microbiology.

## Architectural Design
At its core, the toolkit employs a dispatcher pattern implemented in `Bact-CRISPR-Library.py`. This script uses Python's `argparse` library with subparsers to handle mode-specific commands, ensuring a clean CLI interface. Upon invocation, it constructs a command list via the `build_cmd` function and executes the appropriate sub-script using `subprocess.run`. This approach decouples the user interface from core logic, allowing independent maintenance of sub-scripts.

Sub-scripts share common utilities:
- **GenomeProcessor**: Loads FASTA sequences and parses GFF3 files using Biopython and gffutils, extracting gene features like locus tags, strands, and CDS coordinates.
- **SequenceUtils**: Handles reverse complements, sequence extraction (supporting circular genomes), restriction site checks, and unique barcode generation with GC content and homopolymer filters.
- **SGRNADesigner**: Searches for sgRNAs with mode-specific PAM patterns, avoiding restriction sites.
- **CRISPRDesigner**: Core design logic, including arm searches, oligo assembly, and gene-specific validations.
- **ResultProcessor**: Generates CSV outputs for successes and failures, with detailed columns like Gene_ID, sgRNA_Sequence, and Final_Oligo_for_Synthesis.

For CASTs modes, the multi-purpose `CASTs_designer_v3.py` is invoked with a `--mode` flag (e.g., Knockout_CASTs), adapting PAM searches to CASTs-compatible sequences like AC, GC, CC.

## Mode-Specific Implementations
Each mode targets a distinct genetic manipulation strategy, with tailored parameters and filters.

### Cas9-Based Modes
- **Knockdown_Cas9**: Focuses on inserting sgRNAs upstream of start codons to reduce expression. Searches in configurable ranges (--sgrna_upstream_range), ensures minimum distances (--rha_min_dist_to_atg), and uses homology arms (--HR_len). Designs avoid CDS overlaps and restriction sites.
- **PromoterChange_Cas9**: Replaces promoters by cutting upstream and inserting new sequences. Searches in promoter regions (--promoter_search_size), with optional relaxation (--relax) for adjacent genes. Prioritizes optimal cuts with minimal offsets.
- **Cfusion_Cas9**: Enables C-terminal fusions by inserting at stop codons. Expands searches downstream (--sgrna_search_range, --max_search_expansion) and protects neighboring genes (--strict). Checks for in-frame stop codons in left homology arms.
- **Knockout_Cas9**: Induces frameshifts or deletions in CDS. Targets specific percentages (--ko_search_range), controls deletion sizes (--del_length), and safeguards promoters (--promoter_region, --strict).

### CASTs-Based Modes
- **Knockout_CASTs**: Integrates transposons into CDS for disruption. Targets mid-CDS regions (--target_cds_range) with CASTs-specific PAMs. Designs ensure insertions avoid premature stops and prioritize frameshifts.
- **PromoterChange_CASTs**: Inserts promoters upstream via transposons. Defines insertion windows (--insertion_range_promoter) and checks for upstream gene collisions unless relaxed (--relax). Predicts insertion sites based on sgRNA orientation.

All modes incorporate barcode uniqueness, oligo exemption for cloning sites, and logging for progress (e.g., gene processing every 100 entries).

## Input and Output Specifications
### Inputs
- **FASTA File**: Single contig preferred; uppercased during loading.
- **GFF3 File**: Must include 'gene' and 'CDS' features with 'locus_tag'. Handles multi-exon genes by aggregating CDS coordinates.
- **Synthesis Template**: Plain text with placeholders (e.g., {sgRNA_fwd}, {barcode}). Validated for fixed regions against restriction sites.

### Outputs
- **Success CSV**: Mode-dependent columns, e.g., for Knockdown_Cas9: Gene_ID, Design_Type, sgRNA_Seq, Barcode, Final_Oligo_for_Synthesis.
- **Failed CSV** (e.g., results_failed.csv): Gene_ID, Status, Strand, Relevant_Sequence (CDS for knockouts, upstream for promoter changes).
- **Logging**: Console outputs include genome loading stats, gene counts, and summaries (e.g., "Generated X successful designs, Y failed genes").

Failed designs often stem from insufficient sgRNA candidates or restriction site conflicts; relevant sequences aid manual troubleshooting.

## Performance and Optimization
- **Runtime**: Scales with gene count; ~1-5 seconds per gene on standard hardware. For 5,000 genes, expect 1-2 hours total.
- **Memory**: Low footprint; genome sequences loaded in memory (suitable for bacterial sizes <10 MB).
- **Optimizations**: Sub-scripts use efficient sliding windows for sgRNA searches and sets for barcode uniqueness. Circular genomes handled via sequence doubling.
- **Edge Cases**: Handles short genes (skips if below min arm lengths), overlapping genes (via strict mode), and empty inputs (raises errors).

## Extensions and Customization
- **Adding Modes**: Update subparsers, script_map, and sub_mode_map in the dispatcher. Create new sub-scripts following existing patterns.
- **PAM Customization**: Modify SGRNADesigner in sub-scripts for new patterns.
- **Multi-Threading**: Not native; parallelize by splitting gene lists externally.
- **Integration**: Outputs compatible with oligo synthesis pipelines or CRISPR design software like Benchling.

## License

[To be determined]

## Author

Caizhaohui
