# Bact-CRISPR-Library: Bacterial CRISPR Library Design Tool

`Bact-CRISPR-Library` is a powerful command-line Python tool for automated, high-throughput design of CRISPR-based oligonucleotide libraries for bacterial genomes. It supports multiple CRISPR editing applications, including gene knockouts, promoter replacements, protein C-terminal fusions, gene knockdowns, and INTEGRATE/CASTs-based designs, with intelligent optimization strategies to ensure high-quality and high-success-rate libraries.

## ✨ Features

- **Design Modes**:
  - `Knockout_Cas9`: Gene knockout library with customizable deletion lengths.
  - `PromoterChange_Cas9`: Promoter replacement library for gene regulation.
  - `Cfusion_Cas9`: Protein C-terminal fusion tag library.
  - `Knockdown_Cas9`: Gene knockdown library via regulatory element insertion.
  - `Knockout_CASTs`: Gene knockout library for INTEGRATE/CASTs systems.
  - `PromoterChange_CASTs`: Promoter insertion library for INTEGRATE/CASTs systems.

- **Intelligent sgRNA Design**: Uses a model-based scoring system (MODEL_WEIGHTS) to prioritize optimal sgRNA cutting sites.
- **Flexible Homology Arm Design**: Supports customizable homology arm length ranges with safety checks to protect upstream genes in promoter replacement mode.
- **Circular and Linear Genome Support**: Handles both genome types via the `--genome_type` parameter, ensuring correct boundary processing.
- **Restriction Site Avoidance**: Avoids user-specified restriction sites in variable regions (sgRNA, homology arms, barcode) while exempting a designated cloning site.
- **Automated Barcode Generation**: Generates unique barcodes with GC content (30-70%) and sequence complexity constraints.
- **Template-Based Oligo Synthesis**: Outputs oligonucleotides based on user-provided templates for high flexibility.
- **Comprehensive Output**: Saves detailed design results in CSV format, including sgRNA sequences, homology arms, barcodes, and final oligonucleotides.

## ⚙️ Requirements and Installation

### System Requirements
- **Python Version**: Python 3.8 or higher
- **Dependencies**:
  - `pandas`: For data handling and CSV output
  - `biopython`: For sequence manipulation and FASTA parsing
  - `gffutils`: For GFF3 annotation parsing

### Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/bact-crispr-library.git
   cd bact-crispr-library
   ```

2. Install dependencies:
   ```bash
   pip install pandas biopython gffutils
   ```

3. Prepare input files (FASTA, GFF3, and synthesis template).

## 🚀 Usage

Run the script via the command line with a specified design mode and parameters:

```bash
python bact_crispr_library.py <mode> [common_options] [mode_specific_options]
```

### Common Parameters
These apply to all design modes:

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--input_fna` | **[Required]** Path to the genome FASTA file. | `--input_fna genome.fna` |
| `--input_gff` | **[Required]** Path to the GFF3 annotation file. | `--input_gff annotation.gff` |
| `--output` | **[Required]** Path for the output CSV file. | `--output results.csv` |
| `--synthesis_template` | **[Required]** Path to the oligo template file. | `--synthesis_template template.txt` |
| `--genome_type` | Genome type: `circle` or `linear` (default: `linear`). | `--genome_type circle` |
| `--sgRNA_num` | Number of designs per gene (default: 1). | `--sgRNA_num 2` |
| `--barcode_len` | Barcode length (default: 8). | `--barcode_len 8` |
| `--restriction_site` | Space-separated restriction sites to avoid in variable regions. | `--restriction_site EcoRI HindIII` |
| `--cloning_site` | Optional restriction site allowed in the template for cloning. | `--cloning_site GGTCTC` |

### Mode-Specific Parameters and Examples

#### 1. Knockout_Cas9
Designs small deletions within gene coding regions.

**Parameters**:
- `--pam`: PAM sequence (default: `NGG`).
- `--HR_len`: Homology arm length range (`preferred:min`, default: `50:40`).
- `--del_length`: Deletion length range (`min:max`, default: `50:100`).
- `--promoter_region`: Protected promoter region size (`min:max`, default: `50:150`).
- `--ko_search_range`: CDS percentage range for sgRNA search (`min:max`, default: `5:80`).

**Example**:
```bash
python bact_crispr_library.py Knockout_Cas9 \
  --input_fna ecoli.fna \
  --input_gff ecoli.gff \
  --output ecoli_knockout.csv \
  --genome_type circle \
  --pam NGG \
  --HR_len 46:40 \
  --del_length 20:100 \
  --promoter_region 50:150 \
  --ko_search_range 5:80 \
  --sgRNA_num 2 \
  --restriction_site gaagac ggtctc \
  --cloning_site AAGGTGAGACCAAGGTCTCTCGAC \
  --synthesis_template knockout_template.txt
```

#### 2. PromoterChange_Cas9
Replaces promoters upstream of gene start codons.

**Parameters**:
- `--pam`: PAM sequence (default: `NGG`).
- `--HR_len`: Homology arm length range (`preferred:min`, default: `50:40`).
- `--promoter_search_size`: Upstream region size for sgRNA search (default: 150 bp).

**Example**:
```bash
python bact_crispr_library.py PromoterChange_Cas9 \
  --input_fna ecoli.fna \
  --input_gff ecoli.gff \
  --output ecoli_promoter.csv \
  --genome_type circle \
  --pam NGG \
  --HR_len 50:45 \
  --promoter_search_size 150 \
  --restriction_site gaagac \
  --cloning_site GGTCTCAAGCTT \
  --synthesis_template promoter_template.txt
```

#### 3. Cfusion_Cas9
Inserts tags (e.g., fluorescent proteins) at protein C-termini.

**Parameters**:
- `--pam`: PAM sequence (default: `NGG`).
- `--HR_len`: Homology arm length range (`preferred:min`, default: `50:40`).
- `--sgrna_search_range`: Search range around stop codon (default: 50 bp).

**Example**:
```bash
python bact_crispr_library.py Cfusion_Cas9 \
  --input_fna ecoli.fna \
  --input_gff ecoli.gff \
  --output ecoli_cfusion.csv \
  --genome_type circle \
  --pam NGG \
  --HR_len 50:45 \
  --sgrna_search_range 40 \
  --restriction_site gaagac \
  --cloning_site GGTCTC \
  --synthesis_template cfusion_template.txt
```

#### 4. Knockdown_Cas9
Inserts regulatory elements for gene knockdown.

**Parameters**:
- `--pam`: PAM sequence (default: `NGG`).
- `--HR_len`: Homology arm length range (`preferred:min`, default: `50:40`).
- `--insertion_distance`: Distance range upstream of start codon (`min:max`, default: `20:30`).
- `--sgrna_search_range`: Search range around insertion site (default: 50 bp).

**Example**:
```bash
python bact_crispr_library.py Knockdown_Cas9 \
  --input_fna ecoli.fna \
  --input_gff ecoli.gff \
  --output ecoli_knockdown.csv \
  --genome_type circle \
  --pam NGG \
  --HR_len 50:45 \
  --insertion_distance 20:30 \
  --sgrna_search_range 50 \
  --restriction_site gaagac \
  --cloning_site GGTCTC \
  --synthesis_template knockdown_template.txt
```

#### 5. Knockout_CASTs
Designs knockouts for INTEGRATE/CASTs systems.

**Parameters**:
- `--target_cds_range`: CDS percentage range for sgRNA search (`min:max`, default: `5:80`).

**Example**:
```bash
python bact_crispr_library.py Knockout_CASTs \
  --input_fna ecoli.fna \
  --input_gff ecoli.gff \
  --output ecoli_casts_knockout.csv \
  --genome_type circle \
  --target_cds_range 5:80 \
  --restriction_site gaagac \
  --cloning_site GGTCTC \
  --synthesis_template casts_knockout_template.txt
```

#### 6. PromoterChange_CASTs
Designs promoter insertions for INTEGRATE/CASTs systems.

**Parameters**:
- `--insertion_range_promoter`: Distance range upstream of start codon (`min:max`, default: `25:50`).
- `--sgrna_search_window`: Search window around ideal target (default: 20 bp).

**Example**:
```bash
python bact_crispr_library.py PromoterChange_CASTs \
  --input_fna ecoli.fna \
  --input_gff ecoli.gff \
  --output ecoli_casts_promoter.csv \
  --genome_type circle \
  --insertion_range_promoter 25:50 \
  --sgrna_search_window 20 \
  --restriction_site gaagac \
  --cloning_site GGTCTC \
  --synthesis_template casts_promoter_template.txt
```

## 📄 Input File Formats

### 1. Genome Sequence (`--input_fna`)
Standard FASTA format. For circular genomes, provide a linear sequence, and the tool handles circularity via `--genome_type circle`.

### 2. Gene Annotation (`--input_gff`)
Standard GFF3 format with `gene` and `CDS` features. Each `gene` must have a `locus_tag` attribute as a unique identifier.

### 3. Synthesis Template (`--synthesis_template`)
A text file defining the oligonucleotide structure with placeholders:
- `{sgRNA_fwd}`: sgRNA sequence (20 nt for Cas9, 32 nt for CASTs).
- `{sgRNA_rc}`: Reverse complement of sgRNA.
- `{pam}`: PAM sequence.
- `{pam_rc}`: Reverse complement of PAM.
- `{upstream_arm}`: Upstream homology arm.
- `{downstream_arm}`: Downstream homology arm.
- `{barcode}`: Randomly generated barcode.
- `{exempt_restriction_site}`: Cloning site specified by `--cloning_site`.
- `{insert}`: Insertion sequence (used in Cfusion_Cas9).

**Example Template** (`promoter_template.txt`):
```
GACGACTTGCTATTTCTAGCTCTAAAAC{sgRNA_fwd}GTTAAACCCTATAGTGAGTCGTATTAC{barcode}AGGTGAGACC{exempt_restriction_site}GAC{upstream_arm}{downstream_arm}ACTCGCACTGCTGGTCACTTGGTCTGGA
```

## 📤 Output File Formats

### 1. Main Output (`<output>.csv`)
Contains successful designs with columns:
- `Gene_ID`: Gene locus_tag.
- `Status`: Design status ("Success").
- `sgRNA_Sequence`: Designed sgRNA sequence.
- `Barcode`: Generated barcode.
- `Final_Oligo_for_Synthesis`: Complete oligonucleotide for synthesis.
- `sgRNA_Score`: sgRNA activity score (Cas9 modes).
- `Arm_Length`: Homology arm length (Cas9 modes).
- `Upstream_Arm`, `Downstream_Arm`: Homology arm sequences (Cas9 modes).
- `sgRNA_PAM`: PAM sequence.
- `sgRNA_Strand`: sgRNA strand (forward/reverse).
- `sgRNA_Cut_Site`: Genomic cut site (Cas9 modes).
- `Design_Scenario`: Design strategy details (PromoterChange_Cas9, Knockdown_Cas9).
- `Deletion_Length`, `Deletion_Start`, `Deletion_End`: Deletion details (Knockout_Cas9).
- `Design_Strategy`, `Target_Site`, `Predicted_Insertion_Site`: CASTs-specific fields.

### 2. Failed Genes (`<output>_failed.csv`)
Lists genes that failed design with `Gene_ID` and `Status` ("Failed").

## 📜 License

This project is licensed under the [MIT License](LICENSE).

## 📧 Contact

- **Author**: Zhaohui Cai
- **Email**: cai_zhaohui@163.com
- **GitHub Issues**: Open an issue for questions or support.

## 🛠️ Contributing

Contributions are welcome! To contribute:
1. Fork the repository.
2. Create a feature branch (`git checkout -b feature/your-feature`).
3. Commit changes (`git commit -m 'Add your feature'`).
4. Push to the branch (`git push origin feature/your-feature`).
5. Open a pull request.
