# filename: casts_designer.py

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import gffutils
import tempfile
import os
import argparse
import logging
import random
import sys
from typing import List, Dict, Tuple, Optional, Set
from dataclasses import dataclass
import math

# --- 日志设置 ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- 数据类 ---
@dataclass
class DesignConfig:
    """CASTs模式的统一配置类"""
    sgrna_num: int
    barcode_len: int
    relax: bool
    synthesis_template: Optional[str] = None
    casts_target_start_pct: Optional[float] = None
    casts_target_end_pct: Optional[float] = None
    promoter_insertion_dist_min: Optional[int] = None
    promoter_insertion_dist_max: Optional[int] = None
    sgrna_search_window_ki: Optional[int] = None

@dataclass
class Gene:
    """通用基因信息容器"""
    id: str; start: int; end: int; strand: int
    cds_5prime_start: int; cds_3prime_start: int; cds_3prime_end: int
    cds_min_coord: int; cds_max_coord: int

@dataclass
class CASTsDesignResult:
    """CASTs/INTEGRATE 设计结果容器"""
    gene_id: str
    design_strategy: str
    sgrna_seq: str
    sgrna_strand: str
    pam: str
    target_site: int
    insertion_site: int
    barcode: str
    final_oligo_for_synthesis: str

# --- 通用辅助类 ---
class ConfigParser:
    @staticmethod
    def parse_range_param(param_str: str, param_name: str, sort_values: bool = True) -> Tuple[int, int]:
        try:
            parts = param_str.split(':')
            if len(parts) == 1:
                val = int(parts[0])
                return val, val
            elif len(parts) == 2:
                val1, val2 = int(parts[0]), int(parts[1])
                return sorted((val1, val2)) if sort_values else (val1, val2)
            else:
                raise ValueError("输入必须是单个数字或由冒号分隔的两个数字。")
        except (ValueError, IndexError) as e:
            raise ValueError(f"无效的 {param_name} 格式 '{param_str}'. 错误: {e}") from e

class GenomeProcessor:
    def __init__(self, genome_file: str, gff_file: str):
        self.genome_file = genome_file
        self.gff_file = gff_file
        self.genome_seq = None
        self.genome_len = 0
    
    def load_genome(self) -> None:
        logging.info(f"正在加载基因组: {self.genome_file}...")
        try:
            genome_record = SeqIO.to_dict(SeqIO.parse(self.genome_file, "fasta"))
            if not genome_record: raise ValueError("基因组文件中没有序列")
            chrom_id = list(genome_record.keys())[0]
            self.genome_seq = genome_record[chrom_id].seq.upper()
            self.genome_len = len(self.genome_seq)
            logging.info(f"基因组加载完成: {self.genome_len} bp")
        except Exception as e: logging.error(f"加载基因组失败: {e}"); raise

    def parse_genes(self) -> List[Gene]:
        logging.info(f"正在解析GFF: {self.gff_file}...")
        db_fn = tempfile.NamedTemporaryFile(delete=False).name
        genes = []
        try:
            db = gffutils.create_db(self.gff_file, dbfn=db_fn, force=True, keep_order=True,
                                    merge_strategy='merge', sort_attribute_values=True)
            for gene in db.features_of_type('gene', order_by='start'):
                locus_tag = gene.attributes.get('locus_tag', [None])[0]
                if not locus_tag: continue
                cds_features = list(db.children(gene, featuretype='CDS', order_by='start'))
                if not cds_features: continue
                
                strand = 1 if gene.strand == '+' else -1
                cds_min_coord = min(c.start - 1 for c in cds_features)
                cds_max_coord = max(c.end for c in cds_features)
                
                cds_5prime_start = cds_min_coord if strand == 1 else cds_max_coord
                cds_3prime_start = cds_max_coord - 3 if strand == 1 else cds_min_coord
                cds_3prime_end = cds_max_coord if strand == 1 else cds_min_coord + 3
                
                genes.append(Gene(id=locus_tag, start=gene.start - 1, end=gene.end, strand=strand,
                                  cds_5prime_start=cds_5prime_start, cds_3prime_start=cds_3prime_start,
                                  cds_3prime_end=cds_3prime_end, cds_min_coord=cds_min_coord, cds_max_coord=cds_max_coord))
        finally:
            if os.path.exists(db_fn): os.unlink(db_fn)
        logging.info(f"找到 {len(genes)} 个包含有效CDS的基因")
        return genes

class SequenceUtils:
    @staticmethod
    def get_reverse_complement(seq: str) -> str: return str(Seq(seq).reverse_complement())

    @staticmethod
    def get_sequence(base_genome_len: int, search_genome: Seq, start: int, end: int, genome_type: str) -> str:
        length = end - start
        if length <= 0: return ""
        if genome_type == 'linear':
            start = max(0, start)
            end = min(base_genome_len, end)
            if start >= end: return ""
            return str(search_genome[start:end])
        elif genome_type == 'circle':
            # Using a doubled genome for circularity is a robust way to handle slices over the origin
            search_genome_double = search_genome + search_genome
            return str(search_genome_double[start : start + length])
        else:
            raise ValueError(f"未知的基因组类型: {genome_type}")

    @staticmethod
    def contains_restriction_sites(sequence: str, sites: Optional[List[str]]) -> bool:
        if not sites or not sequence: return False
        seq_upper = sequence.upper()
        for site in sites:
            site_upper = site.upper()
            if site_upper in seq_upper or SequenceUtils.get_reverse_complement(site_upper) in seq_upper:
                return True
        return False

    @staticmethod
    def check_variable_regions_for_sites(final_oligo: str, sgrna_seq: str, barcode: str, restriction_sites: Optional[List[str]]) -> bool:
        if not restriction_sites: return False
        max_site_len = max(len(s) for s in restriction_sites) if restriction_sites else 0
        if max_site_len == 0: return False
        flank_len = max_site_len - 1

        variable_parts = [sgrna_seq, SequenceUtils.get_reverse_complement(sgrna_seq), barcode]
        for part in variable_parts:
            if not part: continue
            pos = final_oligo.find(part)
            if pos != -1:
                start = max(0, pos - flank_len)
                end = min(len(final_oligo), pos + len(part) + flank_len)
                context_to_check = final_oligo[start:end]
                if SequenceUtils.contains_restriction_sites(context_to_check, restriction_sites):
                    logging.debug(f"在 '{part}' 的上下文中发现限制性酶切位点")
                    return True
        return False

    @staticmethod
    def generate_unique_barcode(length: int, existing_barcodes: Set[str], restriction_sites: Optional[List[str]]) -> str:
        max_attempts = 1000
        for _ in range(max_attempts):
            barcode = "".join(random.choices("ATGC", k=length))
            if barcode in existing_barcodes: continue
            gc_count = barcode.count('G') + barcode.count('C')
            min_gc, max_gc = math.ceil(length * 0.3), math.floor(length * 0.7)
            if not (min_gc <= gc_count <= max_gc): continue
            if 'A'*6 in barcode or 'T'*6 in barcode or 'G'*6 in barcode or 'C'*6 in barcode: continue
            if SequenceUtils.contains_restriction_sites(barcode, restriction_sites): continue
            return barcode
        raise RuntimeError(f"在 {max_attempts} 次尝试后未能生成唯一的barcode。")

class SGRNADesigner:
    def __init__(self, config: DesignConfig): 
        self.config = config

    def find_casts_sgrnas_in_region(self, base_genome_len: int, search_genome: Seq, search_start: int, search_end: int, genome_type: str, restriction_sites: Optional[List[str]]) -> List[dict]:
        all_sgrnas, seen = [], set()
        guide_len, pam_len = 32, 2
        
        ## CHANGE: Updated PAM sequences
        pams_fwd = ["AC", "GC", "CC"]
        pams_rev_comp = ["GT", "GC", "GG"]

        buffer = guide_len + pam_len + 1
        search_region_seq = SequenceUtils.get_sequence(base_genome_len, search_genome, search_start - buffer, search_end + buffer, genome_type)
        
        for i in range(search_start, search_end):
            local_i = i - (search_start - buffer)
            
            if local_i + pam_len + guide_len > len(search_region_seq): continue

            # 正链搜索: 5'-[PAM][Spacer]-3'
            window_fwd = search_region_seq[local_i : local_i + pam_len + guide_len]
            current_pam_fwd = window_fwd[:pam_len]
            if current_pam_fwd in pams_fwd:
                guide = window_fwd[pam_len:]
                if guide not in seen and not SequenceUtils.contains_restriction_sites(guide, restriction_sites):
                    all_sgrnas.append({'sgrna_seq': guide, 'sgrna_strand': 'forward', 'pam': current_pam_fwd, 'target_site': i})
                    seen.add(guide)

            # 负链搜索 (在正链上寻找 5'-[Spacer_rc][PAM_rc]-3')
            window_rev = search_region_seq[local_i : local_i + guide_len + pam_len]
            current_pam_rc = window_rev[guide_len:]
            if current_pam_rc in pams_rev_comp:
                guide_rc = window_rev[:guide_len]
                guide = SequenceUtils.get_reverse_complement(guide_rc)
                if guide not in seen and not SequenceUtils.contains_restriction_sites(guide, restriction_sites):
                    original_pam = SequenceUtils.get_reverse_complement(current_pam_rc)
                    # Note: target_site for reverse strand refers to the start of the [guide_rc][pam_rc] block
                    all_sgrnas.append({'sgrna_seq': guide, 'sgrna_strand': 'reverse', 'pam': original_pam, 'target_site': i})
                    seen.add(guide)
        return all_sgrnas

class CRISPRDesigner:
    def __init__(self, config: DesignConfig, genome_seq: Seq, genome_len: int, genome_type: str):
        self.config = config
        self.sgrna_designer = SGRNADesigner(config)
        self.genome_seq = genome_seq
        self.genome_len = genome_len
        self.genome_type = genome_type
        # For circular genomes, create a doubled sequence to handle searches crossing the origin
        self.search_genome = self.genome_seq + self.genome_seq if self.genome_type == 'circle' else self.genome_seq

    def _assemble_final_oligo(self, sgrna_fwd: str, barcode: str) -> str:
        if not self.config.synthesis_template:
            raise ValueError("此模式需要合成模板，但未提供。")
        
        template = self.config.synthesis_template
        replacements = {
            "{sgRNA_fwd}": sgrna_fwd,
            "{sgRNA_rc}": SequenceUtils.get_reverse_complement(sgrna_fwd),
            "{barcode}": barcode
        }
        # Clean up placeholders not used by CASTs design
        for placeholder in ["{pam}", "{pam_rc}", "{upstream_arm}", "{downstream_arm}", "{insert}", "{exempt_restriction_site}"]:
            template = template.replace(placeholder, "")

        for placeholder, value in replacements.items():
            template = template.replace(placeholder, value)
        return template
    
    def design_casts_knockout_for_gene(self, gene: Gene, used_barcodes: Set[str], restriction_sites: Optional[List[str]]) -> List[CASTsDesignResult]:
        cds_len = gene.cds_max_coord - gene.cds_min_coord
        cds_midpoint = gene.cds_min_coord + (cds_len // 2)

        search_start = gene.cds_min_coord + int(cds_len * self.config.casts_target_start_pct)
        search_end = gene.cds_min_coord + int(cds_len * self.config.casts_target_end_pct)
        if search_start >= search_end:
            logging.debug(f"基因 {gene.id}: CDS 长度过短或靶向范围无效。")
            return []
        
        all_sgrnas = self.sgrna_designer.find_casts_sgrnas_in_region(self.genome_len, self.search_genome, search_start, search_end, self.genome_type, restriction_sites)
        if not all_sgrnas: return []

        valid_sgrnas_with_scores = []
        for sgrna in all_sgrnas:
            ## CHANGE: Correct insertion site calculation
            if sgrna['sgrna_strand'] == 'forward':
                # Insertion is 50 bp downstream of the PAM start
                insertion_site = sgrna['target_site'] + 50
            else: # reverse
                # Insertion is 50 bp downstream (rev strand) of the PAM start, which is 50 bp upstream (fwd strand)
                # from the PAM's corresponding position on the fwd strand.
                # `target_site` is the start of [guide_rc][pam_rc] block. PAM on rev strand is at `target_site + 33`.
                insertion_site = (sgrna['target_site'] + 33) - 50 # Simplified to target_site - 17
                
            if gene.cds_min_coord <= insertion_site < gene.cds_max_coord:
                distance_to_midpoint = abs(insertion_site - cds_midpoint)
                sgrna_info = sgrna.copy()
                sgrna_info['insertion_site'] = insertion_site
                sgrna_info['score'] = distance_to_midpoint
                valid_sgrnas_with_scores.append(sgrna_info)
        
        if not valid_sgrnas_with_scores: return []

        valid_sgrnas_with_scores.sort(key=lambda x: x['score'])

        found_designs = []
        # Take the best N sgRNAs based on proximity to CDS center
        sgrnas_to_process = valid_sgrnas_with_scores[:self.config.sgrna_num]

        for sgrna in sgrnas_to_process:
            barcode = SequenceUtils.generate_unique_barcode(self.config.barcode_len, used_barcodes, restriction_sites)
            final_oligo = self._assemble_final_oligo(sgrna_fwd=sgrna['sgrna_seq'], barcode=barcode)
            
            if SequenceUtils.check_variable_regions_for_sites(final_oligo, sgrna['sgrna_seq'], barcode, restriction_sites):
                logging.warning(f"为基因 {gene.id} 的最优sgRNA {sgrna['sgrna_seq']} 生成的oligo存在酶切位点，跳过。")
                continue

            used_barcodes.add(barcode)
            found_designs.append(CASTsDesignResult(
                gene_id=gene.id, design_strategy='knockout', 
                sgrna_seq=sgrna['sgrna_seq'], sgrna_strand=sgrna['sgrna_strand'], 
                pam=sgrna['pam'], target_site=sgrna['target_site'], 
                insertion_site=sgrna['insertion_site'], barcode=barcode, 
                final_oligo_for_synthesis=final_oligo
            ))
            
        return found_designs

    def design_casts_knockin_for_gene(self, gene: Gene, gene_index: int, all_genes: List[Gene], used_barcodes: Set[str], restriction_sites: Optional[List[str]]) -> List[CASTsDesignResult]:
        dist_min = self.config.promoter_insertion_dist_min
        dist_max = self.config.promoter_insertion_dist_max

        # Define the desired insertion region upstream of the gene's 5' start
        if gene.strand == 1: # + strand, upstream has smaller coordinates
            insertion_region_start = gene.cds_5prime_start - dist_max
            insertion_region_end = gene.cds_5prime_start - dist_min
        else: # - strand, upstream has larger coordinates
            insertion_region_start = gene.cds_5prime_start + dist_min
            insertion_region_end = gene.cds_5prime_start + dist_max
            
        # Define a wider search area for sgRNAs around the insertion region
        search_margin = 80
        sgrna_search_start = insertion_region_start - search_margin
        sgrna_search_end = insertion_region_end + search_margin

        all_sgrnas = self.sgrna_designer.find_casts_sgrnas_in_region(
            self.genome_len, self.search_genome, sgrna_search_start, sgrna_search_end, self.genome_type, restriction_sites
        )
        if not all_sgrnas: return []

        valid_designs = []
        for sgrna in all_sgrnas:
            ## CHANGE: Correct insertion site calculation
            if sgrna['sgrna_strand'] == 'forward':
                insertion_site = sgrna['target_site'] + 50
            else: # reverse
                insertion_site = (sgrna['target_site'] + 33) - 50 # Simplified to target_site - 17
            
            # Filter 1: Check if the calculated insertion site is within the desired promoter region
            if not (insertion_region_start <= insertion_site <= insertion_region_end):
                continue

            # Filter 2: If not in relax mode, check for collision with upstream gene's CDS
            if gene_index > 0 and not self.config.relax:
                upstream_gene = all_genes[gene_index - 1]
                # Check only if genes are close (e.g., < 1kb apart)
                if abs(gene.start - upstream_gene.end) < 1000:
                    if upstream_gene.cds_min_coord <= insertion_site < upstream_gene.cds_max_coord:
                        continue

            barcode = SequenceUtils.generate_unique_barcode(self.config.barcode_len, used_barcodes, restriction_sites)
            final_oligo = self._assemble_final_oligo(sgrna_fwd=sgrna['sgrna_seq'], barcode=barcode)

            if SequenceUtils.check_variable_regions_for_sites(final_oligo, sgrna['sgrna_seq'], barcode, restriction_sites):
                continue

            used_barcodes.add(barcode)
            valid_designs.append(CASTsDesignResult(
                gene_id=gene.id, design_strategy='knockin',
                sgrna_seq=sgrna['sgrna_seq'], sgrna_strand=sgrna['sgrna_strand'],
                pam=sgrna['pam'], target_site=sgrna['target_site'],
                insertion_site=insertion_site, barcode=barcode,
                final_oligo_for_synthesis=final_oligo
            ))
            
            if len(valid_designs) >= self.config.sgrna_num:
                break
        
        return valid_designs

class ResultProcessor:
    @staticmethod
    def save_results(designs: List, failed_gene_data: List[Dict], output_file: str, design_type_info: str):
        if designs:
            results = [{
                "Gene_ID": d.gene_id, "Status": "Success", "Design_Strategy": d.design_strategy,
                "sgRNA_Sequence": d.sgrna_seq, "Barcode": d.barcode,
                "Final_Oligo_for_Synthesis": d.final_oligo_for_synthesis, "PAM": d.pam,
                "Strand": d.sgrna_strand, "Target_Site": d.target_site,
                "Predicted_Insertion_Site": d.insertion_site
            } for d in designs]
            df = pd.DataFrame(results)
            cols = ["Gene_ID", "Status", "Design_Strategy", "sgRNA_Sequence", "Barcode", "Final_Oligo_for_Synthesis", "PAM", "Strand", "Target_Site", "Predicted_Insertion_Site"]
            df = df[cols]
            df.to_csv(output_file, index=False)
            logging.info(f"{design_type_info} 结果已保存至: {output_file}")
        else:
            logging.warning(f"没有成功的 {design_type_info} 设计可供保存。")

        if failed_gene_data:
            failed_df = pd.DataFrame(failed_gene_data)
            failed_output_path = os.path.splitext(output_file)[0] + "_failed.csv"
            failed_df.to_csv(failed_output_path, index=False)
            logging.info(f"失败的基因ID及序列信息已保存至: {failed_output_path}")
            
        successful_designs_count = len(designs)
        failed_genes_count = len(failed_gene_data)
        logging.info(f"总结: 生成 {successful_designs_count} 个成功设计, {failed_genes_count} 个基因设计失败。")

# --- 各模式的执行流程函数 ---
def run_casts_knockout_pipeline(args, config, genes, genome_processor):
    designer = CRISPRDesigner(config, genome_processor.genome_seq, genome_processor.genome_len, args.genome_type)
    successful_designs, used_barcodes = [], set()
    designed_gene_ids = set()

    for i, gene in enumerate(genes):
        if (i + 1) % 100 == 0: logging.info(f"正在处理基因 {i+1}/{len(genes)}...")
        results = designer.design_casts_knockout_for_gene(gene, used_barcodes, args.restriction_site)
        if results:
            successful_designs.extend(results)
            designed_gene_ids.add(gene.id)

    # Prepare detailed failed gene data
    failed_genes = [g for g in genes if g.id not in designed_gene_ids]
    failed_gene_data = []
    for g in failed_genes:
        cds_seq = SequenceUtils.get_sequence(designer.genome_len, designer.search_genome, g.cds_min_coord, g.cds_max_coord, designer.genome_type)
        failed_gene_data.append({
            "Gene_ID": g.id, "Status": "Failed",
            "Strand": "+" if g.strand == 1 else "-",
            "Relevant_Sequence": cds_seq
        })

    ResultProcessor.save_results(successful_designs, failed_gene_data, args.output, "Knockout_CASTs")

def run_casts_knockin_pipeline(args, config, genes, genome_processor):
    designer = CRISPRDesigner(config, genome_processor.genome_seq, genome_processor.genome_len, args.genome_type)
    successful_designs, used_barcodes = [], set()
    designed_gene_ids = set()

    for i, gene in enumerate(genes):
        if (i + 1) % 100 == 0: logging.info(f"正在处理基因 {i+1}/{len(genes)}...")
        results = designer.design_casts_knockin_for_gene(gene, i, genes, used_barcodes, args.restriction_site)
        if results:
            successful_designs.extend(results)
            designed_gene_ids.add(gene.id)
    
    # Prepare detailed failed gene data
    failed_genes = [g for g in genes if g.id not in designed_gene_ids]
    failed_gene_data = []
    for g in failed_genes:
        # For knock-in, the relevant sequence is the upstream region
        if g.strand == 1:
            start, end = g.cds_5prime_start - 200, g.cds_5prime_start
        else: # - strand
            start, end = g.cds_5prime_start, g.cds_5prime_start + 200
        
        upstream_seq = SequenceUtils.get_sequence(designer.genome_len, designer.search_genome, start, end, designer.genome_type)
        failed_gene_data.append({
            "Gene_ID": g.id, "Status": "Failed",
            "Strand": "+" if g.strand == 1 else "-",
            "Relevant_Sequence": upstream_seq
        })

    ResultProcessor.save_results(successful_designs, failed_gene_data, args.output, "PromoterChange_CASTs")

def main():
    parser = argparse.ArgumentParser(description="细菌CRISPR-CASTs/INTEGRATE文库设计工具", formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest='mode', required=True, help="可用的设计模式")
    
    base_parser = argparse.ArgumentParser(add_help=False)
    base_parser.add_argument("--input_fna", required=True, help="输入的基因组FASTA文件。")
    base_parser.add_argument("--input_gff", required=True, help="输入的GFF3注释文件。")
    base_parser.add_argument("--output", required=True, help="输出的CSV文件路径。")
    base_parser.add_argument("--genome_type", type=str, default="linear", choices=['linear', 'circle'], help="指定基因组类型 (默认: linear)。")
    base_parser.add_argument("--sgRNA_num", type=int, default=1, help="为每个基因生成的设计方案数量 (默认: 1)。")
    base_parser.add_argument("--barcode_len", type=int, default=8, help="指定生成条形码(barcode)的长度 (默认: 8)。")
    base_parser.add_argument("--relax", action="store_true", help="不考虑设计是否会破坏相邻基因的CDS (仅对PromoterChange模式有效)。")

    oligo_synthesis_parser = argparse.ArgumentParser(add_help=False)
    oligo_synthesis_parser.add_argument("--restriction_site", type=str, nargs='+', help="在可变区中需要避免的限制性酶切位点。")
    oligo_synthesis_parser.add_argument("--synthesis_template", type=str, required=True, help="定义最终oligo结构的文本模板文件路径。")
    
    # --- Mode Definitions ---
    parser_ko_casts = subparsers.add_parser('Knockout_CASTs', parents=[base_parser, oligo_synthesis_parser], help="为INTEGRATE/CASTs系统设计基因敲除sgRNA文库。")
    parser_ko_casts.add_argument("--target_cds_range", type=str, default="5:80", help="在CDS中搜索sgRNA的百分比范围, 格式 '最小:最大'。")

    parser_ki_casts = subparsers.add_parser('PromoterChange_CASTs', parents=[base_parser, oligo_synthesis_parser], help="为INTEGRATE/CASTs系统设计启动子插入sgRNA文库。")
    parser_ki_casts.add_argument("--insertion_range_promoter", type=str, default="25:50", help="插入点距离起始密码子上游的距离范围, 格式 '最小:最大' bp (默认: '25:50')。")
    
    args = parser.parse_args()

    try:
        synthesis_template_content = None
        if hasattr(args, 'synthesis_template') and args.synthesis_template:
            try:
                with open(args.synthesis_template, 'r') as f:
                    synthesis_template_content = f.read().strip()
                logging.info(f"成功从以下路径加载合成模板: {args.synthesis_template}")
            except FileNotFoundError:
                logging.error(f"致命错误: 找不到合成模板文件: {args.synthesis_template}")
                sys.exit(1)
        
        config = DesignConfig(
            sgrna_num=args.sgRNA_num,
            barcode_len=args.barcode_len,
            relax=args.relax,
            synthesis_template=synthesis_template_content,
        )
        
        if args.mode == 'Knockout_CASTs':
            min_pct, max_pct = ConfigParser.parse_range_param(args.target_cds_range, "靶向范围")
            config.casts_target_start_pct, config.casts_target_end_pct = min_pct / 100.0, max_pct / 100.0
        elif args.mode == 'PromoterChange_CASTs':
            config.promoter_insertion_dist_min, config.promoter_insertion_dist_max = ConfigParser.parse_range_param(args.insertion_range_promoter, "启动子插入范围")

        genome_processor = GenomeProcessor(args.input_fna, args.input_gff)
        genome_processor.load_genome()
        genes = genome_processor.parse_genes()

        if args.mode == 'Knockout_CASTs':
            run_casts_knockout_pipeline(args, config, genes, genome_processor)
        elif args.mode == 'PromoterChange_CASTs':
            run_casts_knockin_pipeline(args, config, genes, genome_processor)

    except Exception as e:
        logging.error(f"发生严重错误: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
