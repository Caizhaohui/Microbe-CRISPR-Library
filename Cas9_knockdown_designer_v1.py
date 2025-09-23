#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import math
import random
import tempfile
import logging
import argparse
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Set

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import gffutils

# --- 日志设置 ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- 数据结构 ---
@dataclass
class DesignConfig:
    """Knockdown模式的配置类"""
    arm_search_order: List[int]
    guide_len: int
    sgrna_num: int
    barcode_len: int
    sgrna_upstream_range_min: int
    sgrna_upstream_range_max: int
    rha_min_dist_to_atg: int
    synthesis_template: Optional[str] = None
    cloning_site: Optional[str] = None

@dataclass
class Gene:
    """通用基因信息容器"""
    id: str
    start: int
    end: int
    strand: int
    cds_5prime_start: int
    cds_min_coord: int
    cds_max_coord: int

@dataclass
class KnockinDesignResult:
    """基因敲入设计结果容器"""
    gene_id: str
    design_type: str
    insertion_site: int
    arm_length: int
    sgrna_seq: str
    sgrna_pam: str
    sgrna_strand: str
    upstream_arm: str
    downstream_arm: str
    barcode: str
    final_oligo_for_synthesis: str
    design_scenario: Optional[str] = None
    sgrna_cut_site: Optional[int] = None

# --- 辅助类 ---
class ConfigParser:
    @staticmethod
    def parse_range_param(param_str: str, param_name: str, sort_values: bool = True) -> Tuple[int, int]:
        try:
            parts = param_str.split(':')
            if len(parts) == 1: val = int(parts[0]); return val, val
            elif len(parts) == 2:
                val1, val2 = int(parts[0]), int(parts[1])
                if not sort_values and val1 < val2:
                    raise ValueError(f"在 {param_name} 中，首选长度 ({val1}) 必须大于等于最小长度 ({val2}).")
                return sorted((val1, val2)) if sort_values else (val1, val2)
            raise ValueError("输入必须是单个数字或由冒号分隔的两个数字。")
        except (ValueError, IndexError) as e:
            raise ValueError(f"无效的 {param_name} 格式 '{param_str}'. 错误: {e}") from e

    @staticmethod
    def create_arm_search_order(preferred_hr: int, min_hr: int) -> List[int]:
        return list(range(preferred_hr, min_hr - 1, -1))

class SequenceUtils:
    @staticmethod
    def get_reverse_complement(seq: str) -> str: return str(Seq(seq).reverse_complement())

    @staticmethod
    def get_sequence(search_genome: Seq, start: int, end: int, coord_offset: int = 0) -> str:
        length = end - start
        if length <= 0: return ""
        safe_start = max(0, start + coord_offset)
        safe_end = max(0, end + coord_offset)
        return str(search_genome[safe_start:safe_end])

    @staticmethod
    def contains_restriction_sites(sequence: str, sites: Optional[List[str]]) -> bool:
        if not sites or not sequence: return False
        seq_upper = sequence.upper()
        for site in sites:
            s_upper = site.upper()
            if s_upper in seq_upper or SequenceUtils.get_reverse_complement(s_upper) in seq_upper: return True
        return False

    @staticmethod
    def check_final_oligo_with_exemption(final_oligo: str, restriction_sites: Optional[List[str]], exempt_site: Optional[str]) -> bool:
        if not restriction_sites: return False
        if not exempt_site: return SequenceUtils.contains_restriction_sites(final_oligo, restriction_sites)
        placeholder = 'X' * len(exempt_site)
        masked_oligo = final_oligo.replace(exempt_site, placeholder)
        return SequenceUtils.contains_restriction_sites(masked_oligo, restriction_sites)

    @staticmethod
    def generate_unique_barcode(length: int, existing_barcodes: Set[str], restriction_sites: Optional[List[str]]) -> str:
        max_attempts = 1000
        for _ in range(max_attempts):
            barcode = "".join(random.choices("ATGC", k=length))
            if barcode in existing_barcodes: continue
            gc_count = barcode.count('G') + barcode.count('C')
            if not (math.ceil(length * 0.3) <= gc_count <= math.floor(length * 0.7)): continue
            if any(h in barcode for h in ('A'*6, 'T'*6, 'G'*6, 'C'*6)): continue
            if SequenceUtils.contains_restriction_sites(barcode, restriction_sites): continue
            return barcode
        raise RuntimeError(f"在 {max_attempts} 次尝试后未能生成唯一的barcode。")

class GenomeProcessor:
    def __init__(self, genome_file: str, gff_file: str):
        self.genome_file, self.gff_file = genome_file, gff_file
        self.genome_seq: Optional[Seq] = None
        self.genome_len = 0

    def load_genome(self) -> None:
        logging.info(f"正在加载基因组: {self.genome_file}...")
        try:
            record = next(SeqIO.parse(self.genome_file, "fasta"))
            self.genome_seq, self.genome_len = record.seq.upper(), len(record.seq)
            logging.info(f"基因组加载完成: {self.genome_len} bp")
        except StopIteration: raise ValueError("基因组FASTA文件为空或格式错误。")

    def parse_genes(self) -> List[Gene]:
        logging.info(f"正在解析GFF: {self.gff_file}...")
        db_fn = tempfile.NamedTemporaryFile(delete=False).name
        genes = []
        try:
            db = gffutils.create_db(self.gff_file, dbfn=db_fn, force=True, keep_order=True, merge_strategy='merge')
            for gene_feature in db.features_of_type('gene', order_by='start'):
                locus_tag = gene_feature.attributes.get('locus_tag', [None])[0]
                if not locus_tag: continue
                cds_features = list(db.children(gene_feature, featuretype='CDS', order_by='start'))
                if not cds_features: continue
                cds_starts = [c.start - 1 for c in cds_features]
                cds_ends = [c.end for c in cds_features]
                strand = 1 if gene_feature.strand == '+' else -1
                cds_min_coord, cds_max_coord = min(cds_starts), max(cds_ends)
                cds_5prime_start = cds_min_coord if strand == 1 else cds_max_coord - 3
                genes.append(Gene(id=locus_tag, start=gene_feature.start - 1, end=gene_feature.end, strand=strand,
                                  cds_5prime_start=cds_5prime_start, cds_min_coord=cds_min_coord, cds_max_coord=cds_max_coord))
        finally:
            if os.path.exists(db_fn): os.unlink(db_fn)
        logging.info(f"找到 {len(genes)} 个包含有效CDS的基因")
        return genes

class SGRNADesigner:
    def __init__(self, config: DesignConfig):
        self.config = config

    def match_pam(self, seq: str, pattern: str) -> bool:
        if len(seq) != len(pattern): return False
        return all(p == 'N' or s == p for s, p in zip(seq, pattern))

    def _get_cut_site(self, sgrna_local_start: int, sgrna_strand: str, pam_len: int) -> int:
        cut_offset = -3
        if sgrna_strand == "forward":
            return sgrna_local_start + self.config.guide_len + cut_offset
        else:
            return sgrna_local_start + pam_len - cut_offset

    def find_sgrnas_in_sequence(self, sequence: str, pam: str, restriction_sites: Optional[List[str]]) -> List[dict]:
        all_sgrnas, seen = [], set()
        pam_len, guide_len = len(pam), self.config.guide_len
        pam_fwd_pattern = pam.upper()
        pam_rev_pattern = SequenceUtils.get_reverse_complement(pam_fwd_pattern)

        for i in range(len(sequence) - guide_len - pam_len + 1):
            win_fwd = sequence[i : i + guide_len + pam_len]
            guide = win_fwd[:guide_len]
            pam_fwd = win_fwd[guide_len:]
            if self.match_pam(pam_fwd, pam_fwd_pattern) and guide not in seen and not SequenceUtils.contains_restriction_sites(guide, restriction_sites):
                cut_site = self._get_cut_site(i, "forward", pam_len)
                all_sgrnas.append({'seq': guide, 'pam': pam_fwd, 'strand': 'forward', 'local_start': i, 'cut_site': cut_site})
                seen.add(guide)

            win_rev = sequence[i : i + pam_len + guide_len]
            pam_rc = win_rev[:pam_len]
            guide_rc = win_rev[pam_len:]
            if self.match_pam(pam_rc, pam_rev_pattern):
                guide = SequenceUtils.get_reverse_complement(guide_rc)
                if guide not in seen and not SequenceUtils.contains_restriction_sites(guide, restriction_sites):
                    cut_site = self._get_cut_site(i, "reverse", pam_len)
                    all_sgrnas.append({'seq': guide, 'pam': SequenceUtils.get_reverse_complement(pam_rc), 'strand': 'reverse', 'local_start': i, 'cut_site': cut_site})
                    seen.add(guide)
        
        all_sgrnas.sort(key=lambda x: x['local_start'])
        return all_sgrnas

class CRISPRDesigner:
    def __init__(self, config: DesignConfig, genome_seq: Seq, genome_len: int, genome_type: str):
        self.config = config
        self.sgrna_designer = SGRNADesigner(config)
        self.genome_len = genome_len
        self.search_genome = genome_seq + genome_seq if genome_type == 'circle' else genome_seq
        self.coord_offset = self.genome_len if genome_type == 'circle' else 0

    def _assemble_final_oligo(self, sgrna_fwd: str, pam: str, upstream_arm: str, downstream_arm: str, barcode: str) -> str:
        sgrna_rc, pam_rc = SequenceUtils.get_reverse_complement(sgrna_fwd), SequenceUtils.get_reverse_complement(pam)
        template = self.config.synthesis_template
        if not template: raise ValueError("合成模板未提供。")
        replacements = {"{sgRNA_fwd}": sgrna_fwd, "{sgRNA_rc}": sgrna_rc, "{pam}": pam, "{pam_rc}": pam_rc,
                        "{upstream_arm}": upstream_arm, "{downstream_arm}": downstream_arm, "{barcode}": barcode,
                        "{insert}": "", "{exempt_restriction_site}": self.config.cloning_site or ""}
        for placeholder, value in replacements.items():
            template = template.replace(placeholder, value)
        return template

    def design_knockdown_for_gene(self, gene: Gene, used_barcodes: Set[str], restriction_sites: Optional[List[str]]) -> List[KnockinDesignResult]:
        found_designs = []
        is_forward = gene.strand == 1
        atg_pos = gene.cds_5prime_start

        search_start, search_end = 0, 0
        if is_forward:
            search_start = atg_pos - self.config.sgrna_upstream_range_max
            search_end = atg_pos - self.config.sgrna_upstream_range_min
        else:
            search_start = atg_pos + 3 + self.config.sgrna_upstream_range_min
            search_end = atg_pos + 3 + self.config.sgrna_upstream_range_max

        search_seq = SequenceUtils.get_sequence(self.search_genome, search_start, search_end, self.coord_offset)
        if not search_seq: return []
        
        candidate_sgrnas = self.sgrna_designer.find_sgrnas_in_sequence(search_seq, "NGG", restriction_sites)
        if not candidate_sgrnas: return []

        for sgrna in candidate_sgrnas:
            if len(found_designs) >= self.config.sgrna_num: break
            
            sgrna_start_logical = search_start + sgrna['local_start']
            sgrna_cut_site_logical = search_start + sgrna['cut_site']
            
            if gene.cds_min_coord <= sgrna_cut_site_logical < gene.cds_max_coord:
                continue

            # --- 二级策略循环 ---
            design_found = False
            for strategy in ["split_guide", "replace_target"]:
                if design_found: break
                
                up_arm_end, down_arm_start, insertion_site = 0, 0, 0
                scenario_text = ""

                if strategy == "split_guide":
                    if sgrna['strand'] == 'forward':
                        insertion_site = sgrna_start_logical + 14
                    else: # reverse
                        insertion_site = sgrna_start_logical + 9
                    up_arm_end, down_arm_start = insertion_site, insertion_site
                    scenario_text = f"Split-guide insertion"
                
                else: # strategy == "replace_target"
                    pam_len = len(sgrna['pam'])
                    up_arm_end = sgrna_start_logical
                    down_arm_start = sgrna_start_logical + self.config.guide_len + pam_len
                    insertion_site = up_arm_end
                    scenario_text = "Target replacement (fallback)"
                
                if gene.cds_min_coord <= insertion_site < gene.cds_max_coord: continue

                for arm_len in self.config.arm_search_order:
                    if (is_forward and (atg_pos - down_arm_start) < self.config.rha_min_dist_to_atg) or \
                       (not is_forward and (down_arm_start - (atg_pos + 3)) < self.config.rha_min_dist_to_atg):
                        continue
                    
                    up_arm = SequenceUtils.get_sequence(self.search_genome, up_arm_end - arm_len, up_arm_end, self.coord_offset)
                    down_arm = SequenceUtils.get_sequence(self.search_genome, down_arm_start, down_arm_start + arm_len, self.coord_offset)
                    
                    if len(up_arm) != arm_len or len(down_arm) != arm_len: continue
                    if SequenceUtils.contains_restriction_sites(up_arm, restriction_sites) or SequenceUtils.contains_restriction_sites(down_arm, restriction_sites): continue
                    
                    try:
                        barcode = SequenceUtils.generate_unique_barcode(self.config.barcode_len, used_barcodes, restriction_sites)
                        final_oligo = self._assemble_final_oligo(sgrna['seq'], sgrna['pam'], up_arm, down_arm, barcode)
                        if SequenceUtils.check_final_oligo_with_exemption(final_oligo, restriction_sites, self.config.cloning_site): continue
                        
                        used_barcodes.add(barcode)
                        dist = abs(insertion_site - atg_pos) if is_forward else abs(insertion_site - (atg_pos + 3))
                        found_designs.append(KnockinDesignResult(
                            gene_id=gene.id, design_type='knockdown', insertion_site=insertion_site,
                            arm_length=arm_len, sgrna_seq=sgrna['seq'], sgrna_pam=sgrna['pam'],
                            sgrna_strand=sgrna['strand'], upstream_arm=up_arm, downstream_arm=down_arm,
                            barcode=barcode, final_oligo_for_synthesis=final_oligo,
                            design_scenario=f"{scenario_text} at ATG upstream {dist}bp", sgrna_cut_site=sgrna_cut_site_logical
                        ))
                        design_found = True
                        break # 成功，跳出 arm_len 循环
                    except RuntimeError as e: logging.warning(f"为基因 {gene.id} 生成barcode时出错: {e}")
        return found_designs

# --- 结果处理与主流程 ---
class ResultProcessor:
    @staticmethod
    def save_results(designs, failed_genes, output_file):
        if designs:
            df = pd.DataFrame([d.__dict__ for d in designs])
            if 'sgrna_score' in df.columns: df = df.drop(columns=['sgrna_score'])
            df.rename(columns=lambda c: c.replace('_', ' ').title(), inplace=True)
            df.to_csv(output_file, index=False)
            logging.info(f"Knockdown 结果已保存至: {output_file}")
        else:
            logging.warning("没有成功的 Knockdown 设计可供保存。")
        if failed_genes:
            failed_data = [{"Gene_ID": g.id, "Status": "Failed", "Strand": "+" if g.strand == 1 else "-"} for g in failed_genes]
            failed_df = pd.DataFrame(failed_data)
            failed_output_path = os.path.splitext(output_file)[0] + "_failed.csv"
            failed_df.to_csv(failed_output_path, index=False)
            logging.info(f"失败的基因ID已保存至: {failed_output_path}")
        logging.info(f"总结: 生成 {len(designs)} 个成功设计, {len(failed_genes)} 个基因设计失败。")

def run_knockdown_pipeline(args, config, genes, genome_processor):
    designer = CRISPRDesigner(config, genome_processor.genome_seq, genome_processor.genome_len, args.genome_type)
    successful_designs, used_barcodes = [], set()
    logging.info(f"开始处理 {len(genes)} 个基因...")
    for i, gene in enumerate(genes):
        if (i + 1) % 100 == 0: logging.info(f"正在处理基因 {i+1}/{len(genes)}...")
        results = designer.design_knockdown_for_gene(gene=gene, used_barcodes=used_barcodes, restriction_sites=args.restriction_site)
        if results: successful_designs.extend(results)
    successful_gene_ids = {d.gene_id for d in successful_designs}
    failed_genes = [g for g in genes if g.id not in successful_gene_ids]
    ResultProcessor.save_results(successful_designs, failed_genes, args.output)

def main():
    parser = argparse.ArgumentParser(description="基因下调 (Knockdown) 文库设计工具 (Cas9)", formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("--input_fna", required=True, help="输入的基因组FASTA文件。")
    parser.add_argument("--input_gff", required=True, help="输入的GFF3注释文件。")
    parser.add_argument("--output", required=True, help="输出的CSV文件路径。")
    parser.add_argument("--genome_type", type=str, default="linear", choices=['linear', 'circle'], help="基因组类型。")
    parser.add_argument("--pam", type=str, default="NGG", help="PAM序列。")
    parser.add_argument("--HR_len", type=str, default="50:40", help="同源臂长度范围 '首选:最小'。")
    parser.add_argument("--sgRNA_num", type=int, default=1, help="每个基因的设计数量。")
    parser.add_argument("--barcode_len", type=int, default=10, help="条形码(barcode)的长度。")
    parser.add_argument("--sgrna_upstream_range", type=str, default="20:60", help="sgRNA距离起始密码子上游的搜索范围 '最小:最大' bp。")
    parser.add_argument("--rha_min_dist_to_atg", type=int, default=15, help="插入点与ATG的最小距离。")
    parser.add_argument("--restriction_site", type=str, nargs='+', help="需要避免的限制性酶切位点。")
    parser.add_argument("--synthesis_template", type=str, required=True, help="定义最终oligo结构的模板文件路径。")
    parser.add_argument("--cloning_site", type=str, help="将被豁免检查的克隆位点序列。")
    
    args = parser.parse_args()
    try:
        synthesis_template_content = ""
        if args.synthesis_template:
            with open(args.synthesis_template, 'r') as f: synthesis_template_content = f.read().strip()
            logging.info(f"成功加载合成模板: {args.synthesis_template}")

        preferred_hr, min_hr = ConfigParser.parse_range_param(args.HR_len, "同源臂长度", sort_values=False)
        sgrna_up_min, sgrna_up_max = ConfigParser.parse_range_param(args.sgrna_upstream_range, "sgRNA上游范围")

        config = DesignConfig(
            arm_search_order=ConfigParser.create_arm_search_order(preferred_hr, min_hr),
            guide_len=20,
            sgrna_num=args.sgRNA_num,
            barcode_len=args.barcode_len,
            sgrna_upstream_range_min=sgrna_up_min,
            sgrna_upstream_range_max=sgrna_up_max,
            rha_min_dist_to_atg=args.rha_min_dist_to_atg,
            synthesis_template=synthesis_template_content, 
            cloning_site=args.cloning_site
        )
        
        genome_processor = GenomeProcessor(args.input_fna, args.input_gff)
        genome_processor.load_genome()
        genes = genome_processor.parse_genes()
        
        run_knockdown_pipeline(args, config, genes, genome_processor)

    except Exception as e:
        logging.error(f"发生严重错误: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
