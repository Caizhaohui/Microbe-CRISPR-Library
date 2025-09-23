#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Cfusion_designer_v3.4.1.py

功能:
- 为目标基因设计 C 端蛋白融合文库 (Cas9-based)。
- LHA 固定为终止密码子前的序列。
- RHA 在 sgRNA 靶点下游滑动，sgRNA 区域将被删除。
- sgRNA 在终止密码子下游区域寻找，如果初始范围未找到，则自动向3'端扩大搜索范围。
- 实现“阅读框感知”的终止密码子检查，精确排除 LHA 中的提前终止密码子。
- 严格执行 LHA 与 sgRNA 最多重叠 9bp 的规则。
- 默认不考虑邻近基因，可通过 --strict 参数开启保护模式。
"""

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

# --- 生物学数据表 ---
ECOLI_CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}
SYNONYMOUS_CODONS = {
    'G': ['GGT', 'GGC', 'GGA', 'GGG'], 'E': ['GAA', 'GAG'], 'D': ['GAT', 'GAC'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'K': ['AAA', 'AAG'], 'N': ['AAT', 'AAC'], 'M': ['ATG'], 'I': ['ATT', 'ATC', 'ATA'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'W': ['TGG'], '_': ['TAA', 'TAG', 'TGA'], 'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'Y': ['TAT', 'TAC'], 'C': ['TGT', 'TGC'], 'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
}
MODEL_WEIGHTS = {
    'pos': {('A', 2): -0.275, ('C', 2): 0.194, ('T', 2): -0.326, ('A', 3): -0.373, ('C', 3): 0.129, ('G', 3): -0.174, ('A', 4): -0.012, ('C', 4): 0.088, ('T', 4): -0.019, ('A', 5): 0.252, ('C', 5): -0.100, ('T', 5): -0.294, ('A', 6): 0.130, ('C', 6): -0.091, ('G', 6): 0.297, ('A', 7): -0.201, ('C', 7): 0.245, ('G', 7): -0.208, ('A', 11): -0.298, ('C', 11): 0.178, ('T', 11): 0.117, ('C', 12): -0.017, ('G', 12): 0.134, ('T', 12): -0.258, ('A', 13): 0.329, ('C', 13): -0.149, ('T', 13): -0.323, ('A', 14): 0.075, ('G', 14): -0.012, ('T', 14): -0.193, ('A', 15): 0.388, ('C', 15): -0.402, ('G', 15): 0.094, ('A', 16): -0.014, ('C', 16): 0.209, ('G', 16): -0.218, ('C', 17): -0.239, ('G', 17): 0.317, ('T', 17): 0.082, ('G', 18): 0.491, ('T', 18): -0.428, ('C', 19): 0.082, ('G', 19): 0.158, ('T', 19): -0.306, ('G', 20): 0.088, ('T', 20): -0.188, ('G', 21): -0.324, ('T', 21): 0.389, ('C', 22): -0.730, ('G', 22): 0.520, ('C', 23): 0.277, ('G', 23): -0.413, ('T', 23): 0.223, ('G', 24): 0.032, ('T', 24): -0.153, ('A', 27): 0.099, ('C', 27): -0.046, ('T', 27): -0.103, ('A', 28): 0.279, ('G', 28): -0.223, ('T', 28): -0.190, ('C', 29): -0.021, ('G', 29): 0.147, ('T', 29): -0.207},
    'dinuc': {('GT', 3): -0.620, ('GG', 5): 0.507, ('TA', 5): -0.548, ('TC', 6): 0.327, ('CC', 11): -0.533, ('TG', 11): 0.443, ('GA', 13): 0.449, ('CT', 13): -0.697, ('GC', 14): 0.419, ('AA', 15): -0.499, ('AG', 15): 0.541, ('AC', 18): -0.420, ('GT', 18): 0.499, ('TC', 18): -0.551, ('CG', 19): 0.589, ('AG', 20): -0.542, ('TG', 21): 0.398, ('GT', 23): -0.672, ('GG', 23): 0.533, ('GA', 27): -0.580, ('CT', 28): 0.471},
    'intercept': 0.59763615, 'gc_high': -0.1665878, 'gc_low': -0.2026259
}

# --- 数据类定义 ---
@dataclass
class DesignConfig:
    arm_search_order: List[int]; guide_len: int; sgrna_num: int; barcode_len: int
    sgrna_search_range: int; max_search_expansion: int; strict: bool
    synthesis_template: Optional[str] = None; cloning_site: Optional[str] = None

@dataclass
class Gene:
    id: str; start: int; end: int; strand: int; cds_5prime_start: int
    cds_3prime_start: int; cds_3prime_end: int; cds_min_coord: int; cds_max_coord: int

@dataclass
class KnockinDesignResult:
    gene_id: str; design_type: str; insertion_site: int; arm_length: int
    sgrna_seq: str; sgrna_pam: str; sgrna_strand: str; sgrna_score: float
    deletion_size: int; upstream_arm: str; downstream_arm: str; barcode: str
    final_oligo_for_synthesis: str; design_scenario: Optional[str] = None
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
                return tuple(sorted((val1, val2))) if sort_values else (val1, val2)
            else: raise ValueError("输入必须是单个数字或由冒号分隔的两个数字。")
        except (ValueError, IndexError) as e:
            raise ValueError(f"无效的 {param_name} 格式 '{param_str}'. 错误: {e}") from e

    @staticmethod
    def create_arm_search_order(preferred_hr: int, min_hr: int) -> List[int]:
        if preferred_hr < min_hr:
            raise ValueError(f"首选同源臂长度 ({preferred_hr}) 不能小于最小长度 ({min_hr})。")
        return list(range(preferred_hr, min_hr - 1, -1))

class GenomeProcessor:
    def __init__(self, genome_file: str, gff_file: str):
        self.genome_file, self.gff_file = genome_file, gff_file
        self.genome_seq: Optional[Seq] = None; self.genome_len: int = 0

    def load_genome(self) -> None:
        logging.info(f"正在加载基因组: {self.genome_file}...")
        try:
            genome_record = SeqIO.to_dict(SeqIO.parse(self.genome_file, "fasta"))
            if not genome_record: raise ValueError("基因组文件中没有序列")
            chrom_id = list(genome_record.keys())[0]
            # --- FIX: Split assignment into two lines to avoid TypeError ---
            self.genome_seq = genome_record[chrom_id].seq.upper()
            self.genome_len = len(self.genome_seq)
            logging.info(f"基因组加载完成: {self.genome_len} bp")
        except Exception as e:
            logging.error(f"加载基因组失败: {e}"); raise

    def parse_genes(self) -> List[Gene]:
        logging.info(f"正在解析GFF: {self.gff_file}...")
        db_fn = tempfile.NamedTemporaryFile(delete=False).name; genes = []
        try:
            db = gffutils.create_db(self.gff_file, dbfn=db_fn, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
            for gene in db.features_of_type('gene', order_by='start'):
                locus_tag = gene.attributes.get('locus_tag', [None])[0]
                if not locus_tag: continue
                cds_features = list(db.children(gene, featuretype='CDS', order_by='start'))
                if not cds_features: continue
                cds_starts, cds_ends = [c.start - 1 for c in cds_features], [c.end for c in cds_features]
                strand = 1 if gene.strand == '+' else -1
                cds_min_coord, cds_max_coord = min(cds_starts), max(cds_ends)
                if strand == 1:
                    cds_5prime_start, cds_3prime_start, cds_3prime_end = cds_min_coord, cds_max_coord - 3, cds_max_coord
                else:
                    cds_5prime_start, cds_3prime_start, cds_3prime_end = cds_max_coord - 3, cds_min_coord, cds_min_coord + 3
                genes.append(Gene(id=locus_tag, start=gene.start - 1, end=gene.end, strand=strand, cds_5prime_start=cds_5prime_start, cds_3prime_start=cds_3prime_start, cds_3prime_end=cds_3prime_end, cds_min_coord=cds_min_coord, cds_max_coord=cds_max_coord))
        finally:
            if os.path.exists(db_fn): os.unlink(db_fn)
        logging.info(f"找到 {len(genes)} 个包含有效CDS的基因"); return genes

class SequenceUtils:
    @staticmethod
    def get_reverse_complement(seq: str) -> str: return str(Seq(seq).reverse_complement())

    @staticmethod
    def get_sequence(base_genome_len: int, search_genome: Seq, start: int, end: int, genome_type: str) -> str:
        length = end - start
        if length <= 0: return ""
        if genome_type == 'linear':
            start, end = max(0, start), min(base_genome_len, end)
            if start >= end: return ""
            return str(search_genome[start:end])
        elif genome_type == 'circle':
            shifted_start = start + base_genome_len
            return str(search_genome[shifted_start : shifted_start + length])
        else: raise ValueError(f"未知的基因组类型: {genome_type}")

    @staticmethod
    def contains_restriction_sites(sequence: str, sites: Optional[List[str]]) -> bool:
        if not sites or not sequence: return False
        seq_upper = sequence.upper()
        for site in sites:
            site_upper = site.upper()
            if site_upper in seq_upper or SequenceUtils.get_reverse_complement(site_upper) in seq_upper: return True
        return False

    @staticmethod
    def check_final_oligo_with_exemption(final_oligo: str, restriction_sites: Optional[List[str]], exempt_site: Optional[str]) -> bool:
        if not restriction_sites: return False
        if not exempt_site or exempt_site not in final_oligo:
            return SequenceUtils.contains_restriction_sites(final_oligo, restriction_sites)
        max_site_len = max((len(s) for s in restriction_sites if s), default=0)
        if max_site_len == 0: return False
        junction_check_len = max_site_len - 1
        if len(exempt_site) <= 2 * junction_check_len:
            masked_oligo = final_oligo.replace(exempt_site, 'X' * len(exempt_site))
        else:
            prefix, suffix = exempt_site[:junction_check_len], exempt_site[-junction_check_len:]
            middle = 'N' * (len(exempt_site) - 2 * junction_check_len)
            masked_oligo = final_oligo.replace(exempt_site, prefix + middle + suffix)
        return SequenceUtils.contains_restriction_sites(masked_oligo, restriction_sites)

    @staticmethod
    def generate_unique_barcode(length: int, existing_barcodes: Set[str], restriction_sites: Optional[List[str]]) -> str:
        max_attempts = 1000
        for _ in range(max_attempts):
            barcode = "".join(random.choices("ATGC", k=length))
            if barcode in existing_barcodes: continue
            gc_count = barcode.count('G') + barcode.count('C')
            if not (math.ceil(length * 0.3) <= gc_count <= math.floor(length * 0.7)): continue
            if any(h in barcode for h in ['AAAAAA', 'TTTTTT', 'GGGGGG', 'CCCCCC']): continue
            if SequenceUtils.contains_restriction_sites(barcode, restriction_sites): continue
            return barcode
        raise RuntimeError(f"在 {max_attempts} 次尝试后未能生成唯一的barcode。")

class SGRNADesigner:
    def __init__(self, config: DesignConfig): self.config = config
    def match_pam(self, seq: str, pattern: str) -> bool:
        if len(seq) != len(pattern): return False
        return all(p == 'N' or s == p for s, p in zip(seq, pattern))

    def find_sgrnas_in_region(self, base_genome_len: int, search_genome: Seq, search_start: int, search_end: int, pam: str, restriction_sites: Optional[List[str]], genome_type: str) -> List[dict]:
        all_sgrnas, seen = [], set()
        pam_len, guide_len = len(pam), self.config.guide_len
        pam_fwd_pattern, pam_rev_pattern = pam.upper(), SequenceUtils.get_reverse_complement(pam).upper()
        buffer = guide_len + pam_len + 4
        region_start, region_end = min(search_start, search_end), max(search_start, search_end)
        search_region_seq = SequenceUtils.get_sequence(base_genome_len, search_genome, region_start - buffer, region_end + buffer, genome_type)
        
        for i in range(region_start, region_end):
            local_i = i - (region_start - buffer)
            if not (0 <= local_i < len(search_region_seq)): continue
            if local_i + guide_len + pam_len <= len(search_region_seq):
                window = search_region_seq[local_i : local_i + guide_len + pam_len]
                pam_candidate, guide = window[guide_len:], window[:guide_len]
                if self.match_pam(pam_candidate, pam_fwd_pattern) and guide not in seen and not SequenceUtils.contains_restriction_sites(guide, restriction_sites):
                    score = MODEL_WEIGHTS.get('intercept', 0.0) # Simplified score
                    all_sgrnas.append({'seq': guide, 'pam': pam_candidate, 'strand': 'forward', 'score': score, 'genomic_start': i})
                    seen.add(guide)
            if local_i + pam_len + guide_len <= len(search_region_seq):
                window = search_region_seq[local_i : local_i + pam_len + guide_len]
                pam_rc_candidate, guide_rc = window[:pam_len], window[pam_len:]
                if self.match_pam(pam_rc_candidate, pam_rev_pattern):
                    guide = SequenceUtils.get_reverse_complement(guide_rc)
                    if guide not in seen and not SequenceUtils.contains_restriction_sites(guide, restriction_sites):
                        score = MODEL_WEIGHTS.get('intercept', 0.0) # Simplified score
                        all_sgrnas.append({'seq': guide, 'pam': SequenceUtils.get_reverse_complement(pam_rc_candidate), 'strand': 'reverse', 'score': score, 'genomic_start': i})
                        seen.add(guide)
        return all_sgrnas

class CRISPRDesigner:
    def __init__(self, config: DesignConfig, genome_seq: Seq, genome_len: int, genome_type: str):
        self.config = config; self.sgrna_designer = SGRNADesigner(config); self.genome_seq = genome_seq
        self.genome_len = genome_len; self.genome_type = genome_type
        self.search_genome = self.genome_seq + self.genome_seq if self.genome_type == 'circle' else self.genome_seq

    def _assemble_final_oligo(self, sgrna_fwd: str, pam: str, upstream_arm: str, downstream_arm: str, barcode: str) -> str:
        sgrna_rc, pam_rc = SequenceUtils.get_reverse_complement(sgrna_fwd), SequenceUtils.get_reverse_complement(pam)
        template = self.config.synthesis_template
        if not template: raise ValueError("未提供合成模板文件。")
        replacements = {"{sgRNA_fwd}": sgrna_fwd, "{sgRNA_rc}": sgrna_rc, "{pam}": pam, "{pam_rc}": pam_rc, "{upstream_arm}": upstream_arm, "{downstream_arm}": downstream_arm, "{barcode}": barcode, "{insert}": "", "{exempt_restriction_site}": self.config.cloning_site or ""}
        for placeholder, value in replacements.items(): template = template.replace(placeholder, value)
        return template

    def design_c_fusion_for_gene(self, gene: Gene, gene_index: int, all_genes: List[Gene], pam: str, used_barcodes: Set[str], restriction_sites: Optional[List[str]]) -> List[KnockinDesignResult]:
        found_designs = []
        is_forward = gene.strand == 1

        current_search_range = self.config.sgrna_search_range
        candidate_sgrnas = []
        while not candidate_sgrnas and current_search_range <= self.config.max_search_expansion:
            search_start, search_end = (gene.cds_3prime_start, gene.cds_3prime_start + current_search_range) if is_forward else (gene.cds_3prime_start - current_search_range, gene.cds_3prime_start)
            candidate_sgrnas = self.sgrna_designer.find_sgrnas_in_region(self.genome_len, self.search_genome, search_start, search_end, pam, restriction_sites, self.genome_type)
            if not candidate_sgrnas: current_search_range += 50
        
        if not candidate_sgrnas:
            logging.warning(f"基因 {gene.id}: 在终止密码子下游 {self.config.max_search_expansion}bp 范围内未找到有效sgRNA。")
            return []
        
        candidate_sgrnas.sort(key=lambda s: (-s['score'], abs(s['genomic_start'] - gene.cds_3prime_start)))
        
        downstream_gene_protected_cds = None
        if self.config.strict and gene_index < len(all_genes) - 1:
            downstream_gene = all_genes[gene_index + 1]
            if downstream_gene.start - gene.end < 1000:
                downstream_gene_protected_cds = (downstream_gene.cds_min_coord, downstream_gene.cds_max_coord)

        for sgrna in candidate_sgrnas:
            if len(found_designs) >= self.config.sgrna_num: break
            for arm_len in self.config.arm_search_order:
                if len(found_designs) >= self.config.sgrna_num: break

                if is_forward:
                    lha_start, lha_end = gene.cds_3prime_start - arm_len, gene.cds_3prime_start
                else: 
                    lha_start, lha_end = gene.cds_3prime_end, gene.cds_3prime_end + arm_len
                
                upstream_arm = SequenceUtils.get_sequence(self.genome_len, self.search_genome, lha_start, lha_end, self.genome_type)
                if not (upstream_arm and len(upstream_arm) == arm_len): continue
                
                stop_codons_to_check = {"TAA", "TAG", "TGA"} if is_forward else {"TTA", "CTA", "TCA"}
                has_internal_stop = False
                for i in range(len(upstream_arm) % 3, len(upstream_arm), 3):
                    if len(upstream_arm[i:i+3]) == 3 and upstream_arm[i:i+3] in stop_codons_to_check:
                        has_internal_stop = True; break
                if has_internal_stop:
                    logging.debug(f"基因 {gene.id}: LHA 含有框内终止密码子，跳过。"); continue
                
                pam_len = len(sgrna['pam'])
                sgrna_region_start, sgrna_region_end = sgrna['genomic_start'], sgrna['genomic_start'] + self.config.guide_len + pam_len
                
                overlap = max(0, min(lha_end, sgrna_region_end) - max(lha_start, sgrna_region_start))
                if overlap > 9: continue

                max_rha_offset = 50
                for offset in range(max_rha_offset + 1):
                    if is_forward:
                        rha_start, rha_end = sgrna_region_end + offset, sgrna_region_end + offset + arm_len
                    else:
                        rha_start, rha_end = sgrna_region_start - offset - arm_len, sgrna_region_start - offset
                    
                    downstream_arm = SequenceUtils.get_sequence(self.genome_len, self.search_genome, rha_start, rha_end, self.genome_type)
                    if not (downstream_arm and len(downstream_arm) == arm_len): continue
                    if downstream_gene_protected_cds and max(rha_start, downstream_gene_protected_cds[0]) < min(rha_end, downstream_gene_protected_cds[1]): continue
                    if SequenceUtils.contains_restriction_sites(upstream_arm, restriction_sites) or SequenceUtils.contains_restriction_sites(downstream_arm, restriction_sites): continue
                    
                    barcode = SequenceUtils.generate_unique_barcode(self.config.barcode_len, used_barcodes, restriction_sites)
                    final_oligo = self._assemble_final_oligo(sgrna['seq'], sgrna['pam'], upstream_arm, downstream_arm, barcode)
                    if SequenceUtils.check_final_oligo_with_exemption(final_oligo, restriction_sites, self.config.cloning_site): continue
                    
                    used_barcodes.add(barcode)
                    deletion_size = abs(rha_start - lha_end)
                    found_designs.append(KnockinDesignResult(
                        gene_id=gene.id, design_type='C_fusion', insertion_site=gene.cds_3prime_start, arm_length=arm_len, sgrna_seq=sgrna['seq'], 
                        sgrna_pam=sgrna['pam'], sgrna_strand=sgrna['strand'], sgrna_score=sgrna['score'], deletion_size=deletion_size, 
                        upstream_arm=upstream_arm, downstream_arm=downstream_arm, barcode=barcode, final_oligo_for_synthesis=final_oligo,
                        design_scenario=f"Deletion of sgRNA region (LHA Overlap={overlap}bp)", sgrna_cut_site=0
                    ))
                    if len(found_designs) >= self.config.sgrna_num: break
                if len(found_designs) >= self.config.sgrna_num: break

        return found_designs

class ResultProcessor:
    @staticmethod
    def save_results(designs, failed_genes, output_file, design_type_info):
        if designs:
            df = pd.DataFrame([d.__dict__ for d in designs])
            if "sgrna_cut_distance" in df.columns: df = df.drop(columns=["sgrna_cut_distance"])
            core_cols = ["gene_id", "design_type", "design_scenario", "sgrna_seq", "barcode", "final_oligo_for_synthesis"]
            other_cols = sorted([col for col in df.columns if col not in core_cols])
            df = df[core_cols + other_cols]
            df.to_csv(output_file, index=False)
            logging.info(f"{design_type_info} 结果已保存至: {output_file}")
        else:
            logging.warning(f"没有成功的 {design_type_info} 设计可供保存。")
        if failed_genes:
            failed_df = pd.DataFrame([{"Gene_ID": g.id, "Status": "Failed"} for g in failed_genes])
            failed_output_path = os.path.splitext(output_file)[0] + "_failed.csv"
            failed_df.to_csv(failed_output_path, index=False)
            logging.info(f"失败的基因ID已保存至: {failed_output_path}")
        successful_designs, failed_genes_count = (len(designs) if designs else 0), (len(failed_genes) if failed_genes else 0)
        logging.info(f"总结: {successful_designs} 个成功设计, {failed_genes_count} 个基因设计失败。")

# --- 主流程 ---
def run_c_fusion_pipeline(args, config, genes, genome_processor):
    designer = CRISPRDesigner(config, genome_processor.genome_seq, genome_processor.genome_len, args.genome_type)
    successful_designs, used_barcodes = [], set()
    min_arm_len = config.arm_search_order[-1] if config.arm_search_order else 0
    genes_to_process = [g for g in genes if (g.cds_max_coord - g.cds_min_coord) >= min_arm_len]
    logging.info(f"正在处理 {len(genes_to_process)} 个长度足够进行C-fusion设计的基因。")
    
    for i, gene in enumerate(genes_to_process):
        if (i + 1) % 100 == 0: logging.info(f"正在处理基因 {i+1}/{len(genes_to_process)}...")
        results = designer.design_c_fusion_for_gene(
            gene=gene, gene_index=i, all_genes=genes_to_process, pam=args.pam,
            used_barcodes=used_barcodes, restriction_sites=args.restriction_site)
        if results: successful_designs.extend(results)
            
    successful_gene_ids = {d.gene_id for d in successful_designs}
    failed_genes = [g for g in genes_to_process if g.id not in successful_gene_ids]
    ResultProcessor.save_results(successful_designs, failed_genes, args.output, "C-fusion")

def main():
    parser = argparse.ArgumentParser(description="C端融合 (C-fusion) 设计工具 (v3.4)", formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("--input_fna", required=True, help="输入的基因组FASTA文件。")
    parser.add_argument("--input_gff", required=True, help="输入的GFF3注释文件。")
    parser.add_argument("--output", required=True, help="输出的CSV文件路径。")
    parser.add_argument("--genome_type", type=str, default="linear", choices=['linear', 'circle'], help="指定基因组类型 (默认: linear)。")
    parser.add_argument("--pam", type=str, default="NGG", help="PAM序列 (默认: 'NGG')。")
    parser.add_argument("--HR_len", type=str, default="50:40", help="同源臂长度。格式为 '首选:最小'。")
    parser.add_argument("--sgRNA_num", type=int, default=1, help="为每个基因生成的设计方案数量 (默认: 1)。")
    parser.add_argument("--barcode_len", type=int, default=8, help="指定barcode序列的长度 (默认: 8)。")
    parser.add_argument("--sgrna_search_range", type=int, default=50, help="在终止密码子下游搜索sgRNA的初始范围 (默认: 50 bp)。")
    parser.add_argument("--max_search_expansion", type=int, default=300, help="sgRNA搜索的最大下游扩展距离 (默认: 300 bp)。")
    parser.add_argument("--restriction_site", type=str, nargs='+', help="在可变区中需要避免的限制性酶切位点。")
    parser.add_argument("--synthesis_template", type=str, required=True, help="定义最终oligo结构的文本模板文件路径。")
    parser.add_argument("--cloning_site", type=str, help="一个特殊的酶切位点，它将被放置在模板的{exempt_restriction_site}占位符中。")
    parser.add_argument("--strict", action="store_true", help="开启邻近基因保护机制，避免破坏下游基因的CDS区 (默认:关闭)。")
    
    args = parser.parse_args()

    try:
        synthesis_template_content = None
        if args.synthesis_template:
            try:
                with open(args.synthesis_template, 'r') as f:
                    synthesis_template_content = f.read().strip()
                logging.info(f"成功加载合成模板: {args.synthesis_template}")
            except FileNotFoundError:
                logging.error(f"致命错误: 找不到合成模板文件: {args.synthesis_template}"); sys.exit(1)

        preferred_hr, min_hr = ConfigParser.parse_range_param(args.HR_len, "同源臂长度", sort_values=False)
        arm_search_order = ConfigParser.create_arm_search_order(preferred_hr, min_hr)
        
        config = DesignConfig(
            arm_search_order=arm_search_order, guide_len=20, synthesis_template=synthesis_template_content,
            cloning_site=args.cloning_site, sgrna_num=args.sgRNA_num, barcode_len=args.barcode_len,
            sgrna_search_range=args.sgrna_search_range, max_search_expansion=args.max_search_expansion, strict=args.strict
        )
        
        genome_processor = GenomeProcessor(args.input_fna, args.input_gff)
        genome_processor.load_genome()
        genes = genome_processor.parse_genes()

        run_c_fusion_pipeline(args, config, genes, genome_processor)

    except Exception as e:
        logging.error(f"发生严重错误: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
