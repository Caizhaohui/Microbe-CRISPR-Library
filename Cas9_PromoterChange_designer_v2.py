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
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}
SYNONYMOUS_CODONS = {
    'G': ['GGT', 'GGC', 'GGA', 'GGG'], 'E': ['GAA', 'GAG'], 'D': ['GAT', 'GAC'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'K': ['AAA', 'AAG'], 'N': ['AAT', 'AAC'], 'M': ['ATG'], 'I': ['ATT', 'ATC', 'ATA'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'W': ['TGG'], '_': ['TAA', 'TAG', 'TGA'], 'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'Y': ['TAT', 'TAC'], 'C': ['TGT', 'TGC'], 'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
}
MODEL_WEIGHTS = {
    'pos': {('A', 2): -0.275, ('C', 2): 0.194, ('T', 2): -0.326, ('A', 3): -0.373, ('C', 3): 0.129, ('G', 3): -0.174, ('A', 4): -0.012, ('C', 4): 0.088, ('T', 4): -0.019, ('A', 5): 0.252, ('C', 5): -0.100, ('T', 5): -0.294, ('A', 6): 0.130, ('C', 6): -0.091, ('G', 6): 0.297, ('A', 7): -0.201, ('C', 7): 0.245, ('G', 7): -0.208, ('A', 11): -0.298, ('C', 11): 0.178, ('T', 11): 0.117, ('C', 12): -0.017, ('G', 12): 0.134, ('T', 12): -0.258, ('A', 13): 0.329, ('C', 13): -0.149, ('T', 13): -0.323, ('A', 14): 0.075, ('G', 14): -0.012, ('T', 14): -0.193, ('A', 15): 0.388, ('C', 15): -0.402, ('G', 15): 0.094, ('A', 16): -0.014, ('C', 16): 0.209, ('G', 16): -0.218, ('C', 17): -0.239, ('G', 17): 0.317, ('T', 17): 0.082, ('G', 18): 0.491, ('T', 18): -0.428, ('C', 19): 0.082, ('G', 19): 0.158, ('T', 19): -0.306, ('G', 20): 0.088, ('T', 20): -0.188, ('G', 21): -0.324, ('T', 21): 0.389, ('C', 22): -0.730, ('G', 22): 0.520, ('C', 23): 0.277, ('G', 23): -0.413, ('T', 23): 0.223, ('G', 24): 0.032, ('T', 24): -0.153, ('A', 27): 0.099, ('C', 27): -0.046, ('T', 27): -0.103, ('A', 28): 0.279, ('G', 28): -0.223, ('T', 28): -0.190, ('C', 29): -0.021, ('G', 29): 0.147, ('T', 29): -0.207},
    'dinuc': {('GT', 3): -0.620, ('GG', 5): 0.507, ('TA', 5): -0.548, ('TC', 6): 0.327, ('CC', 11): -0.533, ('TG', 11): 0.443, ('GA', 13): 0.449, ('CT', 13): -0.697, ('GC', 14): 0.419, ('AA', 15): -0.499, ('AG', 15): 0.541, ('AC', 18): -0.420, ('GT', 18): 0.499, ('TC', 18): -0.551, ('CG', 19): 0.589, ('AG', 20): -0.542, ('TG', 21): 0.398, ('GT', 23): -0.672, ('GG', 23): 0.533, ('GA', 27): -0.580, ('CT', 28): 0.471},
    'intercept': 0.59763615, 'gc_high': -0.1665878, 'gc_low': -0.2026259
}

# --- 数据类 ---
@dataclass
class DesignConfig:
    """启动子替换模式的配置类"""
    arm_search_order: List[int]
    guide_len: int
    promoter_search_size: int
    sgrna_num: int
    barcode_len: int
    synthesis_template: Optional[str] = None
    cloning_site: Optional[str] = None

@dataclass
class Gene:
    """通用基因信息容器"""
    id: str; start: int; end: int; strand: int
    cds_5prime_start: int; cds_3prime_start: int; cds_3prime_end: int
    cds_min_coord: int; cds_max_coord: int

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
    sgrna_score: float
    sgrna_cut_distance: int
    upstream_arm: str
    downstream_arm: str
    barcode: str
    final_oligo_for_synthesis: str
    design_scenario: Optional[str] = None
    sgrna_cut_site: Optional[int] = None
    upstream_deletion_len: Optional[int] = None
    design_quality: Optional[str] = None


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
                if sort_values:
                    return sorted((val1, val2))
                else:
                    if val1 < val2:
                        raise ValueError(f"首选长度 ({val1}) 必须大于等于最小长度 ({val2}).")
                    return val1, val2
            else:
                raise ValueError("输入必须是单个数字或由冒号分隔的两个数字。")
        except (ValueError, IndexError) as e:
            raise ValueError(f"无效的 {param_name} 格式 '{param_str}'. 错误: {e}") from e

    @staticmethod
    def create_arm_search_order(preferred_hr: int, min_hr: int) -> List[int]:
        return list(range(preferred_hr, min_hr - 1, -1))

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
                cds_starts = [c.start - 1 for c in cds_features]
                cds_ends = [c.end for c in cds_features]
                strand = 1 if gene.strand == '+' else -1
                cds_min_coord = min(cds_starts)
                cds_max_coord = max(cds_ends)
                cds_5prime_start = cds_min_coord if strand == 1 else cds_max_coord - 3
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
        if length <= 0:
            return ""

        if genome_type == 'linear':
            start = max(0, start)
            end = min(base_genome_len, end)
            if start >= end: return ""
            return str(search_genome[start:end])
        
        elif genome_type == 'circle':
            shifted_start = start + base_genome_len
            return str(search_genome[shifted_start : shifted_start + length])
        else:
            raise ValueError(f"未知的基因组类型: {genome_type}")

    @staticmethod
    def contains_restriction_sites(sequence: str, sites: Optional[List[str]]) -> bool:
        if not sites or not sequence: return False
        seq_upper = sequence.upper()
        for site in sites:
            site_upper = site.upper()
            if site_upper in seq_upper: return True
            if SequenceUtils.get_reverse_complement(site_upper) in seq_upper: return True
        return False

    @staticmethod
    def check_final_oligo_with_exemption(final_oligo: str, restriction_sites: Optional[List[str]], exempt_site: Optional[str]) -> bool:
        if not restriction_sites: return False
        if not exempt_site: return SequenceUtils.contains_restriction_sites(final_oligo, restriction_sites)
        
        max_site_len = 0
        for site in restriction_sites:
            if len(site) > max_site_len: max_site_len = len(site)
        if max_site_len == 0: return False

        junction_check_len = max_site_len - 1
        if len(exempt_site) <= 2 * junction_check_len:
            placeholder = 'X' * len(exempt_site)
            masked_oligo = final_oligo.replace(exempt_site, placeholder)
        else:
            prefix = exempt_site[:junction_check_len]
            suffix = exempt_site[-junction_check_len:]
            middle_len = len(exempt_site) - 2 * junction_check_len
            placeholder = 'X' * middle_len
            masked_exempt_site = prefix + placeholder + suffix
            masked_oligo = final_oligo.replace(exempt_site, masked_exempt_site)
            
        return SequenceUtils.contains_restriction_sites(masked_oligo, restriction_sites)
    
    @staticmethod
    def generate_unique_barcode(length: int, existing_barcodes: Set[str], restriction_sites: Optional[List[str]]) -> str:
        max_attempts = 1000
        for _ in range(max_attempts):
            barcode = "".join(random.choices("ATGC", k=length))
            if barcode in existing_barcodes: continue
            gc_count = barcode.count('G') + barcode.count('C')
            min_gc = math.ceil(length * 0.3)
            max_gc = math.floor(length * 0.7)
            if not (min_gc <= gc_count <= max_gc): continue
            if 'A'*6 in barcode or 'T'*6 in barcode or 'G'*6 in barcode or 'C'*6 in barcode: continue
            if SequenceUtils.contains_restriction_sites(barcode, restriction_sites): continue
            return barcode
        raise RuntimeError(f"在 {max_attempts} 次尝试后未能生成唯一的barcode。")

class SGRNADesigner:
    def __init__(self, config: DesignConfig): 
        self.config = config

    def match_pam(self, seq: str, pattern: str) -> bool:
        if len(seq) != len(pattern): return False
        for i in range(len(pattern)):
            if pattern[i] != 'N' and pattern[i] != seq[i]: return False
        return True

    @staticmethod
    def score_sgrna(sgrna_sequence_30bp: str) -> float:
        if len(sgrna_sequence_30bp) != 30: return 0.0
        guide_seq = sgrna_sequence_30bp[4:24]
        gc_count = guide_seq.count('G') + guide_seq.count('C')
        score = MODEL_WEIGHTS['intercept']
        for i, nt in enumerate(sgrna_sequence_30bp):
            if (nt, i+1) in MODEL_WEIGHTS['pos']: score += MODEL_WEIGHTS['pos'][(nt, i+1)]
        for i in range(29):
            dinuc = sgrna_sequence_30bp[i:i+2]
            if (dinuc, i+1) in MODEL_WEIGHTS['dinuc']: score += MODEL_WEIGHTS['dinuc'][(dinuc, i+1)]
        if gc_count > 10: score += (gc_count - 10) * MODEL_WEIGHTS['gc_high']
        if gc_count < 10: score += (10 - gc_count) * MODEL_WEIGHTS['gc_low']
        return 1 / (1 + math.exp(-score))
    
    def _get_cut_site(self, sgrna_genomic_start: int, sgrna_strand: str, pam_len: int) -> int:
        cut_offset = -3 
        if sgrna_strand == "forward":
            pam_start_pos = sgrna_genomic_start + self.config.guide_len
            return pam_start_pos + cut_offset
        else: # reverse
            pam_start_pos = sgrna_genomic_start
            return pam_start_pos + pam_len + cut_offset
            
    def _calculate_cut_distance(self, cut_site: int, target_site: int) -> int:
        return abs(cut_site - target_site)

    def find_sgrnas_in_region(self, base_genome_len: int, search_genome: Seq, search_start: int, search_end: int, target_site: int, pam: str, restriction_sites: Optional[List[str]], genome_type: str) -> List[dict]:
        all_sgrnas, seen = [], set()
        pam_len, guide_len = len(pam), self.config.guide_len
        pam_fwd_pattern, pam_rev_pattern = pam.upper(), SequenceUtils.get_reverse_complement(pam).upper()

        buffer = guide_len + pam_len + 4
        search_region_seq = SequenceUtils.get_sequence(base_genome_len, search_genome, search_start - buffer, search_end + buffer, genome_type)
        
        for i in range(search_start, search_end):
            local_i = i - (search_start - buffer)

            if local_i + guide_len + pam_len <= len(search_region_seq):
                win_fwd = search_region_seq[local_i : local_i + guide_len + pam_len]
                pam_fwd = str(win_fwd[guide_len:]).upper()
                if self.match_pam(pam_fwd, pam_fwd_pattern):
                    guide = str(win_fwd[:guide_len]).upper()
                    if guide not in seen and not SequenceUtils.contains_restriction_sites(guide, restriction_sites):
                        context = search_region_seq[local_i - 4 : local_i + guide_len + pam_len + 3]
                        if len(context) >= 30:
                            score = self.score_sgrna(context[len(context)-30:])
                            cut_site = self._get_cut_site(i, "forward", pam_len)
                            dist = self._calculate_cut_distance(cut_site, target_site)
                            all_sgrnas.append({'seq': guide, 'pam': pam_fwd, 'strand': 'forward', 'score': score, 'distance': dist, 'cut_site': cut_site, 'genomic_start': i})
                            seen.add(guide)
            
            if local_i + pam_len + guide_len <= len(search_region_seq):
                win_rev = search_region_seq[local_i : local_i + pam_len + guide_len]
                pam_rev_comp = str(win_rev[:pam_len]).upper()
                if self.match_pam(pam_rev_comp, pam_rev_pattern):
                    guide = SequenceUtils.get_reverse_complement(str(win_rev[pam_len:])).upper()
                    if guide not in seen and not SequenceUtils.contains_restriction_sites(guide, restriction_sites):
                        context_plus = search_region_seq[local_i - 3 : local_i + pam_len + guide_len + 4]
                        if len(context_plus) >= 30:
                            context = SequenceUtils.get_reverse_complement(context_plus)
                            score = self.score_sgrna(context[:30])
                            cut_site = self._get_cut_site(i, "reverse", pam_len)
                            dist = self._calculate_cut_distance(cut_site, target_site)
                            all_sgrnas.append({'seq': guide, 'pam': SequenceUtils.get_reverse_complement(pam_rev_comp), 'strand': 'reverse', 'score': score, 'distance': dist, 'cut_site': cut_site, 'genomic_start': i})
                            seen.add(guide)
                            
        all_sgrnas.sort(key=lambda x: (-x['score'], x['distance']))
        return all_sgrnas
        
class CRISPRDesigner:
    def __init__(self, config: DesignConfig, genome_seq: Seq, genome_len: int, genome_type: str):
        self.config = config
        self.sgrna_designer = SGRNADesigner(config)
        self.genome_seq = genome_seq
        self.genome_len = genome_len
        self.genome_type = genome_type
        if self.genome_type == 'circle':
            self.search_genome = self.genome_seq + self.genome_seq
        else:
            self.search_genome = self.genome_seq

    def _assemble_final_oligo(self, sgrna_fwd: str, pam: str, upstream_arm: str, downstream_arm: str, barcode: str) -> str:
        sgrna_rc = SequenceUtils.get_reverse_complement(sgrna_fwd)
        pam_rc = SequenceUtils.get_reverse_complement(pam)
        if self.config.synthesis_template:
            template = self.config.synthesis_template
            replacements = {
                "{sgRNA_fwd}": sgrna_fwd, "{sgRNA_rc}": sgrna_rc, "{pam}": pam, "{pam_rc}": pam_rc,
                "{upstream_arm}": upstream_arm, "{downstream_arm}": downstream_arm, "{barcode}": barcode,
                "{insert}": "",
                "{exempt_restriction_site}": self.config.cloning_site or ""
            }
            for placeholder, value in replacements.items():
                template = template.replace(placeholder, value)
            return template
        else:
            raise ValueError("硬编码的oligo拼接逻辑已被移除，请提供合成模板文件。")

    def _check_coordinate_overlap(self, sgrna: Dict, pam: str, arm_coords: Tuple[int, int, int, int]) -> bool:
        """检查sgRNA靶点区域与同源臂区域是否存在坐标重叠"""
        lha_start, lha_end, rha_start, rha_end = arm_coords
        pam_len = len(pam)
        guide_len = self.config.guide_len
        
        sgrna_target_start, sgrna_target_end = 0, 0
        if sgrna['strand'] == 'forward':
            sgrna_target_start = sgrna['genomic_start']
            sgrna_target_end = sgrna['genomic_start'] + guide_len
        else: # reverse
            sgrna_target_start = sgrna['genomic_start'] + pam_len
            sgrna_target_end = sgrna['genomic_start'] + pam_len + guide_len

        # 检查与上游臂的重叠
        if max(sgrna_target_start, lha_start) < min(sgrna_target_end, lha_end):
            return True
        # 检查与下游臂的重叠
        if max(sgrna_target_start, rha_start) < min(sgrna_target_end, rha_end):
            return True
        
        return False

    def design_promoter_replacement_for_gene(self, gene: Gene, gene_index: int, all_genes: List[Gene], pam: str, used_barcodes: Set[str], restriction_sites: Optional[List[str]], relax_mode: bool = False) -> List[KnockinDesignResult]:
        found_designs = []
        
        # --- 根据基因链方向确定搜索参数 ---
        if gene.strand == 1:
            search_start = gene.cds_5prime_start - self.config.promoter_search_size
            search_end = gene.cds_5prime_start
            anchor_point = gene.cds_5prime_start
        else: # gene.strand == -1
            anchor_point = gene.cds_5prime_start + 3
            search_start = anchor_point
            search_end = anchor_point + self.config.promoter_search_size

        sgrnas_to_process = self.sgrna_designer.find_sgrnas_in_region(
            self.genome_len, self.search_genome, search_start, search_end, anchor_point, pam, restriction_sites, self.genome_type
        )
        logging.info(f"基因 {gene.id} ({'+' if gene.strand == 1 else '-'}): 找到 {len(sgrnas_to_process)} 个候选sgRNA。")

        # --- 两阶段搜索循环 ---
        for pass_type in ["optimal", "tolerant"]:
            if pass_type == "optimal":
                logging.info(f"基因 {gene.id}: 开始严格模式搜索 (无重叠)...")
            else:
                if len(found_designs) >= self.config.sgrna_num: break
                logging.info(f"基因 {gene.id}: 严格模式设计不足，开始宽容模式搜索 (允许重叠)...")

            for sgrna in sgrnas_to_process:
                if gene.cds_min_coord <= sgrna['cut_site'] < gene.cds_max_coord:
                    continue

                design_found_for_this_sgrna = False
                for arm_len in self.config.arm_search_order:
                    offset_range = range(51)
                    for offset in offset_range:
                        # --- 计算同源臂和删除区域坐标 ---
                        lha_start, lha_end, rha_start, rha_end = 0, 0, 0, 0
                        deletion_start, deletion_end = 0, 0
                        
                        if gene.strand == 1:
                            lha_end = sgrna['cut_site'] - offset
                            rha_start = anchor_point
                            if lha_end >= rha_start: continue
                            deletion_start, deletion_end = lha_end, rha_start
                            lha_start, rha_end = lha_end - arm_len, rha_start + arm_len
                        else: # strand == -1
                            lha_end = anchor_point
                            rha_start = sgrna['cut_site'] + offset
                            if rha_start <= lha_end: continue
                            deletion_start, deletion_end = lha_end, rha_start
                            lha_start, rha_end = lha_end - arm_len, rha_start + arm_len

                        # --- 检查邻近基因保护 ---
                        if not relax_mode:
                            if gene.strand == 1 and gene_index > 0:
                                upstream_gene = all_genes[gene_index - 1]
                                if gene.start - upstream_gene.end < 1000:
                                    if max(deletion_start, upstream_gene.cds_min_coord) < min(deletion_end, upstream_gene.cds_max_coord):
                                        continue
                            elif gene.strand == -1 and gene_index < len(all_genes) - 1:
                                downstream_gene = all_genes[gene_index + 1]
                                if downstream_gene.start - gene.end < 1000:
                                    if max(deletion_start, downstream_gene.cds_min_coord) < min(deletion_end, downstream_gene.cds_max_coord):
                                        continue
                        
                        arm_coords = (lha_start, lha_end, rha_start, rha_end)
                        has_overlap = self._check_coordinate_overlap(sgrna, pam, arm_coords)

                        # --- 根据pass_type应用重叠过滤器 ---
                        if pass_type == "optimal" and has_overlap: continue
                        if pass_type == "tolerant" and not has_overlap: continue
                        
                        # --- 获取序列并进行检查 ---
                        upstream_arm = SequenceUtils.get_sequence(self.genome_len, self.search_genome, lha_start, lha_end, self.genome_type)
                        downstream_arm = SequenceUtils.get_sequence(self.genome_len, self.search_genome, rha_start, rha_end, self.genome_type)

                        if not (upstream_arm and downstream_arm and len(upstream_arm) == arm_len and len(downstream_arm) == arm_len and \
                               not SequenceUtils.contains_restriction_sites(upstream_arm, restriction_sites) and \
                               not SequenceUtils.contains_restriction_sites(downstream_arm, restriction_sites)):
                            continue

                        # --- 生成并保存设计 ---
                        barcode = SequenceUtils.generate_unique_barcode(self.config.barcode_len, used_barcodes, restriction_sites)
                        final_oligo = self._assemble_final_oligo(sgrna['seq'], sgrna['pam'], upstream_arm, downstream_arm, barcode)
                        
                        if SequenceUtils.check_final_oligo_with_exemption(final_oligo, restriction_sites, self.config.cloning_site):
                            continue
                        
                        used_barcodes.add(barcode)
                        result = KnockinDesignResult(
                            gene_id=gene.id, design_type='promoter_replace', insertion_site=anchor_point,
                            arm_length=arm_len, sgrna_seq=sgrna['seq'], sgrna_pam=sgrna['pam'],
                            sgrna_strand=sgrna['strand'], sgrna_score=sgrna['score'],
                            sgrna_cut_distance=abs(sgrna['cut_site'] - anchor_point),
                            upstream_arm=upstream_arm, downstream_arm=downstream_arm, barcode=barcode, final_oligo_for_synthesis=final_oligo,
                            sgrna_cut_site=sgrna['cut_site'],
                            upstream_deletion_len=(deletion_end - deletion_start),
                            design_quality="Optimal" if not has_overlap else "Tolerated_Overlap",
                            design_scenario=f"Upstream Cut (Offset={offset})"
                        )
                        found_designs.append(result)
                        if len(found_designs) >= self.config.sgrna_num: design_found_for_this_sgrna = True; break
                    if design_found_for_this_sgrna: break
                if len(found_designs) >= self.config.sgrna_num: break
            if len(found_designs) >= self.config.sgrna_num: break

        if not found_designs: 
            logging.warning(f"基因 {gene.id}: 设计失败。")
        return found_designs[:self.config.sgrna_num]


class ResultProcessor:
    @staticmethod
    def save_results(designs, failed_genes, output_file, design_type_info):
        if not designs:
            logging.warning(f"没有成功的 {design_type_info} 设计可供保存。")
        else:
            results = []
            for d in designs:
                common_dict = {
                    "Gene_ID": d.gene_id, "Status": "Success",
                    "sgRNA_Sequence": d.sgrna_seq.upper(),
                    "Barcode": d.barcode.upper(),
                    "Final_Oligo_for_Synthesis": d.final_oligo_for_synthesis.upper(),
                    "sgRNA_Score": f"{d.sgrna_score:.4f}",
                    "Arm_Length": d.arm_length,
                    "Upstream_Arm": d.upstream_arm.upper(),
                    "Downstream_Arm": d.downstream_arm.upper(),
                    "sgRNA_PAM": d.sgrna_pam.upper(),
                    "sgRNA_Strand": d.sgrna_strand,
                    "sgRNA_Cut_Site": d.sgrna_cut_site,
                    "Design_Type": d.design_type,
                    "Design_Scenario": d.design_scenario,
                    "Design_Quality": d.design_quality
                }
                results.append(common_dict)
            
            df = pd.DataFrame(results)
            core_cols = ["Gene_ID", "Status", "Design_Quality", "sgRNA_Sequence", "Barcode", "Final_Oligo_for_Synthesis"]
            other_cols = sorted([col for col in df.columns if col not in core_cols])
            df = df[core_cols + other_cols]
            df.to_csv(output_file, index=False)
            logging.info(f"{design_type_info} 结果已保存至: {output_file}")

        if failed_genes:
            failed_data = [{"Gene_ID": g.id, "Status": "Failed"} for g in failed_genes]
            failed_df = pd.DataFrame(failed_data)
            failed_output_path = os.path.splitext(output_file)[0] + "_failed.csv"
            failed_df.to_csv(failed_output_path, index=False)
            logging.info(f"失败的基因ID已保存至: {failed_output_path}")

        successful_designs_count = len(designs) if designs else 0
        failed_genes_count = len(failed_genes) if failed_genes else 0
        logging.info(f"总结: {successful_designs_count} 个成功设计, {failed_genes_count} 个基因设计失败。")

def run_promoter_replacement_pipeline(args, config, genes, genome_processor):
    designer = CRISPRDesigner(config, genome_processor.genome_seq, genome_processor.genome_len, args.genome_type)
    successful_designs, used_barcodes = [], set()
    min_arm_len = config.arm_search_order[-1]
    
    genes_to_process = [g for g in genes if g.cds_max_coord - g.cds_min_coord >= min_arm_len]
    logging.info(f"正在处理 {len(genes_to_process)} 个长度足够进行启动子替换设计的基因。")
    
    for i, gene in enumerate(genes_to_process):
        if (i + 1) % 100 == 0: logging.info(f"正在处理基因 {i+1}/{len(genes_to_process)}...")
        
        results = designer.design_promoter_replacement_for_gene(
            gene=gene, 
            gene_index=i, 
            all_genes=genes_to_process, 
            pam=args.pam, 
            used_barcodes=used_barcodes, 
            restriction_sites=args.restriction_site,
            relax_mode=args.relax
        )
        if results: successful_designs.extend(results)
        
    successful_gene_ids = {d.gene_id for d in successful_designs}
    failed_genes = [g for g in genes_to_process if g.id not in successful_gene_ids]
    
    ResultProcessor.save_results(successful_designs, failed_genes, args.output, "Promoter Replacement")

def main():
    parser = argparse.ArgumentParser(description="启动子替换 (Promoter Replacement) 设计工具 (Cas9)", formatter_class=argparse.RawTextHelpFormatter)
    
    # 核心参数
    parser.add_argument("--input_fna", required=True, help="输入的基因组FASTA文件。")
    parser.add_argument("--input_gff", required=True, help="输入的GFF3注释文件。")
    parser.add_argument("--output", required=True, help="输出的CSV文件路径。")
    parser.add_argument("--genome_type", type=str, default="linear", choices=['linear', 'circle'], help="指定基因组类型 (默认: linear)。")
    parser.add_argument("--pam", type=str, default="NGG", help="PAM序列 (默认: 'NGG')。")
    parser.add_argument("--HR_len", type=str, default="50:40", help="同源臂长度。格式为 '首选:最小' (例如, '50:40') 或单个数字。")
    parser.add_argument("--sgRNA_num", type=int, default=1, help="为每个基因生成的设计方案数量 (默认: 1)。")
    parser.add_argument("--promoter_search_size", type=int, default=150, help="在ATG上游搜索sgRNA的区域大小 (bp)。")
    parser.add_argument("--barcode_len", type=int, default=8, help="指定barcode序列的长度 (默认: 8)。")
    
    # 序列筛选参数
    parser.add_argument("--restriction_site", type=str, nargs='+', help="在可变区中需要避免的限制性酶切位点。")
    parser.add_argument("--synthesis_template", type=str, required=True, help="定义最终oligo结构的文本模板文件路径。")
    parser.add_argument("--cloning_site", type=str, help="一个特殊的酶切位点，它将被放置在模板的{exempt_restriction_site}占位符中，并且不会在模板骨架检查中被排除。")
    parser.add_argument("--relax", action="store_true", help="不考虑是否破坏临近基因的CDS区。")
    
    args = parser.parse_args()

    try:
        synthesis_template_content = None
        if args.synthesis_template:
            try:
                with open(args.synthesis_template, 'r') as f:
                    synthesis_template_content = f.read().strip()
                logging.info(f"成功从以下路径加载合成模板: {args.synthesis_template}")
                if args.restriction_site:
                    dummy_filled = synthesis_template_content.replace("{exempt_restriction_site}", "").format(
                        sgRNA_fwd="", sgRNA_rc="", pam="", pam_rc="",
                        upstream_arm="", downstream_arm="", barcode="", insert="",
                    )
                    if SequenceUtils.contains_restriction_sites(dummy_filled, args.restriction_site):
                        logging.error("致命错误: 在您的合成模板的固定骨架中发现了需要排除的限制性酶切位点。")
                        logging.error("请从模板文件中移除该位点，或选择其他需要排除的酶切位点。")
                        sys.exit(1)
            except FileNotFoundError:
                logging.error(f"致命错误: 找不到合成模板文件: {args.synthesis_template}")
                sys.exit(1)
            except KeyError as e:
                logging.error(f"致命错误: 您的合成模板文件包含未知占位符: {e}")
                sys.exit(1)

        preferred_hr, min_hr = ConfigParser.parse_range_param(args.HR_len, "同源臂长度", sort_values=False)
        arm_search_order = ConfigParser.create_arm_search_order(preferred_hr, min_hr)
        
        config = DesignConfig(
            arm_search_order=arm_search_order,
            guide_len=20,
            promoter_search_size=args.promoter_search_size,
            synthesis_template=synthesis_template_content,
            cloning_site=args.cloning_site,
            sgrna_num=args.sgRNA_num,
            barcode_len=args.barcode_len
        )
        
        genome_processor = GenomeProcessor(args.input_fna, args.input_gff)
        genome_processor.load_genome()
        genes = genome_processor.parse_genes()

        run_promoter_replacement_pipeline(args, config, genes, genome_processor)

    except Exception as e:
        logging.error(f"发生严重错误: {e}", exc_info=True)

if __name__ == "__main__":
    main()


