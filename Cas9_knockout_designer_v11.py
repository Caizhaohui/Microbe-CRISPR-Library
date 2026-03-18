# filename: Cas9_knockout_designer_v11.py
# Version 11.0 - 动态间距自适应优化版（基于 V10）
#
# V11 相对 V10 的改进：
#   1. threshold=0 回退时主动最大化设计多样性（核心改进）：
#      - 针对 V10 中 40 个基因（0.4%）因候选池稠密无法达到 100bp 间距的问题
#      - V10 在 threshold=0 回退时，按质量顺序取 top-1 和 top-2，两者切割位点可能相邻
#      - V11 在 threshold=0 回退时，对设计-2 重排候选池：按"与设计-1 最大间距"优先
#        * 设计-1：仍取质量最佳候选（不变）
#        * 设计-2：在全部候选中选间距最大者（同等间距时取质量最佳）
#        * 这保证了即使在无间距约束下，两套设计仍尽可能分散
#      - 对 sgrna_num>=2 的多套设计同理：每套新设计都在已有设计中最大化最小间距
#
# V10 核心特性（继承保留）：
#   1. 删除策略内联：Mt 删除模式通过 --dele_model 参数切换，无猴子补丁
#   2. 穷举非3倍数删除长度：legacy & cut_window 双模式全枚举，避免遗漏最优候选
#   3. 全局候选池 + 多样性排序：每个 sgRNA top_k=5 候选，统一排序后多样性选取
#   4. 线程安全 Barcode 生成：原子化加锁，修复多线程 barcode 重复风险
#   5. CSV 新增列：Design_Index、Cut_Site_Spacing

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
import threading
from concurrent.futures import ThreadPoolExecutor
from typing import List, Dict, Tuple, Optional, Set
from dataclasses import dataclass, field
import math
import time

try:
    from species_manager import SpeciesConfig
    SPECIES_CONFIG_AVAILABLE = True
except ImportError:
    SPECIES_CONFIG_AVAILABLE = False

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# ============================================================
# V11 新增：动态间距计算函数
# ============================================================

def compute_dynamic_spacing(cds_length: int, default_spacing: int = 100) -> int:
    """
    根据 CDS 长度自动计算最优设计间距。
    
    逻辑：
        - CDS_len >= 2 * default_spacing (200bp)：使用标准间距 (100bp)
        - CDS_len < 2 * default_spacing：返回 0（允许完全重叠，保证双设计覆盖率）
    
    Args:
        cds_length: CDS 序列长度（bp）
        default_spacing: 默认标准间距（默认 100bp）
    
    Returns:
        int: 计算得出的最小设计间距（bp）
    """
    if cds_length < 2 * default_spacing:
        return 0  # 短基因：完全放开间距限制
    return default_spacing  # 标准基因：使用默认间距

# ============================================================
# 密码子表与模型权重
# ============================================================

STANDARD_CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

MODEL_WEIGHTS = {
    'pos': {('A', 2): -0.275, ('C', 2): 0.194, ('T', 2): -0.326, ('A', 3): -0.373, ('C', 3): 0.129, ('G', 3): -0.174, ('A', 4): -0.012, ('C', 4): 0.088, ('T', 4): -0.019, ('A', 5): 0.252, ('C', 5): -0.100, ('T', 5): -0.294, ('A', 6): 0.130, ('C', 6): -0.091, ('G', 6): 0.297, ('A', 7): -0.201, ('C', 7): 0.245, ('G', 7): -0.208, ('A', 11): -0.298, ('C', 11): 0.178, ('T', 11): 0.117, ('C', 12): -0.017, ('G', 12): 0.134, ('T', 12): -0.258, ('A', 13): 0.329, ('C', 13): -0.149, ('T', 13): -0.323, ('A', 14): 0.075, ('G', 14): -0.012, ('T', 14): -0.193, ('A', 15): 0.388, ('C', 15): -0.402, ('G', 15): 0.094, ('A', 16): -0.014, ('C', 16): 0.209, ('G', 16): -0.218, ('C', 17): -0.239, ('G', 17): 0.317, ('T', 17): 0.082, ('G', 18): 0.491, ('T', 18): -0.428, ('C', 19): 0.082, ('G', 19): 0.158, ('T', 19): -0.306, ('G', 20): 0.088, ('T', 20): -0.188, ('G', 21): -0.324, ('T', 21): 0.389, ('C', 22): -0.730, ('G', 22): 0.520, ('C', 23): 0.277, ('G', 23): -0.413, ('T', 23): 0.223, ('G', 24): 0.032, ('T', 24): -0.153, ('A', 27): 0.099, ('C', 27): -0.046, ('T', 27): -0.103, ('A', 28): 0.279, ('G', 28): -0.223, ('T', 28): -0.190, ('C', 29): -0.021, ('G', 29): 0.147, ('T', 29): -0.207},
    'dinuc': {('GT', 3): -0.620, ('GG', 5): 0.507, ('TA', 5): -0.548, ('TC', 6): 0.327, ('CC', 11): -0.533, ('TG', 11): 0.443, ('GA', 13): 0.449, ('CT', 13): -0.697, ('GC', 14): 0.419, ('AA', 15): -0.499, ('AG', 15): 0.541, ('AC', 18): -0.420, ('GT', 18): 0.499, ('TC', 18): -0.551, ('CG', 19): 0.589, ('AG', 20): -0.542, ('TG', 21): 0.398, ('GT', 23): -0.672, ('GG', 23): 0.533, ('GA', 27): -0.580, ('CT', 28): 0.471},
    'intercept': 0.59763615, 'gc_high': -0.1665878, 'gc_low': -0.2026259
}

_COMPLEMENT_TABLE = str.maketrans('ATCGatcgNn', 'TAGCtagcNn')

# ============================================================
# 数据类
# ============================================================

@dataclass
class DesignConfig:
    """Knockout_Cas9 统一配置类（V10：新增 deletion_strategy、min_design_spacing）"""
    arm_search_order: List[int]
    guide_len: int
    sgrna_num: int
    barcode_len: int
    strict: bool
    synthesis_template: Optional[str] = None
    cloning_site: Optional[str] = None
    ko_search_start_pct: Optional[float] = None
    ko_search_end_pct: Optional[float] = None
    min_promoter_size: Optional[int] = None
    max_promoter_size: Optional[int] = None
    species: str = "M_thermophila"
    barcode_gc_min: float = 0.30
    barcode_gc_max: float = 0.80
    barcode_repeat_threshold: int = 5
    neighbor_distance_threshold: int = 1500
    max_oligo_length: Optional[int] = None
    pam_len: int = 3
    # V6：双因素删除约束
    del_pct_min: float = 0.0
    del_pct_max: float = 1.0
    del_bp_min: int = 1
    del_bp_max: int = 10000000
    del_has_pct: bool = False
    del_has_bp: bool = False
    # V8：Mt 物种专项参数
    use_mt_pam: bool = False
    deletion_mode: str = "legacy_length"
    cut_window_upstream: int = 20
    cut_window_downstream: int = 100
    # V10 新增：策略模式 + 多样性约束
    deletion_strategy: str = "normal"   # "normal"（V8行为）或 "Mt"（V9行为：链感知窗口+cut_site+3）
    min_design_spacing: int = 100        # 两套设计切割位点的最小间距（bp），0=不限制


@dataclass
class Gene:
    """通用基因信息容器"""
    seqid: str
    id: str
    start: int
    end: int
    strand: int
    cds_5prime_start: int
    cds_min_coord: int
    cds_max_coord: int


@dataclass
class KnockoutDesignResult:
    """基因敲除设计结果容器（V10新增 design_index、cut_site_spacing）"""
    gene_id: str
    deletion_start: int
    deletion_end: int
    deletion_length: int
    strategy: str
    arm_length: int
    upstream_arm: str
    downstream_arm: str
    sgrna_seq: str
    sgrna_pam: str
    sgrna_strand: str
    sgrna_score: float
    sgrna_cut_site: int
    barcode: str
    final_oligo_for_synthesis: str
    design_index: int = 1                    # 本基因第几套设计（1或2）
    cut_site_spacing: Optional[int] = None   # 与设计1切割位点的间距（设计1为None）
    deletion_warning: str = ""               # 删除长度超标警告（空字符串=正常）


# ============================================================
# 辅助类：参数解析
# ============================================================

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
                    return tuple(sorted((val1, val2)))
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

    @staticmethod
    def parse_del_pct_param(param_str: str) -> Tuple[float, float]:
        cleaned = param_str.replace('%', '').strip()
        parts = cleaned.split(':')
        if len(parts) == 1:
            pct_max = float(parts[0]) / 100.0
            pct_min = 0.0
        elif len(parts) == 2:
            pct_min = float(parts[0]) / 100.0
            pct_max = float(parts[1]) / 100.0
        else:
            raise ValueError(f"无效的 --del_length_per 格式: '{param_str}'")
        if not (0.0 <= pct_min <= 1.0 and 0.0 <= pct_max <= 1.0):
            raise ValueError(f"百分比必须在 0-100 范围内: '{param_str}'")
        return min(pct_min, pct_max), max(pct_min, pct_max)

    @staticmethod
    def parse_del_bp_param(param_str: str) -> Tuple[int, int]:
        parts = param_str.strip().split(':')
        if len(parts) == 1:
            bp_max = int(float(parts[0]))
            bp_min = 1
        elif len(parts) == 2:
            bp_min = int(float(parts[0]))
            bp_max = int(float(parts[1]))
        else:
            raise ValueError(f"无效的 --del_length_bp 格式: '{param_str}'")
        return min(bp_min, bp_max), max(bp_min, bp_max)


# ============================================================
# 基因组处理类
# ============================================================

class GenomeProcessor:
    def __init__(self, genome_file: Optional[str] = None, gff_file: Optional[str] = None,
                 gbff_file: Optional[str] = None):
        self.genome_file = genome_file
        self.gff_file = gff_file
        self.gbff_file = gbff_file
        self.genome_seqs: Dict[str, Seq] = {}
        self.genome_lens: Dict[str, int] = {}
        self.genome_seq = None
        self.genome_len = 0

    def load_genome(self) -> None:
        source_file = self.gbff_file if self.gbff_file else self.genome_file
        logging.info(f"正在加载基因组: {source_file}...")
        try:
            if self.gbff_file:
                genome_record = SeqIO.to_dict(SeqIO.parse(self.gbff_file, "genbank"))
            else:
                genome_record = SeqIO.to_dict(SeqIO.parse(self.genome_file, "fasta"))
            if not genome_record:
                raise ValueError("基因组文件中没有序列")
            self.genome_seqs = {rid: rec.seq.upper() for rid, rec in genome_record.items()}
            self.genome_lens = {rid: len(seq) for rid, seq in self.genome_seqs.items()}
            chrom_id = next(iter(self.genome_seqs.keys()))
            self.genome_seq = self.genome_seqs[chrom_id]
            self.genome_len = self.genome_lens[chrom_id]
            total_bp = sum(self.genome_lens.values())
            logging.info(f"基因组加载完成: {len(self.genome_seqs)} 条序列, 总长度 {total_bp} bp")
        except Exception as e:
            logging.error(f"加载基因组失败: {e}")
            raise

    def parse_genes(self) -> List[Gene]:
        if self.gbff_file:
            return self._parse_genes_from_gbff()
        logging.info(f"正在解析GFF: {self.gff_file}...")
        db_fn = tempfile.NamedTemporaryFile(delete=False).name
        genes = []
        try:
            db = gffutils.create_db(self.gff_file, dbfn=db_fn, force=True, keep_order=True,
                                    merge_strategy='merge', sort_attribute_values=True)
            missing_seqid_count = 0
            for gene in db.features_of_type('gene', order_by='start'):
                locus_tag = gene.attributes.get('locus_tag', [None])[0]
                if not locus_tag:
                    continue
                if gene.seqid not in self.genome_seqs:
                    missing_seqid_count += 1
                    continue
                cds_features = list(db.children(gene, featuretype='CDS', order_by='start'))
                if not cds_features:
                    continue
                cds_starts = [c.start - 1 for c in cds_features]
                cds_ends = [c.end for c in cds_features]
                strand = 1 if gene.strand == '+' else -1
                cds_min_coord = min(cds_starts)
                cds_max_coord = max(cds_ends)
                cds_5prime_start = cds_min_coord if strand == 1 else cds_max_coord - 3
                genes.append(Gene(seqid=gene.seqid, id=locus_tag, start=gene.start - 1, end=gene.end,
                                  strand=strand, cds_5prime_start=cds_5prime_start,
                                  cds_min_coord=cds_min_coord, cds_max_coord=cds_max_coord))
        finally:
            try:
                if os.path.exists(db_fn):
                    time.sleep(0.1)
                    os.unlink(db_fn)
            except (PermissionError, OSError) as e:
                logging.warning(f"无法删除临时文件 {db_fn}: {e}")
        genes.sort(key=lambda g: (g.seqid, g.start))
        if missing_seqid_count > 0:
            logging.warning(f"有 {missing_seqid_count} 个基因所在seqid不在FASTA中，已跳过")
        logging.info(f"找到 {len(genes)} 个包含有效CDS的基因")
        return genes

    def _parse_genes_from_gbff(self) -> List[Gene]:
        logging.info(f"正在解析GBFF注释: {self.gbff_file}...")
        genes: List[Gene] = []
        missing_locus_tag_count = 0
        for record in SeqIO.parse(self.gbff_file, "genbank"):
            feature_map: Dict[str, Dict] = {}
            for feat in record.features:
                if feat.type not in {"gene", "CDS"}:
                    continue
                qualifiers = feat.qualifiers
                locus_tag = (qualifiers.get("locus_tag", [None])[0]
                             or qualifiers.get("gene", [None])[0]
                             or qualifiers.get("old_locus_tag", [None])[0])
                if not locus_tag:
                    missing_locus_tag_count += 1
                    continue
                entry = feature_map.setdefault(locus_tag, {"gene": None, "cds_starts": [], "cds_ends": [], "strand": None})
                strand = -1 if feat.location.strand == -1 else 1
                if entry["strand"] is None:
                    entry["strand"] = strand
                start = int(feat.location.start)
                end = int(feat.location.end)
                if feat.type == "gene":
                    entry["gene"] = (start, end, strand)
                else:
                    entry["cds_starts"].append(start)
                    entry["cds_ends"].append(end)
            for locus_tag, entry in feature_map.items():
                cds_starts = entry["cds_starts"]
                cds_ends = entry["cds_ends"]
                if not cds_starts:
                    continue
                gene_info = entry["gene"]
                if gene_info is not None:
                    gene_start, gene_end, strand = gene_info
                else:
                    gene_start, gene_end = min(cds_starts), max(cds_ends)
                    strand = entry["strand"] if entry["strand"] is not None else 1
                cds_min_coord = min(cds_starts)
                cds_max_coord = max(cds_ends)
                cds_5prime_start = cds_min_coord if strand == 1 else cds_max_coord - 3
                genes.append(Gene(seqid=record.id, id=locus_tag, start=gene_start, end=gene_end,
                                  strand=strand, cds_5prime_start=cds_5prime_start,
                                  cds_min_coord=cds_min_coord, cds_max_coord=cds_max_coord))
        genes.sort(key=lambda g: (g.seqid, g.start))
        if missing_locus_tag_count > 0:
            logging.warning(f"GBFF中有 {missing_locus_tag_count} 个特征缺少locus_tag，已跳过")
        logging.info(f"从GBFF找到 {len(genes)} 个包含有效CDS的基因")
        return genes


# ============================================================
# 序列工具类
# ============================================================

class SequenceUtils:
    @staticmethod
    def get_reverse_complement(seq: str) -> str:
        return seq.translate(_COMPLEMENT_TABLE)[::-1]

    @staticmethod
    def get_sequence(base_genome_len: int, search_genome, start: int, end: int, genome_type: str) -> str:
        length = end - start
        if length <= 0:
            return ""
        if genome_type == 'linear':
            start = max(0, start)
            end = min(base_genome_len, end)
            if start >= end:
                return ""
            return str(search_genome[start:end])
        elif genome_type == 'circle':
            shifted_start = start + base_genome_len
            return str(search_genome[shifted_start: shifted_start + length])
        else:
            raise ValueError(f"未知的基因组类型: {genome_type}")

    @staticmethod
    def contains_restriction_sites(sequence: str, sites: Optional[List[str]]) -> bool:
        if not sites or not sequence:
            return False
        seq_upper = sequence.upper()
        for site in sites:
            site_upper = site.upper()
            if site_upper in seq_upper:
                return True
            if SequenceUtils.get_reverse_complement(site_upper) in seq_upper:
                return True
        return False

    _KNOWN_PLACEHOLDERS = [
        '{sgRNA_fwd}', '{sgRNA_rc}', '{pam}', '{pam_rc}',
        '{upstream_arm}', '{downstream_arm}', '{barcode}',
        '{exempt_restriction_site}', '{insert}'
    ]

    @staticmethod
    def parse_template_segments(template: str) -> List[Tuple[str, str]]:
        segments: List[Tuple[str, str]] = []
        pos = 0
        while pos < len(template):
            next_ph_pos = len(template)
            next_ph = None
            for ph in SequenceUtils._KNOWN_PLACEHOLDERS:
                idx = template.find(ph, pos)
                if idx != -1 and idx < next_ph_pos:
                    next_ph_pos = idx
                    next_ph = ph
            if next_ph is None:
                segments.append(('fixed', template[pos:]))
                break
            if next_ph_pos > pos:
                segments.append(('fixed', template[pos:next_ph_pos]))
            segments.append(('variable', next_ph))
            pos = next_ph_pos + len(next_ph)
        return segments

    @staticmethod
    def check_junctions_for_restriction_sites(final_oligo: str, template_segments: List[Tuple[str, str]],
                                               variable_lengths: Dict[str, int],
                                               restriction_sites: Optional[List[str]]) -> bool:
        if not restriction_sites or not template_segments:
            return False
        max_site_len = max(len(s) for s in restriction_sites)
        w = max_site_len - 1
        seg_bounds: List[Tuple[str, int, int]] = []
        pos = 0
        for seg_type, seg_val in template_segments:
            seg_len = len(seg_val) if seg_type == 'fixed' else variable_lengths.get(seg_val, 0)
            seg_bounds.append((seg_type, pos, pos + seg_len))
            pos += seg_len
        for i in range(len(seg_bounds) - 1):
            type_a, _, end_a = seg_bounds[i]
            type_b, start_b, _ = seg_bounds[i + 1]
            if type_a == type_b:
                continue
            window_start = max(0, end_a - w)
            window_end = min(len(final_oligo), start_b + w)
            junction_seq = final_oligo[window_start:window_end]
            if SequenceUtils.contains_restriction_sites(junction_seq, restriction_sites):
                return True
        return False

    @staticmethod
    def generate_unique_barcode(length: int, existing_barcodes: Set[str],
                                restriction_sites: Optional[List[str]],
                                gc_min: float = 0.35, gc_max: float = 0.65,
                                repeat_threshold: int = 5) -> str:
        """非线程安全版（兼容保留，内部单线程调用时使用）"""
        max_attempts = 1000
        for _ in range(max_attempts):
            barcode = "".join(random.choices("ATGC", k=length))
            if barcode in existing_barcodes:
                continue
            gc_count = barcode.count('G') + barcode.count('C')
            min_gc = math.ceil(length * gc_min)
            max_gc = math.floor(length * gc_max)
            if not (min_gc <= gc_count <= max_gc):
                continue
            if any(nt * (repeat_threshold + 1) in barcode for nt in 'ATGC'):
                continue
            if SequenceUtils.contains_restriction_sites(barcode, restriction_sites):
                continue
            return barcode
        raise RuntimeError(f"在 {max_attempts} 次尝试后未能生成唯一的barcode。")

    @staticmethod
    def generate_unique_barcode_atomic(length: int, used_barcodes: Set[str],
                                       lock: threading.Lock,
                                       restriction_sites: Optional[List[str]],
                                       gc_min: float = 0.30, gc_max: float = 0.80,
                                       repeat_threshold: int = 5) -> str:
        """
        V10：线程安全原子化 barcode 生成。
        策略：先检查本地约束（无锁），通过后加锁原子性地检查全局唯一性并写入。
        修复 V8 多线程下跨基因 barcode 重复的竞争条件。
        """
        max_attempts = 2000
        for _ in range(max_attempts):
            barcode = "".join(random.choices("ATGC", k=length))
            # 本地约束（无需加锁，速度快）
            gc_count = barcode.count('G') + barcode.count('C')
            min_gc = math.ceil(length * gc_min)
            max_gc = math.floor(length * gc_max)
            if not (min_gc <= gc_count <= max_gc):
                continue
            if any(nt * (repeat_threshold + 1) in barcode for nt in 'ATGC'):
                continue
            if SequenceUtils.contains_restriction_sites(barcode, restriction_sites):
                continue
            # 原子加锁：check-and-add，确保全局唯一
            with lock:
                if barcode not in used_barcodes:
                    used_barcodes.add(barcode)
                    return barcode
            # 被其他线程抢占，继续生成新候选
        raise RuntimeError(f"在 {max_attempts} 次尝试后未能生成唯一的barcode。")


# ============================================================
# sgRNA 设计类
# ============================================================

class SGRNADesigner:
    def __init__(self, config: DesignConfig):
        self.config = config

    def match_pam(self, seq: str, pattern: str) -> bool:
        if len(seq) != len(pattern):
            return False
        for i in range(len(pattern)):
            if pattern[i] != 'N' and pattern[i] != seq[i]:
                return False
        return True

    @staticmethod
    def score_sgrna(sgrna_sequence_30bp: str) -> float:
        if len(sgrna_sequence_30bp) != 30:
            return 0.0
        guide_seq = sgrna_sequence_30bp[4:24]
        gc_count = guide_seq.count('G') + guide_seq.count('C')
        score = MODEL_WEIGHTS['intercept']
        for i, nt in enumerate(sgrna_sequence_30bp):
            if (nt, i + 1) in MODEL_WEIGHTS['pos']:
                score += MODEL_WEIGHTS['pos'][(nt, i + 1)]
        for i in range(29):
            dinuc = sgrna_sequence_30bp[i:i + 2]
            if (dinuc, i + 1) in MODEL_WEIGHTS['dinuc']:
                score += MODEL_WEIGHTS['dinuc'][(dinuc, i + 1)]
        if gc_count > 10:
            score += (gc_count - 10) * MODEL_WEIGHTS['gc_high']
        if gc_count < 10:
            score += (10 - gc_count) * MODEL_WEIGHTS['gc_low']
        return 1 / (1 + math.exp(-score))

    def _get_cut_site(self, sgrna_genomic_start: int, sgrna_strand: str, pam_len: int) -> int:
        """
        计算 Cas9 切割位点。
        normal 策略：PAM 上游 -3 bp（正/反链统一）
        Mt 策略（V9逻辑）：正链同 normal；反向链改为 PAM 右侧 +3 bp（基因组正链坐标），
                         校正反链 cut_site 在基因组坐标上的偏差（相比 normal 差 +6 bp）
        """
        if sgrna_strand == "forward":
            pam_start_pos = sgrna_genomic_start + self.config.guide_len
            return pam_start_pos - 3
        else:  # reverse
            pam_start_pos = sgrna_genomic_start
            if self.config.deletion_strategy == "Mt":
                # V9 Mt 校正：反向链 cut_site 在 PAM 右侧 +3
                return pam_start_pos + pam_len + 3
            else:
                return pam_start_pos + pam_len - 3

    def find_sgrnas_in_region(self, base_genome_len: int, search_genome,
                               search_start: int, search_end: int, pam: str,
                               restriction_sites: Optional[List[str]], genome_type: str) -> List[dict]:
        all_sgrnas, seen = [], set()
        pam_len, guide_len = len(pam), self.config.guide_len
        buffer = guide_len + pam_len + 4
        region_start = search_start - buffer
        search_region_seq = SequenceUtils.get_sequence(base_genome_len, search_genome,
                                                        region_start, search_end + buffer, genome_type)
        region_str = search_region_seq if isinstance(search_region_seq, str) else str(search_region_seq)
        region_len = len(region_str)
        local_search_start = buffer
        local_search_end = buffer + (search_end - search_start)

        _contains_rs = SequenceUtils.contains_restriction_sites
        _rc = SequenceUtils.get_reverse_complement
        _score = self.score_sgrna
        _cut = self._get_cut_site

        if pam.upper() == 'NGG':
            _mt_pam = self.config.use_mt_pam
            p = 0
            while True:
                p = region_str.find('GG', p)
                if p < 0:
                    break
                local_i = p - guide_len - 1
                if local_search_start <= local_i < local_search_end and local_i >= 0 and local_i + guide_len + pam_len <= region_len:
                    if _mt_pam and p > 0 and region_str[p - 1].upper() == 'G':
                        p += 1
                        continue
                    guide = region_str[local_i: local_i + guide_len]
                    if guide not in seen and not _contains_rs(guide, restriction_sites):
                        ctx = region_str[local_i - 4: local_i + guide_len + pam_len + 3]
                        if len(ctx) >= 30:
                            pam_seq = region_str[local_i + guide_len: local_i + guide_len + pam_len]
                            genomic_i = local_i + search_start - buffer
                            all_sgrnas.append({'seq': guide, 'pam': pam_seq, 'strand': 'forward',
                                               'score': _score(ctx[len(ctx) - 30:]),
                                               'cut_site': _cut(genomic_i, "forward", pam_len)})
                            seen.add(guide)
                p += 1
            p = 0
            while True:
                p = region_str.find('CC', p)
                if p < 0:
                    break
                local_i = p
                if local_search_start <= local_i < local_search_end and local_i + pam_len + guide_len <= region_len:
                    if _mt_pam and local_i + 2 < region_len and region_str[local_i + 2].upper() == 'C':
                        p += 1
                        continue
                    guide = _rc(region_str[local_i + pam_len: local_i + pam_len + guide_len])
                    if guide not in seen and not _contains_rs(guide, restriction_sites):
                        ctx_plus = region_str[local_i - 3: local_i + pam_len + guide_len + 4]
                        if len(ctx_plus) >= 30:
                            pam_rc = region_str[local_i: local_i + pam_len]
                            genomic_i = local_i + search_start - buffer
                            all_sgrnas.append({'seq': guide, 'pam': _rc(pam_rc), 'strand': 'reverse',
                                               'score': _score(_rc(ctx_plus)[:30]),
                                               'cut_site': _cut(genomic_i, "reverse", pam_len)})
                            seen.add(guide)
                p += 1
        else:
            pam_fwd_pattern = pam.upper()
            pam_rev_pattern = _rc(pam).upper()
            for i in range(search_start, search_end):
                local_i = i - region_start
                if local_i + guide_len + pam_len <= region_len:
                    pam_fwd = region_str[local_i + guide_len: local_i + guide_len + pam_len]
                    if self.match_pam(pam_fwd, pam_fwd_pattern):
                        guide = region_str[local_i: local_i + guide_len]
                        if guide not in seen and not _contains_rs(guide, restriction_sites):
                            ctx = region_str[local_i - 4: local_i + guide_len + pam_len + 3]
                            if len(ctx) >= 30:
                                all_sgrnas.append({'seq': guide, 'pam': pam_fwd, 'strand': 'forward',
                                                   'score': _score(ctx[len(ctx) - 30:]),
                                                   'cut_site': _cut(i, "forward", pam_len)})
                                seen.add(guide)
                if local_i + pam_len + guide_len <= region_len:
                    pam_rev_comp = region_str[local_i: local_i + pam_len]
                    if self.match_pam(pam_rev_comp, pam_rev_pattern):
                        guide = _rc(region_str[local_i + pam_len: local_i + pam_len + guide_len])
                        if guide not in seen and not _contains_rs(guide, restriction_sites):
                            ctx_plus = region_str[local_i - 3: local_i + pam_len + guide_len + 4]
                            if len(ctx_plus) >= 30:
                                all_sgrnas.append({'seq': guide, 'pam': _rc(pam_rev_comp), 'strand': 'reverse',
                                                   'score': _score(_rc(ctx_plus)[:30]),
                                                   'cut_site': _cut(i, "reverse", pam_len)})
                                seen.add(guide)
        return all_sgrnas


# ============================================================
# CRISPR 设计核心类
# ============================================================

class CRISPRDesigner:
    def __init__(self, config: DesignConfig, genome_seqs: Dict[str, Seq],
                 genome_lens: Dict[str, int], genome_type: str):
        self.config = config
        self.sgrna_designer = SGRNADesigner(config)
        self.genome_seqs = genome_seqs
        self.genome_lens = genome_lens
        self.genome_type = genome_type
        self.search_genomes: Dict[str, str] = {}
        for seqid, seq in self.genome_seqs.items():
            s = str(seq)
            self.search_genomes[seqid] = s + s if self.genome_type == 'circle' else s
        self.template_fixed_length = self._calculate_template_fixed_length()
        self.arm_lengths_for_assembly = self._build_arm_length_candidates()
        self.template_segments = SequenceUtils.parse_template_segments(
            self.config.synthesis_template
        ) if self.config.synthesis_template else []

    def _calculate_template_fixed_length(self) -> int:
        if not self.config.synthesis_template:
            return 0
        template = self.config.synthesis_template
        placeholder_lengths = {
            "{sgRNA_fwd}": self.config.guide_len,
            "{sgRNA_rc}": self.config.guide_len,
            "{pam}": self.config.pam_len,
            "{pam_rc}": self.config.pam_len,
            "{barcode}": self.config.barcode_len,
            "{exempt_restriction_site}": len(self.config.cloning_site) if self.config.cloning_site else 0,
            "{insert}": 0,
            "{upstream_arm}": 0,
            "{downstream_arm}": 0,
        }
        fixed_len = len(template)
        for placeholder, replacement_len in placeholder_lengths.items():
            count = template.count(placeholder)
            if count:
                fixed_len -= count * len(placeholder)
                fixed_len += count * replacement_len
        return fixed_len

    def _build_arm_length_candidates(self) -> List[int]:
        base_arms = self.config.arm_search_order
        min_arm = min(base_arms) if base_arms else 0
        if not self.config.max_oligo_length:
            return base_arms
        total_arm_budget = self.config.max_oligo_length - self.template_fixed_length
        exact_arm = total_arm_budget // 2
        if total_arm_budget % 2 != 0:
            logging.info(f"oligo余量为奇数({total_arm_budget})，同源臂取 {exact_arm} bp")
        if exact_arm < min_arm:
            logging.error(f"max_oligo_length={self.config.max_oligo_length} 过小；同源臂 {exact_arm} bp < 最小值 {min_arm} bp")
            sys.exit(1)
        logging.info(f"精确同源臂长度: {exact_arm} bp（固定部分 {self.template_fixed_length} bp + 2×{exact_arm} bp = {self.template_fixed_length + 2 * exact_arm} bp）")
        return [exact_arm]

    def _get_seq_context(self, seqid: str) -> Tuple[int, str]:
        if seqid not in self.genome_lens or seqid not in self.search_genomes:
            raise ValueError(f"seqid '{seqid}' 不在已加载的序列中")
        return self.genome_lens[seqid], self.search_genomes[seqid]

    def _assemble_final_oligo(self, sgrna_fwd: str, pam: str, upstream_arm: str,
                               downstream_arm: str, barcode: str) -> str:
        sgrna_rc = SequenceUtils.get_reverse_complement(sgrna_fwd)
        pam_rc = SequenceUtils.get_reverse_complement(pam) if pam else ""
        if self.config.synthesis_template:
            template = self.config.synthesis_template
            replacements = {
                "{sgRNA_fwd}": sgrna_fwd, "{sgRNA_rc}": sgrna_rc,
                "{pam}": pam, "{pam_rc}": pam_rc,
                "{upstream_arm}": upstream_arm, "{downstream_arm}": downstream_arm,
                "{barcode}": barcode,
                "{exempt_restriction_site}": self.config.cloning_site or ""
            }
            template = template.replace("{insert}", "")
            for placeholder, value in replacements.items():
                template = template.replace(placeholder, value)
            return template
        else:
            raise ValueError("此模式需要合成模板，但未提供。")

    def _compute_valid_del_range(self, cut_site: int, cds_min: int, cds_max: int,
                                  min_del: int, max_del: int) -> Tuple[int, int]:
        if cut_site < cds_min or cut_site > cds_max:
            return (min_del, min_del - 1)
        max_left = cut_site - cds_min
        max_right = cds_max - cut_site
        phys_max_del = max_left + max_right
        effective_max = min(max_del, phys_max_del)
        return (min_del, effective_max)

    def _check_constraints_fast(self, del_start: int, del_end: int, cds_5prime: int,
                                 strand: int, protected_zones: List[Tuple[int, int]]) -> bool:
        if protected_zones:
            for z_start, z_end in protected_zones:
                if max(del_start, z_start) < min(del_end, z_end):
                    return False
        return True

    def _rank_candidate(self, del_start: int, del_end: int, del_len: int,
                        cds_5prime: int, strand: int) -> Tuple:
        """legacy 模式排序键：起始密码子 > 移码 > 5'端位置（小值优先）"""
        if strand == 1:
            is_start_deleted = (del_start <= cds_5prime < del_end)
            position_key = del_start
        else:
            is_start_deleted = (del_start < cds_5prime + 3 <= del_end)
            position_key = -del_end
        is_frameshift = del_len % 3 != 0
        return (not is_start_deleted, not is_frameshift, position_key)

    def _rank_candidate_cut_window(self, del_start: int, del_end: int, del_len: int,
                                    cds_5prime: int, strand: int,
                                    target_window_size: int) -> Tuple:
        """cut_window 模式排序键：起始密码子 > 移码 > 窗口完整度（小值优先）"""
        if strand == 1:
            is_start_deleted = (del_start <= cds_5prime < del_end)
        else:
            is_start_deleted = (del_start < cds_5prime + 3 <= del_end)
        is_frameshift = del_len % 3 != 0
        window_deviation = abs(del_len - target_window_size)
        return (not is_start_deleted, not is_frameshift, window_deviation)

    def _generate_cut_window_candidates(self, sgrna: Dict, cds_min: int, cds_max: int,
                                         strand: int, cds_5prime: int,
                                         protected_zones: List[Tuple[int, int]],
                                         top_k: int = 5,
                                         downstream_override: Optional[int] = None) -> List[Dict]:
        """
        V10 cut_window 模式：穷举窗口内所有有效删除长度，返回评分最优的 top_k 个候选。

        Mt 策略（V9逻辑）：链感知窗口方向
          - 正链：窗口 [cut-upstream, cut+downstream]，删除区间左对齐
          - 反向链：窗口 [cut-downstream, cut+upstream]（翻转），删除区间右对齐
        normal 策略（V8逻辑）：始终左对齐窗口

        downstream_override: 若指定，覆盖 config.cut_window_downstream（用于自适应窗口扩展）。
        """
        cut_site = sgrna['cut_site']
        upstream = self.config.cut_window_upstream
        downstream = downstream_override if downstream_override is not None else self.config.cut_window_downstream
        max_del_len = upstream + downstream

        if self.config.deletion_strategy == "Mt":
            # V9 链感知：PAM 方向决定窗口朝向
            pam_dir = 1 if sgrna.get('strand') == 'forward' else -1
            upstream_pos = cut_site - pam_dir * upstream
            downstream_pos = cut_site + pam_dir * downstream
            ds = max(min(upstream_pos, downstream_pos), cds_min)
            de = min(max(upstream_pos, downstream_pos), cds_max)
        else:
            # normal：固定左对齐
            pam_dir = 1
            ds = max(cut_site - upstream, cds_min)
            de = min(cut_site + downstream, cds_max)

        # cut_site 必须在窗口内（保证编辑发生在cut_site处）
        if not (ds <= cut_site < de):
            return []

        max_available = min(de - ds, max_del_len)
        if max_available <= 0:
            return []

        target_window_size = max_del_len
        ranked = []

        # V10：穷举所有有效长度（1 到 max_available）
        for dl in range(1, max_available + 1):
            if self.config.deletion_strategy == "Mt" and pam_dir == -1:
                # 反向链：右对齐
                del_end = de
                del_start = del_end - dl
            else:
                # 正链或 normal：左对齐
                del_start = ds
                del_end = del_start + dl

            if del_start < cds_min or del_end > cds_max:
                continue
            if not (del_start <= cut_site < del_end):
                continue
            if not self._check_constraints_fast(del_start, del_end, cds_5prime, strand, protected_zones):
                continue

            rank = self._rank_candidate_cut_window(del_start, del_end, dl, cds_5prime, strand, target_window_size)
            ranked.append((rank, {'del_start': del_start, 'del_end': del_end, 'del_len': dl}))

        ranked.sort(key=lambda x: x[0])
        return [r[1] for r in ranked[:top_k]]

    def _generate_legacy_candidates(self, sgrna: Dict, cds_min: int, cds_max: int,
                                     min_del: int, max_del: int, strand: int, cds_5prime: int,
                                     protected_zones: List[Tuple[int, int]],
                                     top_k: int = 5) -> List[Dict]:
        """
        V10 legacy 模式：穷举有效范围内所有删除长度，以 cut_site 为中心居中定位，
        返回评分最优的 top_k 个候选。替代 V8 的 ~20点分层采样。
        """
        cut_site = sgrna['cut_site']
        min_eff, max_eff = self._compute_valid_del_range(cut_site, cds_min, cds_max, min_del, max_del)
        if max_eff < min_eff:
            return []

        ranked = []
        for del_len in range(min_eff, max_eff + 1):
            # V6 兼容：以 cut_site 为中心居中定位
            del_half = del_len >> 1
            del_start = cut_site - del_half
            del_end = del_start + del_len

            if del_start < cds_min or del_end > cds_max:
                continue
            if not self._check_constraints_fast(del_start, del_end, cds_5prime, strand, protected_zones):
                continue

            rank = self._rank_candidate(del_start, del_end, del_len, cds_5prime, strand)
            ranked.append((rank, {'del_start': del_start, 'del_end': del_end, 'del_len': del_len}))

        ranked.sort(key=lambda x: x[0])
        return [r[1] for r in ranked[:top_k]]

    def _generate_candidate_list(self, sgrna: Dict, cds_min: int, cds_max: int,
                                  min_del: int, max_del: int, strand: int, cds_5prime: int,
                                  protected_zones: List[Tuple[int, int]],
                                  top_k: int = 5) -> List[Dict]:
        """调度：cut_window 或 legacy 模式，返回 top_k 候选列表"""
        if self.config.deletion_mode == 'cut_window':
            return self._generate_cut_window_candidates(sgrna, cds_min, cds_max, strand, cds_5prime,
                                                         protected_zones, top_k)
        else:
            return self._generate_legacy_candidates(sgrna, cds_min, cds_max, min_del, max_del,
                                                     strand, cds_5prime, protected_zones, top_k)

    def design_knockout_for_gene(self, gene: Gene, gene_index: int, all_genes: List[Gene],
                                  pam: str, used_barcodes: Set[str], barcode_lock: threading.Lock,
                                  restriction_sites: Optional[List[str]]) -> List[KnockoutDesignResult]:
        """
        单基因敲除设计（V10）。
        关键改进：全局候选池 + 多样性感知选取（min_design_spacing），梯度回退保证覆盖率。
        """
        base_genome_len, search_genome = self._get_seq_context(gene.seqid)

        # ---- 构建邻近基因保护区 ----
        protected_zones = []
        if self.config.strict:
            neighbor_threshold = self.config.neighbor_distance_threshold
            for offset in [-1, 1]:
                ni = gene_index + offset
                if 0 <= ni < len(all_genes):
                    neighbor = all_genes[ni]
                    dist = (gene.start - neighbor.end) if offset == -1 else (neighbor.start - gene.end)
                    if neighbor.seqid == gene.seqid and dist < neighbor_threshold:
                        protected_zones.append((neighbor.cds_min_coord, neighbor.cds_max_coord))
                        if neighbor.strand == 1:
                            p_start = neighbor.cds_5prime_start - self.config.max_promoter_size
                            p_end = neighbor.cds_5prime_start - self.config.min_promoter_size
                        else:
                            p_start = neighbor.cds_5prime_start + self.config.min_promoter_size
                            p_end = neighbor.cds_5prime_start + self.config.max_promoter_size
                        protected_zones.append((p_start, p_end))

        cds_len = gene.cds_max_coord - gene.cds_min_coord
        search_start = gene.cds_min_coord + int(cds_len * self.config.ko_search_start_pct)
        search_end = gene.cds_min_coord + int(cds_len * self.config.ko_search_end_pct)
        if search_start >= search_end:
            return []

        all_sgrnas = self.sgrna_designer.find_sgrnas_in_region(
            base_genome_len, search_genome, search_start, search_end, pam, restriction_sites, self.genome_type
        )
        if not all_sgrnas:
            return []

        # ---- 双因素删除长度约束 ----
        pct_min_bp = max(1, int(cds_len * self.config.del_pct_min))
        pct_max_bp = int(cds_len * self.config.del_pct_max)
        if self.config.del_has_bp and self.config.del_has_pct:
            min_del = max(self.config.del_bp_min, pct_min_bp)
            max_del = min(self.config.del_bp_max, pct_max_bp)
            if max_del < min_del and self.config.del_bp_min > pct_max_bp:
                min_del = pct_min_bp
        elif self.config.del_has_bp:
            min_del = self.config.del_bp_min
            max_del = self.config.del_bp_max
        else:
            min_del = pct_min_bp
            max_del = pct_max_bp
        if max_del < min_del:
            return []

        strand = gene.strand
        cds_5prime = gene.cds_5prime_start
        cds_min, cds_max = gene.cds_min_coord, gene.cds_max_coord

        # ---- 阶段一：构建全局候选池 ----
        # 每个 sgRNA 贡献 top_k=5 个最优候选（穷举非3倍数长度后排序截取）
        all_candidates = []
        for sgrna in all_sgrnas:
            cands = self._generate_candidate_list(
                sgrna, cds_min, cds_max, min_del, max_del, strand, cds_5prime, protected_zones, top_k=5
            )
            for cand in cands:
                cand['sgrna'] = sgrna
                # 计算全局排序键
                del_s, del_e, del_l = cand['del_start'], cand['del_end'], cand['del_len']
                if self.config.deletion_mode == 'cut_window':
                    if strand == 1:
                        is_start_del = (del_s <= cds_5prime < del_e)
                        cut_dist = sgrna['cut_site'] - cds_5prime
                    else:
                        is_start_del = (del_s < cds_5prime + 3 <= del_e)
                        cut_dist = cds_5prime - sgrna['cut_site']
                    is_frameshift = del_l % 3 != 0
                    window_dev = abs(del_l - (self.config.cut_window_upstream + self.config.cut_window_downstream))
                    # V11 修复：补充 not is_frameshift 优先级，确保 frameshift 优先于 in-frame 选取
                    cand['sort_key'] = (not is_start_del, not is_frameshift, cut_dist, -sgrna['score'], window_dev)
                else:
                    cand['sort_key'] = self._rank_candidate(del_s, del_e, del_l, cds_5prime, strand) + (-sgrna['score'],)
                all_candidates.append(cand)

        if not all_candidates:
            logging.warning(f"基因 {gene.id}: 找到sgRNA，但未能生成有效删除候选。")
            return []

        all_candidates.sort(key=lambda x: x['sort_key'])

        # ---- 阶段二：多样性感知选取 + 同源臂/oligo构建 ----
        # V11 核心改进：根据 CDS 长度动态调整设计间距。
        # 长基因保留标准间距，短基因直接切换到 0bp 以优先保证双设计覆盖率。
        # 同时保留 threshold=0 时按最大切割位点间距重排，以及 frameshift 优先池策略。
        min_spacing = compute_dynamic_spacing(cds_len, self.config.min_design_spacing)
        if min_spacing > 0:
            spacing_thresholds = [min_spacing, max(1, min_spacing // 2), 0]
        else:
            spacing_thresholds = [0]

        found_designs: List[KnockoutDesignResult] = []

        def _build_one_design(cand_list, existing_cut_sites, sp_thr):
            """从 cand_list 中找第一个可构建完整设计的候选，返回 KnockoutDesignResult 或 None。"""
            for cand in cand_list:
                cut_site = cand['sgrna']['cut_site']
                if existing_cut_sites and sp_thr > 0:
                    if any(abs(cut_site - cs) < sp_thr for cs in existing_cut_sites):
                        continue

                del_start, del_end = cand['del_start'], cand['del_end']
                upstream_arm = downstream_arm = None
                arm_len = 0
                for try_arm in self.arm_lengths_for_assembly:
                    up = SequenceUtils.get_sequence(base_genome_len, search_genome,
                                                     del_start - try_arm, del_start, self.genome_type)
                    dn = SequenceUtils.get_sequence(base_genome_len, search_genome,
                                                     del_end, del_end + try_arm, self.genome_type)
                    if not (up and dn and len(up) == try_arm and len(dn) == try_arm):
                        continue
                    if (SequenceUtils.contains_restriction_sites(up, restriction_sites) or
                            SequenceUtils.contains_restriction_sites(dn, restriction_sites)):
                        continue
                    upstream_arm, downstream_arm, arm_len = up, dn, try_arm
                    break
                if not upstream_arm:
                    continue

                try:
                    barcode = SequenceUtils.generate_unique_barcode_atomic(
                        self.config.barcode_len, used_barcodes, barcode_lock,
                        restriction_sites,
                        gc_min=self.config.barcode_gc_min,
                        gc_max=self.config.barcode_gc_max,
                        repeat_threshold=self.config.barcode_repeat_threshold
                    )
                except RuntimeError:
                    continue

                sgrna = cand['sgrna']
                final_oligo = self._assemble_final_oligo(sgrna['seq'], sgrna['pam'],
                                                          upstream_arm, downstream_arm, barcode)
                sgrna_rc_seq = SequenceUtils.get_reverse_complement(sgrna['seq'])
                pam_rc_seq = SequenceUtils.get_reverse_complement(sgrna['pam']) if sgrna['pam'] else ''
                variable_lengths = {
                    '{sgRNA_fwd}': len(sgrna['seq']),
                    '{sgRNA_rc}': len(sgrna_rc_seq),
                    '{pam}': len(sgrna['pam']),
                    '{pam_rc}': len(pam_rc_seq),
                    '{upstream_arm}': len(upstream_arm),
                    '{downstream_arm}': len(downstream_arm),
                    '{barcode}': len(barcode),
                    '{exempt_restriction_site}': len(self.config.cloning_site) if self.config.cloning_site else 0,
                    '{insert}': 0,
                }
                if SequenceUtils.check_junctions_for_restriction_sites(
                        final_oligo, self.template_segments, variable_lengths, restriction_sites):
                    with barcode_lock:
                        used_barcodes.discard(barcode)
                    continue

                design_index = len(existing_cut_sites) + 1
                spacing_val = abs(cut_site - existing_cut_sites[0]) if existing_cut_sites else None
                strategy = "frameshift" if cand['del_len'] % 3 != 0 else "in_frame_del"
                return KnockoutDesignResult(
                    gene_id=gene.id,
                    deletion_start=del_start, deletion_end=del_end,
                    deletion_length=cand['del_len'], strategy=strategy,
                    arm_length=arm_len, upstream_arm=upstream_arm, downstream_arm=downstream_arm,
                    sgrna_seq=sgrna['seq'], sgrna_pam=sgrna['pam'],
                    sgrna_strand=sgrna['strand'], sgrna_score=sgrna['score'],
                    sgrna_cut_site=cut_site, barcode=barcode,
                    final_oligo_for_synthesis=final_oligo,
                    design_index=design_index, cut_site_spacing=spacing_val
                )
            return None

        def _try_select_from_pool(pool: List) -> List[KnockoutDesignResult]:
            """
            在给定候选池内按间距梯度选取 sgrna_num 套设计。
            返回已选设计列表（< sgrna_num 表示候选池不足，调用方负责回收 barcode）。
            """
            for spacing_threshold in spacing_thresholds:
                local_found: List[KnockoutDesignResult] = []
                local_cuts: List[int] = []

                for _ in range(self.config.sgrna_num):
                    # threshold=0 且已有设计时，重排以最大化切割位点间距多样性
                    if spacing_threshold == 0 and local_cuts:
                        ex = local_cuts[:]
                        cands_this_round = sorted(
                            pool,
                            key=lambda c, _ex=ex: (
                                -min(abs(c['sgrna']['cut_site'] - cs) for cs in _ex),
                                c['sort_key']
                            )
                        )
                    else:
                        cands_this_round = pool

                    result = _build_one_design(cands_this_round, local_cuts, spacing_threshold)
                    if result:
                        local_found.append(result)
                        local_cuts.append(result.sgrna_cut_site)
                    else:
                        break

                if len(local_found) >= self.config.sgrna_num:
                    return local_found

                # 本轮不足 → 回收 barcode，降低间距重试
                if spacing_threshold > 0:
                    with barcode_lock:
                        for r in local_found:
                            used_barcodes.discard(r.barcode)
                    local_found = []
                    local_cuts = []

            # 全部梯度耗尽，返回最终（可能不足）的结果
            return local_found

        # ---- Frameshift 优先策略 ----
        # 先尝试 frameshift-only 候选池；若不足 sgrna_num 套，回收 barcode 后用全部候选重试
        frameshift_pool = [c for c in all_candidates if c['del_len'] % 3 != 0]

        if frameshift_pool:
            found_designs = _try_select_from_pool(frameshift_pool)
            if len(found_designs) < self.config.sgrna_num:
                # frameshift-only 不足，回收 barcode 并用全部候选重试
                with barcode_lock:
                    for r in found_designs:
                        used_barcodes.discard(r.barcode)
                found_designs = _try_select_from_pool(all_candidates)
        else:
            # 无 frameshift 候选（极端情况）：直接用全部候选
            found_designs = _try_select_from_pool(all_candidates)

        # ---- Mt cut_window 自适应窗口扩展 ----
        # 在 Mt+cut_window 模式下，若标准窗口设计不足，步长10bp逐步扩展下游窗口直至完成双套设计
        standard_max_del = self.config.cut_window_upstream + self.config.cut_window_downstream
        actual_dn_used = self.config.cut_window_downstream
        if (len(found_designs) < self.config.sgrna_num and
                self.config.deletion_mode == 'cut_window' and
                self.config.deletion_strategy == 'Mt'):
            relaxed_dn = self.config.cut_window_downstream + 10
            max_relaxed_dn = self.config.cut_window_downstream + 300  # 上限扩展 300bp
            while len(found_designs) < self.config.sgrna_num and relaxed_dn <= max_relaxed_dn:
                # 回收本轮 barcode，从头尝试
                with barcode_lock:
                    for r in found_designs:
                        used_barcodes.discard(r.barcode)
                found_designs = []
                # 用扩展下游窗口重建候选池
                expanded_candidates = []
                for sgrna in all_sgrnas:
                    cands = self._generate_cut_window_candidates(
                        sgrna, cds_min, cds_max, strand, cds_5prime, protected_zones,
                        top_k=5, downstream_override=relaxed_dn
                    )
                    for cand in cands:
                        cand['sgrna'] = sgrna
                        del_s, del_e, del_l = cand['del_start'], cand['del_end'], cand['del_len']
                        if strand == 1:
                            is_start_del = (del_s <= cds_5prime < del_e)
                            cut_dist = sgrna['cut_site'] - cds_5prime
                        else:
                            is_start_del = (del_s < cds_5prime + 3 <= del_e)
                            cut_dist = cds_5prime - sgrna['cut_site']
                        is_frameshift = del_l % 3 != 0
                        window_dev = abs(del_l - (self.config.cut_window_upstream + relaxed_dn))
                        cand['sort_key'] = (not is_start_del, not is_frameshift, cut_dist,
                                            -sgrna['score'], window_dev)
                        expanded_candidates.append(cand)
                if expanded_candidates:
                    expanded_candidates.sort(key=lambda x: x['sort_key'])
                    fs_pool = [c for c in expanded_candidates if c['del_len'] % 3 != 0]
                    if fs_pool:
                        found_designs = _try_select_from_pool(fs_pool)
                        if len(found_designs) < self.config.sgrna_num:
                            with barcode_lock:
                                for r in found_designs:
                                    used_barcodes.discard(r.barcode)
                            found_designs = _try_select_from_pool(expanded_candidates)
                    else:
                        found_designs = _try_select_from_pool(expanded_candidates)
                if len(found_designs) >= self.config.sgrna_num:
                    actual_dn_used = relaxed_dn
                    logging.info(f"基因 {gene.id}: cut_window 下游扩展至 {relaxed_dn}bp "
                                 f"（{self.config.cut_window_upstream}:{relaxed_dn}="
                                 f"{self.config.cut_window_upstream + relaxed_dn}bp）后完成设计")
                    break
                relaxed_dn += 10

        # 为使用扩展窗口且删除长度超过标准上限的设计添加警告
        for design in found_designs:
            if design.deletion_length > standard_max_del:
                design.deletion_warning = (
                    f"扩展窗口：删除长度 {design.deletion_length}bp 超过标准 "
                    f"{standard_max_del}bp"
                    f"（窗口 {self.config.cut_window_upstream}:{actual_dn_used}"
                    f"={self.config.cut_window_upstream + actual_dn_used}bp）"
                )

        if found_designs and len(found_designs) < self.config.sgrna_num:
            actual_spacing = found_designs[-1].cut_site_spacing if len(found_designs) > 1 else None
            if actual_spacing is not None:
                logging.info(f"基因 {gene.id}: 梯度回退后2套设计实际间距={actual_spacing}bp（目标={min_spacing}bp）")
            else:
                logging.info(f"基因 {gene.id}: 仅找到1套设计（候选池已耗尽）")

        return found_designs


# ============================================================
# 结果处理类
# ============================================================

class ResultProcessor:
    @staticmethod
    def save_results(designs: List, failed_gene_data: List[Dict], output_file: str):
        design_type_info = "Knockout_Cas9_V11"
        if designs:
            df = pd.DataFrame([d.__dict__ for d in designs])
            df['Status'] = "Success"
            # V11：保留 design_index 和 cut_site_spacing 列。
            cols_in_order = [
                'gene_id', 'Status', 'design_index', 'cut_site_spacing',
                'strategy', 'deletion_warning', 'sgrna_seq', 'sgrna_score', 'barcode',
                'final_oligo_for_synthesis', 'arm_length', 'upstream_arm', 'downstream_arm',
                'sgrna_pam', 'sgrna_strand', 'sgrna_cut_site',
                'deletion_length', 'deletion_start', 'deletion_end'
            ]
            df = df.reindex(columns=cols_in_order).rename(columns=lambda c: c.replace('_', ' ').title())
            df.to_csv(output_file, index=False)
            logging.info(f"{design_type_info} 结果已保存至: {output_file}")
        else:
            logging.warning(f"没有成功的 {design_type_info} 设计可供保存。")

        if failed_gene_data:
            failed_df = pd.DataFrame(failed_gene_data)
            failed_output_path = os.path.splitext(output_file)[0] + "_failed.csv"
            failed_df.to_csv(failed_output_path, index=False)
            logging.info(f"失败基因信息已保存至: {failed_output_path}")

        logging.info(f"总结: 生成 {len(designs) if designs else 0} 个成功设计, "
                     f"{len(failed_gene_data)} 个基因设计失败。")


# ============================================================
# 主流程
# ============================================================

def run_knockout_pipeline(args, config, genes, genome_processor, num_workers=8):
    """
    V11：多线程并行 + 原子化 barcode 管理。
    关键修复：所有线程共享 used_barcodes + barcode_lock，传入 design_knockout_for_gene，
    消除 V8 中每线程独立 barcode 集合导致的跨基因重复风险。
    """
    designer = CRISPRDesigner(config, genome_processor.genome_seqs,
                               genome_processor.genome_lens, args.genome_type)
    successful_designs: List[KnockoutDesignResult] = []
    used_barcodes: Set[str] = set()
    barcode_lock = threading.Lock()
    gene_design_counts: Dict[str, int] = {}
    t_start = time.time()

    def design_single_gene(idx_gene_tuple):
        i, gene = idx_gene_tuple
        if gene.cds_5prime_start is None:
            return (gene.id, [])
        # V11：传入共享 used_barcodes 和 barcode_lock（原子化barcode生成）
        results = designer.design_knockout_for_gene(
            gene, i, genes, args.pam, used_barcodes, barcode_lock, args.restriction_site
        )
        return (gene.id, results)

    logging.info(f"启用多线程加速 (num_workers={num_workers})")

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures_to_gene = {
            executor.submit(design_single_gene, (i, gene)): gene.id
            for i, gene in enumerate(genes)
        }
        completed_count = 0
        for future in futures_to_gene:
            gene_id, results = future.result()
            gene_design_counts[gene_id] = len(results)
            if results:
                successful_designs.extend(results)
            completed_count += 1
            if completed_count % 500 == 0:
                elapsed = time.time() - t_start
                logging.info(f"正在处理基因 {completed_count}/{len(genes)}... 已用时 {elapsed:.1f}s")

    elapsed_total = time.time() - t_start
    logging.info(f"全部基因设计完成，总耗时: {elapsed_total:.1f}s")

    # 统计两套设计的多样性指标
    genes_with_2 = sum(1 for c in gene_design_counts.values() if c >= 2)
    spacings = [r.cut_site_spacing for r in successful_designs if r.cut_site_spacing is not None]
    if spacings:
        avg_sp = sum(spacings) / len(spacings)
        min_sp = min(spacings)
        logging.info(f"设计多样性统计: 双套设计基因数={genes_with_2}, "
                     f"切割位点间距 平均={avg_sp:.1f}bp 最小={min_sp}bp")

    # 失败/不完整基因处理
    failed_genes = [g for g in genes if gene_design_counts.get(g.id, 0) == 0]
    partial_genes = [g for g in genes if 0 < gene_design_counts.get(g.id, 0) < config.sgrna_num]
    failed_gene_data = []
    for g in failed_genes + partial_genes:
        base_genome_len, search_genome = designer._get_seq_context(g.seqid)
        cds_seq = SequenceUtils.get_sequence(base_genome_len, search_genome,
                                              g.cds_min_coord, g.cds_max_coord, designer.genome_type)
        status = "Failed" if gene_design_counts.get(g.id, 0) == 0 else "Partial"
        failed_gene_data.append({
            "Gene_ID": g.id, "Status": status,
            "Designed_Count": gene_design_counts.get(g.id, 0),
            "Required_Count": config.sgrna_num,
            "SeqID": g.seqid,
            "Strand": "+" if g.strand == 1 else "-",
            "Relevant_Sequence": cds_seq
        })

    if partial_genes:
        logging.warning(f"有 {len(partial_genes)} 个基因未达到目标设计数 {config.sgrna_num}（已标记为Partial）。")

    ResultProcessor.save_results(successful_designs, failed_gene_data, args.output)


def main():
    parser = argparse.ArgumentParser(
        description="基因敲除文库设计工具 (Cas9) V11 - 动态间距自适应版",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--input_fna", help="输入的基因组FASTA文件（与--input_gff配对）。")
    parser.add_argument("--input_gff", help="输入的GFF3注释文件（与--input_fna配对）。")
    parser.add_argument("--input_gbff", help="输入的GBFF文件（同时包含序列和注释）。")
    parser.add_argument("--output", required=True, help="输出的CSV文件路径。")
    parser.add_argument("--genome_type", type=str, default="linear", choices=['linear', 'circle'])
    parser.add_argument("--species", type=str, default="M_thermophila",
                        choices=['M_thermophila', 'E_coli'])
    parser.add_argument("--sgRNA_num", type=int, default=2, help="每个基因生成的设计方案数量（默认2）。")
    parser.add_argument("--barcode_len", type=int, default=10)
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--promoter_region", type=str, default="200:500")
    parser.add_argument("--HR_len", type=str, default="65")
    parser.add_argument("--del_length_per", type=str, default=None,
                        help="按CDS比例限制删除长度，如 '20%%:80%%'")
    parser.add_argument("--del_length_bp", type=str, default=None,
                        help="按绝对碱基数限制删除长度，如 '500:2000'")
    parser.add_argument("--restriction_site", type=str, nargs='+')
    parser.add_argument("--synthesis_template", type=str, required=True)
    parser.add_argument("--cloning_site", type=str)
    parser.add_argument("--pam", type=str, default="NGG")
    parser.add_argument("--ko_search_range", type=str, default="5:80")
    parser.add_argument("--max_oligo_length", type=int)
    parser.add_argument("--num_workers", type=int, default=8)

    # V8 兼容参数
    parser.add_argument("--deletion_mode", type=str, default="auto",
                        choices=['auto', 'cut_window', 'legacy_length'],
                        help="删除区间模式（auto=物种感知，Mt物种默认cut_window）")
    parser.add_argument("--cut_window", type=str, default="20:100",
                        help="cut_window模式窗口大小'上游:下游'，默认'20:100'")
    parser.add_argument("--no_mt_pam", action="store_true",
                        help="关闭Mt最优PAM过滤（默认M_thermophila自动启用AGG/TGG/CGG过滤）")

        # V10/V11 参数
    parser.add_argument("--dele_model", type=str, default="normal", choices=["normal", "Mt"],
                        help="删除策略（V10内联替代V9猴子补丁）：\n"
                             "  normal：V8原始行为（cut_site均使用-3偏移，窗口左对齐）\n"
                             "  Mt：V9 M.thermophila优化行为（反链cut_site+3校正，链感知窗口方向）\n"
                             "  注：M_thermophila物种自动选择Mt（可被此参数覆盖）")
    parser.add_argument("--min_design_spacing", type=int, default=100,
                        help="V11 动态间距的基准值（bp，默认100）。\n"
                            "  CDS 长度 < 2×该值时自动切换为 0bp，以保证短基因覆盖率。\n"
                            "  长基因若无法满足，自动梯度回退至一半间距 → 0bp。")

    args = parser.parse_args()

    try:
        if SPECIES_CONFIG_AVAILABLE:
            try:
                species_cfg = SpeciesConfig(args.species)
                min_p, max_p = ConfigParser.parse_range_param(args.promoter_region, "启动子保护区域")
                preferred_hr, min_hr = ConfigParser.parse_range_param(args.HR_len, "同源臂最短长度", sort_values=False)
            except Exception as e:
                logging.warning(f"加载物种配置失败，使用命令行参数: {e}")
                min_p, max_p = ConfigParser.parse_range_param(args.promoter_region, "启动子保护区域")
                preferred_hr, min_hr = ConfigParser.parse_range_param(args.HR_len, "同源臂最短长度", sort_values=False)
        else:
            min_p, max_p = ConfigParser.parse_range_param(args.promoter_region, "启动子保护区域")
            preferred_hr, min_hr = ConfigParser.parse_range_param(args.HR_len, "同源臂最短长度", sort_values=False)

        # 加载合成模板
        synthesis_template_content = None
        if args.synthesis_template:
            try:
                with open(args.synthesis_template, 'r') as f:
                    synthesis_template_content = f.read().strip()
                logging.info(f"成功加载合成模板: {args.synthesis_template}")
            except FileNotFoundError:
                logging.error(f"找不到合成模板文件: {args.synthesis_template}")
                sys.exit(1)

        arm_search_order = ConfigParser.create_arm_search_order(preferred_hr, min_hr)

        config = DesignConfig(
            arm_search_order=arm_search_order,
            guide_len=20,
            synthesis_template=synthesis_template_content,
            cloning_site=args.cloning_site,
            sgrna_num=args.sgRNA_num,
            barcode_len=args.barcode_len,
            strict=args.strict,
            min_promoter_size=min_p,
            max_promoter_size=max_p,
            species=args.species,
            barcode_gc_min=0.30,
            barcode_gc_max=0.80,
            barcode_repeat_threshold=5,
            neighbor_distance_threshold=1500,
            max_oligo_length=args.max_oligo_length,
            pam_len=len(args.pam),
            min_design_spacing=args.min_design_spacing
        )

        # 解析 deletion_mode
        actual_deletion_mode = args.deletion_mode
        if actual_deletion_mode == 'auto':
            actual_deletion_mode = 'cut_window' if args.species == 'M_thermophila' else 'legacy_length'
        config.deletion_mode = actual_deletion_mode

        if actual_deletion_mode == 'cut_window':
            cw_parts = args.cut_window.strip().split(':')
            if len(cw_parts) != 2:
                parser.error("--cut_window 格式错误，请使用'上游:下游'，例如 '20:100'")
            config.cut_window_upstream = int(cw_parts[0])
            config.cut_window_downstream = int(cw_parts[1])

        # V10：解析 dele_model（策略模式，替代V9猴子补丁）
        # 规则：明确指定 --dele_model 时优先采用；否则 M_thermophila 物种自动选 Mt
        if args.dele_model != "normal":
            config.deletion_strategy = args.dele_model
        else:
            config.deletion_strategy = "Mt" if args.species == "M_thermophila" else "normal"

        # Mt PAM 过滤
        config.use_mt_pam = (args.species == 'M_thermophila' and not args.no_mt_pam)

        # 双因素约束解析
        if actual_deletion_mode == 'legacy_length' and args.del_length_per is None and args.del_length_bp is None:
            parser.error("legacy_length模式必须至少指定 --del_length_per 或 --del_length_bp")
        if args.del_length_per is not None:
            config.del_pct_min, config.del_pct_max = ConfigParser.parse_del_pct_param(args.del_length_per)
            config.del_has_pct = True
        if args.del_length_bp is not None:
            config.del_bp_min, config.del_bp_max = ConfigParser.parse_del_bp_param(args.del_length_bp)
            config.del_has_bp = True

        min_pct, max_pct = ConfigParser.parse_range_param(args.ko_search_range, "敲除搜索范围", sort_values=True)
        config.ko_search_start_pct, config.ko_search_end_pct = min_pct / 100.0, max_pct / 100.0

        if actual_deletion_mode == 'cut_window':
            del_len_desc = (f"cut_window: cut_site-{config.cut_window_upstream}bp 到 "
                            f"cut_site+{config.cut_window_downstream}bp（共{config.cut_window_upstream + config.cut_window_downstream}bp）")
        else:
            del_desc_parts = []
            if config.del_has_pct:
                del_desc_parts.append(f"CDS的 {config.del_pct_min * 100:.0f}%-{config.del_pct_max * 100:.0f}%")
            if config.del_has_bp:
                bp_max_str = str(config.del_bp_max) if config.del_bp_max < 1000000 else "无上限"
                del_desc_parts.append(f"{config.del_bp_min}-{bp_max_str} bp")
            del_len_desc = "交集约束: " + " & ".join(del_desc_parts) if len(del_desc_parts) > 1 else del_desc_parts[0]

        logging.info("配置概览:")
        logging.info(f"  版本: V11（动态间距自适应版）")
        logging.info(f"  物种: {config.species}")
        logging.info(f"  删除策略: {config.deletion_strategy}（原V9 --dele_model）")
        logging.info(f"  删除模式: {actual_deletion_mode}")
        logging.info(f"  删除长度/窗口: {del_len_desc}")
        logging.info(f"  Mt最优PAM过滤: {'启用 (AGG/TGG/CGG)' if config.use_mt_pam else '关闭 (所有NGG)'}")
        if config.min_design_spacing > 0:
            logging.info(
                f"  V11动态间距: 基准 {config.min_design_spacing} bp"
                f"（CDS < {2 * config.min_design_spacing} bp 自动切换为 0bp；"
                f"长基因回退至 {max(1, config.min_design_spacing // 2)}bp → 0bp）"
            )
        else:
            logging.info("  V11动态间距: 已禁用，所有基因均使用 0bp")
        logging.info(f"  每sgRNA候选数: top_k=5（穷举非3倍数长度后截取）")
        logging.info(f"  Barcode GC范围: {config.barcode_gc_min:.0%}-{config.barcode_gc_max:.0%}")
        logging.info(f"  同源臂最短长度: {min(config.arm_search_order)} bp")
        if config.max_oligo_length:
            logging.info(f"  目标oligo长度: {config.max_oligo_length} bp")

        use_gbff_mode = bool(args.input_gbff)
        if use_gbff_mode:
            if args.input_fna or args.input_gff:
                raise ValueError("使用--input_gbff时，不应同时提供--input_fna/--input_gff。")
            genome_processor = GenomeProcessor(gbff_file=args.input_gbff)
        else:
            if not args.input_fna or not args.input_gff:
                raise ValueError("未使用--input_gbff时，必须同时提供--input_fna和--input_gff。")
            genome_processor = GenomeProcessor(genome_file=args.input_fna, gff_file=args.input_gff)

        genome_processor.load_genome()
        genes = genome_processor.parse_genes()
        run_knockout_pipeline(args, config, genes, genome_processor, num_workers=args.num_workers)

    except Exception as e:
        logging.error(f"发生严重错误: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
