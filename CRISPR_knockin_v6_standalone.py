#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CRISPR_knockin_v6_standalone.py — Fully self-contained N_start / C_stop dual-mode knockin library design
=========================================================================================================
Author : Zhaohui Cai with GitHub Copilot Assistant
Version: 6.0-standalone
Date   : 2026-03-19

This is a standalone version of CRISPR_knockin_v6.py that does NOT depend on any
external knockin_J23119RBS_V*.py files. All logic from V36→V37→V38→V39→V40→V41→V42→
V43→V44→V2→V3→V4→V5→V6 is inlined into this single file.

New in V6 (mutation policy + biological filtering refinement):
        * PAM inside the current target CDS: still try synonymous mutation first.
        * PAM inside the neighboring CDS: also try synonymous mutation first.
        * Conservative amino-acid substitution is now a strict fallback only when
            synonymous mutation has no valid solution.
        * Conservative fallback candidates are ranked by minimal codon edit count,
            so 1-nt solutions are preferred over 2-nt solutions when both break the PAM.
        * Rank2 diversity now runs only after biological filtering, so candidates with
            no amino-acid change are exhausted before conservative-fallback designs are used.
        * The "Optimal" label remains distance-based, while amino-acid cost is reported
            separately as [AA_Change:Conservative] instead of being hidden inside quality labels.

New in V5 (C_stop overlap redesign — preserve fusion + neighboring CDS):
    * Fixed the remaining overlap bug in C_stop when the downstream CDS overlaps the
        target stop codon or extends further upstream into the target CDS.
    * Root cause in V4: shared bases were appended to LHA. That preserved the neighbor
        start but also reintroduced the target stop codon before the fusion payload.
    * V5 moves the rescued overlap segment to the LEFT of RHA instead:
        - LHA now stays upstream of the target stop codon, so C-terminal fusion remains open.
        - RHA begins with a rescued downstream-CDS prefix, then continues with true right-arm
            homology starting AFTER the stop codon.
        - overlap>3 now also copies the neighbor CDS bases that lie upstream of the stop codon.
    * Restriction-site fixing and HA trimming now treat this rescued RHA prefix as
        non-homology sequence, so effective RHA length checks stay biologically correct.

New in V4 (C_stop RHA overlap fix — protect neighboring gene CDS):
    * Fixed RHA start position for P2/P3 strategies when adjacent genes share bases.
    * Root cause: hardcoded rha_seq[3:] skipped shared bases when overlap_bases > 0
        (e.g., TGATG pattern where 'A' is both stop codon end and downstream ATG start).
    * overlap>=4: [LargeOverlap:Nbp] tag added to Info column for manual review.

New in V3 (C_stop Rank2 diversity improvement):
  * Option C + Option A diversity algorithm for Rank2 selection in C_stop mode.
  * New CLI parameter: --rank2_sim_max (default 50). Controls same-strategy HA
    similarity upper bound for Rank2. When Rank1 strategy == Rank2 strategy, only
    accept Rank2 if max(SIM_LHA, SIM_RHA) < rank2_sim_max.
  * Priority order for Rank2:
    1. Different strategy class (P1 vs P2/P3, etc.) — Option C.
    2. Same strategy, but max(SIM_LHA, SIM_RHA) < rank2_sim_max — Option A.
    3. Final fallback: least similar valid candidate regardless of threshold.
  * Diversity source annotated in Info column:
    [Diverse:StrategyDiff] / [Diverse:SeqDiff,SIM_max=XX%] / [FallbackDiverse:SIM_max=XX%]

Features (unchanged from v2):
  * N_start mode: identical to V46 — sgRNA/HDR around the start codon.
  * C_stop  mode: sgRNA/HDR around the stop codon (C-terminal fusion).
  * sgRNA genome-wide specificity filter (off-target index).
  * _compute_ha_budget() fix (V45): ha_budget=120.
  * Mutation-aware balanced HA trimming (V42).
  * RHA budget-fill padding (V43 Stage 4.5).
  * Deterministic barcode generation with similarity flagging (V44).
  * BsaI/BbsI strict exclusion (V39).
  * Flexible HA lengths with max-oligo-length budget enforcement (V40).
  * RHA start codon gate (V41).
  * CCN downstream P3B only (V41).
  * Dual design per gene (V38).
  * Tiered mutation logic: Level 1 (synonymous) + Level 2 (conservative) (V36).

Usage:
  python CRISPR_knockin_v2_standalone.py --model N_start [args ...]
  python CRISPR_knockin_v2_standalone.py --model C_stop \\
      --template Cfusion_library_oligo_template.fasta [args ...]
"""

import argparse
import os
import random
import re
import sys
from collections import Counter
from enum import Enum

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


# ===========================================================================
# Constants
# ===========================================================================

class Colors:
    HEADER  = '\033[95m'
    OKBLUE  = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL    = '\033[91m'
    ENDC    = '\033[0m'


# Standard Genetic Code
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

AA_TO_CODONS = {}
for _codon, _aa in CODON_TABLE.items():
    AA_TO_CODONS.setdefault(_aa, []).append(_codon)

# Conservative Amino Acid Groups — biochemical 4-class classification
# Nonpolar hydrophobic: G A V L I P M F W
# Polar uncharged:      S T Y N Q C
# Basic (+):            K R H
# Acidic (-):           D E
_NONPOLAR  = ['G', 'A', 'V', 'L', 'I', 'P', 'M', 'F', 'W']
_POLAR     = ['S', 'T', 'Y', 'N', 'Q', 'C']
_BASIC     = ['K', 'R', 'H']
_ACIDIC    = ['D', 'E']
CONSERVATIVE_GROUPS = {
    aa: [x for x in grp if x != aa]
    for grp in (_NONPOLAR, _POLAR, _BASIC, _ACIDIC)
    for aa in grp
}

# Valid prokaryotic start codons (E. coli genetic code 11)
VALID_START_CODONS = frozenset({'ATG', 'GTG', 'TTG', 'CTG', 'ATT', 'ATC', 'ATA'})

STOP_CODONS = frozenset({'TAA', 'TAG', 'TGA'})

_BALANCE_TOL = 5  # max allowed |eff_lha - eff_rha| imbalance (bp)

# Template backbone region with 2 intentional BsaI Golden Gate sites
TEMPLATE_BACKBONE_BSAI = "GGTGTGAGACCATAGGTCTCTCAAG"

# All forbidden RE recognition sequences (forward strand; RC computed at runtime)
FORBIDDEN_RE_SEQS = ['GGTCTC', 'GAAGAC']


# ===========================================================================
# Standardized Error Codes (V37)
# ===========================================================================

class DesignError(Enum):
    P1_CLEAN_DEL_TOO_LONG          = "P1_Clean_Deletion_Too_Long"
    P1_CLEAN_LHA_RHA_INVALID       = "P1_Clean_LHA_or_RHA_Invalid"
    P1_PARTIAL_SEED_DEL_TOO_LONG   = "P1_Partial_Seed_Deletion_Too_Long"
    P1_PARTIAL_SEED_LHA_RHA_INVALID= "P1_Partial_Seed_LHA_or_RHA_Invalid"
    P1_PARTIAL_PAM_DEL_TOO_LONG    = "P1_Partial_PAM_Deletion_Too_Long"
    P1_PARTIAL_PAM_LHA_RHA_INVALID = "P1_Partial_PAM_LHA_or_RHA_Invalid"
    P2_RHA_MUT_FAILED              = "P2_RHA_Mutation_Failed"
    P2_LHA_MUT_SHIFT_OOR           = "P2_LHA_Mut_Shift_Out_Of_Range"
    P2_LHA_MUT_FAILED              = "P2_LHA_Mutation_Failed"
    P2_LHA_RHA_INVALID             = "P2_LHA_RHA_Invalid_Or_Too_Long"
    P2_REQUIRED_LEN_EXCEEDS        = "P2_Required_Len_Exceeds"
    P2_RHA_GC_OR_TOO_SHORT         = "P2_RHA_GC_Or_Too_Short"
    P3A_MUT_FAILED                 = "P3A_Mutation_Failed"
    P3A_LEN_RHA_INVALID            = "P3A_Length_Or_RHA_GC_Invalid"
    P3B_REQ_LEN_EXCEEDS            = "P3B_Required_Len_Exceeds"
    P3B_MUT_FAILED                 = "P3B_Mutation_Failed"
    P3B_LHA_MISSING                = "P3B_LHA_Missing"
    P3B_RHA_GC_INVALID             = "P3B_RHA_GC_Invalid"
    P3B_MUT_IDX_OOR                = "P3B_Mut_Idx_Out_Of_Range"
    SGRNA_GC_INVALID               = "SgRNA_GC_Invalid"
    TARGET_WINDOW_OOB              = "Target_Window_Out_Of_Bounds"
    NO_PAM_FOUND                   = "No_PAM_Found"
    FILTER_MIN_EFFECTIVE_HA        = "Filter_Min_Effective_HA"
    FILTER_LENGTH_HARD             = "Filter_Length_Hard_Limit"


# ===========================================================================
# Constraint Validator (V37)
# ===========================================================================

class ConstraintValidator:
    def __init__(self, args):
        self.args = args

    def calc_gc(self, seq):
        if not seq:
            return 0.0
        seq = seq.upper()
        return round((seq.count('G') + seq.count('C')) / len(seq) * 100, 2)

    def is_gc_valid(self, seq, is_sgrna=False):
        gc = self.calc_gc(seq)
        if is_sgrna:
            return self.args.min_gc_sgrna <= gc <= self.args.max_gc_sgrna
        return gc >= self.args.min_gc_ha

    def check_recut_risk(self, final_oligo, spacer_seq, pam_seq):
        final_oligo = final_oligo.upper()
        spacer_seq  = spacer_seq.upper()
        pam_seq     = pam_seq.upper()
        target_fwd  = spacer_seq + pam_seq
        target_rc   = str(Seq(target_fwd).reverse_complement())
        if target_fwd in final_oligo:
            return True, "Risk: Intact Fwd Target"
        if target_rc in final_oligo:
            return True, "Risk: Intact Rev Target"
        seed_len = self.args.seed_region
        if len(spacer_seq) >= seed_len:
            seed_seq = spacer_seq[-seed_len:] + pam_seq
            seed_rc  = str(Seq(seed_seq).reverse_complement())
            if seed_seq in final_oligo or seed_rc in final_oligo:
                return True, f"Risk: Intact Seed+PAM ({seed_len}bp)"
        return False, "Safe"


# ===========================================================================
# Helper functions (V44)
# ===========================================================================

def _seq_sim_right(s1, s2):
    """Right-aligned similarity: compare from the 3' end backward."""
    if not s1 or not s2:
        return 0.0
    n = min(len(s1), len(s2))
    match = sum(1 for i in range(1, n + 1) if s1[-i] == s2[-i])
    return match / n


def _seq_sim_left(s1, s2):
    """Left-aligned similarity: compare from the 5' end forward."""
    if not s1 or not s2:
        return 0.0
    n = min(len(s1), len(s2))
    match = sum(1 for i in range(n) if s1[i] == s2[i])
    return match / n


def _merge_mutation_class(*classes):
    if any(cls == 'Conservative' for cls in classes):
        return 'Conservative'
    if any(cls == 'Silent' for cls in classes):
        return 'Silent'
    if any(cls == 'UTR' for cls in classes):
        return 'UTR'
    return 'None'


# ===========================================================================
# KnockinDesigner — fully self-contained
# ===========================================================================

class KnockinDesigner:
    """N_start (V46-equivalent) + C_stop dual-mode knockin designer.

    All logic from V36→V44 + V2 is inlined. No external knockin_*.py imports.
    """

    # ------------------------------------------------------------------
    # Initialisation (V36 base + V37 modules + V40 budget + V44 deterministic barcode)
    # ------------------------------------------------------------------

    def __init__(self, args):
        self.args = args
        self.model = getattr(args, 'model', 'N_start')
        self.offtarget_index = {}

        # V44: deterministic barcode seed (must be before barcode pool generation)
        seed = getattr(args, 'barcode_seed', 42)
        random.seed(seed)

        # V37 modules
        self.validator = ConstraintValidator(args)

        # V36 base init
        self.payload_seq  = self._load_payload()
        self.template_seq = self._load_template()
        self.results      = []
        self.stats = {
            "total": 0, "p1_clean": 0, "p1_flexible": 0,
            "p2_bridge_nomut": 0, "p2_bridge_mut": 0,
            "p3a_lha": 0, "p3b_rha": 0,
            "failed": 0, "error": 0,
            "optimal_selected": 0, "suboptimal_selected": 0,
            "aa_change_selected": 0,
            "sgrna_offtarget": 0,
        }
        self.current_record_seq = None

        # V37 granular failure reasons
        self.fail_reason_sequences = {reason.value: [] for reason in DesignError}

        # Barcode library
        print(f"{Colors.OKGREEN}[INFO] Generating restriction-site-free barcode library...{Colors.ENDC}")
        self.barcode_library = self._generate_barcode_library(
            length=self.args.barcode_len,
            size=10000,
            restriction_sites=self.args.restriction_site,
        )
        # V44: sorted for determinism
        self.barcode_library = sorted(self.barcode_library)
        print(f"{Colors.OKBLUE}  Generated {len(self.barcode_library)} clean barcodes "
              f"of {self.args.barcode_len} bp{Colors.ENDC}")
        self.barcode_index = 0

        # V40: HA budget
        self.ha_budget = self._compute_ha_budget()
        print(
            f"{Colors.OKBLUE}  Oligo cap  : {args.max_oligo_len} bp{Colors.ENDC}\n"
            f"{Colors.OKBLUE}  HA budget  : {self.ha_budget} bp  "
            f"(each arm: {args.min_effective_ha}~{self.ha_budget - args.min_effective_ha} bp){Colors.ENDC}"
        )

    # ------------------------------------------------------------------
    # Payload & template loading (V36)
    # ------------------------------------------------------------------

    def _load_payload(self):
        if os.path.isfile(self.args.payload):
            try:
                record = next(SeqIO.parse(self.args.payload, "fasta"))
                return str(record.seq).strip().upper()
            except Exception:
                with open(self.args.payload) as f:
                    return f.read().strip().upper()
        return self.args.payload.strip().upper()

    def _load_template(self):
        content = ""
        if os.path.isfile(self.args.template):
            try:
                with open(self.args.template, 'r') as f:
                    content = f.read()
            except Exception as e:
                print(f"{Colors.FAIL}[Error] Failed to load template: {e}{Colors.ENDC}")
                sys.exit(1)
        else:
            content = self.args.template
        content = content.replace('\r\n', '\n')
        if content.strip().startswith(">"):
            lines = content.split('\n')
            clean_lines = [line.strip() for line in lines
                           if line.strip() and not line.strip().startswith(">")]
            return "".join(clean_lines)
        return content.strip()

    # ------------------------------------------------------------------
    # Barcode generation (V36 + V44 sorted)
    # ------------------------------------------------------------------

    def _generate_barcode_library(self, length, size, restriction_sites):
        barcodes   = set()
        attempts   = 0
        max_attempts = size * 100
        while len(barcodes) < size and attempts < max_attempts:
            attempts += 1
            candidate = ''.join(random.choices('ATCG', k=length))
            candidate_upper = candidate.upper()
            candidate_rc    = str(Seq(candidate).reverse_complement()).upper()
            has_site = False
            for site in restriction_sites:
                site_upper = site.upper()
                site_rc    = str(Seq(site_upper).reverse_complement())
                if (site_upper in candidate_upper or site_rc in candidate_upper
                        or site_upper in candidate_rc or site_rc in candidate_rc):
                    has_site = True
                    break
            if not has_site:
                if re.search(r'(.)\1{5,}', candidate):
                    continue
                barcodes.add(candidate)
        return list(barcodes)

    def generate_barcode(self, length):
        barcode = self.barcode_library[self.barcode_index % len(self.barcode_library)]
        self.barcode_index += 1
        return barcode

    # ------------------------------------------------------------------
    # GC helpers (V36)
    # ------------------------------------------------------------------

    def _calc_gc(self, seq):
        if not seq:
            return 0.0
        seq = seq.upper()
        return round((seq.count('G') + seq.count('C')) / len(seq) * 100, 2)

    def _is_gc_valid(self, seq, is_sgrna=False):
        return self.validator.is_gc_valid(seq, is_sgrna)

    # ------------------------------------------------------------------
    # HA budget (V45 fix)
    # ------------------------------------------------------------------

    def _compute_ha_budget(self):
        """Return max total bp available for LHA + RHA combined.

        The new template uses {sgRNA_PAM_fwd} (23 bp) and {sgRNA_PAM_rc}
        (23 bp).  The V44/V40 parent falls back to max_oligo_len-(169+barcode)
        which gives 121 instead of correct 120.
        """
        try:
            probe = self.template_seq.format(
                sgRNA_PAM_fwd='A' * 23,
                sgRNA_PAM_rc='T' * 23,
                sgRNA_fwd='A' * 20,
                sgRNA_rc='T' * 20,
                pam='NNN',
                pam_rc='NNN',
                LHA='',
                RHA='',
                barcode='A' * self.args.barcode_len,
                payload='',
                J23118_RBS='',
            )
            return self.args.max_oligo_len - len(probe)
        except Exception:
            return self.args.max_oligo_len - (170 + self.args.barcode_len)

    # ------------------------------------------------------------------
    # Sequence extraction — N_start (V36) + start codon capture (V41)
    # ------------------------------------------------------------------

    def get_upstream_and_rha(self, record, feature):
        if self.model == 'N_start':
            return self._get_upstream_rha_nstart(record, feature)
        return self._get_upstream_rha_cstop(record, feature)

    def _get_upstream_rha_nstart(self, record, feature):
        """N_start: extract upstream + RHA, capture start codon (V41)."""
        start  = int(feature.location.start)
        end    = int(feature.location.end)
        strand = feature.location.strand
        full_seq = record.seq
        seq_len  = len(full_seq)
        buffer_len = 500
        scan_len   = self.args.search_window + self.args.lha_len + self.args.rha_len + 200

        if strand == 1:
            req_len = max(scan_len, self.args.lha_len + buffer_len)
            p_start = max(0, start - req_len)
            upstream_raw = full_seq[p_start : start]
            rha_raw      = full_seq[start : min(seq_len, start + self.args.rha_len + buffer_len)]
        else:
            req_len    = max(scan_len, self.args.lha_len + buffer_len)
            p_rha_start = max(0, end - (self.args.rha_len + buffer_len))
            phys_rha    = full_seq[p_rha_start : end]
            rha_raw     = phys_rha.reverse_complement()
            p_up_end    = min(seq_len, end + req_len)
            phys_up     = full_seq[end : p_up_end]
            upstream_raw = phys_up.reverse_complement()

        # V41: capture start codon
        self._current_start_codon = self._get_start_codon(record, feature)
        return str(upstream_raw).upper(), str(rha_raw).upper()

    def _get_upstream_rha_cstop(self, record, feature):
        """C_stop: extract upstream + RHA around stop codon."""
        start    = int(feature.location.start)
        end      = int(feature.location.end)
        strand   = feature.location.strand
        full_seq = record.seq
        seq_len  = len(full_seq)
        rha_buffer = self.args.rha_len + 500
        req_up     = self.args.search_window + self.args.lha_len + 400

        if strand == 1:
            stop_start   = end - 3
            phy_up_s     = max(0, stop_start - req_up)
            upstream_raw = str(full_seq[phy_up_s : stop_start]).upper()
            rha_raw      = str(full_seq[stop_start :
                                        min(seq_len, stop_start + rha_buffer)]).upper()
        else:
            phy_up_e     = min(end, start + 3 + req_up)
            upstream_raw = str(full_seq[start + 3 : phy_up_e].reverse_complement()).upper()
            phy_rha_e    = start + 3
            phy_rha_s    = max(0, start + 3 - rha_buffer)
            rha_raw      = str(full_seq[phy_rha_s : phy_rha_e].reverse_complement()).upper()

        return upstream_raw, rha_raw

    def _build_cstop_rha_with_overlap_rescue(self, upstream_seq, rha_seq,
                                             overlap_bases, max_total_len):
        """Build C_stop RHA as rescued neighbor-CDS prefix + true downstream homology."""
        stop_len            = min(3, len(rha_seq))
        pre_stop_overlap    = min(len(upstream_seq), max(0, overlap_bases - stop_len))
        stop_shared_start   = max(0, stop_len - overlap_bases)
        upstream_rescue     = (upstream_seq[-pre_stop_overlap:]
                               if pre_stop_overlap > 0 else "")
        stop_rescue         = rha_seq[stop_shared_start:stop_len]
        rescue_seq          = upstream_rescue + stop_rescue
        homology_start      = stop_len
        max_homology_len    = max(0, max_total_len - len(rescue_seq))
        homology_seq        = rha_seq[homology_start:homology_start + max_homology_len]
        return rescue_seq, rescue_seq + homology_seq, len(rescue_seq), homology_seq

    def _map_cstop_full_idx_to_rha(self, full_idx, junction_idx, overlap_bases):
        """Map index in upstream+rha space to V5 rescued RHA index; return None if deleted."""
        pre_stop_overlap  = max(0, overlap_bases - 3)
        stop_shared_start = max(0, 3 - overlap_bases)
        rescue_len        = pre_stop_overlap + max(0, 3 - stop_shared_start)

        if full_idx < junction_idx - pre_stop_overlap:
            return None
        if full_idx < junction_idx:
            return full_idx - (junction_idx - pre_stop_overlap)
        if full_idx < junction_idx + stop_shared_start:
            return None
        if full_idx < junction_idx + 3:
            return pre_stop_overlap + (full_idx - (junction_idx + stop_shared_start))
        return rescue_len + (full_idx - (junction_idx + 3))

    # ------------------------------------------------------------------
    # Start codon helpers (V41)
    # ------------------------------------------------------------------

    def _get_start_codon(self, record, feature):
        strand = feature.location.strand
        seq    = record.seq
        if strand == 1:
            sc = str(seq[feature.location.start : feature.location.start + 3]).upper()
        else:
            sc = str(seq[feature.location.end - 3 : feature.location.end]
                     .reverse_complement()).upper()
        return sc

    # ------------------------------------------------------------------
    # CDS validation (V36)
    # ------------------------------------------------------------------

    def _validate_cds_boundary(self, record, feature):
        qualifiers   = feature.qualifiers
        protein_id   = qualifiers.get('protein_id', [None])[0]
        translation  = qualifiers.get('translation', [None])[0]
        if not protein_id or not translation:
            return False, "Missing protein_id/translation"
        translation = translation.replace(" ", "").strip()
        cds_seq = feature.extract(record.seq)
        if len(cds_seq) % 3 != 0:
            return False, "CDS length not multiple of 3"
        match_found = False
        try:
            translated_cds = str(Seq(cds_seq).translate(table=11, cds=True))
            if translated_cds.endswith('*'):
                translated_cds = translated_cds[:-1]
            if translated_cds == translation:
                match_found = True
        except Exception:
            pass
        if not match_found:
            translated = str(Seq(cds_seq).translate(table=11))
            if translated.endswith('*'):
                translated = translated[:-1]
            if translated != translation:
                return False, "Translation mismatch"
        strand = feature.location.strand
        if strand == 1:
            start_codon_pos = int(feature.location.start)
            stop_codon_pos  = int(feature.location.end) - 3
        else:
            start_codon_pos = int(feature.location.end) - 3
            stop_codon_pos  = int(feature.location.start)
        qualifiers['_cds_valid']       = [True]
        qualifiers['_start_codon_pos'] = [str(start_codon_pos)]
        qualifiers['_stop_codon_pos']  = [str(stop_codon_pos)]
        qualifiers['_cds_strand']      = [str(strand)]
        return True, "OK"

    def _validate_cds_stop_codon(self, record, feature):
        """Verify stop codon is TAA/TAG/TGA."""
        strand = feature.location.strand
        seq    = record.seq
        if strand == 1:
            stop_seq = str(seq[int(feature.location.end) - 3 :
                               int(feature.location.end)]).upper()
        else:
            stop_seq = str(seq[int(feature.location.start) :
                               int(feature.location.start) + 3
                               ].reverse_complement()).upper()
        if stop_seq not in STOP_CODONS:
            return False, f"Non-standard stop codon: {stop_seq}"
        return True, stop_seq

    # ------------------------------------------------------------------
    # Codon resolution (V36)
    # ------------------------------------------------------------------

    def _resolve_codon_at_genomic(self, genomic_pos, cds_feature):
        if cds_feature is None:
            return None
        valid_flag = cds_feature.qualifiers.get('_cds_valid', [False])[0]
        if not valid_flag:
            return None
        seq       = self.current_record_seq
        cds_start = int(cds_feature.location.start)
        cds_end   = int(cds_feature.location.end)
        strand    = cds_feature.location.strand
        if strand == 1:
            offset = genomic_pos - cds_start
            if offset < 0:
                return None
            frame       = offset % 3
            codon_start = genomic_pos - frame
            codon       = str(seq[codon_start : codon_start + 3]).upper()
            return codon, (codon_start, codon_start + 2), strand, frame
        else:
            cds_end_incl = cds_end - 1
            offset = cds_end_incl - genomic_pos
            if offset < 0:
                return None
            frame  = offset % 3
            wstart = genomic_pos - (2 - frame)
            wend   = genomic_pos + frame
            segment = str(seq[wstart : wend + 1]).upper()
            codon   = str(Seq(segment).reverse_complement())
            return codon, (wstart, wend), strand, frame

    # ------------------------------------------------------------------
    # Mutation engine (V36/V37)
    # ------------------------------------------------------------------

    def perform_precise_mutation(self, sequence, pam_start_idx, pam_type,
                                 ref_feature, target_strand, genomic_junction,
                                 junction_idx_in_seq, ref_cds_label="Upstream"):
        sequence = sequence.upper()
        seq_list = list(sequence)

        if pam_start_idx < 0 or pam_start_idx + 3 > len(sequence):
            return sequence, False, "PAM out of bounds", -1, "Unknown", 'None'

        original_pam = sequence[pam_start_idx : pam_start_idx + 3]
        pam_indices  = ([pam_start_idx + 1, pam_start_idx + 2]
                        if pam_type == 'NGG'
                        else [pam_start_idx, pam_start_idx + 1])

        def to_genomic(idx):
            offset = idx - junction_idx_in_seq
            if target_strand == 1:
                return (genomic_junction - (-offset)) if offset < 0 else (genomic_junction + offset)
            else:
                return genomic_junction - offset - 1  # FIX: unified formula for LHA & RHA

        # === Level 1: Silent Mutation ===
        for mut_idx in pam_indices:
            if mut_idx >= len(sequence):
                continue
            gpos      = to_genomic(mut_idx)
            loc_label = "LHA" if mut_idx < junction_idx_in_seq else "RHA"

            in_coding = False
            if ref_feature:
                ref_start = int(ref_feature.location.start)
                ref_end   = int(ref_feature.location.end)
                if ref_start <= gpos < ref_end:
                    in_coding = True

            # Case A: Non-coding (UTR) → Any mutation is "Silent"
            if not in_coding:
                base     = seq_list[mut_idx]
                possible = [b for b in 'ATCG' if b != base]
                for new_base in possible:
                    tmp = seq_list[:]
                    tmp[mut_idx] = new_base
                    cand_pam = "".join(tmp[pam_start_idx : pam_start_idx + 3])
                    if ((pam_type == 'NGG' and cand_pam[1:] != "GG")
                            or (pam_type == 'CCN' and cand_pam[:2] != "CC")):
                        seq_list[mut_idx] = new_base
                        return ("".join(seq_list), True,
                                f"Level 1 (UTR): {original_pam}->{cand_pam}",
                            mut_idx, loc_label, 'UTR')
                continue

            # Case B: Coding Region → Try Synonymous
            codon_info = self._resolve_codon_at_genomic(gpos, ref_feature)
            if not codon_info:
                continue
            current_codon_seq, (c_start, c_end), c_strand, index_in_codon = codon_info
            aa = CODON_TABLE.get(current_codon_seq)
            if not aa:
                continue
            syn_codons = AA_TO_CODONS.get(aa, [])
            for syn_codon in syn_codons:
                if syn_codon == current_codon_seq:
                    continue
                if c_strand == 1:
                    bases_on_top = list(syn_codon)
                else:
                    bases_on_top = list(str(Seq(syn_codon).reverse_complement()))
                tmp_seq = seq_list[:]
                valid_mapping = True
                genomic_positions = [c_start, c_start + 1, c_start + 2]
                for k, g_p in enumerate(genomic_positions):
                    base_to_write = bases_on_top[k]
                    if target_strand == -1:
                        base_to_write = str(Seq(base_to_write).reverse_complement())
                    if target_strand == 1:
                        offset = g_p - genomic_junction
                    else:
                        offset = genomic_junction - g_p - 1  # FIX: off-by-one for neg strand
                    idx_in_seq = junction_idx_in_seq + offset
                    if 0 <= idx_in_seq < len(sequence):
                        if base_to_write != seq_list[idx_in_seq]:
                            if idx_in_seq not in pam_indices:
                                valid_mapping = False
                                break
                        tmp_seq[idx_in_seq] = base_to_write
                    else:
                        valid_mapping = False
                        break
                if not valid_mapping:
                    continue
                cand_pam   = "".join(tmp_seq[pam_start_idx : pam_start_idx + 3])
                is_broken  = ((pam_type == 'NGG' and cand_pam[1:] != "GG")
                              or (pam_type == 'CCN' and cand_pam[:2] != "CC"))
                if is_broken:
                    return ("".join(tmp_seq), True,
                            f"Level 1 (Silent): {current_codon_seq}({aa})->{syn_codon} | PAM:{original_pam}->{cand_pam}",
                            mut_idx, loc_label, 'Silent')

        # === Level 2: Conservative Missense Fallback ===
        conservative_hits = []
        for mut_idx in pam_indices:
            if mut_idx >= len(sequence):
                continue
            gpos      = to_genomic(mut_idx)
            loc_label = "LHA" if mut_idx < junction_idx_in_seq else "RHA"
            if not ref_feature:
                continue
            ref_start = int(ref_feature.location.start)
            ref_end   = int(ref_feature.location.end)
            if not (ref_start <= gpos < ref_end):
                continue
            codon_info = self._resolve_codon_at_genomic(gpos, ref_feature)
            if not codon_info:
                continue
            current_codon_seq, (c_start, c_end), c_strand, _ = codon_info
            aa = CODON_TABLE.get(current_codon_seq)
            if not aa:
                continue
            similar_aas = CONSERVATIVE_GROUPS.get(aa, [])
            for target_aa in similar_aas:
                target_codons = AA_TO_CODONS.get(target_aa, [])
                for cand_codon in target_codons:
                    if cand_codon == current_codon_seq:
                        continue
                    if c_strand == 1:
                        bases_on_top = list(cand_codon)
                    else:
                        bases_on_top = list(str(Seq(cand_codon).reverse_complement()))
                    tmp_seq = seq_list[:]
                    valid_mapping = True
                    genomic_positions = [c_start, c_start + 1, c_start + 2]
                    for k, g_p in enumerate(genomic_positions):
                        base_to_write = bases_on_top[k]
                        if target_strand == -1:
                            base_to_write = str(Seq(base_to_write).reverse_complement())
                        if target_strand == 1:
                            offset = g_p - genomic_junction
                        else:
                            offset = genomic_junction - g_p - 1  # FIX: off-by-one for neg strand
                        idx_in_seq = junction_idx_in_seq + offset
                        if 0 <= idx_in_seq < len(sequence):
                            if base_to_write != seq_list[idx_in_seq]:
                                if idx_in_seq not in pam_indices:
                                    valid_mapping = False
                                    break
                            tmp_seq[idx_in_seq] = base_to_write
                        else:
                            valid_mapping = False
                            break
                    if not valid_mapping:
                        continue
                    cand_pam  = "".join(tmp_seq[pam_start_idx : pam_start_idx + 3])
                    is_broken = ((pam_type == 'NGG' and cand_pam[1:] != "GG")
                                 or (pam_type == 'CCN' and cand_pam[:2] != "CC"))
                    if is_broken:
                        edit_count = sum(1 for a, b in zip(current_codon_seq, cand_codon) if a != b)
                        conservative_hits.append({
                            'edit_count': edit_count,
                            'mut_idx': mut_idx,
                            'loc_label': loc_label,
                            'cand_seq': "".join(tmp_seq),
                            'cand_codon': cand_codon,
                            'target_aa': target_aa,
                            'cand_pam': cand_pam,
                            'aa': aa,
                            'current_codon_seq': current_codon_seq,
                        })

        if conservative_hits:
            conservative_hits.sort(key=lambda item: (item['edit_count'], item['mut_idx'], item['cand_codon']))
            best_hit = conservative_hits[0]
            return (
                best_hit['cand_seq'], True,
                f"Level 2 (Conservative): {best_hit['current_codon_seq']}({best_hit['aa']})->{best_hit['cand_codon']}({best_hit['target_aa']}) | PAM:{original_pam}->{best_hit['cand_pam']}",
                best_hit['mut_idx'], best_hit['loc_label'], 'Conservative'
            )

        return sequence, False, "Mutation Failed (No Valid Candidates)", -1, "Unknown", 'None'

    def _candidate_mutation_class(self, cand):
        return cand.get('Mutation_Class', 'None')

    def _candidate_has_aa_change(self, cand):
        return self._candidate_mutation_class(cand) == 'Conservative'

    def _candidate_bio_rank(self, cand):
        return 1 if self._candidate_has_aa_change(cand) else 0

    def _quality_status(self, cand):
        return ("(Optimal)"
                if cand.get('Mut_Dist', 0) <= self.args.max_mut_dist
                else "(Suboptimal [Mut_Dist])")

    def _aa_cost_tag(self, cand):
        if self._candidate_has_aa_change(cand):
            return "[AA_Change:Conservative]"
        return ""

    # ------------------------------------------------------------------
    # RE helpers (V36 + V39)
    # ------------------------------------------------------------------

    def _map_idx_to_genomic(self, idx, junction_idx, genomic_junction, strand):
        if strand == 1:
            return genomic_junction - (junction_idx - idx)
        return genomic_junction + (junction_idx - idx) - 1  # FIX: off-by-one for neg strand

    def _map_genomic_to_idx(self, gpos, junction_idx, genomic_junction, strand):
        if strand == 1:
            return junction_idx - (genomic_junction - gpos)
        return junction_idx - (gpos - genomic_junction) - 1  # FIX: off-by-one for neg strand

    def _silent_fix_restriction_site(self, seq, site, ref_feature, target_strand,
                                      genomic_junction, junction_idx):
        if not ref_feature:
            return seq, False, "No ref feature"
        seq      = seq.upper()
        seq_list = list(seq)
        site_u   = site.upper()
        site_rc  = str(Seq(site_u).reverse_complement())
        site_len = len(site_u)

        for pattern, _ in [(site_u, False), (site_rc, True)]:
            start = 0
            while True:
                pos = seq.find(pattern, start)
                if pos == -1:
                    break
                start = pos + 1
                for offset in range(site_len):
                    idx  = pos + offset
                    gpos = self._map_idx_to_genomic(idx, junction_idx, genomic_junction, target_strand)
                    if not (int(ref_feature.location.start) <= gpos < int(ref_feature.location.end)):
                        continue
                    codon_info = self._resolve_codon_at_genomic(gpos, ref_feature)
                    if not codon_info:
                        continue
                    current_codon_seq, (c_start, c_end), c_strand, _ = codon_info
                    aa = CODON_TABLE.get(current_codon_seq)
                    if not aa:
                        continue
                    for syn_codon in AA_TO_CODONS.get(aa, []):
                        if syn_codon == current_codon_seq:
                            continue
                        if c_strand == 1:
                            bases_on_top = list(syn_codon)
                        else:
                            bases_on_top = list(str(Seq(syn_codon).reverse_complement()))
                        tmp = seq_list[:]
                        genomic_positions = [c_start, c_start + 1, c_start + 2]
                        valid = True
                        for k, g_p in enumerate(genomic_positions):
                            base_to_write = bases_on_top[k]
                            if target_strand == -1:
                                base_to_write = str(Seq(base_to_write).reverse_complement())
                            idx_in_seq = self._map_genomic_to_idx(
                                g_p, junction_idx, genomic_junction, target_strand)
                            if 0 <= idx_in_seq < len(tmp):
                                tmp[idx_in_seq] = base_to_write
                            else:
                                valid = False
                                break
                        if not valid:
                            continue
                        candidate_seq = "".join(tmp)
                        if site_u not in candidate_seq and site_rc not in candidate_seq:
                            return candidate_seq, True, f"Silent Mut {site_u}"
        return seq, False, "No Silent Fix"

    def _sanitize_arm(self, seq, arm_label, ref_feature, target_strand,
                       genomic_junction, junction_idx):
        updated_seq = seq
        logs    = []
        changed = False
        for site in self.args.restriction_site:
            updated_seq, ok, log = self._silent_fix_restriction_site(
                updated_seq, site, ref_feature, target_strand, genomic_junction, junction_idx)
            if ok:
                changed = True
                logs.append(f"{arm_label}:{log}")
        return updated_seq, logs, changed

    @staticmethod
    def _seq_contains_re(seq):
        """Return True if seq contains any BsaI or BbsI site (fwd or RC)."""
        if not seq:
            return False
        seq_upper = seq.upper()
        for site_fwd in FORBIDDEN_RE_SEQS:
            site_rc = str(Seq(site_fwd).reverse_complement())
            if site_fwd in seq_upper or site_rc in seq_upper:
                return True
        return False

    def _find_variable_re_violations(self, final_oligo):
        """Scan oligo for BsaI/BbsI sites outside the template backbone."""
        oligo_upper    = final_oligo.upper()
        backbone_start = oligo_upper.find(TEMPLATE_BACKBONE_BSAI)
        backbone_end   = (backbone_start + len(TEMPLATE_BACKBONE_BSAI)
                          if backbone_start >= 0 else -1)
        violations = []
        for site_fwd in FORBIDDEN_RE_SEQS:
            site_rc = str(Seq(site_fwd).reverse_complement())
            for pattern in [site_fwd, site_rc]:
                start = 0
                while True:
                    pos = oligo_upper.find(pattern, start)
                    if pos == -1:
                        break
                    if (backbone_start >= 0
                            and backbone_start <= pos
                            and pos + len(pattern) <= backbone_end):
                        start = pos + 1
                        continue
                    violations.append((pattern, pos))
                    start = pos + 1
        return violations

    # ------------------------------------------------------------------
    # Recut risk & restriction site analysis (V36 + V39)
    # ------------------------------------------------------------------

    def check_recut_risk(self, final_oligo, spacer_seq, pam_seq):
        return self.validator.check_recut_risk(final_oligo, spacer_seq, pam_seq)

    def _find_site_positions(self, seq, site):
        if not seq or not site:
            return []
        positions  = []
        seq_upper  = seq.upper()
        site_upper = site.upper()
        site_rc    = str(Seq(site_upper).reverse_complement())
        start = 0
        while True:
            pos = seq_upper.find(site_upper, start)
            if pos == -1:
                break
            positions.append((pos, '+'))
            start = pos + 1
        start = 0
        while True:
            pos = seq_upper.find(site_rc, start)
            if pos == -1:
                break
            positions.append((pos, '-'))
            start = pos + 1
        return positions

    def _identify_junction_boundaries(self, final_oligo, design_data, site_positions):
        oligo_upper = final_oligo.upper()
        lha_seq     = design_data['LHA']
        rha_seq     = design_data['RHA']
        sgrna_seq   = design_data['sgRNA']
        barcode     = design_data['Barcode']
        payload     = self.payload_seq

        components = []
        lha_pos     = oligo_upper.find(lha_seq.upper())
        rha_pos     = oligo_upper.find(rha_seq.upper())
        sgrna_fwd   = oligo_upper.find(sgrna_seq.upper())
        payload_pos = oligo_upper.find(payload.upper())
        barcode_pos = oligo_upper.find(barcode.upper())

        if lha_pos >= 0:     components.append(('LHA',     lha_pos,     lha_pos + len(lha_seq)))
        if rha_pos >= 0:     components.append(('RHA',     rha_pos,     rha_pos + len(rha_seq)))
        if sgrna_fwd >= 0:   components.append(('sgRNA',   sgrna_fwd,   sgrna_fwd + len(sgrna_seq)))
        if payload_pos >= 0: components.append(('Payload', payload_pos, payload_pos + len(payload)))
        if barcode_pos >= 0: components.append(('Barcode', barcode_pos, barcode_pos + len(barcode)))
        components.sort(key=lambda x: x[1])

        junction_sources = []
        site_len = 6
        for site_pos, strand in site_positions:
            site_end   = site_pos + site_len
            spanning   = [n for n, s, e in components if site_pos < e and site_end > s]
            if len(spanning) == 0:
                left_comp  = None
                right_comp = None
                for n, s, e in components:
                    if e <= site_pos:
                        left_comp = n
                    if s >= site_end and right_comp is None:
                        right_comp = n
                if left_comp and right_comp:
                    junction_sources.append(f"[{left_comp}|{right_comp}]")
                elif left_comp:
                    junction_sources.append(f"[{left_comp}|Backbone]")
                elif right_comp:
                    junction_sources.append(f"[Backbone|{right_comp}]")
                else:
                    junction_sources.append("[Backbone]")
            elif len(spanning) > 1:
                junction_sources.append(f"[{'|'.join(spanning)}]")
        return junction_sources

    def check_restriction_sites(self, design_data, final_oligo):
        if not self.args.restriction_site:
            return "", "None", ""
        found_sites_details = []
        found_sites_list    = []
        final_oligo_upper   = final_oligo.upper()
        restriction_entries = []

        def count_sites(seq, site):
            return len(self._find_site_positions(seq, site)) if seq else 0

        for site in self.args.restriction_site:
            site_u      = site.upper()
            total_count = count_sites(final_oligo_upper, site)
            if total_count == 0:
                continue
            lha_count     = len(self._find_site_positions(design_data['LHA'], site))
            rha_count     = len(self._find_site_positions(design_data['RHA'], site))
            sgrna_count   = len(self._find_site_positions(design_data['sgRNA'], site))
            payload_count = len(self._find_site_positions(self.payload_seq, site))
            barcode_count = len(self._find_site_positions(design_data['Barcode'], site))

            parts_total    = lha_count + rha_count + sgrna_count + payload_count + barcode_count
            junction_count = total_count - parts_total

            junction_details = ""
            if junction_count > 0:
                all_sp        = self._find_site_positions(final_oligo, site)
                component_sp  = (self._find_site_positions(design_data['LHA'], site)
                                 + self._find_site_positions(design_data['RHA'], site)
                                 + self._find_site_positions(design_data['sgRNA'], site)
                                 + self._find_site_positions(self.payload_seq, site)
                                 + self._find_site_positions(design_data['Barcode'], site))
                junction_sp   = [p for p in all_sp if p not in component_sp]
                j_sources     = self._identify_junction_boundaries(final_oligo, design_data, junction_sp)
                if j_sources:
                    junction_details = ",".join(
                        [f"Junction:{idx+1}{src}" for idx, src in enumerate(j_sources)])

            comp_parts = []
            if lha_count > 0:     comp_parts.append(f"LHA:{lha_count}")
            if rha_count > 0:     comp_parts.append(f"RHA:{rha_count}")
            if sgrna_count > 0:   comp_parts.append(f"sgRNA:{sgrna_count}")
            if payload_count > 0: comp_parts.append(f"Payload:{payload_count}")
            if barcode_count > 0: comp_parts.append(f"Barcode:{barcode_count}")
            if junction_count > 0:
                comp_parts.append(junction_details if junction_details else f"Junction:{junction_count}")

            if comp_parts:
                restriction_entries.append(f"{site_u}({','.join(comp_parts)})")
            found_sites_list.append(site_u)

            locs = []
            if lha_count > 0:                 locs.append("LHA")
            if rha_count > 0:                 locs.append("RHA")
            if sgrna_count > 0:               locs.append("sgRNA")
            if payload_count > 0:             locs.append("Payload")
            if barcode_count > 0:             locs.append("Barcode")
            if total_count > parts_total:     locs.append("Backbone/Junction")
            found_sites_details.append(f"{site_u}: {', '.join(locs) if locs else 'Unknown'} ({total_count})")

        details_str = " | ".join(found_sites_details) if found_sites_details else ""
        types_str   = ";".join(found_sites_list) if found_sites_list else "None"
        re_str      = ",".join(restriction_entries) if restriction_entries else ""
        return details_str, types_str, re_str

    # ------------------------------------------------------------------
    # LHA optimisation (V36)
    # ------------------------------------------------------------------

    def _optimize_lha_window(self, upstream_seq, initial_end_idx,
                              min_required_len, max_allowed_len=None):
        if max_allowed_len is None:
            max_allowed_len = self.args.lha_len
        if max_allowed_len <= 0:
            return None, 0
        start_idx   = max(0, initial_end_idx - max_allowed_len)
        current_lha = upstream_seq[start_idx : initial_end_idx]
        if len(current_lha) < min_required_len:
            return None, 0
        if self._is_gc_valid(current_lha):
            return current_lha, 0
        max_cut = len(current_lha) - min_required_len
        for cut_len in range(1, max_cut + 1):
            cand_lha = current_lha[cut_len:]
            if self._is_gc_valid(cand_lha):
                return cand_lha, cut_len
        return None, 0

    # ------------------------------------------------------------------
    # Oligo assembly (V2 override with new template placeholders + V39 suppression)
    # ------------------------------------------------------------------

    def assemble_final_oligo(self, design_data):
        barcode   = self.generate_barcode(self.args.barcode_len)
        sgrna_seq = design_data['sgRNA']
        pam_seq   = design_data['PAM']
        sgrna_rc  = str(Seq(sgrna_seq).reverse_complement())
        pam_rc    = str(Seq(pam_seq).reverse_complement())

        if pam_seq[:2] == 'CC':
            sgRNA_PAM_fwd = sgrna_seq + str(Seq(pam_seq).reverse_complement())
        else:
            sgRNA_PAM_fwd = sgrna_seq + pam_seq
        sgRNA_PAM_rc = str(Seq(sgRNA_PAM_fwd).reverse_complement())

        context = {
            "sgRNA_PAM_fwd": sgRNA_PAM_fwd,
            "sgRNA_PAM_rc":  sgRNA_PAM_rc,
            "sgRNA_fwd": sgrna_seq,  "sgRNA_rc": sgrna_rc,
            "pam": pam_seq,          "pam_rc":   pam_rc,
            "LHA": design_data['LHA'], "RHA": design_data['RHA'],
            "barcode": barcode,
            "payload": self.payload_seq, "J23118_RBS": self.payload_seq,
        }
        try:
            final_oligo = self.template_seq.format(**context)
            is_risky, risk_msg = self.check_recut_risk(final_oligo, sgrna_seq, pam_seq)
            if "Mut" in design_data['Strategy'] or "Priority1" in design_data['Strategy']:
                is_risky = False

            # V39: suppress template-inherent RC risk
            if risk_msg == "Risk: Intact Rev Target":
                is_risky = False
                risk_msg = "Safe"

            design_data['Barcode'] = barcode
            re_details, re_types, restriction_sites = self.check_restriction_sites(
                design_data, final_oligo)

            # V39: clear if all sites are within backbone
            if restriction_sites and final_oligo:
                true_violations = self._find_variable_re_violations(final_oligo)
                if not true_violations:
                    restriction_sites = ""

            return (final_oligo, barcode, is_risky, risk_msg,
                    re_details, re_types, restriction_sites)
        except KeyError as e:
            return (f"Error: Missing placeholder {e}", barcode,
                    True, "Template Error", "", "", "")

    # ------------------------------------------------------------------
    # Mutation-aware balanced HA trim (V42)
    # ------------------------------------------------------------------

    def _mutation_aware_trim(self, lha, rha, lha_mut_idx, rha_mut_idx,
                             rha_prefix_len=0):
        budget = self.ha_budget
        min_ha = self.args.min_effective_ha
        lha_len = len(lha)
        rha_len = len(rha)

        E_L = (lha_mut_idx + 1)       if lha_mut_idx is not None else lha_len
        protected_rha_start = rha_prefix_len
        if rha_mut_idx is not None and rha_mut_idx >= rha_prefix_len:
            protected_rha_start = rha_mut_idx
        E_R = rha_len - protected_rha_start

        if E_L < min_ha or E_R < min_ha:
            return None

        max_trim_L = E_L - min_ha
        max_trim_R = E_R - min_ha

        excess = lha_len + rha_len - budget
        if excess <= 0:
            return lha, rha, lha_mut_idx, E_L, E_R

        if max_trim_L + max_trim_R < excess:
            return None

        D      = E_L - E_R
        trim_L = (excess + D) // 2
        trim_R = excess - trim_L

        trim_L = max(0, min(trim_L, max_trim_L))
        trim_R = max(0, min(trim_R, max_trim_R))

        deficit = excess - trim_L - trim_R
        if deficit > 0:
            eff_L_cur = E_L - trim_L
            eff_R_cur = E_R - trim_R
            if eff_L_cur >= eff_R_cur:
                add = min(deficit, max_trim_L - trim_L)
                trim_L += add
                deficit -= add
                if deficit > 0:
                    add = min(deficit, max_trim_R - trim_R)
                    trim_R += add
                    deficit -= add
            else:
                add = min(deficit, max_trim_R - trim_R)
                trim_R += add
                deficit -= add
                if deficit > 0:
                    add = min(deficit, max_trim_L - trim_L)
                    trim_L += add
                    deficit -= add
            if deficit > 0:
                return None

        new_lha     = lha[trim_L:]
        new_rha     = rha[:rha_len - trim_R]
        new_lha_mut = (lha_mut_idx - trim_L) if lha_mut_idx is not None else None
        eff_lha     = E_L - trim_L
        eff_rha     = E_R - trim_R
        return new_lha, new_rha, new_lha_mut, eff_lha, eff_rha

    # ------------------------------------------------------------------
    # Off-target index (V46 / V2)
    # ------------------------------------------------------------------

    def _build_offtarget_index(self, genome_str):
        g      = genome_str.upper()
        g_circ = g + g[:22]
        n      = len(g_circ)
        counts = Counter()
        for i in range(1, n - 22):
            dinuc = g_circ[i:i + 2]
            if dinuc == 'GG':
                spacer = g_circ[i - 21:i - 1]
                if len(spacer) == 20:
                    counts[spacer] += 1
            elif dinuc == 'CC':
                t_seq = g_circ[i + 3:i + 23]
                if len(t_seq) == 20:
                    spacer = str(Seq(t_seq).reverse_complement())
                    counts[spacer] += 1
        return counts

    def _build_and_print_index(self, genome_str):
        print(f"{Colors.OKBLUE}[INFO] Building genome-wide spacer index "
              f"({len(genome_str):,} bp)...{Colors.ENDC}")
        self.offtarget_index = self._build_offtarget_index(genome_str)
        unique_spacers  = len(self.offtarget_index)
        total_hits      = sum(self.offtarget_index.values())
        unique_targets  = sum(1 for v in self.offtarget_index.values() if v == 1)
        print(f"{Colors.OKBLUE}  {unique_spacers:,} unique spacers | "
              f"{total_hits:,} total hits | "
              f"{unique_targets:,} genome-unique (count=1){Colors.ENDC}")
        print(f"{Colors.OKBLUE}  Filter threshold: count > {self.args.max_offtargets + 1} → rejected  "
              f"(--max_offtargets {self.args.max_offtargets}){Colors.ENDC}")

    # ------------------------------------------------------------------
    # design_target() — dispatch by mode
    # ------------------------------------------------------------------

    def design_target(self, upstream_seq, rha_seq, safe_limit, overlap,
                      upstream_feat, current_feat, strand, genomic_junction):
        if self.model == 'N_start':
            return self._design_target_nstart(
                upstream_seq, rha_seq, safe_limit, overlap,
                upstream_feat, current_feat, strand, genomic_junction)
        return self._design_target_cstop(
            upstream_seq, rha_seq, safe_limit, overlap,
            upstream_feat, current_feat, strand, genomic_junction)

    # ------------------------------------------------------------------
    # N_start design_target (V41 — fixed CCN logic)
    # ------------------------------------------------------------------

    def _design_target_nstart(self, upstream_seq, rha_seq, safe_upstream_limit,
                               overlap_bases, upstream_feature, current_feature,
                               target_strand, genomic_junction):
        candidates   = []
        junction_idx = len(upstream_seq)
        reason_counts = Counter() if self.args.report_fail_reasons else None
        full = upstream_seq + rha_seq
        start = max(0, junction_idx - self.args.search_window)
        end   = junction_idx + self.args.search_window

        def bump(reason_enum):
            if reason_counts is not None:
                reason_counts[reason_enum.value] += 1

        pam_found = False
        for i in range(start, end):
            if i + 23 > len(full):
                break
            ptype = None
            if full[i:i+2] == "GG":
                ptype = 'NGG'
            elif full[i:i+2] == "CC":
                ptype = 'CCN'
            if ptype is None:
                continue

            pam_found = True
            if ptype == 'NGG':
                pam_s    = i - 1
                sp_e     = pam_s
                sp_s     = sp_e - 20
                cut      = pam_s - 3
                sp_seq   = full[sp_s:sp_e]
                pam_seq  = full[pam_s:pam_s+3]
                t_s, t_e = sp_s, pam_s + 3
            else:
                pam_s    = i
                sp_s     = pam_s + 3
                sp_e     = sp_s + 20
                cut      = pam_s + 6
                t_seq    = full[sp_s:sp_e]
                sp_seq   = str(Seq(t_seq).reverse_complement())
                pam_seq  = full[pam_s:pam_s+3]
                t_s, t_e = pam_s, sp_e

            if t_s < 0 or t_e > len(full):
                bump(DesignError.TARGET_WINDOW_OOB)
                continue
            dist = abs(cut - junction_idx)

            if self.args.target_gene:
                print(f"[DEBUG] PAM at {pam_s} ({ptype}), Dist={dist}, Seq={pam_seq}")

            if not self.validator.is_gc_valid(sp_seq, is_sgrna=True):
                if self.args.target_gene:
                    print(f"[DEBUG] -> GC Invalid: {self.validator.calc_gc(sp_seq)}")
                bump(DesignError.SGRNA_GC_INVALID)
                continue

            cut_offset = cut - junction_idx

            # ====== Priority 1: Deletion Strategy (target fully upstream) ======
            if t_e <= junction_idx:
                # P1_Clean
                del_len_clean = junction_idx - t_s
                if del_len_clean <= safe_upstream_limit and del_len_clean <= self.args.max_deletion:
                    req_lha_len = self.args.min_effective_ha
                    lha_fin, shift = self._optimize_lha_window(full, t_s, req_lha_len)
                    std_rha = rha_seq[:self.args.rha_len]
                    if lha_fin and self.validator.is_gc_valid(std_rha):
                        note_local = f"UTR Deletion: {del_len_clean}bp"
                        candidates.append({
                            'sgRNA': sp_seq, 'PAM': pam_seq, 'LHA': lha_fin, 'RHA': std_rha,
                            'Distance': dist, 'Mut_Dist': 0, 'Strategy': 'Priority1_Clean',
                            'Note_Extra': note_local
                        })
                    else:
                        bump(DesignError.P1_CLEAN_LHA_RHA_INVALID)
                else:
                    bump(DesignError.P1_CLEAN_DEL_TOO_LONG)

                # P1_Partial
                candidates_p1_partial = []

                if ptype == 'NGG':
                    crit_start_seed = max(0, pam_s - self.args.seed_region)
                    del_len_seed    = junction_idx - crit_start_seed
                    crit_start_pam  = pam_s
                    del_len_pam     = junction_idx - crit_start_pam

                    if del_len_seed <= self.args.max_deletion and del_len_seed <= safe_upstream_limit:
                        lha_fin_seed, shift = self._optimize_lha_window(
                            full, crit_start_seed, self.args.min_effective_ha)
                        std_rha = rha_seq[:self.args.rha_len]
                        if lha_fin_seed and self.validator.is_gc_valid(std_rha):
                            retained_len = crit_start_seed - t_s
                            candidates_p1_partial.append({
                                'sgRNA': sp_seq, 'PAM': pam_seq, 'LHA': lha_fin_seed, 'RHA': std_rha,
                                'Distance': dist, 'Mut_Dist': 0, 'Strategy': 'Priority1_Partial_Seed',
                                'Note_Extra': f"Partial Deletion (Seed+PAM): -{del_len_seed}bp | Retained Spacer: {retained_len}bp"
                            })
                        else:
                            bump(DesignError.P1_PARTIAL_SEED_LHA_RHA_INVALID)
                    else:
                        bump(DesignError.P1_PARTIAL_SEED_DEL_TOO_LONG)

                    if del_len_pam <= self.args.max_deletion and del_len_pam <= safe_upstream_limit:
                        if not candidates_p1_partial:
                            lha_fin_pam, shift = self._optimize_lha_window(
                                full, crit_start_pam, self.args.min_effective_ha)
                            std_rha = rha_seq[:self.args.rha_len]
                            if lha_fin_pam and self.validator.is_gc_valid(std_rha):
                                retained_len = crit_start_pam - t_s
                                candidates_p1_partial.append({
                                    'sgRNA': sp_seq, 'PAM': pam_seq, 'LHA': lha_fin_pam, 'RHA': std_rha,
                                    'Distance': dist, 'Mut_Dist': 0, 'Strategy': 'Priority1_Partial_PAM',
                                    'Note_Extra': f"Partial Deletion (PAM Only): -{del_len_pam}bp | Retained Spacer: {retained_len}bp"
                                })
                            else:
                                bump(DesignError.P1_PARTIAL_PAM_LHA_RHA_INVALID)
                    else:
                        bump(DesignError.P1_PARTIAL_PAM_DEL_TOO_LONG)

                elif ptype == 'CCN':
                    # V41 FIX (A): guard del_len_seed > 0
                    crit_end_seed = min(len(full), pam_s + 3 + self.args.seed_region)
                    del_len_seed  = crit_end_seed - junction_idx
                    crit_end_pam  = pam_s + 3
                    del_len_pam   = crit_end_pam - junction_idx

                    if del_len_seed > 0 and del_len_seed <= self.args.max_deletion:
                        lha_fin_seed, shift = self._optimize_lha_window(
                            upstream_seq, junction_idx, self.args.min_effective_ha)
                        std_rha = rha_seq[:self.args.rha_len]
                        if lha_fin_seed and std_rha and self.validator.is_gc_valid(std_rha):
                            retained_len = t_e - crit_end_seed
                            candidates_p1_partial.append({
                                'sgRNA': sp_seq, 'PAM': pam_seq, 'LHA': lha_fin_seed, 'RHA': std_rha,
                                'Distance': dist, 'Mut_Dist': 0, 'Strategy': 'Priority1_Partial_Seed',
                                'Note_Extra': f"Partial Deletion (PAM+Seed): -{del_len_seed}bp | Retained Spacer: {retained_len}bp"
                            })
                        else:
                            bump(DesignError.P1_PARTIAL_SEED_LHA_RHA_INVALID)
                    else:
                        bump(DesignError.P1_PARTIAL_SEED_DEL_TOO_LONG)

                    if del_len_pam > 0 and del_len_pam <= self.args.max_deletion:
                        if not candidates_p1_partial:
                            lha_fin_pam, shift = self._optimize_lha_window(
                                upstream_seq, junction_idx, self.args.min_effective_ha)
                            std_rha = rha_seq[:self.args.rha_len]
                            if lha_fin_pam and std_rha and self.validator.is_gc_valid(std_rha):
                                retained_len = t_e - crit_end_pam
                                candidates_p1_partial.append({
                                    'sgRNA': sp_seq, 'PAM': pam_seq, 'LHA': lha_fin_pam, 'RHA': std_rha,
                                    'Distance': dist, 'Mut_Dist': 0, 'Strategy': 'Priority1_Partial_PAM',
                                    'Note_Extra': f"Partial Deletion (PAM Only): -{del_len_pam}bp | Retained Spacer: {retained_len}bp"
                                })
                            else:
                                bump(DesignError.P1_PARTIAL_PAM_LHA_RHA_INVALID)
                    else:
                        bump(DesignError.P1_PARTIAL_PAM_DEL_TOO_LONG)

                candidates.extend(candidates_p1_partial)

            # ====== V41 FIX (B): CCN Downstream — P3B Only ======
            elif ptype == 'CCN' and t_s >= junction_idx:
                req_rha_base = self.args.min_effective_ha + cut_offset
                req_lha_base = self.args.min_effective_ha

                if req_rha_base <= self.args.rha_len and req_lha_base <= self.args.lha_len:
                    pam_rel = pam_s - junction_idx

                    mut_rha, ok, log, mut_idx, _, mut_class = self.perform_precise_mutation(
                        rha_seq, pam_rel, ptype, current_feature, target_strand,
                        genomic_junction, 0, ref_cds_label="Current"
                    )
                    if ok:
                        grab_len = min(self.args.rha_len, len(mut_rha))
                        if mut_idx >= grab_len:
                            bump(DesignError.P3B_MUT_IDX_OOR)
                        else:
                            rha_extracted = mut_rha[:grab_len]
                            lha_base, _ = self._optimize_lha_window(
                                upstream_seq, junction_idx,
                                self.args.min_effective_ha, self.args.lha_len)
                            rescue_seq    = rha_seq[:overlap_bases]
                            lha_extracted = (lha_base + rescue_seq) if lha_base else None
                            if lha_extracted and self.validator.is_gc_valid(rha_extracted):
                                candidates.append({
                                    'sgRNA': sp_seq, 'PAM': pam_seq,
                                    'LHA': lha_extracted, 'RHA': rha_extracted,
                                    'Distance': dist, 'Mut_Dist': dist,
                                    'Strategy': 'Priority3_Mut_RHA',
                                    'Note_Extra': log, 'RHA_Mut_Index': mut_idx,
                                    'Mutation_Class': mut_class,
                                    'RHA_Buffer': mut_rha,
                                })
                            else:
                                bump(DesignError.P3B_LHA_MISSING if not lha_extracted
                                     else DesignError.P3B_RHA_GC_INVALID)
                    else:
                        bump(DesignError.P3B_MUT_FAILED)
                else:
                    bump(DesignError.P3B_REQ_LEN_EXCEEDS)

            # ====== Priority 2: Bridge Strategy (target spans junction) ======
            elif t_s < junction_idx < t_e:
                req_lha_base = self.args.min_effective_ha + (abs(cut_offset) if cut_offset < 0 else 0)
                req_rha_base = self.args.min_effective_ha + (cut_offset if cut_offset > 0 else 0)

                if req_lha_base <= self.args.lha_len and req_rha_base <= self.args.rha_len:
                    grab_lha_len  = min(self.args.lha_len, junction_idx)
                    raw_lha_part  = upstream_seq[junction_idx - grab_lha_len : junction_idx]
                    rescue_seq    = rha_seq[:overlap_bases]
                    lha_cand_full = raw_lha_part + rescue_seq

                    grab_rha_len = min(self.args.rha_len, len(rha_seq))
                    raw_rha      = rha_seq[:grab_rha_len]

                    eff_rha_check = len(raw_rha) - (cut_offset if cut_offset > 0 else 0)
                    if eff_rha_check >= self.args.min_effective_ha and self.validator.is_gc_valid(raw_rha):
                        if ptype == 'NGG':
                            physical_target = full[sp_s : pam_s + 3]
                        else:
                            physical_target = full[pam_s : sp_e]

                        seed_len = self.args.seed_region
                        seed_seq = (physical_target[-(seed_len + 3):]
                                    if ptype == 'NGG'
                                    else physical_target[:(3 + seed_len)])

                        risk_in_lha = ((physical_target in lha_cand_full)
                                       or (seed_seq in lha_cand_full))
                        risk_in_rha = ((physical_target in raw_rha)
                                       or (seed_seq in raw_rha))

                        lha_final = rha_final = None
                        lha_mut_idx = rha_mut_idx = None
                        note     = ""
                        mut_dist = 0
                        is_mut_strategy = False
                        mutation_class = 'None'

                        # 1. RHA Mutation
                        if risk_in_rha:
                            pam_in_rha_idx = pam_s - junction_idx
                            mut_rha, ok, log, mut_idx, loc_lbl, mut_class = self.perform_precise_mutation(
                                raw_rha, pam_in_rha_idx, ptype, current_feature,
                                target_strand, genomic_junction, 0, ref_cds_label="Current")
                            if ok:
                                rha_final   = mut_rha
                                note       += f"RHA Mut: {log} [{loc_lbl}] | "
                                is_mut_strategy = True
                                rha_mut_idx = mut_idx
                                mutation_class = _merge_mutation_class(mutation_class, mut_class)
                            else:
                                rha_final = None
                                bump(DesignError.P2_RHA_MUT_FAILED)
                        else:
                            rha_final = raw_rha

                        # 2. LHA Mutation
                        if rha_final and risk_in_lha:
                            start_of_lha_cand = junction_idx - len(raw_lha_part)
                            pam_in_lha_idx    = pam_s - start_of_lha_cand
                            mut_lha, ok, log, mut_idx, loc_lbl, mut_class = self.perform_precise_mutation(
                                lha_cand_full, pam_in_lha_idx, ptype, upstream_feature,
                                target_strand, genomic_junction, len(raw_lha_part),
                                ref_cds_label="Upstream")
                            if ok:
                                req_opt          = req_lha_base + overlap_bases
                                lha_final, shift = self._optimize_lha_window(mut_lha, len(mut_lha), req_opt)
                                if lha_final is not None:
                                    lha_mut_idx = mut_idx - shift
                                    if lha_mut_idx < 0:
                                        lha_final = None
                                        bump(DesignError.P2_LHA_MUT_SHIFT_OOR)
                                note    += f"LHA Mut: {log} [{loc_lbl}]"
                                is_mut_strategy = True
                                mut_dist = len(mut_lha) - mut_idx
                                mutation_class = _merge_mutation_class(mutation_class, mut_class)
                            else:
                                lha_final = None
                                bump(DesignError.P2_LHA_MUT_FAILED)
                        elif rha_final:
                            req_opt          = req_lha_base + overlap_bases
                            lha_final, shift = self._optimize_lha_window(
                                lha_cand_full, len(lha_cand_full), req_opt)
                            if shift > 0:
                                note += f"LHA GC Shift -{shift}bp"

                        if lha_final and rha_final and len(lha_final) <= self.args.lha_len:
                            stgy       = 'Priority2_Bridge_Mut' if is_mut_strategy else 'Priority2_Bridge'
                            note_clean = note.strip(" | ")
                            cand = {
                                'sgRNA': sp_seq, 'PAM': pam_seq,
                                'LHA': lha_final, 'RHA': rha_final,
                                'Distance': dist, 'Strategy': stgy,
                                'Note_Extra': note_clean, 'Mut_Dist': mut_dist,
                                'Mutation_Class': mutation_class,
                            }
                            if is_mut_strategy:
                                if 'LHA Mut' in note and lha_mut_idx is not None:
                                    cand['LHA_Mut_Index'] = lha_mut_idx
                                if 'RHA Mut' in note and rha_mut_idx is not None:
                                    cand['RHA_Mut_Index'] = rha_mut_idx
                            candidates.append(cand)
                        else:
                            bump(DesignError.P2_LHA_RHA_INVALID)
                    else:
                        bump(DesignError.P2_RHA_GC_OR_TOO_SHORT)
                else:
                    bump(DesignError.P2_REQUIRED_LEN_EXCEEDS)

            # ====== Priority 3B: RHA Mutation (NGG target fully downstream) ======
            elif t_s >= junction_idx:
                req_rha_base = self.args.min_effective_ha + cut_offset
                req_lha_base = self.args.min_effective_ha

                if req_rha_base <= self.args.rha_len and req_lha_base <= self.args.lha_len:
                    pam_rel = pam_s - junction_idx
                    mut_rha, ok, log, mut_idx, _, mut_class = self.perform_precise_mutation(
                        rha_seq, pam_rel, ptype, current_feature, target_strand,
                        genomic_junction, 0, ref_cds_label="Current")
                    if ok:
                        grab_len = min(self.args.rha_len, len(mut_rha))
                        if mut_idx >= grab_len:
                            bump(DesignError.P3B_MUT_IDX_OOR)
                            continue
                        rha_extracted = mut_rha[:grab_len]
                        lha_base, _ = self._optimize_lha_window(
                            upstream_seq, junction_idx,
                            self.args.min_effective_ha, self.args.lha_len)
                        rescue_seq    = rha_seq[:overlap_bases]
                        lha_extracted = (lha_base + rescue_seq) if lha_base else None
                        if lha_extracted and self.validator.is_gc_valid(rha_extracted):
                            candidates.append({
                                'sgRNA': sp_seq, 'PAM': pam_seq,
                                'LHA': lha_extracted, 'RHA': rha_extracted,
                                'Distance': dist, 'Mut_Dist': dist,
                                'Strategy': 'Priority3_Mut_RHA',
                                'Note_Extra': log, 'RHA_Mut_Index': mut_idx,
                                'Mutation_Class': mut_class,
                                'RHA_Buffer': mut_rha,
                            })
                        else:
                            bump(DesignError.P3B_LHA_MISSING if not lha_extracted
                                 else DesignError.P3B_RHA_GC_INVALID)
                    else:
                        bump(DesignError.P3B_MUT_FAILED)
                else:
                    bump(DesignError.P3B_REQ_LEN_EXCEEDS)

        # --- Post-filter (V41) ---
        filtered = []
        for cand in candidates:
            lha_len_c = len(cand['LHA']) if cand['LHA'] else 0
            rha_len_c = len(cand['RHA']) if cand['RHA'] else 0
            eff_lha = lha_len_c
            eff_rha = rha_len_c
            lha_mi  = cand.get('LHA_Mut_Index')
            rha_mi  = cand.get('RHA_Mut_Index')
            if lha_mi is not None:
                eff_lha = lha_mi + 1
            if rha_mi is not None:
                eff_rha = rha_len_c - rha_mi
            if eff_lha < self.args.min_effective_ha or eff_rha < self.args.min_effective_ha:
                bump(DesignError.FILTER_MIN_EFFECTIVE_HA)
                continue
            if lha_len_c > (self.args.lha_len + 20) or rha_len_c > (self.args.rha_len + 20):
                bump(DesignError.FILTER_LENGTH_HARD)
                continue
            filtered.append(cand)

        if not pam_found:
            bump(DesignError.NO_PAM_FOUND)

        self._last_design_reasons = reason_counts if reason_counts is not None else None
        return sorted(filtered, key=lambda x: x['Distance'])

    # ------------------------------------------------------------------
    # C_stop design_target
    # ------------------------------------------------------------------

    def _design_target_cstop(self, upstream_seq, rha_seq,
                              safe_downstream_limit, overlap_bases,
                              neighbor_feature, current_feature,
                              target_strand, genomic_junction):
        candidates   = []
        junction_idx = len(upstream_seq)
        full         = upstream_seq + rha_seq
        search_s     = max(0, junction_idx - self.args.search_window)
        search_e     = junction_idx + self.args.search_window

        rha_gj = (int(current_feature.location.end)
                  if target_strand == 1
                  else int(current_feature.location.start))

        for i in range(search_s, search_e):
            if i + 23 > len(full):
                break
            ptype = None
            if full[i:i+2] == "GG":
                ptype = 'NGG'
            elif full[i:i+2] == "CC":
                ptype = 'CCN'
            if ptype is None:
                continue

            if ptype == 'NGG':
                pam_s    = i - 1
                sp_e     = pam_s
                sp_s     = sp_e - 20
                cut      = pam_s - 3
                sp_seq   = full[sp_s:sp_e]
                pam_seq  = full[pam_s:pam_s+3]
                t_s, t_e = sp_s, pam_s + 3
            else:
                pam_s    = i
                sp_s     = pam_s + 3
                sp_e     = sp_s + 20
                cut      = pam_s + 6
                t_seq    = full[sp_s:sp_e]
                sp_seq   = str(Seq(t_seq).reverse_complement())
                pam_seq  = full[pam_s:pam_s+3]
                t_s, t_e = pam_s, sp_e

            if t_s < 0 or t_e > len(full):
                continue
            dist       = abs(cut - junction_idx)
            cut_offset = cut - junction_idx

            if self.args.target_gene:
                print(f"[DEBUG C_stop] pam={pam_s} ({ptype}) dist={dist} "
                      f"t_s={t_s} t_e={t_e} junc={junction_idx}")

            if not self.validator.is_gc_valid(sp_seq, is_sgrna=True):
                continue

            # P1 — target fully downstream
            if t_s >= junction_idx:
                del_len = t_e - junction_idx
                if (del_len >= 3
                        and del_len <= safe_downstream_limit
                        and del_len <= self.args.max_deletion):
                    lha_fin, _shift = self._optimize_lha_window(
                        upstream_seq, junction_idx,
                        self.args.min_effective_ha, self.args.lha_len)
                    (rha_rescue, _rha_total, rha_prefix_len,
                     _rha_homology) = self._build_cstop_rha_with_overlap_rescue(
                        upstream_seq, rha_seq, overlap_bases, self.args.rha_len)
                    homology_room = max(0, self.args.rha_len - rha_prefix_len)
                    rha_core      = rha_seq[del_len : del_len + homology_room]
                    rha_extracted = rha_rescue + rha_core
                    if (lha_fin and rha_extracted
                            and len(rha_core) >= self.args.min_effective_ha
                            and self.validator.is_gc_valid(rha_extracted)):
                        note_extra = f"CStop Del: {del_len}bp (incl. stop codon)"
                        if overlap_bases > 0:
                            note_extra += f" | NeighborCDS_Rescue:{overlap_bases}bp"
                        candidates.append({
                            'sgRNA':      sp_seq, 'PAM': pam_seq,
                            'LHA':        lha_fin, 'RHA': rha_extracted,
                            'Distance':   dist, 'Mut_Dist': 0,
                            'Strategy':   'CStop_P1_Del_Downstream',
                            'Note_Extra': note_extra,
                            'Mutation_Class': 'None',
                            'RHA_Prefix_Len': rha_prefix_len,
                            'RHA_Junction_Idx': rha_prefix_len,
                        })

            # P2 — target spans junction
            elif t_s < junction_idx < t_e:
                req_lha = (self.args.min_effective_ha
                           + (abs(cut_offset) if cut_offset < 0 else 0))
                rha_penalty = max(0, cut_offset - 3)
                req_rha     = self.args.min_effective_ha + rha_penalty
                (rha_rescue, _rha_total, rha_prefix_len,
                 _rha_homology) = self._build_cstop_rha_with_overlap_rescue(
                    upstream_seq, rha_seq, overlap_bases, self.args.rha_len)
                if (req_lha > self.args.lha_len
                        or req_rha + rha_prefix_len > self.args.rha_len):
                    continue

                grab_lha_len  = min(self.args.lha_len, junction_idx)
                raw_lha_part  = upstream_seq[junction_idx - grab_lha_len : junction_idx]
                lha_cand_full = raw_lha_part

                homology_room = max(0, self.args.rha_len - rha_prefix_len)
                raw_rha_core  = rha_seq[3 : 3 + homology_room]
                raw_rha       = rha_rescue + raw_rha_core

                eff_rha_check = len(raw_rha_core) - rha_penalty
                if (eff_rha_check < self.args.min_effective_ha
                        or not self.validator.is_gc_valid(raw_rha)):
                    continue

                if ptype == 'NGG':
                    physical_target = full[sp_s : pam_s + 3]
                    seed_seq        = physical_target[-(self.args.seed_region + 3):]
                else:
                    physical_target = full[pam_s : sp_e]
                    seed_seq        = physical_target[:(3 + self.args.seed_region)]

                risk_in_lha = ((physical_target in lha_cand_full)
                               or (seed_seq in lha_cand_full))
                risk_in_rha = ((physical_target in raw_rha)
                               or (seed_seq in raw_rha))

                lha_final = rha_final = None
                lha_mut_idx = rha_mut_idx = None
                note     = ""
                mut_dist = 0
                is_mut   = False
                mutation_class = 'None'
                if overlap_bases > 0:
                    note += f"NeighborCDS_Rescue:{overlap_bases}bp | "

                if risk_in_rha:
                    pam_in_rha = self._map_cstop_full_idx_to_rha(
                        pam_s, junction_idx, overlap_bases)
                    if 0 <= pam_in_rha < len(raw_rha):
                        mut_rha, ok, log, mut_idx, loc_lbl, mut_class = self.perform_precise_mutation(
                            raw_rha, pam_in_rha, ptype,
                            neighbor_feature, target_strand, rha_gj, rha_prefix_len,
                            ref_cds_label="Neighbor")
                        if ok:
                            rha_final   = mut_rha
                            note       += f"RHA Mut: {log} [{loc_lbl}] | "
                            is_mut      = True
                            rha_mut_idx = mut_idx
                            mutation_class = _merge_mutation_class(mutation_class, mut_class)
                        else:
                            rha_final = None
                    else:
                        rha_final = raw_rha
                else:
                    rha_final = raw_rha

                if rha_final and risk_in_lha:
                    start_of_lha = junction_idx - grab_lha_len
                    pam_in_lha   = pam_s - start_of_lha
                    if 0 <= pam_in_lha < len(lha_cand_full):
                        mut_lha, ok, log, mut_idx, loc_lbl, mut_class = self.perform_precise_mutation(
                            lha_cand_full, pam_in_lha, ptype,
                            current_feature, target_strand,
                            genomic_junction, len(raw_lha_part),
                            ref_cds_label="Current")
                        if ok:
                            req_opt          = req_lha
                            lha_final, shift = self._optimize_lha_window(
                                mut_lha, len(mut_lha), req_opt)
                            if lha_final is not None:
                                lha_mut_idx = mut_idx - shift
                                if lha_mut_idx < 0:
                                    lha_final = None
                            note    += f"LHA Mut: {log} [{loc_lbl}]"
                            is_mut   = True
                            mut_dist = len(mut_lha) - mut_idx
                            mutation_class = _merge_mutation_class(mutation_class, mut_class)
                        else:
                            lha_final = None
                    else:
                        lha_final = None
                elif rha_final:
                    req_opt          = req_lha
                    lha_final, shift = self._optimize_lha_window(
                        lha_cand_full, len(lha_cand_full), req_opt)
                    if shift > 0:
                        note += f"LHA GC Shift -{shift}bp"

                if lha_final and rha_final and len(lha_final) <= self.args.lha_len:
                    stgy       = 'CStop_P2_Bridge_Mut' if is_mut else 'CStop_P2_Bridge'
                    if overlap_bases >= 4:
                        note += f" [LargeOverlap:{overlap_bases}bp]"
                    note_clean = note.strip(" | ")
                    cand = {
                        'sgRNA': sp_seq, 'PAM': pam_seq,
                        'LHA': lha_final, 'RHA': rha_final,
                        'Distance': dist, 'Strategy': stgy,
                        'Note_Extra': note_clean, 'Mut_Dist': mut_dist,
                        'Mutation_Class': mutation_class,
                        'RHA_Prefix_Len': rha_prefix_len,
                        'RHA_Junction_Idx': rha_prefix_len,
                    }
                    if is_mut:
                        if 'LHA Mut' in note and lha_mut_idx is not None:
                            cand['LHA_Mut_Index'] = lha_mut_idx
                        if 'RHA Mut' in note and rha_mut_idx is not None:
                            cand['RHA_Mut_Index'] = rha_mut_idx
                    candidates.append(cand)

            # P3 — target fully upstream (coding before stop codon)
            elif t_e <= junction_idx:
                req_lha = (self.args.min_effective_ha
                           + (abs(cut_offset) if cut_offset < 0 else 0))
                req_rha = self.args.min_effective_ha
                (rha_rescue, rha_extracted, rha_prefix_len,
                 rha_homology) = self._build_cstop_rha_with_overlap_rescue(
                    upstream_seq, rha_seq, overlap_bases, self.args.rha_len)
                if (req_lha > self.args.lha_len
                        or req_rha + rha_prefix_len > self.args.rha_len):
                    continue

                mut_upstream, ok, log, mut_idx, _, mut_class = self.perform_precise_mutation(
                    upstream_seq, pam_s, ptype,
                    current_feature, target_strand,
                    genomic_junction, junction_idx,
                    ref_cds_label="Current")
                if not ok:
                    continue

                lha_base, _gc_shift = self._optimize_lha_window(
                    mut_upstream, len(mut_upstream),
                    self.args.min_effective_ha, self.args.lha_len)
                if lha_base is None:
                    continue

                total_trim  = len(mut_upstream) - len(lha_base)
                new_mut_idx = mut_idx - total_trim
                if new_mut_idx < 0:
                    continue

                lha_extracted = lha_base

                if (lha_extracted and rha_extracted
                        and len(rha_homology) >= self.args.min_effective_ha
                        and self.validator.is_gc_valid(rha_extracted)):
                    if overlap_bases > 0:
                        log += f" | NeighborCDS_Rescue:{overlap_bases}bp"
                    if overlap_bases >= 4:
                        log += f" [LargeOverlap:{overlap_bases}bp]"
                    candidates.append({
                        'sgRNA': sp_seq, 'PAM': pam_seq,
                        'LHA': lha_extracted, 'RHA': rha_extracted,
                        'Distance': dist, 'Mut_Dist': dist,
                        'Strategy': 'CStop_P3_Mut_LHA',
                        'Note_Extra': log, 'LHA_Mut_Index': new_mut_idx,
                        'Mutation_Class': mut_class,
                        'RHA_Prefix_Len': rha_prefix_len,
                        'RHA_Junction_Idx': rha_prefix_len,
                    })

        return candidates

    # ------------------------------------------------------------------
    # _select_top_designs() — dispatch by mode
    # ------------------------------------------------------------------

    def _select_top_designs(self, candidates, ctx, feature, genomic_junction,
                            num_designs=2):
        if self.model == 'N_start':
            return self._select_top_designs_nstart(
                candidates, ctx, feature, genomic_junction, num_designs)
        return self._select_top_designs_cstop(
            candidates, ctx, feature, genomic_junction, num_designs)

    # ------------------------------------------------------------------
    # N_start _select_top_designs (V43 stages 0-8 + V44 sim + V46 off-target)
    # ------------------------------------------------------------------

    def _select_top_designs_nstart(self, candidates, ctx, feature,
                                    genomic_junction, num_designs=2):
        budget = self.ha_budget

        # ── Stage 0a: Off-target specificity filter (V46) ──
        threshold = self.args.max_offtargets + 1
        filtered_ot = []
        for cand in candidates:
            sgrna    = cand.get('sgRNA', '').upper()
            ot_count = self.offtarget_index.get(sgrna, 0)
            if ot_count > threshold:
                self.stats['sgrna_offtarget'] = self.stats.get('sgrna_offtarget', 0) + 1
                continue
            cand_a = dict(cand)
            cand_a['OffTarget_Count'] = ot_count
            filtered_ot.append(cand_a)

        # ── Stage 0b: Start codon gate (V41) ──
        expected_sc = getattr(self, '_current_start_codon', None)
        if expected_sc and expected_sc in VALID_START_CODONS:
            filtered_ot = [c for c in filtered_ot
                           if c.get('RHA', '')[:3].upper() == expected_sc]

        def sort_key(c):
            is_safe     = c['Mut_Dist'] <= self.args.max_mut_dist
            strat       = c['Strategy']
            strat_score = 3 if 'Priority1' in strat else (2 if 'Priority2' in strat else 1)
            return (is_safe, strat_score, -c['Distance'])

        used_sgrnas   = set()
        accepted_list = []

        for cand in sorted(filtered_ot, key=sort_key, reverse=True):
            if len(accepted_list) >= num_designs:
                break
            sgrna = cand.get('sgRNA', '')
            if sgrna in used_sgrnas:
                continue

            # ── Stage 1: Hard-reject sgRNA containing BsaI/BbsI ──
            if self._seq_contains_re(sgrna):
                continue

            cand_use = dict(cand)

            # ── Stage 2: Sanitize LHA/RHA ──
            if self.args.restriction_site:
                lha_seq      = cand_use.get('LHA')
                rha_seq_cand = cand_use.get('RHA')
                if lha_seq:
                    lha_fixed, lha_logs, lha_changed = self._sanitize_arm(
                        lha_seq, 'LHA', ctx['up_feat'],
                        feature.location.strand, genomic_junction, len(lha_seq))
                else:
                    lha_fixed, lha_logs, lha_changed = lha_seq, [], False
                if rha_seq_cand:
                    rha_fixed, rha_logs, rha_changed = self._sanitize_arm(
                        rha_seq_cand, 'RHA', feature,
                        feature.location.strand, genomic_junction, 0)
                else:
                    rha_fixed, rha_logs, rha_changed = rha_seq_cand, [], False
                if lha_changed or rha_changed:
                    if lha_fixed and not self._is_gc_valid(lha_fixed):
                        lha_fixed = lha_seq
                        lha_logs.append("LHA:FixReverted_GC")
                    if rha_fixed and not self._is_gc_valid(rha_fixed):
                        rha_fixed = rha_seq_cand
                        rha_logs.append("RHA:FixReverted_GC")
                    cand_use['LHA'] = lha_fixed
                    cand_use['RHA'] = rha_fixed
                    logs = lha_logs + rha_logs
                    if logs:
                        extra = cand_use.get('Note_Extra', '')
                        cand_use['Note_Extra'] = ((extra + " | " + ";".join(logs))
                                                   if extra else ";".join(logs))

            # ── Stage 3: Hard-reject if LHA/RHA still contains RE ──
            if self._seq_contains_re(cand_use.get('LHA', '')):
                continue
            if self._seq_contains_re(cand_use.get('RHA', '')):
                continue

            # ── Stage 4: Mutation-aware balanced HA trim (V42) ──
            lha_raw   = cand_use.get('LHA', '')
            rha_raw   = cand_use.get('RHA', '')
            lha_mut_i = cand_use.get('LHA_Mut_Index')
            rha_mut_i = cand_use.get('RHA_Mut_Index')
            trim_result = self._mutation_aware_trim(lha_raw, rha_raw, lha_mut_i, rha_mut_i)
            if trim_result is None:
                continue
            lha_t, rha_t, new_lha_mut_i, eff_lha, eff_rha = trim_result

            # ── Stage 4.5: RHA budget-fill padding (V43) ──
            rha_padded = False
            rha_buf    = cand_use.get('RHA_Buffer')
            if rha_buf and len(lha_t) + len(rha_t) < budget:
                desired = budget - len(lha_t)
                if desired <= len(rha_buf):
                    if rha_mut_i is None or rha_mut_i < desired:
                        ext_rha = rha_buf[:desired]
                        if self.validator.is_gc_valid(ext_rha):
                            rha_t      = ext_rha
                            eff_rha    = ((len(rha_t) - rha_mut_i)
                                          if rha_mut_i is not None else len(rha_t))
                            rha_padded = True

            note_parts = []
            if len(lha_t) != len(lha_raw) or len(rha_t) != len(rha_raw):
                note_parts.append(
                    f"Trimmed:LHA{len(lha_raw)}->{len(lha_t)},"
                    f"RHA{len(rha_raw)}->{len(rha_t)}")
            elif rha_padded:
                note_parts.append(
                    f"Padded:LHA{len(lha_raw)},"
                    f"RHA{len(cand_use.get('RHA',''))}->{len(rha_t)}")
            if lha_mut_i is not None or rha_mut_i is not None:
                note_parts.append(f"Eff:LHA={eff_lha}bp,RHA={eff_rha}bp")
            if note_parts:
                addition = " | ".join(note_parts)
                extra    = cand_use.get('Note_Extra', '')
                cand_use['Note_Extra'] = ((extra + " | " + addition) if extra else addition)

            cand_use['LHA'] = lha_t
            cand_use['RHA'] = rha_t
            if new_lha_mut_i is not None:
                cand_use['LHA_Mut_Index'] = new_lha_mut_i

            # ── Stage 5: Assemble final oligo ──
            (final_oligo, barcode, risky, msg,
             re_details, re_types, restriction_sites) = self.assemble_final_oligo(cand_use)

            # ── Stage 6: Full-oligo non-backbone RE violation scan ──
            if self._find_variable_re_violations(final_oligo):
                continue

            # ── Stage 7: Hard total-length check ──
            if len(final_oligo) > self.args.max_oligo_len:
                continue

            # ── Stage 8: Double-RE-type filter ──
            types_set = {t for t in re_types.split(';') if t and t != "None"}
            if len(types_set) >= 2:
                continue

            used_sgrnas.add(sgrna)
            accepted_list.append(
                (cand_use, final_oligo, barcode, risky, msg,
                 re_details, re_types, restriction_sites))

        # V44: similarity flagging
        if len(accepted_list) >= 2:
            lha1 = accepted_list[0][0].get('LHA', '')
            rha1 = accepted_list[0][0].get('RHA', '')
            lha2 = accepted_list[1][0].get('LHA', '')
            rha2 = accepted_list[1][0].get('RHA', '')
            sim_lha   = round(_seq_sim_right(lha1, lha2) * 100)
            sim_rha   = round(_seq_sim_left(rha1, rha2) * 100)
            threshold = getattr(self.args, 'sim_threshold', 85)
            if sim_lha >= threshold and sim_rha >= threshold:
                tag = f"[SIM:LHA={sim_lha}%,RHA={sim_rha}%]"
                for idx in range(2):
                    c     = accepted_list[idx][0]
                    extra = c.get('Note_Extra', '')
                    c['Note_Extra'] = (extra + ' ' + tag).strip() if extra else tag

        return accepted_list

    # ------------------------------------------------------------------
    # C_stop _select_top_designs (with off-target + stages 1-8)
    # ------------------------------------------------------------------

    def _select_top_designs_cstop(self, candidates, ctx, feature,
                                   genomic_junction, num_designs=2):
        """V3: Option-C (strategy diversity) + Option-A (sim threshold) + fallback."""
        budget        = self.ha_budget
        rha_gj        = (int(feature.location.end)
                          if feature.location.strand == 1
                          else int(feature.location.start))
        rank2_sim_max = getattr(self.args, 'rank2_sim_max', 50)

        def sort_key(c):
            bio_rank   = self._candidate_bio_rank(c)
            dist_rank  = 0 if c['Mut_Dist'] <= self.args.max_mut_dist else 1
            strat      = c['Strategy']
            strat_rank = 0 if 'P1' in strat else (1 if 'P2' in strat else 2)
            return (bio_rank, dist_rank, strat_rank, -c['Distance'])

        def strat_class(strat):
            """Return coarse strategy class: 'P1', 'P2', or 'P3'."""
            if 'P1' in strat:
                return 'P1'
            if 'P2' in strat:
                return 'P2'
            return 'P3'

        # ── Stage 0: Off-target filter ──
        ot_threshold = self.args.max_offtargets + 1
        filtered_candidates = []
        for cand in candidates:
            sgrna    = cand.get('sgRNA', '').upper()
            ot_count = self.offtarget_index.get(sgrna, 0)
            if ot_count > ot_threshold:
                self.stats['sgrna_offtarget'] = self.stats.get('sgrna_offtarget', 0) + 1
                continue
            cand_a = dict(cand)
            cand_a['OffTarget_Count'] = ot_count
            filtered_candidates.append(cand_a)

        sorted_cands = sorted(filtered_candidates, key=sort_key)
        used_sgrnas  = set()

        # ── Nested helper: run Stages 1-8 on a candidate (non-destructive to dict) ──
        # Returns validated tuple or None on rejection. Does NOT mutate used_sgrnas.
        def try_validate(cand):
            sgrna = cand.get('sgRNA', '')
            if sgrna in used_sgrnas:
                return None
            # Stage 1
            if self._seq_contains_re(sgrna):
                return None

            cand_use = dict(cand)

            # Stage 2: Sanitize LHA / RHA
            if self.args.restriction_site:
                lha_seq      = cand_use.get('LHA')
                rha_seq_cand = cand_use.get('RHA')
                rha_junc_idx = cand_use.get('RHA_Junction_Idx', 0)
                if lha_seq:
                    lha_fixed, lha_logs, lha_changed = self._sanitize_arm(
                        lha_seq, 'LHA', feature,
                        feature.location.strand, genomic_junction, len(lha_seq))
                else:
                    lha_fixed, lha_logs, lha_changed = lha_seq, [], False
                if rha_seq_cand:
                    rha_fixed, rha_logs, rha_changed = self._sanitize_arm(
                        rha_seq_cand, 'RHA', ctx['neighbor_feat'],
                        feature.location.strand, rha_gj, rha_junc_idx)
                else:
                    rha_fixed, rha_logs, rha_changed = rha_seq_cand, [], False
                if lha_changed or rha_changed:
                    if lha_fixed and not self._is_gc_valid(lha_fixed):
                        lha_fixed = lha_seq
                        lha_logs.append("LHA:FixReverted_GC")
                    if rha_fixed and not self._is_gc_valid(rha_fixed):
                        rha_fixed = rha_seq_cand
                        rha_logs.append("RHA:FixReverted_GC")
                    cand_use['LHA'] = lha_fixed
                    cand_use['RHA'] = rha_fixed
                    logs = lha_logs + rha_logs
                    if logs:
                        extra = cand_use.get('Note_Extra', '')
                        cand_use['Note_Extra'] = ((extra + " | " + ";".join(logs))
                                                   if extra else ";".join(logs))

            # Stage 3
            if self._seq_contains_re(cand_use.get('LHA', '')):
                return None
            if self._seq_contains_re(cand_use.get('RHA', '')):
                return None

            # Stage 4: Mutation-aware HA trim
            lha_raw   = cand_use.get('LHA', '')
            rha_raw   = cand_use.get('RHA', '')
            lha_mut_i = cand_use.get('LHA_Mut_Index')
            rha_mut_i = cand_use.get('RHA_Mut_Index')
            rha_prefix_len = cand_use.get('RHA_Prefix_Len', 0)
            trim_result = self._mutation_aware_trim(
                lha_raw, rha_raw, lha_mut_i, rha_mut_i, rha_prefix_len)
            if trim_result is None:
                return None
            lha_t, rha_t, new_lha_mut_i, eff_lha, eff_rha = trim_result

            # Stage 4.5: RHA budget-fill
            rha_padded = False
            rha_buf    = cand_use.get('RHA_Buffer')
            if rha_buf and len(lha_t) + len(rha_t) < budget:
                desired = budget - len(lha_t)
                if desired <= len(rha_buf):
                    if rha_mut_i is None or rha_mut_i < desired:
                        ext_rha = rha_buf[:desired]
                        if self.validator.is_gc_valid(ext_rha):
                            rha_t      = ext_rha
                            eff_rha    = ((len(rha_t) - rha_mut_i)
                                          if rha_mut_i is not None else len(rha_t))
                            rha_padded = True

            note_parts = []
            if len(lha_t) != len(lha_raw) or len(rha_t) != len(rha_raw):
                note_parts.append(
                    f"Trimmed:LHA{len(lha_raw)}->{len(lha_t)},"
                    f"RHA{len(rha_raw)}->{len(rha_t)}")
            elif rha_padded:
                note_parts.append(
                    f"Padded:LHA{len(lha_raw)},"
                    f"RHA{len(cand_use.get('RHA',''))}->{len(rha_t)}")
            if lha_mut_i is not None or rha_mut_i is not None:
                note_parts.append(f"Eff:LHA={eff_lha}bp,RHA={eff_rha}bp")
            if note_parts:
                addition = " | ".join(note_parts)
                extra    = cand_use.get('Note_Extra', '')
                cand_use['Note_Extra'] = ((extra + " | " + addition) if extra else addition)

            cand_use['LHA'] = lha_t
            cand_use['RHA'] = rha_t
            if new_lha_mut_i is not None:
                cand_use['LHA_Mut_Index'] = new_lha_mut_i

            # Stage 5
            (final_oligo, barcode, risky, msg,
             re_details, re_types, restriction_sites) = self.assemble_final_oligo(cand_use)

            # Stage 6
            if self._find_variable_re_violations(final_oligo):
                return None
            # Stage 7
            if len(final_oligo) > self.args.max_oligo_len:
                return None
            # Stage 8
            types_set = {t for t in re_types.split(';') if t and t != "None"}
            if len(types_set) >= 2:
                return None

            return (cand_use, final_oligo, barcode, risky, msg,
                    re_details, re_types, restriction_sites)

        # ── Find Rank1: best valid candidate from sorted pool ──
        rank1_result = None
        rank1_idx    = -1
        for idx, cand in enumerate(sorted_cands):
            result = try_validate(cand)
            if result is not None:
                used_sgrnas.add(cand.get('sgRNA', ''))
                rank1_result = result
                rank1_idx    = idx
                break

        if rank1_result is None:
            return []  # No valid design at all

        accepted_list = [rank1_result]
        if num_designs == 1:
            return accepted_list

        # ── Find Rank2 with biological filtering first, then diversity logic ──
        rank1_cand       = sorted_cands[rank1_idx]
        rank1_sc         = strat_class(rank1_cand.get('Strategy', ''))
        rank1_lha        = rank1_result[0].get('LHA', '')
        rank1_rha        = rank1_result[0].get('RHA', '')
        remaining_cands  = sorted_cands[rank1_idx + 1:]

        rank2_result  = None
        diversity_tag = ''

        validated_remaining = []
        for cand in remaining_cands:
            result = try_validate(cand)
            if result is None:
                continue
            validated_remaining.append((cand, result, self._candidate_bio_rank(result[0])))

        for bio_rank in sorted({entry[2] for entry in validated_remaining}):
            tier_entries = [entry for entry in validated_remaining if entry[2] == bio_rank]
            diff_strat_pool = [entry for entry in tier_entries
                               if strat_class(entry[0].get('Strategy', '')) != rank1_sc]
            same_strat_pool = [entry for entry in tier_entries
                               if strat_class(entry[0].get('Strategy', '')) == rank1_sc]

            for cand, result, _bio in diff_strat_pool:
                used_sgrnas.add(cand.get('sgRNA', ''))
                rank2_result  = result
                diversity_tag = '[Diverse:StrategyDiff]'
                break

            if rank2_result is not None:
                break

            best_fallback = None
            for cand, result, _bio in same_strat_pool:
                lha2    = result[0].get('LHA', '')
                rha2    = result[0].get('RHA', '')
                sim_lha = round(_seq_sim_right(rank1_lha, lha2) * 100)
                sim_rha = round(_seq_sim_left(rank1_rha,  rha2) * 100)
                sim_max = max(sim_lha, sim_rha)
                if best_fallback is None or sim_max < best_fallback[0]:
                    best_fallback = (sim_max, cand, result)
                if sim_max < rank2_sim_max:
                    used_sgrnas.add(cand.get('sgRNA', ''))
                    rank2_result  = result
                    diversity_tag = f'[Diverse:SeqDiff,SIM_max={sim_max}%]'
                    break

            if rank2_result is not None:
                break

            if best_fallback is not None:
                sim_max, fb_cand, fb_result = best_fallback
                used_sgrnas.add(fb_cand.get('sgRNA', ''))
                rank2_result  = fb_result
                diversity_tag = f'[FallbackDiverse:SIM_max={sim_max}%]'
                break

        if rank2_result is not None:
            # Append diversity tag to Rank2's Note_Extra
            if diversity_tag:
                c_dict = rank2_result[0]
                extra  = c_dict.get('Note_Extra', '')
                c_dict['Note_Extra'] = ((extra + ' ' + diversity_tag).strip()
                                        if extra else diversity_tag)
            accepted_list.append(rank2_result)

        return accepted_list

    # ------------------------------------------------------------------
    # run() — dispatch by mode, with off-target index building
    # ------------------------------------------------------------------

    def run(self):
        if self.model == 'N_start':
            self._run_nstart()
        else:
            self._run_cstop()

    def _run_nstart(self):
        """N_start mode: build off-target index, run V38-style pipeline, inject OffTarget_Count."""
        print(f"{Colors.HEADER}[INFO] N_start mode — parsing GenBank for off-target index...{Colors.ENDC}")
        try:
            records = list(SeqIO.parse(self.args.gbff, "genbank"))
            if not records:
                raise ValueError("No records found in GenBank file.")
        except Exception as e:
            print(f"{Colors.FAIL}[Error] {e}{Colors.ENDC}")
            sys.exit(1)

        genome_str = "".join(str(r.seq) for r in records)
        self._build_and_print_index(genome_str)

        num_designs = getattr(self.args, 'num_designs', 2)

        for record in records:
            self.current_record_seq = record.seq
            cds_features = [f for f in record.features if f.type == "CDS"]
            cds_features.sort(key=lambda f: f.location.start)

            # Build inter-gene context
            gene_context = {}
            for i, f in enumerate(cds_features):
                locus  = f.qualifiers.get('locus_tag', f.qualifiers.get('gene', ['Unknown']))[0]
                strand = f.location.strand
                safe, ov, up_feat = 100, 0, None
                if strand == 1:
                    if i > 0:
                        prev    = cds_features[i - 1]
                        up_feat = prev
                        d = int(f.location.start) - int(prev.location.end)
                        safe = max(0, d)
                        ov   = max(0, -d)
                else:
                    if i < len(cds_features) - 1:
                        nxt     = cds_features[i + 1]
                        up_feat = nxt
                        d = int(nxt.location.start) - int(f.location.end)
                        safe = max(0, d)
                        ov   = max(0, -d)
                gene_context[locus] = {'safe': safe, 'overlap': ov, 'up_feat': up_feat}

            for feature in cds_features:
                qualifiers = feature.qualifiers
                gene_id    = qualifiers.get('locus_tag',
                              qualifiers.get('gene', ['Unknown']))[0]
                if self.args.target_gene and gene_id != self.args.target_gene:
                    continue

                ok, msg = self._validate_cds_boundary(record, feature)
                if not ok:
                    self.results.append({
                        'Design_Rank': 1, 'Gene_ID': gene_id,
                        'Note': f"Skipped: CDS boundary error ({msg})",
                    })
                    self.stats["failed"] += 1
                    continue

                product = qualifiers.get('product', [''])[0]
                self.stats["total"] += 1
                ctx = gene_context.get(gene_id, {'safe': 100, 'overlap': 0, 'up_feat': None})
                genomic_junction = (int(feature.location.start)
                                    if feature.location.strand == 1
                                    else int(feature.location.end))

                try:
                    up, rha = self.get_upstream_and_rha(record, feature)
                    candidates = self.design_target(
                        up, rha, ctx['safe'], ctx['overlap'],
                        ctx['up_feat'], feature, feature.location.strand, genomic_junction)
                except Exception as exc:
                    self.stats["error"] += 1
                    print(f"{Colors.FAIL}[Error] {gene_id}: {exc}{Colors.ENDC}")
                    continue

                if not candidates:
                    self.results.append({
                        'Design_Rank': 1, 'Gene_ID': gene_id, 'Note': 'Failed',
                    })
                    self.stats["failed"] += 1
                    continue

                accepted_list = self._select_top_designs(
                    candidates, ctx, feature, genomic_junction, num_designs=num_designs)

                if not accepted_list:
                    self.results.append({
                        'Design_Rank': 1, 'Gene_ID': gene_id,
                        'Note': 'Failed (Double Restriction Sites)',
                    })
                    self.stats["failed"] += 1
                    continue

                # Stats for rank-1
                best_cand = accepted_list[0][0]
                strat_r1  = best_cand['Strategy']
                if best_cand['Mut_Dist'] <= self.args.max_mut_dist:
                    self.stats["optimal_selected"] += 1
                else:
                    self.stats["suboptimal_selected"] += 1
                if self._candidate_has_aa_change(best_cand):
                    self.stats["aa_change_selected"] += 1

                if   strat_r1 == 'Priority1_Clean':           self.stats["p1_clean"] += 1
                elif 'Priority1_Partial' in strat_r1:         self.stats["p1_flexible"] += 1
                elif strat_r1 == 'Priority2_Bridge':          self.stats["p2_bridge_nomut"] += 1
                elif strat_r1 == 'Priority2_Bridge_Mut':      self.stats["p2_bridge_mut"] += 1
                elif strat_r1 == 'Priority3_Mut_LHA':         self.stats["p3a_lha"] += 1
                elif strat_r1 == 'Priority3_Mut_RHA':         self.stats["p3b_rha"] += 1

                # Append rows
                for rank, (cand_use, final_oligo, barcode, risky, msg,
                            re_details, re_types, restriction_sites) \
                        in enumerate(accepted_list, 1):
                    strat  = cand_use['Strategy']
                    status = self._quality_status(cand_use)
                    aa_tag = self._aa_cost_tag(cand_use)
                    note_extra = " | ".join(
                        part for part in [cand_use.get('Note_Extra', '').strip(' | '), aa_tag, status] if part
                    )
                    if risky:
                        print(f"{Colors.FAIL}   [RISK] {gene_id:<15} Rank{rank} | {msg}{Colors.ENDC}")
                    self.results.append({
                        'Design_Rank':       rank,
                        'Gene_ID':           gene_id,
                        'Product':           product,
                        'Strategy':          strat,
                        'Info':              note_extra,
                        'sgRNA':             cand_use['sgRNA'],
                        'sgRNA_GC':          self._calc_gc(cand_use['sgRNA']),
                        'PAM':               cand_use['PAM'],
                        'LHA':               cand_use['LHA'],
                        'LHA_Length':        len(cand_use['LHA']),
                        'LHA_GC':            self._calc_gc(cand_use['LHA']),
                        'RHA':               cand_use['RHA'],
                        'RHA_Length':        len(cand_use['RHA']),
                        'RHA_GC':            self._calc_gc(cand_use['RHA']),
                        'Restriction_Sites': restriction_sites,
                        'Mut_Dist':          cand_use['Mut_Dist'],
                        'Barcode':           barcode,
                        'Final_Oligo':       final_oligo,
                        'Note':              f"Success ({strat})",
                    })

        # Save CSV
        df = pd.DataFrame(self.results)
        cols = [
            'Design_Rank', 'Gene_ID', 'Product', 'Strategy', 'Info',
            'sgRNA', 'sgRNA_GC', 'PAM',
            'LHA', 'LHA_Length', 'LHA_GC',
            'RHA', 'RHA_Length', 'RHA_GC',
            'Restriction_Sites', 'Mut_Dist', 'Barcode', 'Final_Oligo', 'Note',
        ]
        final_cols = [c for c in cols if c in df.columns]
        df[final_cols].to_csv(self.args.output, index=False)

        # Summary
        total_success = self.stats['optimal_selected'] + self.stats['suboptimal_selected']
        total_p1 = self.stats['p1_clean'] + self.stats['p1_flexible']
        total_p2 = self.stats['p2_bridge_nomut'] + self.stats['p2_bridge_mut']
        total_p3 = self.stats['p3a_lha'] + self.stats['p3b_rha']
        rank2_count = (len(df[df['Design_Rank'] == 2]) if 'Design_Rank' in df.columns else 0)

        print(f"\n{Colors.HEADER}{'=' * 70}")
        print(f" Summary Statistics (N_start mode)")
        print(f"{'=' * 70}{Colors.ENDC}")
        print(f" Total Genes Processed  : {self.stats['total']}")
        if self.stats['total']:
            print(f" Genes with ≥1 Design   : {total_success} ({total_success / self.stats['total'] * 100:.1f}%)")
            print(f" Genes with 2nd Design  : {rank2_count}")
            print(f" Failed Designs         : {self.stats['failed']} ({self.stats['failed'] / self.stats['total'] * 100:.1f}%)")
        print(f"\n{Colors.OKBLUE} Strategy Distribution (Rank-1):{Colors.ENDC}")
        if total_success > 0:
            print(f" ├─ Priority 1 (Deletion)      : {total_p1:4d} ({total_p1 / total_success * 100:5.1f}%)")
            print(f" │  ├─ P1_Clean (Full Delete)  : {self.stats['p1_clean']:4d} ({self.stats['p1_clean'] / total_success * 100:5.1f}%)")
            print(f" │  └─ P1_Partial (Opt)        : {self.stats['p1_flexible']:4d} ({self.stats['p1_flexible'] / total_success * 100:5.1f}%)")
            print(f" ├─ Priority 2 (Bridge)        : {total_p2:4d} ({total_p2 / total_success * 100:5.1f}%)")
            print(f" │  ├─ P2_Bridge (No Mut)      : {self.stats['p2_bridge_nomut']:4d} ({self.stats['p2_bridge_nomut'] / total_success * 100:5.1f}%)")
            print(f" │  └─ P2_Bridge_Mut           : {self.stats['p2_bridge_mut']:4d} ({self.stats['p2_bridge_mut'] / total_success * 100:5.1f}%)")
            print(f" └─ Priority 3 (Mutation)      : {total_p3:4d} ({total_p3 / total_success * 100:5.1f}%)")
            print(f"    ├─ P3A_Mut_LHA             : {self.stats['p3a_lha']:4d} ({self.stats['p3a_lha'] / total_success * 100:5.1f}%)")
            print(f"    └─ P3B_Mut_RHA             : {self.stats['p3b_rha']:4d} ({self.stats['p3b_rha'] / total_success * 100:5.1f}%)")
            opt = self.stats['optimal_selected']
            sub = self.stats['suboptimal_selected']
            aa_change = self.stats['aa_change_selected']
            print(f"\n{Colors.OKGREEN} Distance Metrics:{Colors.ENDC}")
            print(f" ├─ Optimal  (Mut_Dist ≤ {self.args.max_mut_dist}bp) : {opt:4d} ({opt / total_success * 100:5.1f}%)")
            print(f" └─ Suboptimal (Mut_Dist > {self.args.max_mut_dist}bp): {sub:4d} ({sub / total_success * 100:5.1f}%)")
            print(f"{Colors.OKGREEN} Amino-acid Impact:{Colors.ENDC}")
            print(f" └─ Rank-1 conservative fallback designs            : {aa_change:4d} ({aa_change / total_success * 100:5.1f}%)")

        # Inject OffTarget_Count column
        self._inject_offtarget_column()
        self._print_offtarget_summary()

    def _run_cstop(self):
        """Full C_stop run with off-target index."""
        print(f"{Colors.HEADER}[INFO] C_stop mode — parsing GenBank...{Colors.ENDC}")
        try:
            records = list(SeqIO.parse(self.args.gbff, "genbank"))
            if not records:
                print(f"{Colors.FAIL}[Error] No records found in {self.args.gbff}{Colors.ENDC}")
                sys.exit(1)
        except Exception as exc:
            print(f"{Colors.FAIL}[Error] Failed to parse GenBank: {exc}{Colors.ENDC}")
            sys.exit(1)

        genome_str = "".join(str(r.seq) for r in records)
        self._build_and_print_index(genome_str)

        num_designs = getattr(self.args, 'num_designs', 2)

        for record in records:
            self.current_record_seq = record.seq
            cds_features = [f for f in record.features if f.type == "CDS"]
            cds_features.sort(key=lambda f: f.location.start)

            gene_context = {}
            for i, f in enumerate(cds_features):
                locus  = f.qualifiers.get('locus_tag',
                          f.qualifiers.get('gene', ['Unknown']))[0]
                strand = f.location.strand
                safe, ov, neighbor = 100, 0, None
                if strand == 1:
                    if i < len(cds_features) - 1:
                        nxt      = cds_features[i + 1]
                        neighbor = nxt
                        d        = int(nxt.location.start) - int(f.location.end)
                        safe = max(0, d) + 3
                        ov   = max(0, -d)
                else:
                    if i > 0:
                        prev     = cds_features[i - 1]
                        neighbor = prev
                        d        = int(f.location.start) - int(prev.location.end)
                        safe = max(0, d) + 3
                        ov   = max(0, -d)
                gene_context[locus] = {
                    'safe': safe, 'overlap': ov,
                    'neighbor_feat': neighbor, 'up_feat': None,
                }

            for feature in cds_features:
                qualifiers = feature.qualifiers
                gene_id    = qualifiers.get('locus_tag',
                              qualifiers.get('gene', ['Unknown']))[0]
                if self.args.target_gene and gene_id != self.args.target_gene:
                    continue

                ok, msg = self._validate_cds_boundary(record, feature)
                if not ok:
                    self.results.append({
                        'Design_Rank': 1, 'Gene_ID': gene_id,
                        'Note': f"Skipped: CDS boundary error ({msg})",
                    })
                    self.stats["failed"] += 1
                    continue

                ok, stop_msg = self._validate_cds_stop_codon(record, feature)
                if not ok:
                    self.results.append({
                        'Design_Rank': 1, 'Gene_ID': gene_id,
                        'Note': f"Skipped: {stop_msg}",
                    })
                    self.stats["failed"] += 1
                    continue

                product = qualifiers.get('product', [''])[0]
                self.stats["total"] += 1
                ctx    = gene_context.get(gene_id,
                          {'safe': 100, 'overlap': 0,
                           'neighbor_feat': None, 'up_feat': None})
                strand = feature.location.strand
                if strand == 1:
                    genomic_junction = int(feature.location.end) - 3
                else:
                    genomic_junction = int(feature.location.start) + 3

                self._current_start_codon = None

                try:
                    up, rha = self.get_upstream_and_rha(record, feature)
                    candidates = self.design_target(
                        up, rha, ctx['safe'], ctx['overlap'],
                        ctx['neighbor_feat'], feature, strand, genomic_junction)
                except Exception as exc:
                    self.stats["error"] += 1
                    print(f"{Colors.FAIL}[Error] {gene_id}: {exc}{Colors.ENDC}")
                    if self.args.target_gene:
                        import traceback; traceback.print_exc()
                    continue

                if not candidates:
                    self.results.append({
                        'Design_Rank': 1, 'Gene_ID': gene_id, 'Note': 'Failed',
                    })
                    self.stats["failed"] += 1
                    continue

                accepted_list = self._select_top_designs(
                    candidates, ctx, feature, genomic_junction, num_designs=num_designs)

                if not accepted_list:
                    self.results.append({
                        'Design_Rank': 1, 'Gene_ID': gene_id,
                        'Note': 'Failed (Double Restriction Sites)',
                    })
                    self.stats["failed"] += 1
                    continue

                best_cand = accepted_list[0][0]
                strat_r1  = best_cand['Strategy']
                if best_cand['Mut_Dist'] <= self.args.max_mut_dist:
                    self.stats["optimal_selected"] += 1
                else:
                    self.stats["suboptimal_selected"] += 1
                if self._candidate_has_aa_change(best_cand):
                    self.stats["aa_change_selected"] += 1

                if   strat_r1 == 'CStop_P1_Del_Downstream': self.stats["p1_clean"] += 1
                elif strat_r1 in ('CStop_P2_Bridge',
                                  'CStop_P2_Bridge_Mut'):   self.stats["p2_bridge_nomut"] += 1
                elif strat_r1 == 'CStop_P3_Mut_LHA':        self.stats["p3b_rha"] += 1

                # V44 similarity flagging
                if len(accepted_list) >= 2:
                    lha1 = accepted_list[0][0].get('LHA', '')
                    rha1 = accepted_list[0][0].get('RHA', '')
                    lha2 = accepted_list[1][0].get('LHA', '')
                    rha2 = accepted_list[1][0].get('RHA', '')
                    sim_lha   = round(_seq_sim_right(lha1, lha2) * 100)
                    sim_rha   = round(_seq_sim_left(rha1, rha2) * 100)
                    thresh    = getattr(self.args, 'sim_threshold', 85)
                    if sim_lha >= thresh and sim_rha >= thresh:
                        tag = f"[SIM:LHA={sim_lha}%,RHA={sim_rha}%]"
                        for idx in range(2):
                            c     = accepted_list[idx][0]
                            extra = c.get('Note_Extra', '')
                            c['Note_Extra'] = (extra + ' ' + tag).strip() if extra else tag

                for rank, (cand_use, final_oligo, barcode, risky, msg,
                            re_details, re_types, restriction_sites) \
                        in enumerate(accepted_list, 1):
                    strat  = cand_use['Strategy']
                    status = self._quality_status(cand_use)
                    aa_tag = self._aa_cost_tag(cand_use)
                    note_extra = " | ".join(
                        part for part in [cand_use.get('Note_Extra', '').strip(' | '), aa_tag, status] if part
                    )
                    if risky:
                        print(f"{Colors.FAIL}   [RISK] {gene_id:<15} Rank{rank} | {msg}{Colors.ENDC}")
                    self.results.append({
                        'Design_Rank':       rank,
                        'Gene_ID':           gene_id,
                        'Product':           product,
                        'Strategy':          strat,
                        'Info':              note_extra,
                        'sgRNA':             cand_use['sgRNA'],
                        'sgRNA_GC':          self._calc_gc(cand_use['sgRNA']),
                        'PAM':               cand_use['PAM'],
                        'LHA':               cand_use['LHA'],
                        'LHA_Length':        len(cand_use['LHA']),
                        'LHA_GC':            self._calc_gc(cand_use['LHA']),
                        'RHA':               cand_use['RHA'],
                        'RHA_Length':        len(cand_use['RHA']),
                        'RHA_GC':            self._calc_gc(cand_use['RHA']),
                        'Restriction_Sites': restriction_sites,
                        'Mut_Dist':          cand_use['Mut_Dist'],
                        'Barcode':           barcode,
                        'Final_Oligo':       final_oligo,
                        'Note':              f"Success ({strat})",
                    })

        # Save CSV
        df = pd.DataFrame(self.results)
        cols = [
            'Design_Rank', 'Gene_ID', 'Product', 'Strategy', 'Info',
            'sgRNA', 'sgRNA_GC', 'PAM',
            'LHA', 'LHA_Length', 'LHA_GC',
            'RHA', 'RHA_Length', 'RHA_GC',
            'Restriction_Sites', 'Mut_Dist', 'Barcode', 'Final_Oligo', 'Note',
        ]
        final_cols = [c for c in cols if c in df.columns]
        df[final_cols].to_csv(self.args.output, index=False)

        self._inject_offtarget_column()

        total_success = self.stats['optimal_selected'] + self.stats['suboptimal_selected']
        rank2_count   = (len(df[df['Design_Rank'] == 2]) if 'Design_Rank' in df.columns else 0)

        print(f"\n{Colors.HEADER}{'=' * 70}")
        print(f" Summary Statistics (C_stop mode)")
        print(f"{'=' * 70}{Colors.ENDC}")
        print(f" Total Genes Processed  : {self.stats['total']}")
        if self.stats['total']:
            print(f" Genes with ≥1 Design   : {total_success}"
                  f" ({total_success / self.stats['total'] * 100:.1f}%)")
            print(f" Genes with 2nd Design  : {rank2_count}")
            print(f" Failed Designs         : {self.stats['failed']}"
                  f" ({self.stats['failed'] / self.stats['total'] * 100:.1f}%)")
        print(f"\n{Colors.OKBLUE} Strategy Distribution (Rank-1):{Colors.ENDC}")
        if total_success:
            p1 = self.stats['p1_clean']
            p2 = self.stats['p2_bridge_nomut']
            p3 = self.stats['p3b_rha']
            print(f" ├─ CStop_P1_Del_Downstream : {p1:4d} ({p1/total_success*100:5.1f}%)")
            print(f" ├─ CStop_P2_Bridge[_Mut]   : {p2:4d} ({p2/total_success*100:5.1f}%)")
            print(f" └─ CStop_P3_Mut_LHA        : {p3:4d} ({p3/total_success*100:5.1f}%)")
            opt = self.stats['optimal_selected']
            sub = self.stats['suboptimal_selected']
            aa_change = self.stats['aa_change_selected']
            print(f"\n{Colors.OKGREEN} Distance Metrics:{Colors.ENDC}")
            print(f" ├─ Optimal  (Mut_Dist ≤ {self.args.max_mut_dist}bp) : {opt:4d} ({opt/total_success*100:5.1f}%)")
            print(f" └─ Suboptimal (Mut_Dist > {self.args.max_mut_dist}bp): {sub:4d} ({sub/total_success*100:5.1f}%)")
            print(f"{Colors.OKGREEN} Amino-acid Impact:{Colors.ENDC}")
            print(f" └─ Rank-1 conservative fallback designs            : {aa_change:4d} ({aa_change/total_success*100:5.1f}%)")
        else:
            print(" No successful designs to display.")

        self._print_offtarget_summary()

    # ------------------------------------------------------------------
    # Post-processing helpers
    # ------------------------------------------------------------------

    def _inject_offtarget_column(self):
        try:
            df = pd.read_csv(self.args.output)

            def lookup_ot(row):
                sgrna = str(row.get('sgRNA', '')).upper().strip()
                if not sgrna or sgrna == 'NAN':
                    return 0
                return self.offtarget_index.get(sgrna, 0)

            df['sgRNA_OffTarget_Count'] = df.apply(lookup_ot, axis=1)
            cols = list(df.columns)
            if 'PAM' in cols:
                pam_idx = cols.index('PAM')
                cols.remove('sgRNA_OffTarget_Count')
                cols.insert(pam_idx + 1, 'sgRNA_OffTarget_Count')
                df = df[cols]
            df.to_csv(self.args.output, index=False)
            print(f"{Colors.OKGREEN}[V2] sgRNA_OffTarget_Count column injected into "
                  f"{self.args.output}{Colors.ENDC}")
        except Exception as e:
            print(f"{Colors.WARNING}[V2] Could not inject OffTarget_Count: {e}{Colors.ENDC}")

    def _print_offtarget_summary(self):
        ot_rejected = self.stats.get('sgrna_offtarget', 0)
        if ot_rejected:
            print(f"{Colors.WARNING}  [V2] Off-target filter rejected "
                  f"{ot_rejected} sgRNA candidates "
                  f"(max_offtargets={self.args.max_offtargets}){Colors.ENDC}")
        else:
            print(f"{Colors.OKBLUE}  [V2] Off-target filter: 0 candidates rejected "
                  f"(all unique or below threshold){Colors.ENDC}")


# ===========================================================================
# CLI
# ===========================================================================

def get_args():
    parser = argparse.ArgumentParser(
        description=(
            "CRISPR Knockin Designer V6 (Standalone) — N_start (V46-equivalent) "
            "or C_stop (stop-codon-targeted) mode with sgRNA specificity filter"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument('--model', choices=['N_start', 'C_stop'], default='N_start',
                        help="Insertion mode.")
    parser.add_argument("--gbff",    default="MG1655_genomic.gbff",    help="GenBank file.")
    parser.add_argument("--payload", default="J23119_RBS",             help="Payload sequence or file.")
    parser.add_argument("--template",
                        default="Knockin_J23100RBS_library_oligo_template.fasta",
                        help="Oligo template FASTA.")
    parser.add_argument("--output",  default="knockin_output.csv",     help="Output CSV file.")
    parser.add_argument("--target_gene", help="Design only this locus tag (debug mode).")

    parser.add_argument("--num_designs",   type=int,   default=2)

    parser.add_argument("--lha_len",         type=int,   default=55)
    parser.add_argument("--rha_len",         type=int,   default=55)
    parser.add_argument("--min_effective_ha",type=int,   default=40)
    parser.add_argument("--max_oligo_len",   type=int,   default=300)
    parser.add_argument("--min_gc_ha",       type=float, default=15.0)
    parser.add_argument("--min_gc_sgrna",    type=float, default=20.0)
    parser.add_argument("--max_gc_sgrna",    type=float, default=80.0)
    parser.add_argument("--search_window",   type=int,   default=50)
    parser.add_argument("--max_deletion",    type=int,   default=150)
    parser.add_argument("--seed_region",     type=int,   default=12)
    parser.add_argument("--max_mut_dist",    type=int,   default=12)

    parser.add_argument("--barcode_len",  type=int, default=10)
    parser.add_argument("--restriction_site", action='append',
                        default=['GGTCTC', 'GAAGAC'],
                        help="RE sequences to exclude.")

    parser.add_argument("--barcode_seed",  type=int, default=42,
                        help="RNG seed for deterministic barcode generation.")
    parser.add_argument("--sim_threshold", type=int, default=85,
                        help="LHA AND RHA identity %% threshold for SIM flagging.")
    parser.add_argument("--rank2_sim_max", type=int, default=50,
                        help="[C_stop V3] Max allowed max(SIM_LHA,SIM_RHA) %% for same-strategy "
                             "Rank2 (Option A). Candidates exceeding this trigger final fallback.")

    parser.add_argument("--max_offtargets", type=int, default=0,
                        help="Max allowed off-target hits (0 = genome-unique only).")

    parser.add_argument("--report_fail_reasons",   action='store_true')
    parser.add_argument("--export_fail_sequences", action='store_true')

    return parser.parse_args()


if __name__ == "__main__":
    args     = get_args()
    designer = KnockinDesigner(args)
    designer.run()
    print(
        f"\n{Colors.OKGREEN}[Done] Results saved to {args.output}{Colors.ENDC}\n"
        f"{Colors.OKBLUE}  model={args.model} | barcode_seed={args.barcode_seed} | "
        f"sim_threshold={args.sim_threshold}% | rank2_sim_max={args.rank2_sim_max}% | "
        f"max_offtargets={args.max_offtargets} | "
        f"lha={args.lha_len}bp | rha={args.rha_len}bp | "
        f"max_oligo={args.max_oligo_len}bp{Colors.ENDC}"
    )
