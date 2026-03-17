#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CRISPR_knockin.py — N_start / C_stop dual-mode knockin library design
======================================================================
Author : Zhaohui Cai with GitHub Copilot Assistant
Version: 1.0
Date   : 2026-03-18

Modes
-----
  N_start (default): identical to V44 — sgRNA/HDR around the start codon;
                     insert payload (J23119 RBS) immediately BEFORE the ATG.

  C_stop           : sgRNA/HDR around the stop codon; insert C-terminal
                     fusion tag immediately BEFORE the stop codon so that
                     the payload replaces the stop codon in the construct.

Key design principles for C_stop
---------------------------------
  * junction = physical start of stop codon
      Strand +1: int(feature.location.end) - 3
      Strand -1: int(feature.location.start) + 3
  * upstream_seq = last N bp of coding EXCLUDING stop codon (anchored at junc)
  * rha_seq[0:3]  = stop codon (preserved in full sequence for PAM search)
  * rha_seq[3:]   = 3'UTR / downstream intergenic
  * LHA ends AT the junction (right before stop codon)
  * RHA starts AFTER the stop codon (rha_seq[3:] slice)
  * safe_downstream_limit = dist(gene.end, downstream.start) + 3 bp
  * _sanitize_arm LHA ref = current feature (coding gene)
  * _sanitize_arm RHA ref = downstream neighbour feature

Strategies for C_stop
----------------------
  CStop_P1_Del_Downstream  - target fully in stop codon or 3'UTR:
                              delete [junction, t_e), RHA starts at t_e
  CStop_P2_Bridge          - target spans stop codon junction (no mut needed)
  CStop_P2_Bridge_Mut      - target spans stop codon junction (LHA/RHA mutated)
  CStop_P3_Mut_LHA         - target fully in coding before stop codon:
                              mutate LHA synonymously / conservatively

Silent mutation strategy applies equally to LHA and RHA for both modes.

Usage
-----
  python CRISPR_knockin.py --model N_start [V44 args ...]
  python CRISPR_knockin.py --model C_stop \\
      --template Cfusion_library_oligo_template.fasta [args ...]
"""

import argparse
import random
import sys
from collections import Counter

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from knockin_J23119RBS_V44 import (
    GenBankLibraryDesigner as V44Designer,
    _seq_sim_right,
    _seq_sim_left,
    Colors,
    VALID_START_CODONS,
    _BALANCE_TOL,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

STOP_CODONS = frozenset({'TAA', 'TAG', 'TGA'})

# ---------------------------------------------------------------------------
# KnockinDesigner
# ---------------------------------------------------------------------------

class KnockinDesigner(V44Designer):
    """N_start (V44 passthrough) + C_stop dual-mode knockin designer.

    For N_start: every method delegates to super(), giving identical output
    to V44 when the same arguments are supplied.

    For C_stop: overrides get_upstream_and_rha(), run(), design_target(),
    and _select_top_designs() to implement stop-codon-targeted insertion.
    """

    # ------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------

    def __init__(self, args):
        self.model = getattr(args, 'model', 'N_start')
        super().__init__(args)          # seeds RNG, builds barcode pool, computes ha_budget

    # ------------------------------------------------------------------
    # Extra C_stop validation: stop codon must be a standard stop codon
    # ------------------------------------------------------------------

    def _validate_cds_stop_codon(self, record, feature):
        """Verify stop codon is TAA/TAG/TGA.  Call AFTER _validate_cds_boundary."""
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
    # Sequence extraction — dispatch by mode
    # ------------------------------------------------------------------

    def get_upstream_and_rha(self, record, feature):
        if self.model == 'N_start':
            return super().get_upstream_and_rha(record, feature)
        return self._get_upstream_rha_cstop(record, feature)

    def _get_upstream_rha_cstop(self, record, feature):
        """Extract (upstream_seq, rha_seq) for C_stop mode.

        upstream_seq        : last req_len bp of coding EXCLUDING stop codon.
        rha_seq[0:3]        : stop codon (used for PAM search only).
        rha_seq[3:]         : 3'UTR / downstream intergenic.
        """
        start    = int(feature.location.start)
        end      = int(feature.location.end)
        strand   = feature.location.strand
        full_seq = record.seq
        seq_len  = len(full_seq)

        rha_buffer = self.args.rha_len + 500
        req_up     = self.args.search_window + self.args.lha_len + 400

        if strand == 1:
            stop_start = end - 3          # first physical base of stop codon
            phy_up_s   = max(0, stop_start - req_up)
            upstream_raw = str(full_seq[phy_up_s : stop_start]).upper()
            rha_raw      = str(full_seq[stop_start :
                                        min(seq_len, stop_start + rha_buffer)]).upper()
        else:
            # Stop codon occupies physical [start, start+3) on minus strand.
            # Coding direction: right-to-left physically.
            # "Last req_up coding bases before the stop codon" in display direction
            # = the req_up bases physically CLOSEST to [start+3, …), i.e.
            #   physical [start+3, start+3+req_up).
            # (Taking bases from the high-coordinate / start-codon end instead
            #  would place the junction wildly outside the extracted window.)
            phy_up_e     = min(end, start + 3 + req_up)
            upstream_raw = str(full_seq[start + 3 : phy_up_e].reverse_complement()).upper()

            # rha_seq = RC( seq[start-rha_buffer : start+3] )
            #         = stop_codon (RC of seq[start:start+3])  +  3'UTR downstream
            phy_rha_e = start + 3
            phy_rha_s = max(0, start + 3 - rha_buffer)
            rha_raw   = str(full_seq[phy_rha_s : phy_rha_e].reverse_complement()).upper()

        return upstream_raw, rha_raw

    # ------------------------------------------------------------------
    # run() — dispatch by mode
    # ------------------------------------------------------------------

    def run(self):
        if self.model == 'N_start':
            return super().run()
        self._run_cstop()

    def _run_cstop(self):
        """Full C_stop run.  Mirror of V38 run() with stop-codon junction logic."""
        print(f"{Colors.HEADER}[INFO] C_stop mode — parsing GenBank...{Colors.ENDC}")
        try:
            records = list(SeqIO.parse(self.args.gbff, "genbank"))
            if not records:
                print(f"{Colors.FAIL}[Error] No records found in {self.args.gbff}{Colors.ENDC}")
                sys.exit(1)
        except Exception as exc:
            print(f"{Colors.FAIL}[Error] Failed to parse GenBank: {exc}{Colors.ENDC}")
            sys.exit(1)

        num_designs = getattr(self.args, 'num_designs', 2)

        for record in records:
            self.current_record_seq = record.seq
            cds_features = [f for f in record.features if f.type == "CDS"]
            cds_features.sort(key=lambda f: f.location.start)

            # Build DOWNSTREAM gene context (mirror of N_start's upstream context)
            gene_context = {}
            for i, f in enumerate(cds_features):
                locus  = f.qualifiers.get('locus_tag',
                          f.qualifiers.get('gene', ['Unknown']))[0]
                strand = f.location.strand
                safe, ov, neighbor = 100, 0, None

                if strand == 1:
                    # downstream for Strand+1 = next CDS in sorted order
                    if i < len(cds_features) - 1:
                        nxt      = cds_features[i + 1]
                        neighbor = nxt
                        d        = int(nxt.location.start) - int(f.location.end)
                        # +3 because the junction is at stop_start = gene.end-3,
                        # not at gene.end.  Deletion in C_stop counts from junction.
                        safe = max(0, d) + 3
                        ov   = max(0, -d)
                else:
                    # downstream for Strand-1 = previous CDS (lower index)
                    if i > 0:
                        prev     = cds_features[i - 1]
                        neighbor = prev
                        d        = int(f.location.start) - int(prev.location.end)
                        safe = max(0, d) + 3
                        ov   = max(0, -d)

                gene_context[locus] = {
                    'safe':          safe,
                    'overlap':       ov,
                    'neighbor_feat': neighbor,
                    # Maintain up_feat key for compatibility with parent stage-2
                    # (not used by _select_top_designs_cstop)
                    'up_feat':       None,
                }

            # Per-gene design loop
            for feature in cds_features:
                qualifiers = feature.qualifiers
                gene_id    = qualifiers.get('locus_tag',
                              qualifiers.get('gene', ['Unknown']))[0]

                if self.args.target_gene and gene_id != self.args.target_gene:
                    continue

                # CDS validation
                ok, msg = self._validate_cds_boundary(record, feature)
                if not ok:
                    self.results.append({
                        'Design_Rank': 1, 'Gene_ID': gene_id,
                        'Note': f"Skipped: CDS boundary error ({msg})",
                    })
                    self.stats["failed"] += 1
                    continue

                # Extra stop-codon check
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

                # C_stop genomic_junction = first physical base of stop codon
                if strand == 1:
                    genomic_junction = int(feature.location.end) - 3
                else:
                    genomic_junction = int(feature.location.start) + 3

                # Disable start-codon gate (Stage 0) — not applicable for C_stop
                self._current_start_codon = None

                try:
                    up, rha = self.get_upstream_and_rha(record, feature)
                    candidates = self.design_target(
                        up, rha, ctx['safe'], ctx['overlap'],
                        ctx['neighbor_feat'], feature, strand, genomic_junction,
                    )
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
                    candidates, ctx, feature, genomic_junction,
                    num_designs=num_designs,
                )

                if not accepted_list:
                    self.results.append({
                        'Design_Rank': 1, 'Gene_ID': gene_id,
                        'Note': 'Failed (Double Restriction Sites)',
                    })
                    self.stats["failed"] += 1
                    continue

                # Stats for rank-1
                best_cand   = accepted_list[0][0]
                strat_r1    = best_cand['Strategy']
                is_optimal  = best_cand['Mut_Dist'] <= self.args.max_mut_dist
                if is_optimal:
                    self.stats["optimal_selected"] += 1
                else:
                    self.stats["suboptimal_selected"] += 1

                if   strat_r1 == 'CStop_P1_Del_Downstream': self.stats["p1_clean"] += 1
                elif strat_r1 in ('CStop_P2_Bridge',
                                  'CStop_P2_Bridge_Mut'):   self.stats["p2_bridge_nomut"] += 1
                elif strat_r1 == 'CStop_P3_Mut_LHA':        self.stats["p3b_rha"] += 1

                # Similarity flagging (same AND-condition logic as V44)
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
                            cand  = accepted_list[idx][0]
                            extra = cand.get('Note_Extra', '')
                            cand['Note_Extra'] = (extra + ' ' + tag).strip() if extra else tag

                # Record results
                for rank, (cand_use, final_oligo, barcode, risky, msg,
                            re_details, re_types, restriction_sites) \
                        in enumerate(accepted_list, 1):
                    strat  = cand_use['Strategy']
                    status = ("(Optimal)"
                              if cand_use['Mut_Dist'] <= self.args.max_mut_dist
                              else "(Suboptimal [Mut_Dist])")
                    note_extra = (cand_use.get('Note_Extra', '') + f" | {status}").strip()

                    if risky:
                        print(f"{Colors.FAIL}   [RISK] {gene_id:<15} Rank{rank} | {msg}{Colors.ENDC}")

                    self.results.append({
                        'Design_Rank':        rank,
                        'Gene_ID':            gene_id,
                        'Product':            product,
                        'Strategy':           strat,
                        'Info':               note_extra,
                        'sgRNA':              cand_use['sgRNA'],
                        'sgRNA_GC':           self._calc_gc(cand_use['sgRNA']),
                        'PAM':                cand_use['PAM'],
                        'LHA':                cand_use['LHA'],
                        'LHA_Length':         len(cand_use['LHA']),
                        'LHA_GC':             self._calc_gc(cand_use['LHA']),
                        'RHA':                cand_use['RHA'],
                        'RHA_Length':         len(cand_use['RHA']),
                        'RHA_GC':             self._calc_gc(cand_use['RHA']),
                        'Restriction_Sites':  restriction_sites,
                        'Mut_Dist':           cand_use['Mut_Dist'],
                        'Barcode':            barcode,
                        'Final_Oligo':        final_oligo,
                        'Note':               f"Success ({strat})",
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
        rank2_count   = (len(df[df['Design_Rank'] == 2])
                         if 'Design_Rank' in df.columns else 0)

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
            print(f"\n{Colors.OKGREEN} Quality Metrics:{Colors.ENDC}")
            print(f" ├─ Optimal  (Mut_Dist ≤ {self.args.max_mut_dist}bp) : {opt:4d} ({opt/total_success*100:5.1f}%)")
            print(f" └─ Suboptimal (Mut_Dist > {self.args.max_mut_dist}bp): {sub:4d} ({sub/total_success*100:5.1f}%)")
        else:
            print(" No successful designs to display.")

    # ------------------------------------------------------------------
    # design_target() — dispatch by mode
    # ------------------------------------------------------------------

    def design_target(self, upstream_seq, rha_seq, safe_limit, overlap,
                      upstream_feat, current_feat, strand, genomic_junction):
        if self.model == 'N_start':
            return super().design_target(
                upstream_seq, rha_seq, safe_limit, overlap,
                upstream_feat, current_feat, strand, genomic_junction,
            )
        return self._design_target_cstop(
            upstream_seq, rha_seq, safe_limit, overlap,
            upstream_feat, current_feat, strand, genomic_junction,
        )

    def _design_target_cstop(self, upstream_seq, rha_seq,
                              safe_downstream_limit, overlap_bases,
                              neighbor_feature, current_feature,
                              target_strand, genomic_junction):
        """C_stop sgRNA/HDR design around the stop codon.

        Parameters
        ----------
        upstream_seq          : coding sequence before stop codon (ends at junction)
        rha_seq               : rha_seq[0:3] = stop codon, rha_seq[3:] = 3'UTR
        safe_downstream_limit : max deletion length permitted toward downstream gene
        overlap_bases         : coding overlap with downstream neighbour
        neighbor_feature      : downstream CDS feature (ref for RHA mutation)
        current_feature       : current CDS feature   (ref for LHA mutation)
        target_strand         : +1 or -1
        genomic_junction      : physical position of stop codon start
        """
        candidates   = []
        junction_idx = len(upstream_seq)
        full         = upstream_seq + rha_seq

        search_s = max(0, junction_idx - self.args.search_window)
        search_e = junction_idx + self.args.search_window

        # RHA genomic junction (where rha_seq[3] maps to; = gene boundary)
        # Used for mutation coordinate mapping in RHA.
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
            else:                                    # CCN
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

            # ──────────────────────────────────────────────────────────
            # P1 — target fully in stop codon or 3'UTR (t_s >= junction)
            # Deletion: [junction_idx, t_e) = stop codon + UTR up to t_e
            # LHA anchored to coding; RHA starts after the deletion
            # ──────────────────────────────────────────────────────────
            if t_s >= junction_idx:
                del_len = t_e - junction_idx
                if del_len <= safe_downstream_limit and del_len <= self.args.max_deletion:
                    lha_fin, _shift = self._optimize_lha_window(
                        upstream_seq, junction_idx,
                        self.args.min_effective_ha, self.args.lha_len,
                    )
                    rha_extracted = full[t_e : t_e + self.args.rha_len]
                    if (lha_fin and rha_extracted
                            and self.validator.is_gc_valid(rha_extracted)):
                        candidates.append({
                            'sgRNA':      sp_seq,
                            'PAM':        pam_seq,
                            'LHA':        lha_fin,
                            'RHA':        rha_extracted,
                            'Distance':   dist,
                            'Mut_Dist':   0,
                            'Strategy':   'CStop_P1_Del_Downstream',
                            'Note_Extra': f"CStop Del: {del_len}bp (incl. stop codon)",
                        })

            # ──────────────────────────────────────────────────────────
            # P2 — target spans the stop codon junction
            # LHA = last lha_len bp of coding + rescue (3'UTR overlap)
            # RHA = rha_seq[3:] (skip stop codon)
            # ──────────────────────────────────────────────────────────
            elif t_s < junction_idx < t_e:
                req_lha = (self.args.min_effective_ha
                           + (abs(cut_offset) if cut_offset < 0 else 0))
                req_rha = (self.args.min_effective_ha
                           + (cut_offset if cut_offset > 0 else 0))

                if req_lha > self.args.lha_len or req_rha > self.args.rha_len:
                    continue

                grab_lha_len  = min(self.args.lha_len, junction_idx)
                raw_lha_part  = upstream_seq[junction_idx - grab_lha_len : junction_idx]
                rescue_seq    = rha_seq[3 : 3 + overlap_bases]  # 3'UTR rescue
                lha_cand_full = raw_lha_part + rescue_seq

                grab_rha_len  = min(self.args.rha_len, max(0, len(rha_seq) - 3))
                raw_rha       = rha_seq[3 : 3 + grab_rha_len]   # skip stop codon

                eff_rha_check = len(raw_rha) - (cut_offset if cut_offset > 0 else 0)
                if (eff_rha_check < self.args.min_effective_ha
                        or not self.validator.is_gc_valid(raw_rha)):
                    continue

                # Determine PAM/seed re-cut risk in each arm
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

                lha_final    = rha_final    = None
                lha_mut_idx  = rha_mut_idx  = None
                note         = ""
                mut_dist     = 0
                is_mut       = False

                # 1. Try RHA mutation if needed (PAM in 3'UTR — ref=neighbor)
                if risk_in_rha:
                    pam_in_rha = pam_s - junction_idx - 3   # relative to raw_rha
                    if 0 <= pam_in_rha < len(raw_rha):
                        mut_rha, ok, log, mut_idx, loc_lbl = self.perform_precise_mutation(
                            raw_rha, pam_in_rha, ptype,
                            neighbor_feature, target_strand, rha_gj, 0,
                            ref_cds_label="Neighbor",
                        )
                        if ok:
                            rha_final   = mut_rha
                            note       += f"RHA Mut: {log} [{loc_lbl}] | "
                            is_mut      = True
                            rha_mut_idx = mut_idx
                        else:
                            rha_final = None           # cannot fix RHA → skip
                    else:
                        # PAM falls in the stop-codon region (deleted) — no RHA risk
                        rha_final = raw_rha
                else:
                    rha_final = raw_rha

                # 2. Try LHA mutation if needed (PAM in coding — ref=current)
                if rha_final and risk_in_lha:
                    start_of_lha  = junction_idx - grab_lha_len
                    pam_in_lha    = pam_s - start_of_lha
                    if 0 <= pam_in_lha < len(lha_cand_full):
                        mut_lha, ok, log, mut_idx, loc_lbl = self.perform_precise_mutation(
                            lha_cand_full, pam_in_lha, ptype,
                            current_feature, target_strand,
                            genomic_junction, len(raw_lha_part),
                            ref_cds_label="Current",
                        )
                        if ok:
                            req_opt   = req_lha + overlap_bases
                            lha_final, shift = self._optimize_lha_window(
                                mut_lha, len(mut_lha), req_opt,
                            )
                            if lha_final is not None:
                                lha_mut_idx = mut_idx - shift
                                if lha_mut_idx < 0:
                                    lha_final = None    # mutation trimmed away
                            note    += f"LHA Mut: {log} [{loc_lbl}]"
                            is_mut   = True
                            mut_dist = len(mut_lha) - mut_idx
                        else:
                            lha_final = None
                    else:
                        lha_final = None
                elif rha_final:
                    req_opt   = req_lha + overlap_bases
                    lha_final, shift = self._optimize_lha_window(
                        lha_cand_full, len(lha_cand_full), req_opt,
                    )
                    if shift > 0:
                        note += f"LHA GC Shift -{shift}bp"

                if lha_final and rha_final and len(lha_final) <= self.args.lha_len:
                    stgy       = 'CStop_P2_Bridge_Mut' if is_mut else 'CStop_P2_Bridge'
                    note_clean = note.strip(" | ")
                    cand = {
                        'sgRNA':      sp_seq,
                        'PAM':        pam_seq,
                        'LHA':        lha_final,
                        'RHA':        rha_final,
                        'Distance':   dist,
                        'Strategy':   stgy,
                        'Note_Extra': note_clean,
                        'Mut_Dist':   mut_dist,
                    }
                    if is_mut:
                        if 'LHA Mut' in note and lha_mut_idx is not None:
                            cand['LHA_Mut_Index'] = lha_mut_idx
                        if 'RHA Mut' in note and rha_mut_idx is not None:
                            cand['RHA_Mut_Index'] = rha_mut_idx
                    candidates.append(cand)

            # ──────────────────────────────────────────────────────────
            # P3 — target fully in coding before stop codon (t_e <= junction)
            # Mutate LHA (coding region) to destroy PAM
            # RHA = rha_seq[3:3+rha_len] — fixed, skips stop codon
            # ──────────────────────────────────────────────────────────
            elif t_e <= junction_idx:
                req_lha = (self.args.min_effective_ha
                           + (abs(cut_offset) if cut_offset < 0 else 0))
                req_rha = self.args.min_effective_ha

                if req_lha > self.args.lha_len or req_rha > self.args.rha_len:
                    continue

                # Mutate upstream_seq (coding) at pam_s
                mut_upstream, ok, log, mut_idx, _ = self.perform_precise_mutation(
                    upstream_seq, pam_s, ptype,
                    current_feature, target_strand,
                    genomic_junction, junction_idx,
                    ref_cds_label="Current",
                )
                if not ok:
                    continue

                # Extract LHA from the right end of the mutated coding sequence.
                # total_trim = initial left-trim (to fit lha_len) + GC shift.
                # Must use total_trim, NOT just _gc_shift, to adjust the mutation index.
                lha_base, _gc_shift = self._optimize_lha_window(
                    mut_upstream, len(mut_upstream),
                    self.args.min_effective_ha, self.args.lha_len,
                )
                if lha_base is None:
                    continue

                total_trim  = len(mut_upstream) - len(lha_base)
                new_mut_idx = mut_idx - total_trim
                if new_mut_idx < 0:
                    continue        # mutation was trimmed away from LHA window

                rescue_seq    = rha_seq[3 : 3 + overlap_bases]
                lha_extracted = lha_base + rescue_seq
                rha_extracted = rha_seq[3 : 3 + self.args.rha_len]

                if (lha_extracted and rha_extracted
                        and self.validator.is_gc_valid(rha_extracted)):
                    candidates.append({
                        'sgRNA':          sp_seq,
                        'PAM':            pam_seq,
                        'LHA':            lha_extracted,
                        'RHA':            rha_extracted,
                        'Distance':       dist,
                        'Mut_Dist':       dist,
                        'Strategy':       'CStop_P3_Mut_LHA',
                        'Note_Extra':     log,
                        'LHA_Mut_Index':  new_mut_idx,
                    })

        return candidates

    # ------------------------------------------------------------------
    # _select_top_designs() — dispatch by mode
    # ------------------------------------------------------------------

    def _select_top_designs(self, candidates, ctx, feature, genomic_junction,
                            num_designs=2):
        if self.model == 'N_start':
            return super()._select_top_designs(
                candidates, ctx, feature, genomic_junction, num_designs,
            )
        return self._select_top_designs_cstop(
            candidates, ctx, feature, genomic_junction, num_designs,
        )

    def _select_top_designs_cstop(self, candidates, ctx, feature, genomic_junction,
                                   num_designs=2):
        """C_stop pipeline — V43 stages with corrected arm sanitise references.

        Differences from V43 _select_top_designs:
          Stage 0  SKIPPED  — RHA is in 3'UTR, not at start codon.
          Stage 2 LHA ref   = current feature (coding gene).
          Stage 2 RHA ref   = ctx['neighbor_feat'] (downstream CDS).
          Stage 2 RHA gjunc = gene.end (Strand+1) or gene.start (Strand-1).
        All other stages (1–8) are identical to V43.
        """
        budget = self.ha_budget

        # RHA genomic junction for sanitise coordinate mapping
        rha_gj = (int(feature.location.end)
                  if feature.location.strand == 1
                  else int(feature.location.start))

        def sort_key(c):
            is_safe     = c['Mut_Dist'] <= self.args.max_mut_dist
            strat       = c['Strategy']
            strat_score = 3 if 'P1' in strat else (2 if 'P2' in strat else 1)
            return (is_safe, strat_score, -c['Distance'])

        used_sgrnas  = set()
        accepted_list = []

        for cand in sorted(candidates, key=sort_key, reverse=True):
            if len(accepted_list) >= num_designs:
                break

            sgrna = cand.get('sgRNA', '')
            if sgrna in used_sgrnas:
                continue

            # ── Stage 1: Hard-reject sgRNA containing BsaI/BbsI ──
            if self._seq_contains_re(sgrna):
                continue

            cand_use = dict(cand)

            # ── Stage 2: Sanitise LHA (ref=current gene) / RHA (ref=neighbor) ──
            if self.args.restriction_site:
                lha_seq      = cand_use.get('LHA')
                rha_seq_cand = cand_use.get('RHA')

                if lha_seq:
                    lha_fixed, lha_logs, lha_changed = self._sanitize_arm(
                        lha_seq, 'LHA',
                        feature,                        # C_stop: LHA ref = current gene
                        feature.location.strand, genomic_junction, len(lha_seq),
                    )
                else:
                    lha_fixed, lha_logs, lha_changed = lha_seq, [], False

                if rha_seq_cand:
                    rha_fixed, rha_logs, rha_changed = self._sanitize_arm(
                        rha_seq_cand, 'RHA',
                        ctx['neighbor_feat'],           # C_stop: RHA ref = downstream CDS
                        feature.location.strand, rha_gj, 0,
                    )
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
                        cand_use['Note_Extra'] = (
                            (extra + " | " + ";".join(logs)) if extra
                            else ";".join(logs)
                        )

            # ── Stage 3: Hard-reject if LHA/RHA still contains RE ──
            if self._seq_contains_re(cand_use.get('LHA', '')):
                continue
            if self._seq_contains_re(cand_use.get('RHA', '')):
                continue

            # ── Stage 4: Mutation-aware balanced HA trim ──
            lha_raw   = cand_use.get('LHA', '')
            rha_raw   = cand_use.get('RHA', '')
            lha_mut_i = cand_use.get('LHA_Mut_Index')
            rha_mut_i = cand_use.get('RHA_Mut_Index')

            trim_result = self._mutation_aware_trim(
                lha_raw, rha_raw, lha_mut_i, rha_mut_i,
            )
            if trim_result is None:
                continue

            lha_t, rha_t, new_lha_mut_i, eff_lha, eff_rha = trim_result

            # ── Stage 4.5: Budget-fill RHA from buffer (if present) ──
            rha_padded = False
            rha_buf    = cand_use.get('RHA_Buffer')
            if rha_buf and len(lha_t) + len(rha_t) < budget:
                desired   = budget - len(lha_t)
                if desired <= len(rha_buf):
                    if rha_mut_i is None or rha_mut_i < desired:
                        ext_rha = rha_buf[:desired]
                        if self.validator.is_gc_valid(ext_rha):
                            rha_t      = ext_rha
                            eff_rha    = ((len(rha_t) - rha_mut_i)
                                          if rha_mut_i is not None else len(rha_t))
                            rha_padded = True

            # Build note fragments
            note_parts = []
            if len(lha_t) != len(lha_raw) or len(rha_t) != len(rha_raw):
                note_parts.append(
                    f"Trimmed:LHA{len(lha_raw)}->{len(lha_t)},"
                    f"RHA{len(rha_raw)}->{len(rha_t)}"
                )
            elif rha_padded:
                note_parts.append(
                    f"Padded:LHA{len(lha_raw)},"
                    f"RHA{len(cand_use.get('RHA',''))}->{len(rha_t)}"
                )
            if lha_mut_i is not None or rha_mut_i is not None:
                note_parts.append(f"Eff:LHA={eff_lha}bp,RHA={eff_rha}bp")

            if note_parts:
                addition = " | ".join(note_parts)
                extra    = cand_use.get('Note_Extra', '')
                cand_use['Note_Extra'] = (extra + " | " + addition) if extra else addition

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
                 re_details, re_types, restriction_sites)
            )

        return accepted_list


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def get_args():
    parser = argparse.ArgumentParser(
        description=(
            "CRISPR Knockin Designer — N_start (V44-identical) "
            "or C_stop (stop-codon-targeted) mode"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument('--model', choices=['N_start', 'C_stop'], default='N_start',
                        help="Insertion mode: N_start (before ATG) or C_stop (before stop codon).")
    parser.add_argument("--gbff",    default="MG1655_genomic.gbff",    help="GenBank file.")
    parser.add_argument("--payload", default="J23119_RBS",             help="Payload sequence or file.")
    parser.add_argument("--template",
                        default="Knockin_J23100RBS_library_oligo_template.fasta",
                        help=("Oligo template FASTA.  "
                              "For C_stop use Cfusion_library_oligo_template.fasta."))
    parser.add_argument("--output",  default="knockin_output.csv",     help="Output CSV file.")
    parser.add_argument("--target_gene", help="Design only this locus tag (debug mode).")

    parser.add_argument("--num_designs",   type=int,   default=2)

    # Design constraints
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
                        help="RE sequences to exclude (may be specified multiple times).")

    # V44 extras
    parser.add_argument("--barcode_seed",  type=int, default=42,
                        help="RNG seed for deterministic barcode generation.")
    parser.add_argument("--sim_threshold", type=int, default=85,
                        help="LHA AND RHA identity %% threshold for SIM flagging.")

    parser.add_argument("--report_fail_reasons",   action='store_true')
    parser.add_argument("--export_fail_sequences", action='store_true')

    return parser.parse_args()


if __name__ == "__main__":
    args    = get_args()
    designer = KnockinDesigner(args)
    designer.run()
    print(
        f"\n{Colors.OKGREEN}[Done] Results saved to {args.output}{Colors.ENDC}\n"
        f"{Colors.OKBLUE}  model={args.model} | barcode_seed={args.barcode_seed} | "
        f"sim_threshold={args.sim_threshold}% | "
        f"lha={args.lha_len}bp | rha={args.rha_len}bp | "
        f"max_oligo={args.max_oligo_len}bp{Colors.ENDC}"
    )
