import argparse
import sys

import Cas9_knockout_designer_v8 as v8


def _v9_get_cut_site_mt(self, sgrna_genomic_start: int, sgrna_strand: str, pam_len: int) -> int:
    if sgrna_strand == "forward":
        pam_start_pos = sgrna_genomic_start + self.config.guide_len
        return pam_start_pos - 3
    pam_start_pos = sgrna_genomic_start
    return pam_start_pos + pam_len + 3


def _v9_generate_candidate_cut_window_mt(self, sgrna, cds_min, cds_max, strand, cds_5prime, protected_zones):
    cut_site = sgrna['cut_site']
    upstream = self.config.cut_window_upstream
    downstream = self.config.cut_window_downstream
    max_del_len = upstream + downstream

    pam_dir = 1 if sgrna.get('strand') == 'forward' else -1
    upstream_pos = cut_site - pam_dir * upstream
    downstream_pos = cut_site + pam_dir * downstream
    ds = max(min(upstream_pos, downstream_pos), cds_min)
    de = min(max(upstream_pos, downstream_pos), cds_max)

    if not (ds <= cut_site < de):
        return None

    max_available = min(de - ds, max_del_len)
    if max_available <= 0:
        return None

    dl = max_available
    while dl > 0 and dl % 3 == 0:
        dl -= 1
    if dl <= 0:
        return None

    if pam_dir == 1:
        del_start = ds
        del_end = del_start + dl
    else:
        del_end = de
        del_start = del_end - dl

    if del_start < cds_min or del_end > cds_max:
        return None
    if not (del_start <= cut_site < del_end):
        return None
    if not self._check_constraints_fast(del_start, del_end, cds_5prime, strand, protected_zones):
        return None

    return {'del_start': del_start, 'del_end': del_end, 'del_len': dl}


def _parse_v9_front_args():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--dele_model", type=str, default="normal", choices=["normal", "Mt"])
    parser.add_argument("--species", type=str, default="M_thermophila")
    args, _ = parser.parse_known_args()
    return args


def _apply_v9_mode(argv):
    front_args = _parse_v9_front_args()
    model = front_args.dele_model

    cleaned = []
    skip_next = False
    for i, token in enumerate(argv):
        if skip_next:
            skip_next = False
            continue
        if token == "--dele_model":
            skip_next = True
            continue
        if token.startswith("--dele_model="):
            continue
        cleaned.append(token)

    if model == "Mt":
        v8.SGRNADesigner._get_cut_site = _v9_get_cut_site_mt
        v8.CRISPRDesigner._generate_candidate_cut_window = _v9_generate_candidate_cut_window_mt
        cleaned.extend(["--deletion_mode", "cut_window"])
        print("[V9审计] 当前采用 Mt 删除策略")
    else:
        cleaned.extend(["--deletion_mode", "legacy_length"])
        print("[V9审计] 当前采用 normal 删除策略")

    return cleaned


if __name__ == "__main__":
    sys.argv = _apply_v9_mode(sys.argv)
    v8.main()
