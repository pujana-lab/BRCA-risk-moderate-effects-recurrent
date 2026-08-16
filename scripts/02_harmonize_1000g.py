#!/usr/bin/env python3
"""Harmonize filtered variants to a GRCh37 1000G EUR BIM reference."""
from __future__ import annotations

import argparse
from collections import defaultdict
import pandas as pd

from common import allele_key, is_indel_alleles, is_palindromic, normalize_chr, read_table, read_yaml, write_table


def read_bim(path: str) -> pd.DataFrame:
    bim = pd.read_csv(path, sep=r"\s+", header=None,
                      names=["BIM_CHR", "BIM_SNP", "CM", "BIM_BP", "BIM_A1", "BIM_A2"],
                      dtype={0: str, 1: str, 4: str, 5: str})
    bim["BIM_CHR"] = normalize_chr(bim["BIM_CHR"])
    bim["BIM_BP"] = pd.to_numeric(bim["BIM_BP"], errors="coerce").astype("Int64")
    bim["BIM_A1"] = bim["BIM_A1"].str.upper()
    bim["BIM_A2"] = bim["BIM_A2"].str.upper()
    bim["PAIR_KEY"] = [allele_key(a, b) for a, b in zip(bim.BIM_A1, bim.BIM_A2)]
    return bim


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--bim", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--qc-report", required=False)
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    hcfg = cfg["harmonization"]
    df = read_table(args.input)
    bim = read_bim(args.bim)

    # Mark reference positions with >1 distinct allele pair as multiallelic.
    pos_pairs = bim.groupby(["BIM_CHR", "BIM_BP"])["PAIR_KEY"].nunique()
    multiallelic_pos = set(pos_pairs[pos_pairs > 1].index)

    # Match by chromosome, position and unordered allele pair.
    idx = defaultdict(list)
    for i, row in bim.iterrows():
        idx[(row.BIM_CHR, int(row.BIM_BP), row.PAIR_KEY)].append(i)

    records = []
    excluded = []
    for _, row in df.iterrows():
        chrom = str(row.CHR)
        bp = int(row.BP)
        a1, a2 = str(row.A1).upper(), str(row.A2).upper()
        pal = is_palindromic(a1, a2)
        indel = is_indel_alleles(a1, a2)
        multi = (chrom, bp) in multiallelic_pos
        hits = idx.get((chrom, bp, allele_key(a1, a2)), [])

        reasons = []
        if hcfg.get("remove_indels", True) and indel:
            reasons.append("INDEL")
        if hcfg.get("remove_multiallelic_positions", True) and multi:
            reasons.append("MULTIALLELIC_POSITION")
        if hcfg.get("require_unique_bim_match", True) and len(hits) != 1:
            reasons.append("NON_UNIQUE_BIM_MATCH" if hits else "NO_BIM_MATCH")

        if reasons:
            excluded.append({**row.to_dict(), "QC_EXCLUSION_REASON": ";".join(reasons)})
            continue
        hit = bim.loc[hits[0]]
        same = a1 == hit.BIM_A1 and a2 == hit.BIM_A2
        swapped = a1 == hit.BIM_A2 and a2 == hit.BIM_A1
        if not same and not swapped:
            excluded.append({**row.to_dict(), "QC_EXCLUSION_REASON": "ALLELE_MISMATCH"})
            continue
        if swapped and not hcfg.get("allow_swapped_allele_order", True):
            excluded.append({**row.to_dict(), "QC_EXCLUSION_REASON": "SWAPPED_ORDER_DISALLOWED"})
            continue

        rec = row.to_dict()
        rec.update({
            "BIM_SNP": hit.BIM_SNP,
            "BIM_CHR": hit.BIM_CHR,
            "BIM_BP": int(hit.BIM_BP),
            "BIM_A1": hit.BIM_A1,
            "BIM_A2": hit.BIM_A2,
            "MATCH_METHOD": "POSITION_ALLELES" if same else "POSITION_ALLELES_SWAPPED",
            "ALLELE_ORIENTATION": "SAME" if same else "SWAPPED",
            "MATCH_STATUS": "MATCHED_UNIQUE",
            "PALINDROMIC": pal,
            "IS_INDEL": indel,
            "MULTIALLELIC_POSITION": multi,
            "DUPLICATE_BIM_STATUS": "NONE",
            "QC_EXCLUSION_REASON": "",
        })
        records.append(rec)

    out = pd.DataFrame(records)
    write_table(out, args.output)
    if args.qc_report:
        write_table(pd.DataFrame(excluded), args.qc_report)
    print(f"Harmonized: {len(out):,}; excluded: {len(excluded):,}")


if __name__ == "__main__":
    main()
