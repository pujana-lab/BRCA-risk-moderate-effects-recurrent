#!/usr/bin/env python3
"""Parse PLINK .clumped output and assign every retained SNP to an index clump."""
from __future__ import annotations

import argparse
import re
import pandas as pd
from common import read_table, write_table


def parse_sp2(value: str) -> list[str]:
    if value is None or pd.isna(value):
        return []
    v = str(value).strip()
    if not v or v == "NONE":
        return []
    out = []
    for token in v.split(","):
        token = token.strip()
        # PLINK commonly writes rs123(1) style annotations.
        token = re.sub(r"\([^)]*\)$", "", token)
        if token and token != "NONE":
            out.append(token)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--clumped", required=True)
    ap.add_argument("--variants", required=True, help="Retained annotated SNP table")
    ap.add_argument("--gwas", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    variants = read_table(args.variants)
    cl = pd.read_csv(args.clumped, sep=r"\s+", dtype=str)
    if "SNP" not in cl.columns:
        raise SystemExit("PLINK .clumped file has no SNP column")

    assignment = {}
    for row in cl.itertuples(index=False):
        idx = str(row.SNP)
        assignment.setdefault(idx, idx)
        sp2 = parse_sp2(getattr(row, "SP2", None))
        for snp in sp2:
            if snp in assignment and assignment[snp] != idx:
                raise SystemExit(f"SNP {snp} assigned to multiple clumps; check overlap settings")
            assignment[snp] = idx

    variants["CLUMP_INDEX_SNP"] = variants["BIM_SNP"].map(assignment)
    missing = variants["CLUMP_INDEX_SNP"].isna()
    if missing.any():
        ex = variants.loc[missing, "BIM_SNP"].head(10).tolist()
        raise SystemExit(
            f"{missing.sum()} retained variants were not present in PLINK clump output. Examples: {ex}. "
            "Check p1/p2, reference membership, or PLINK output."
        )
    variants["GWAS"] = args.gwas
    write_table(variants, args.output)
    print(f"Assigned {len(variants):,} SNPs to {variants.CLUMP_INDEX_SNP.nunique():,} clumps")


if __name__ == "__main__":
    main()
