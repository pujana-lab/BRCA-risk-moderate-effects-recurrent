#!/usr/bin/env python3
"""Convert SNP clumps into gene × GWAS clump counts for one mapping mode."""
from __future__ import annotations

import argparse
import pandas as pd
from common import read_table, split_genes, write_table

MAP_COLUMNS = {
    "body": "GENE_BODY_SYMBOLS",
    "10kb": "GENES_WITHIN_10KB",
    "100kb": "GENES_WITHIN_100KB",
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--mapping", choices=sorted(MAP_COLUMNS), required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    df = read_table(args.input)
    col = MAP_COLUMNS[args.mapping]
    required = {"GWAS", "BIM_SNP", "CLUMP_INDEX_SNP", "ABS_BETA", col}
    if not required.issubset(df.columns):
        raise SystemExit(f"Missing columns: {sorted(required - set(df.columns))}")

    rows = []
    for r in df.itertuples(index=False):
        genes = split_genes(getattr(r, col))
        for gene in genes:
            rows.append({
                "GENE": gene,
                "GWAS": r.GWAS,
                "SNP": r.BIM_SNP,
                "CLUMP_INDEX_SNP": r.CLUMP_INDEX_SNP,
                "ABS_BETA": float(r.ABS_BETA),
            })
    long = pd.DataFrame(rows)
    if long.empty:
        write_table(pd.DataFrame(columns=["GENE","GWAS","N_CLUMPS","MAX_ABS_BETA","N_MAPPED_SNPS"]), args.output)
        return

    out = (long.groupby(["GENE", "GWAS"], as_index=False)
           .agg(N_CLUMPS=("CLUMP_INDEX_SNP", "nunique"),
                MAX_ABS_BETA=("ABS_BETA", "max"),
                N_MAPPED_SNPS=("SNP", "nunique")))
    write_table(out, args.output)
    print(f"Summarized {out.GENE.nunique():,} genes in {out.GWAS.nunique()} GWAS")


if __name__ == "__main__":
    main()
