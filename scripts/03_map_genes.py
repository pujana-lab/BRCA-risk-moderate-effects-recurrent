#!/usr/bin/env python3
"""Annotate harmonized SNPs to gene body, ±10 kb and ±100 kb."""
from __future__ import annotations

import argparse
import numpy as np
import pandas as pd

from common import normalize_chr, read_table, write_table


def map_window(snps: pd.DataFrame, genes: pd.DataFrame, window: int):
    syms, entrez, counts = [], [], []
    by_chr = {c: g.reset_index(drop=True) for c, g in genes.groupby("CHR", sort=False)}
    for row in snps.itertuples(index=False):
        g = by_chr.get(str(row.CHR))
        if g is None:
            hits = g
        else:
            pos = int(row.BP)
            hits = g[(g.START - window <= pos) & (g.END + window >= pos)]
        if hits is None or hits.empty:
            counts.append(0); syms.append(""); entrez.append("")
        else:
            h = hits.sort_values(["START", "END", "SYMBOL"]).drop_duplicates("SYMBOL")
            counts.append(len(h))
            syms.append(";".join(h.SYMBOL.astype(str)))
            entrez.append(";".join(h.ENTREZ_ID.astype(str)))
    return counts, syms, entrez


def nearest_gene(snps: pd.DataFrame, genes: pd.DataFrame):
    syms, entrez, dist = [], [], []
    by_chr = {c: g.reset_index(drop=True) for c, g in genes.groupby("CHR", sort=False)}
    for row in snps.itertuples(index=False):
        g = by_chr.get(str(row.CHR))
        if g is None or g.empty:
            syms.append(""); entrez.append(""); dist.append(np.nan); continue
        pos = int(row.BP)
        d = np.where(pos < g.START, g.START - pos, np.where(pos > g.END, pos - g.END, 0))
        m = np.min(d)
        h = g[d == m].sort_values("SYMBOL").drop_duplicates("SYMBOL")
        syms.append(";".join(h.SYMBOL.astype(str)))
        entrez.append(";".join(h.ENTREZ_ID.astype(str)))
        dist.append(int(m))
    return syms, entrez, dist


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    snps = read_table(args.input)
    genes = read_table(args.genes)
    required = {"CHR", "START", "END", "SYMBOL", "ENTREZ_ID"}
    if not required.issubset(genes.columns):
        raise SystemExit(f"Gene reference must contain {sorted(required)}")
    genes["CHR"] = normalize_chr(genes["CHR"])
    genes["START"] = pd.to_numeric(genes["START"], errors="raise").astype(int)
    genes["END"] = pd.to_numeric(genes["END"], errors="raise").astype(int)

    for window, suffix in [(0, "GENE_BODY"), (10000, "GENES_WITHIN_10KB"), (100000, "GENES_WITHIN_100KB")]:
        n, s, e = map_window(snps, genes, window)
        if window == 0:
            snps["N_GENE_BODY"] = n
            snps["GENE_BODY_SYMBOLS"] = s
            snps["GENE_BODY_ENTREZ_IDS"] = e
        else:
            snps[f"N_{suffix}"] = n
            snps[suffix] = s
            snps["ENTREZ_" + suffix.replace("GENES_", "")] = e

    s, e, d = nearest_gene(snps, genes)
    snps["NEAREST_GENE_SYMBOLS"] = s
    snps["NEAREST_GENE_ENTREZ_IDS"] = e
    snps["NEAREST_GENE_DISTANCE_BP"] = d
    write_table(snps, args.output)
    print(f"Annotated {len(snps):,} variants")


if __name__ == "__main__":
    main()
