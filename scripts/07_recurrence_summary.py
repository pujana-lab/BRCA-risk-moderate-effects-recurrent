#!/usr/bin/env python3
"""Combine gene × GWAS summaries into one point per recurrent gene."""
from __future__ import annotations

import argparse
import pandas as pd
from common import read_table, write_table

CANONICAL_GWAS = ["BRCA1", "BRCA2", "BRCA1_TNBC"]
DISPLAY = {"BRCA1": "BRCA1", "BRCA2": "BRCA2", "BRCA1_TNBC": "BRCA1–TNBC"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", required=True)
    ap.add_argument("--min-gwas", type=int, default=2)
    ap.add_argument("--min-clumps", type=int, default=1)
    ap.add_argument("--output", required=True)
    ap.add_argument("--qualifying-gwas-output", required=False)
    args = ap.parse_args()

    df = pd.concat([read_table(x) for x in args.inputs], ignore_index=True)
    df = df.loc[df.N_CLUMPS >= args.min_clumps].copy()

    q = df.groupby("GENE")["GWAS"].nunique()
    genes = set(q[q >= args.min_gwas].index)
    qdf = df[df.GENE.isin(genes)].copy()

    rows = []
    for gene, g in qdf.groupby("GENE"):
        present = [x for x in CANONICAL_GWAS if x in set(g.GWAS)]
        # Preserve unknown labels deterministically if present.
        present += sorted(set(g.GWAS) - set(CANONICAL_GWAS))
        rows.append({
            "GENE": gene,
            "N_QUALIFYING_GWAS": g.GWAS.nunique(),
            "TOTAL_CLUMPS": int(g.N_CLUMPS.sum()),
            "MAX_ABS_BETA": float(g.MAX_ABS_BETA.max()),
            "GWAS_COMBINATION": " + ".join(DISPLAY.get(x, x) for x in present),
            "GWAS_KEYS": ";".join(present),
        })
    out = pd.DataFrame(rows).sort_values(
        ["N_QUALIFYING_GWAS", "TOTAL_CLUMPS", "MAX_ABS_BETA", "GENE"],
        ascending=[False, False, False, True]
    ) if rows else pd.DataFrame(columns=["GENE","N_QUALIFYING_GWAS","TOTAL_CLUMPS","MAX_ABS_BETA","GWAS_COMBINATION","GWAS_KEYS"])
    write_table(out, args.output)
    if args.qualifying_gwas_output:
        write_table(qdf.sort_values(["GENE", "GWAS"]), args.qualifying_gwas_output)
    print(f"Recurrent genes: {len(out):,}")


if __name__ == "__main__":
    main()
