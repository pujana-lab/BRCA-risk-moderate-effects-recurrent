#!/usr/bin/env python3
"""Add prior-risk and optional archived top-recurrence display categories."""
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd
from common import read_table, write_table


def load_gene_list(path: str | None) -> set[str]:
    if not path:
        return set()
    p = Path(path)
    if not p.exists():
        return set()
    genes = set()
    for line in p.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            genes.add(line.split("\t")[0])
    return genes


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--prior-risk", required=False,
                    help="TSV containing a GENE column; multiple sources/categories allowed")
    ap.add_argument("--top10-reference", required=False,
                    help="Archived validation list; not a reconstructed ranking formula")
    ap.add_argument("--use-archived-top10-for-display", action="store_true")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    df = read_table(args.input)
    prior = set()
    if args.prior_risk and Path(args.prior_risk).exists():
        pr = read_table(args.prior_risk)
        if "GENE" not in pr:
            raise SystemExit("Prior-risk TSV must contain GENE")
        prior = set(pr.GENE.dropna().astype(str))
    top = load_gene_list(args.top10_reference) if args.use_archived_top10_for_display else set()

    df["IS_PRIOR_RISK"] = df.GENE.isin(prior)
    df["IS_TOP10_RECURRENCE"] = df.GENE.isin(top)
    # Match the publication legend precedence: prior-risk is orange even if also top-ranked.
    df["GENE_CATEGORY"] = "Other qualifying recurrent genes"
    df.loc[df.IS_TOP10_RECURRENCE, "GENE_CATEGORY"] = "Top 10 recurrence-ranked genes"
    df.loc[df.IS_PRIOR_RISK, "GENE_CATEGORY"] = "Prior-risk-set genes"
    write_table(df, args.output)


if __name__ == "__main__":
    main()
