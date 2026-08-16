#!/usr/bin/env python3
"""Prepare PLINK clump input and apply the unresolved palindromic policy explicitly."""
from __future__ import annotations

import argparse
import pandas as pd
from common import read_table, read_yaml, write_table


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--retained-output", required=True,
                    help="Annotated SNP table actually passed to clumping")
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    policy = cfg["harmonization"].get("exclude_palindromic_before_clumping")
    if policy is None:
        raise SystemExit(
            "exclude_palindromic_before_clumping is unresolved. Set it explicitly to true or false in config."
        )

    df = read_table(args.input)
    if policy:
        df = df.loc[~df["PALINDROMIC"].astype(bool)].copy()
    if df["BIM_SNP"].duplicated().any():
        raise SystemExit("BIM_SNP is not unique in clump input; resolve duplicates before clumping.")

    cl = df[["BIM_SNP", "P"]].rename(columns={"BIM_SNP": "SNP"})
    write_table(cl, args.output)
    write_table(df, args.retained_output)
    print(f"Prepared {len(cl):,} variants for PLINK clumping (palindromic excluded={policy})")


if __name__ == "__main__":
    main()
