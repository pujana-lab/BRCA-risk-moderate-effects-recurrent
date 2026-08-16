#!/usr/bin/env python3
"""Generate a parameter-grid manifest for reconstructing unknown PLINK clumping settings.

This script intentionally does not claim any candidate is the original setting. It writes
candidate combinations that can be run through the pipeline and scored against archived
strict recurrent-gene fingerprints.
"""
from __future__ import annotations

import argparse
import itertools
import pandas as pd


def parse_list(x, cast):
    return [cast(v.strip()) for v in x.split(",") if v.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r2", default="0.1,0.2,0.3,0.5")
    ap.add_argument("--kb", default="250,500,1000")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    rows = [{"r2": r, "kb": kb} for r, kb in itertools.product(parse_list(args.r2, float), parse_list(args.kb, int))]
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(rows)} candidate parameter combinations. These are search candidates, not recovered original values.")


if __name__ == "__main__":
    main()
