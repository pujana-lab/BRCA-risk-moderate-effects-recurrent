#!/usr/bin/env python3
"""Compare reconstructed recurrent genes with an archived gene-list fingerprint."""
from __future__ import annotations

import argparse
from pathlib import Path
from common import read_table


def read_list(path: str) -> set[str]:
    return {x.strip() for x in Path(path).read_text().splitlines() if x.strip() and not x.startswith("#")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--observed", required=True, help="Recurrence summary TSV with GENE")
    ap.add_argument("--expected", required=True)
    args = ap.parse_args()
    obs = set(read_table(args.observed).GENE.astype(str))
    exp = read_list(args.expected)
    print(f"expected={len(exp)} observed={len(obs)} overlap={len(exp & obs)}")
    print("missing_expected:", ", ".join(sorted(exp - obs)) or "none")
    print("unexpected_observed:", ", ".join(sorted(obs - exp)) or "none")
    if obs != exp:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
