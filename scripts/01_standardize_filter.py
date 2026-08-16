#!/usr/bin/env python3
"""Standardize one GWAS and apply the confirmed candidate-variant filters."""
from __future__ import annotations

import argparse
import numpy as np
import pandas as pd

from common import read_yaml, read_table, write_table, normalize_chr


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--gwas", required=True, help="GWAS key from config, e.g. BRCA1")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    gcfg = cfg["gwas"][args.gwas]
    cols = gcfg["columns"]
    df = read_table(gcfg["input"])

    rename = {}
    canonical = {
        "id": "ORIGINAL_ID", "chr": "CHR", "bp": "BP", "a1": "A1", "a2": "A2",
        "beta": "BETA", "se": "SE", "p": "P", "eaf": "EAF", "maf": "MAF", "info": "INFO"
    }
    for key, canon in canonical.items():
        src = cols.get(key)
        if src is not None:
            rename[src] = canon
    df = df.rename(columns=rename)

    required = ["CHR", "BP", "A1", "A2", "BETA", "SE", "P"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing required standardized columns: {missing}")

    if "ORIGINAL_ID" not in df:
        df["ORIGINAL_ID"] = (
            df["CHR"].astype(str) + "_" + df["BP"].astype(str) + "_" +
            df["A1"].astype(str) + "_" + df["A2"].astype(str)
        )

    df["CHR"] = normalize_chr(df["CHR"])
    for c in ["BP", "BETA", "SE", "P"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    if "EAF" in df:
        df["EAF"] = pd.to_numeric(df["EAF"], errors="coerce")
    if "MAF" in df:
        df["MAF"] = pd.to_numeric(df["MAF"], errors="coerce")
    elif "EAF" in df:
        df["MAF"] = np.minimum(df["EAF"], 1.0 - df["EAF"])
    else:
        raise SystemExit("MAF is required, or EAF must be available so MAF can be derived.")

    df["A1"] = df["A1"].astype(str).str.upper()
    df["A2"] = df["A2"].astype(str).str.upper()
    df["ABS_BETA"] = df["BETA"].abs()
    df["Z"] = df["BETA"] / df["SE"]
    df["VALID_SE"] = df["SE"].notna() & (df["SE"] > 0)
    df["VALID_P"] = df["P"].notna() & df["P"].between(0, 1, inclusive="both")
    df["VALID_MAF"] = df["MAF"].notna() & df["MAF"].between(0, 0.5, inclusive="both")

    fcfg = cfg["variant_filter"]
    keep = (
        df["VALID_SE"] & df["VALID_P"] & df["VALID_MAF"] &
        (df["P"] <= float(fcfg["p_max"])) &
        (df["ABS_BETA"] >= float(fcfg["abs_beta_min"])) &
        (df["MAF"] >= float(fcfg["maf_min"]))
    )

    info_min = gcfg.get("info_min", fcfg.get("default_info_min"))
    if info_min is not None:
        if "INFO" not in df:
            raise SystemExit("INFO threshold requested but no INFO column is configured.")
        df["INFO"] = pd.to_numeric(df["INFO"], errors="coerce")
        keep &= df["INFO"].ge(float(info_min))

    out = df.loc[keep].copy()
    out.insert(0, "PHENOTYPE", args.gwas)
    write_table(out, args.output)
    print(f"{args.gwas}: retained {len(out):,}/{len(df):,} variants")


if __name__ == "__main__":
    main()
