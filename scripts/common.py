#!/usr/bin/env python3
from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterable

import pandas as pd
import yaml


def read_yaml(path: str | Path) -> dict:
    with open(path, "rt", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def read_table(path: str | Path, **kwargs) -> pd.DataFrame:
    path = Path(path)
    return pd.read_csv(path, sep="\t", low_memory=False, **kwargs)


def write_table(df: pd.DataFrame, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def normalize_chr(series: pd.Series) -> pd.Series:
    s = series.astype(str).str.replace(r"^chr", "", regex=True)
    return s.replace({"23": "X", "24": "Y", "25": "XY", "26": "MT", "M": "MT"})


def is_indel_alleles(a1: str, a2: str) -> bool:
    return len(str(a1)) != 1 or len(str(a2)) != 1


def is_palindromic(a1: str, a2: str) -> bool:
    pair = {str(a1).upper(), str(a2).upper()}
    return pair in ({"A", "T"}, {"C", "G"})


def allele_key(a1: str, a2: str) -> tuple[str, str]:
    return tuple(sorted((str(a1).upper(), str(a2).upper())))


def split_genes(value) -> list[str]:
    if pd.isna(value) or value == "":
        return []
    return [x.strip() for x in str(value).split(";") if x.strip()]
