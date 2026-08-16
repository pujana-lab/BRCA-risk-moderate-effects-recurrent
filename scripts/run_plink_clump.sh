#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "Usage: $0 CONFIG_YAML CLUMP_INPUT_TSV OUT_PREFIX GWAS_LABEL" >&2
  exit 2
fi
CONFIG=$1
ASSOC=$2
OUT=$3
LABEL=$4

read_cfg() {
  python - "$CONFIG" "$1" <<'PY'
import sys, yaml
cfg=yaml.safe_load(open(sys.argv[1]))
x=cfg
for part in sys.argv[2].split('.'):
    x=x[part]
print("" if x is None else str(x))
PY
}

PLINK=$(read_cfg ld_clumping.plink_executable)
REF=$(read_cfg ld_clumping.reference_prefix)
P1=$(read_cfg ld_clumping.p1)
P2=$(read_cfg ld_clumping.p2)
R2=$(read_cfg ld_clumping.r2)
KB=$(read_cfg ld_clumping.kb)
OVERLAP=$(read_cfg ld_clumping.allow_overlap)

if [[ -z "$R2" || -z "$KB" ]]; then
  echo "ERROR: ld_clumping.r2 and ld_clumping.kb must be explicitly set; original values are unresolved." >&2
  exit 3
fi

mkdir -p "$(dirname "$OUT")"
cmd=("$PLINK" --bfile "$REF" --clump "$ASSOC" --clump-snp-field SNP --clump-field P \
     --clump-p1 "$P1" --clump-p2 "$P2" --clump-r2 "$R2" --clump-kb "$KB" --out "$OUT")
if [[ "$OVERLAP" == "True" || "$OVERLAP" == "true" ]]; then
  cmd+=(--clump-allow-overlap)
fi
printf 'Running %q ' "${cmd[@]}"; echo
"${cmd[@]}"
