#!/usr/bin/env python3
"""Generate the one-point-per-gene cross-GWAS recurrent-gene publication plot."""
from __future__ import annotations

import argparse
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd

from common import read_table

# These visual settings are chosen to match the archived figure structure.
# Exact rendering can differ by matplotlib/font version.
SHAPES = {
    "BRCA1 + BRCA2 + BRCA1–TNBC": "o",
    "BRCA1 + BRCA1–TNBC": "s",
    "BRCA2 + BRCA1–TNBC": "^",
    "BRCA1 + BRCA2": "D",
}
COLORS = {
    "Other qualifying recurrent genes": "#bdbdbd",
    "Prior-risk-set genes": "#f28e0b",
    "Top 10 recurrence-ranked genes": "#31ad43",
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--window-kb", type=int, default=100)
    ap.add_argument("--label-beta-gt", type=float, default=2.0)
    ap.add_argument("--label-clumps-gt", type=int, default=4)
    ap.add_argument("--svg", required=True)
    ap.add_argument("--png", required=False)
    args = ap.parse_args()

    df = read_table(args.input)
    if "GENE_CATEGORY" not in df:
        df["GENE_CATEGORY"] = "Other qualifying recurrent genes"

    fig, ax = plt.subplots(figsize=(12.8, 9.6))
    for _, r in df.iterrows():
        marker = SHAPES.get(r.GWAS_COMBINATION, "o")
        color = COLORS.get(r.GENE_CATEGORY, COLORS["Other qualifying recurrent genes"])
        ax.scatter(r.TOTAL_CLUMPS, r.MAX_ABS_BETA, s=70, marker=marker,
                   facecolor=color, edgecolor="black", linewidth=0.7, zorder=3)

    should_label = (
        (df.GENE_CATEGORY != "Other qualifying recurrent genes") |
        (df.MAX_ABS_BETA > args.label_beta_gt) |
        (df.TOTAL_CLUMPS > args.label_clumps_gt)
    )
    labels = df.loc[should_label].copy()
    # Deterministic lightweight offsets; manual fine-tuning can be added in an external label-offset TSV.
    offsets = [(8, 8), (8, -12), (-40, 8), (8, 14), (-45, -12), (12, -18)]
    for i, (_, r) in enumerate(labels.iterrows()):
        dx, dy = offsets[i % len(offsets)]
        bold = r.GENE_CATEGORY == "Top 10 recurrence-ranked genes"
        ax.annotate(r.GENE, (r.TOTAL_CLUMPS, r.MAX_ABS_BETA), xytext=(dx, dy),
                    textcoords="offset points", fontsize=8.5,
                    fontweight="bold" if bold else "normal",
                    arrowprops=dict(arrowstyle="-", lw=0.6, color="black"))

    ax.set_title(f"Cross-GWAS recurrent genes at ±{args.window_kb} kb", fontsize=18, pad=12)
    ax.set_xlabel(f"Total mapped LD-independent clumps across qualifying GWAS (±{args.window_kb} kb)", fontsize=12)
    ax.set_ylabel("Maximum |BETA| across qualifying GWAS", fontsize=12)
    ax.set_xlim(left=max(0.9, df.TOTAL_CLUMPS.min() - 0.9) if len(df) else 0.9)
    ax.set_ylim(bottom=0)

    shape_handles = [Line2D([0],[0], marker=m, color="black", markerfacecolor="white",
                            linestyle="None", markersize=8, label=k) for k,m in SHAPES.items()]
    leg1 = ax.legend(handles=shape_handles, title="GWAS combination (shape)",
                     loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=True)
    ax.add_artist(leg1)
    color_order = ["Other qualifying recurrent genes", "Prior-risk-set genes", "Top 10 recurrence-ranked genes"]
    color_handles = [Line2D([0],[0], marker="o", color="black", markerfacecolor=COLORS[k],
                            linestyle="None", markersize=8, label=k) for k in color_order]
    ax.legend(handles=color_handles, title="Gene category (color)",
              loc="upper left", bbox_to_anchor=(1.01, 0.64), frameon=True)

    fig.text(0.08, 0.02,
             f"Labels include all prior-risk-set genes, all displayed top 10 recurrence-ranked genes, "
             f"and every gene with maximum |BETA| > {args.label_beta_gt:g} or total LD clumps > {args.label_clumps_gt}.",
             fontsize=8.5)
    fig.subplots_adjust(right=0.75, bottom=0.09)

    svg = Path(args.svg); svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(svg, format="svg", bbox_inches="tight")
    if args.png:
        png = Path(args.png); png.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(png, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
