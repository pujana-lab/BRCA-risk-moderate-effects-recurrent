# BRCA cross-GWAS recurrent-gene pipeline

A **data-free, reproducible workflow** for reconstructing the BRCA1/BRCA2/TNBC cross-GWAS recurrent-gene analysis and publication figure.

The repository is designed to contain code, configuration, documentation, and validation fingerprints only. Original GWAS summary statistics, 1000 Genomes genotype files, and external annotation datasets stay outside version control.

## What the pipeline does

```text
3 GWAS summary-statistic files
        │
        ├── standardize columns
        ├── compute |BETA| and Z = BETA/SE
        └── filter P <= 1e-3, |BETA| >= 0.25, MAF >= 0.001
                         │
                         ▼
          1000G EUR / GRCh37 harmonization
          - unique position + allele-pair match
          - same OR swapped BIM allele order accepted
          - indels excluded
          - multiallelic/non-unique positions excluded
                         │
                         ▼
             gene mapping: body / ±10 kb / ±100 kb
                         │
                         ▼
              [explicit palindromic policy]
                         │
                         ▼
              PLINK LD clumping with 1000G EUR
                         │
                         ▼
              gene × GWAS independent-clump counts
                         │
                         ▼
               cross-GWAS recurrent-gene summary
                         │
                         ▼
               prior-risk annotation + plotting
```

## Confirmed thresholds

```yaml
p_max: 1.0e-3
abs_beta_min: 0.25
maf_min: 0.001
```

The recovered BRCA1–BCAC TNBC intermediate table contains **4,394 variants**, with max P = 0.001, min |BETA| = 0.2501, and min MAF = 0.001, matching the original candidate count.

## Important QC behavior

The original harmonization accepted both **same** and **swapped** BIM allele order. The recovered TNBC table contains 4,356 `POSITION_ALLELES_SWAPPED` matches and 38 same-order matches.

Indels, multiallelic positions, and non-unique reference matches were excluded.

**Palindromic A/T and C/G handling is deliberately unresolved.** Palindromic variants are still present in the harmonized/gene-annotated intermediate file, so any exclusion occurred later. The config requires the user to choose `true` or `false` before clumping rather than silently guessing.

## LD-clumping parameters: not yet recovered

The exact original PLINK `--clump-r2` and `--clump-kb` settings are currently unknown. They are therefore `null` in `config/config.example.yaml`, and `run_plink_clump.sh` refuses to run until they are explicitly supplied.

Use the archived strict recurrent-gene fingerprints in `reference/` to reconstruct the settings rather than assuming PLINK defaults.

## Two recurrence profiles

### Publication-style figure

A gene qualifies if it occurs in at least **2 GWAS** with at least **1 LD-independent clump in each qualifying GWAS**.

The plotted gene-level point uses:

- x = sum of independent clumps across qualifying GWAS;
- y = maximum |BETA| across qualifying GWAS;
- shape = GWAS combination;
- color = recurrent/prior-risk/top-recurrence category.

### Earlier strict analysis

A gene qualifies if it occurs in >1 GWAS and has **>=2 LD-independent clumps in every counted GWAS**.

Archived strict fingerprints:

- ±10 kb: KCNAB1, TRPM3
- ±100 kb: COX7A2L, EML4, KCNAB1, KCNG3, MTA3, RCN1, STK39, TMEM106B, TRPM3, TTC28, VWDE

## Repository setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
cp config/config.example.yaml config/config.yaml
```

Place external inputs under `data/` (ignored by git), then edit `config/config.yaml` with the real paths and raw-column mappings.

A PLINK 1.9-compatible executable must be available as `plink` or configured under `ld_clumping.plink_executable`.

## Example: one GWAS branch

```bash
python scripts/01_standardize_filter.py \
  --config config/config.yaml --gwas BRCA1 \
  --output results/BRCA1/01_filtered.tsv

python scripts/02_harmonize_1000g.py \
  --config config/config.yaml \
  --input results/BRCA1/01_filtered.tsv \
  --bim data/reference/g1000_eur.bim \
  --output results/BRCA1/02_g1000_clean.tsv \
  --qc-report results/BRCA1/02_excluded.tsv

python scripts/03_map_genes.py \
  --input results/BRCA1/02_g1000_clean.tsv \
  --genes data/reference/genes_grch37.tsv \
  --output results/BRCA1/03_gene_annotated.tsv

python scripts/04_prepare_clump_input.py \
  --config config/config.yaml \
  --input results/BRCA1/03_gene_annotated.tsv \
  --output results/BRCA1/04_clump_input.tsv \
  --retained-output results/BRCA1/04_clump_variants.tsv

scripts/run_plink_clump.sh \
  config/config.yaml results/BRCA1/04_clump_input.tsv \
  results/BRCA1/05_plink BRCA1

python scripts/05_parse_plink_clumps.py \
  --clumped results/BRCA1/05_plink.clumped \
  --variants results/BRCA1/04_clump_variants.tsv \
  --gwas BRCA1 \
  --output results/BRCA1/05_clump_assignments.tsv
```

Repeat for BRCA2 and BRCA1_TNBC, then summarize the desired mapping window:

```bash
python scripts/06_summarize_gene_clumps.py --input results/BRCA1/05_clump_assignments.tsv --mapping 100kb --output results/BRCA1/06_gene_clumps_100kb.tsv
python scripts/06_summarize_gene_clumps.py --input results/BRCA2/05_clump_assignments.tsv --mapping 100kb --output results/BRCA2/06_gene_clumps_100kb.tsv
python scripts/06_summarize_gene_clumps.py --input results/BRCA1_TNBC/05_clump_assignments.tsv --mapping 100kb --output results/BRCA1_TNBC/06_gene_clumps_100kb.tsv

python scripts/07_recurrence_summary.py \
  --inputs results/BRCA1/06_gene_clumps_100kb.tsv results/BRCA2/06_gene_clumps_100kb.tsv results/BRCA1_TNBC/06_gene_clumps_100kb.tsv \
  --min-gwas 2 --min-clumps 1 \
  --output results/recurrent_100kb.tsv \
  --qualifying-gwas-output results/recurrent_100kb_by_gwas.tsv
```

Annotate and plot:

```bash
python scripts/08_annotate_gene_categories.py \
  --input results/recurrent_100kb.tsv \
  --prior-risk data/annotations/prior_risk_union.tsv \
  --top10-reference reference/archived_top10_100kb.txt \
  --use-archived-top10-for-display \
  --output results/recurrent_100kb_annotated.tsv

python scripts/09_plot_publication.py \
  --input results/recurrent_100kb_annotated.tsv \
  --window-kb 100 \
  --svg results/figures/crossGWAS_recurrent_genes_100kb.svg \
  --png results/figures/crossGWAS_recurrent_genes_100kb.png
```

The `--use-archived-top10-for-display` flag is only for visual validation while the original top-10 ranking rule remains unrecovered.

## Reconstructing the LD parameters

Generate a candidate grid:

```bash
python scripts/search_clump_parameters.py \
  --r2 0.1,0.2,0.3,0.5 \
  --kb 250,500,1000 \
  --output results/clump_parameter_grid.tsv
```

For each candidate setting, rerun clumping and the **strict** recurrence profile (`--min-clumps 2`), then compare against:

```bash
python scripts/10_compare_fingerprint.py \
  --observed results/strict_recurrent_100kb.tsv \
  --expected reference/expected_strict_recurrent_genes_100kb.txt
```

Do not label a parameter combination “original” merely because it is a PLINK default; require concordance with archived outputs.

## Hi-C extension

The repository reserves an optional third SNP→gene mapping mode for the later normal-human-breast Hi-C analysis. It is intentionally separated from the ±10/±100 kb code because the attached publication figure is positional and because the exact final GSE261230 SNP→gene table should be supplied externally.

## Reproducibility status

**Recovered:** candidate thresholds, Z calculation, GRCh37/1000G EUR harmonization logic, swapped-allele acceptance, indel/multiallelic/non-unique exclusions, positional mapping windows, strict and relaxed recurrence definitions, gene-level plot axes, shape logic, category logic, and label criteria.

**Still to recover:** palindromic SNP policy before clumping, exact PLINK r²/window parameters, and exact top-10 recurrence ranking formula.
