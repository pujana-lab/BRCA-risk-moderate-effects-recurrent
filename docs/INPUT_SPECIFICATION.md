# Input specification

## 1. GWAS files

Raw source files may use different headers. Map each source into these canonical fields in `config/config.yaml`:

| Canonical field | Required | Meaning |
|---|---:|---|
| ID | yes | Original variant identifier |
| CHR | yes | Chromosome, GRCh37 |
| BP | yes | 1-based base-pair position |
| A1 | yes | Effect allele in the GWAS |
| A2 | yes | Other allele |
| BETA | yes | Log-effect estimate |
| SE | yes | Standard error of BETA |
| P | yes | Association P value |
| EAF | recommended | Effect-allele frequency |
| MAF | yes or derivable | Minor-allele frequency |
| INFO | optional | Imputation quality when applicable |

`01_standardize_filter.py` derives MAF from EAF when MAF is not supplied.

## 2. 1000 Genomes EUR PLINK reference

The repository expects an external GRCh37 EUR reference prefix such as:

```text
data/reference/g1000_eur.bed
data/reference/g1000_eur.bim
data/reference/g1000_eur.fam
```

The `.bim` file is used for harmonization and the full PLINK reference is used for LD clumping.

## 3. Gene coordinate reference

A tab-separated file with:

```text
CHR  START  END  SYMBOL  ENTREZ_ID
```

Coordinates must be GRCh37 and use 1-based inclusive gene boundaries for this implementation.

## 4. Prior-risk annotation

A data-free repository should not bundle proprietary/restricted source tables. Build an external union table with at least:

```text
GENE  SOURCE  CATEGORY
```

Multiple rows per gene are allowed. `CATEGORY` may contain labels such as causal, established, possible, new, GWAS-catalog, or literature-derived. The plotting script collapses any listed gene to the `prior-risk-set` display category unless a different policy is configured.

## 5. Optional human normal-breast Hi-C map

The positional publication figure can be reproduced without Hi-C. To support the later third mapping mode, supply a precomputed table such as:

```text
CHR  BP  SNP  GENE  SAMPLE  CONTACT_SCORE
```

or adapt the mapping step to the final GSE261230-derived format.
