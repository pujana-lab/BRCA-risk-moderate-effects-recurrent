# Known and unresolved parameters

## Confirmed

- Genome build: GRCh37.
- LD population/reference: 1000 Genomes EUR.
- Candidate threshold: P <= 1e-3.
- Candidate threshold: |BETA| >= 0.25.
- MAF threshold: >= 0.001.
- Compute Z = BETA / SE.
- Accept same or swapped BIM allele ordering.
- Exclude indels.
- Exclude multiallelic reference positions.
- Exclude non-unique/duplicate BIM matches.
- Positional gene maps: gene body, ±10 kb, ±100 kb.
- Earlier strict recurrence: >=2 GWAS and >=2 clumps in each counted GWAS.
- Later/publication-style recurrence: >=2 GWAS and >=1 clump in each counted GWAS.

## Unresolved and intentionally NOT guessed

### Palindromic A/T and C/G SNPs
They are present in the harmonized/gene-annotated TNBC intermediate file. If they were excluded, the exclusion happened later. `exclude_palindromic_before_clumping` therefore defaults to `null` and the pipeline refuses to silently decide.

### PLINK LD-clumping parameters
The original `--clump-r2` and `--clump-kb` values have not been recovered. They are `null` in the example config.

### Top-10 recurrence ranking rule
The final figure identifies ten green genes, but the exact score/tie-break rule has not yet been recovered. Their archived identities are stored as a validation fingerprint only.

## Suggested reconstruction strategy for clumping

Run a parameter grid against the three original cleaned candidate files and compare the resulting gene-level clump counts with archived fingerprints. Candidate values should be explicitly recorded in a reconstruction notebook/log; do not silently choose a PLINK default and call it original.
