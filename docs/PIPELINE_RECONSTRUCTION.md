# Reconstructed analysis pipeline

This repository reconstructs the pipeline used to derive the cross-GWAS recurrent-gene figure from three BRCA-related GWAS summary-statistic datasets. No original GWAS, genotype reference panel, or restricted source data are stored in the repository.

## Confirmed workflow

1. Import and standardize three GWAS summary-statistic files:
   - CIMBA BRCA1 carriers
   - CIMBA BRCA2 carriers
   - BRCA1 carrier + BCAC TNBC meta-analysis
2. Compute `ABS_BETA = abs(BETA)` and `Z = BETA / SE`.
3. Retain candidate variants with:
   - `P <= 1e-3`
   - `ABS_BETA >= 0.25`
   - `MAF >= 0.001`
4. Harmonize against the 1000 Genomes EUR GRCh37 BIM/reference:
   - require a unique position + allele-pair match;
   - accept direct or swapped BIM allele order;
   - exclude indels;
   - exclude multiallelic positions;
   - exclude duplicate/non-unique BIM matches.
5. Annotate candidate SNPs to genes using gene body, ±10 kb, and ±100 kb mappings.
6. Perform LD clumping using the same 1000 Genomes EUR reference.
7. Convert SNP clumps to gene × GWAS independent-clump counts.
8. Identify recurrent genes across at least two GWAS.
9. Integrate prior-risk annotations.
10. Generate the publication-style figure with:
    - x = total mapped LD-independent clumps across qualifying GWAS;
    - y = maximum absolute beta across qualifying GWAS;
    - point shape = GWAS combination;
    - point color = gene category;
    - labels for prior-risk genes, top recurrence-ranked genes, genes with max |BETA| > 2, or genes with total clumps > 4.

## Important distinction: two recurrence definitions existed

### Publication-style recurrent set
A gene qualifies when it is represented in at least two GWAS and has at least one independent clump in every counted GWAS.

### Earlier strict recurrent set
A gene qualifies only when it is represented in more than one GWAS and has at least two independent clumps in every counted GWAS.

The strict archived set is retained under `reference/` as a regression fingerprint.

## Harmonization evidence recovered from the TNBC intermediate file

The archived file `BRCA1_BCAC_TNBC_META_P_le_1e-3_g1000_clean_gene_annotated.tsv` contains 4,394 rows, matching the remembered final TNBC candidate count. In that intermediate file:

- maximum P = 0.001;
- minimum ABS_BETA = 0.2501;
- minimum MAF = 0.001;
- all 4,394 rows are unique 1000G matches;
- all retained rows are non-indels and non-multiallelic;
- 4,356 rows use swapped BIM allele order and 38 use the same allele order;
- 677 palindromic A/T or C/G SNPs are still present.

Therefore direct and swapped allele ordering were both valid. Palindromic exclusion, if used, occurred after this intermediate table and before/within the clumping stage.

## Still unresolved

The exact original PLINK clumping values for `--clump-r2` and `--clump-kb` are not preserved in the currently recovered material. They must be reconstructed before claiming bit-for-bit reproduction of the original clump counts.

The exact rule that generated the green “top 10 recurrence-ranked genes” in the final ±100 kb figure is also not yet recovered. The archived labels are included only as a validation fingerprint, not as proof of the ranking formula.
