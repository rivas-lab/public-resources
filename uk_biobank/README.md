# UK Biobank

This directory contains resources for the UK Biobank.

## `bims_combined.vep.tsv`

This file contains some annotations from VEP for the autosomal variants in the
plink `bim` files generated from the 500k UK Biobank release. 

## `variant_filter_table.tsv`
This file has statistics for variants on the arrays plus information for
variant filtering.

### Statistics and annotation
The missingness rate, allele frequency, and Hardy-Weinberg p-values are for the
500k release (though these values were calculated on the subset of individuals
we are using for GWAS). I took these values from the plink output files. 

The column `ld_indep` indicates whether a variant is part of the LD-pruned set
of variants. LD pruning was performed using plink with "--indep 50 5 2".

There are columns indicating whether a variant is specific to either array. 

### Filter columns
There are also columns with filtering information. Some filters have been
applied to all variants while others have only been applied to a subset (often
the subset that we have annotated as PTVs). For some columns, 1 in a column
means that the variant did not pass that filter and should be removed. The
final column contains the number of filters a variant did not pass. So if a
variant has a zero in the final column, it passed all filters. The current
filters are

- Missingness < 1% (calculated on an array-specific basis for array-specific
  variants)
- HWE p > 1e-7
- MCPI = Manual Cluster Plot Inspection. These are variants that we have looked
  at the cluster plots for and see problems. We haven't looked at all variants.

The "gaf" column indicates variants whose allele frequencies differ
substantially from gnomAD NFE. This was performed for a subset of PTVS. We
matched PTVs to PTVs annotated in gnomAD (gnomad.exomes.r2.0.1.sites.vcf.gz)
based on genomic position, reference, and alternate alleles and compared the
allele frequencies in the UKB and gnomAD by (1) performing a Fisher’s exact
test using the minor allele counts from the 337,208 UKB subjects and the minor
allele counts from gnomAD and (2) calculating the log odds ratio of observing
the minor allele in the UKB versus gnomAD. We stratified PTVs by minor allele
frequency into the following three bins: (0.01%, 0.1%], (0.1%, 1%], (1%, 50%].
For bin (0.01%, 0.1%], we removed PTVs with Fisher p < 1e-7 and an absolute log
odds ratio greater than 3. For bin (0.1%, 1%], we removed PTVs with Fisher p <
1e-7 and an absolute log odds ratio greater than 2. For bin (1%, 100%], we
removed PTVs with Fisher p < 1e-7 and an absolute log odds ratio greater than
1. The values in this column are

- NA: The variant was not considered for this filter.
- PASS: The variant passed this filter.
- FAIL: The variant failed this filter.

The "mgi" column stands for manual gnomAD inspection. For this column, we 
manually reviewed a subset of PTVs that we could not match to the gnomAD exome
data (gnomad.exomes.r2.0.1.sites.vcf.gz) on the gnomAD browser to determine
whether they were likely to accurately type a PTV in gnomAD. This column has a
few different values:

- NA: We didn't manually inspect this variant.
- PASS: The PTV was present on the gnomAD browser but was not included in the
  exome data and the allele frequencies between UKB and gnomAD matched
reasonably.
- NOT_PTV: The UKB array likely typed a non-PTV.
- FAIL: We marked variants as failed if the variant was not present on the
  browser, the variant was multiallelic, the allele frequencies differ
drastically between UKB and gnomAD, or other signs that the array genotypes may
not be accurate relative to the annotation.

## UK Biobank summary statistics

This directory contains summary statistics for UK Biobank data using only directly genotyped variants.

## `column_descripitions`

This file contains column descriptions. 

## `variant_filter_table.tsv`
This file has statistics for variants on the arrays plus information for
variant filtering.
