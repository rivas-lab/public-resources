# Digical phenotyping in UK Biobank

In UK Biobank, multiple data fields are available for case definition of diseases. We mapped the disease names recorded in the including hospital in-patient records (ICD10 codes, [UK Biobank Data Coding 19](http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=19)) and self-reported disease status ([UK Biobank Data Coding 6](http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=6)). In our digital phenotyping manuscript[1], we showed their concordance using genetic correlation analysis.

Here, we provide the mapping file: [`UKB_coding_mapping_table.tsv`](UKB_coding_mapping_table.tsv). The table file has the following 6 columns:

- `GBE_ID`: our phenotype code used in [Globel Biobank Engine](http://gbe.stanford.edu/)
- `UKB_coding6`: disease codes in [UK Biobank Data Coding 6](http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=6)
- `UKB_coding19`: disease codes in [UK Biobank Data Coding 19](http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=19)
- `GBE_short_name`: the phenotype name
- `GBE_NAME`: the full phenotype name without abbreviation
- `URL`: the link to phenotype source breakdown page

Note, `UKB_coding6` and `UKB_coding19` may have multiple values, deliminated with `,` symbol.

When you find it useful, please consider citing our digital phenotyping manuscript[1].

## Reference

1. C. DeBoever, Y. Tanigawa, M. Aguirre, G. McInnes, A. Lavertu, M. A. Rivas, Assessing Digital Phenotyping to Enhance Genetic Studies of Human Diseases. The American Journal of Human Genetics. 106, 611-622 (2020). [doi:10.1016/j.ajhg.2020.03.007](https://doi.org/10.1016/j.ajhg.2020.03.007)
