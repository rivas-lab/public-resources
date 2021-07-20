# Digical phenotyping in UK Biobank

In UK Biobank, multiple data fields are available for case definition of diseases. To capture the maximum number of cases in the genetic analysis, we combined information from multiple data fields to define the disease case status.

## Non-cancer disease traits

We mapped the disease names recorded in the including hospital in-patient records (ICD10 codes, [UK Biobank Data Coding 19](http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=19)) and self-reported disease status ([UK Biobank Data Coding 6](http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=6)). In our digital phenotyping manuscript[1], we showed their concordance using genetic correlation analysis.

Here, we provide the mapping file: [`UKB_coding_mapping_table.tsv`](UKB_coding_mapping_table.tsv). The table file has the following 6 columns:

- `GBE_ID`: our phenotype code used in [Globel Biobank Engine](http://gbe.stanford.edu/)
- `UKB_coding6`: disease codes (self-reported disease status) in [UK Biobank Data Coding 6](http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=6)
- `UKB_coding19`: disease codes (ICD10 codes) in [UK Biobank Data Coding 19](http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=19)
- `GBE_short_name`: the phenotype name
- `GBE_NAME`: the full phenotype name without abbreviation
- `URL`: the link to phenotype source breakdown page

Note, `UKB_coding6` and `UKB_coding19` may have multiple values, deliminated with `,` symbol. For example, myasthenia gravis has duplicated codes (1260 and 1437) in [UK Biobank Data Coding 6](http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=6).

When you find it useful to use this resource, please consider citing our digital phenotyping paper[1].

## Cancer

Similary, we prepared a manually curated mapping table between the self-reported cancer ([UK Biobank Data Coding 3](http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=3) available in [UK Biobank Data Field 21001](https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20001)) and ICD-10 code. Here, we have the mapping file ([`self_report_ICD10_mapping_treeRespect.tsv`](self_report_ICD10_mapping_treeRespect.tsv)) that contains the following 3 columns:

- meaning: types of cancer
- self-reported coding: disease codes (self-reported cancer) in [UK Biobank Data Coding 3](http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=3)
- ICD-10 Codes: disease codes (ICD10 codes) in [UK Biobank Data Coding 19](http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=19)

When you find it useful to use this resource, please consider citing our digital phenotyping paper[1].

## Reference

1. C. DeBoever, Y. Tanigawa, M. Aguirre, G. McInnes, A. Lavertu, M. A. Rivas, Assessing Digital Phenotyping to Enhance Genetic Studies of Human Diseases. The American Journal of Human Genetics. 106, 611-622 (2020). [doi:10.1016/j.ajhg.2020.03.007](https://doi.org/10.1016/j.ajhg.2020.03.007)
