# DeGAs (Decomposition of Genetic Associations)

In this DeGAs directory, we deposit analysis scripts used in the following publication:

Y. Tanigawa*, J. Li*, et al., Components of genetic associations across 2,138 phenotypes in the UK Biobank highlight adipocyte biology. Nature Communications (2019). [https://doi.org/10.1038/s41467-019-11953-9](https://doi.org/10.1038/s41467-019-11953-9)

The article PDF file (via the ShardIt link): [https://rdcu.be/bQrf2](https://rdcu.be/bQrf2)


## Directory

`notebook`: This directory has notebooks that were used to generate the main and supplementary figures. 
`R_code`: This directory has R script used for truncated singular value decomposition (TSVD)
`src`: This directory has bash and Python scripts for data pre-processing as well as common functions used in the notebook. 

## Dataset availability

Most of the analysis notebooks take Python numpy (npz) data as the input. The datasets used in our publication is available at figshare:

Tanigawa, Yosuke; Rivas, Manuel (2019): Decomposed matrices used for the analysis described in 'Components of genetic associations across 2,138 phenotypes in the UK Biobank highlight adipocyte biology'. figshare. Dataset.
[https://doi.org/10.35092/yhjc.9202247.v1](https://doi.org/10.35092/yhjc.9202247.v1)

## Interactive web application
We also provide an interactive web application as a part of Global Biobank Engine: [https://gbe.stanford.edu/degas](https://gbe.stanford.edu/degas)

## Abstract

Population-based biobanks with genomic and dense phenotype data provide opportunities for generating effective therapeutic hypotheses and understanding the genomic role in disease predisposition. To characterize latent components of genetic associations, we apply truncated singular value decomposition (DeGAs) to matrices of summary statistics derived from genome-wide association analyses across 2,138 phenotypes measured in 337,199 White British individuals in the UK Biobank study. We systematically identify key components of genetic associations and the contributions of variants, genes, and phenotypes to each component. As an illustration of the utility of the approach to inform downstream experiments, we report putative loss of function variants, rs114285050 (GPR151) and rs150090666 (PDE3B), that substantially contribute to obesity-related traits and experimentally demonstrate the role of these genes in adipocyte biology. Our approach to dissect components of genetic associations across the human phenome will accelerate biomedical hypothesis generation by providing insights on previously unexplored latent structures.


## Reference

### Manuscript

Y. Tanigawa*, J. Li*, et al., Components of genetic associations across 2,138 phenotypes in the UK Biobank highlight adipocyte biology. Nature Communications (2019). [https://doi.org/10.1038/s41467-019-11953-9](https://doi.org/10.1038/s41467-019-11953-9)

The article PDF file (via the ShardIt link): [https://rdcu.be/bQrf2](https://rdcu.be/bQrf2)


### Dataset

Tanigawa, Yosuke; Rivas, Manuel (2019): Decomposed matrices used for the analysis described in 'Components of genetic associations across 2,138 phenotypes in the UK Biobank highlight adipocyte biology'. figshare. Dataset.
[https://doi.org/10.35092/yhjc.9202247.v1](https://doi.org/10.35092/yhjc.9202247.v1)


### Interactive web application

The DeGAs app in the Global Biobank Engine: [https://gbe.stanford.edu/degas](https://gbe.stanford.edu/degas)


### Video presentation of the publication and the tutorial of the web application

DeGAs App -- Global Biobank Engine: [https://www.youtube.com/watch?v=-ATVvLnFr9k](https://www.youtube.com/watch?v=-ATVvLnFr9k)
