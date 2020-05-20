# In silico analysis of alternative splicing on drug-target gene interactions

## Background
This repository contains scripts used to reproduce study published in [Ji, Mishra, Davuluri, Scientific Reports, 2020](https://www.nature.com/articles/s41598-019-56894-x). The aim is to perform a global survey and computational analysis of protein isoforms/splice variants that may affect drug-gene interactions.

## Prerequisites
* R (≥ 3.4.3)
* Bioconductor
* Data from DGIDB, Ensembl, BioLIP, TCGA, GTEx

## Databases
* The Drug Gene Interaction Database (DGIDB) (http://dgidb.org/) - Interactions TSV, Genes TSV, Drugs TSV
* Ensembl GRCh38.p12/GENCODE v30 BioMart (http://useast.ensembl.org/index.html) - Gene, transcript and protein annotations
* BioLIP protein function database (https://zhanglab.ccmb.med.umich.edu/BioLiP/) - Non-redundant set
* DrugBank, CHEMBL and PDB - for ID mapping and conversion
* TCGA, GTEx expression data from UCSC Xena Browser (https://xena.ucsc.edu/)

## File desciption
* isoform_drug_targets.R - major R notebook containing the analysis workflow (summary statistics, multiple alignments)
* expression_tcga_gtex.R - RNA-Seq expression analysis of splice variants as shown in Fig. 7 and 8
* helpers.R - script containing some auxiliary functions to be used in the analysis.

## License
[The MIT License](LICENSE.md)
