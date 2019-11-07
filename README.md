# TITANS

## Introduction

TITANS (TrIo-based Transcriptome-wide AssociatioN Study) is a statistical framework to conduct TWAS in proband-parent trios. The preprint is available at ???.


### Prerequisites

The software is developed and tested in Linux environment. Since TITANS is a statistical analysis based on phased trio genotyped data, we recommand performing phasing and imputation using [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!). You can use [bcftools](http://samtools.github.io/bcftools/bcftools.html) and the following codes to perform post-imputation quality control

```bash
$ bcftools filter -i 'INFO/MAF>0.01 && INFO/R2>0.8' $INPUT -o $OUTPUT
```
where `$INPUT` is the [vcf.gz](https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/) file processed by Michigan Imputation Server. The `bcftools` will output a QCed vcf file. 

**Note:** For vcf.gz or vcf files not from Michigan Imputation Server, they need fields `INFO/MAF` and `INFO/R2` to perform quality control.

The following R packages are required for association tests:

* data.table

* dplyr

* optparse

* reshape2

* survival

* tidyverse


### Getting started

#### 1. Clone the git repository

```bash
$ git clone https://github.com/qlu-lab/TITANS.git
cd TITANS
```

#### 2. Download and upzip the precomputed prediction model (99Mb before unzipping)

```bash
$ wget -O GEmodel.zip -L https://uwmadison.box.com/shared/static/qriv1whlpoxzr0dkbeqmqvaot8yw5g6f.zip
$ unzip GEmodel.zip
```

After upzipping, there will be directories containing gene expression (GE) imputation models from 12 brain tissues

| Tissue | Assay | # Genes | Study |
|------|-----|------------|----------------------------------------|
| Brain_Anterior_cingulate_cortex_BA24 | RNA-seq | 10807 | GTEx |
| Brain_Caudate_basal_ganglia | RNA-seq | 11653 | GTEx |
| Brain_Cerebellar_Hemisphere | RNA-seq | 10886 | GTEx |
| Brain_Cerebellum | RNA-seq | 11231 | GTEx |
| Brain_Cortex | RNA-seq | 11469 | GTEx |
| Brain_Frontal_Cortex_BA9 | RNA-seq | 11546 | GTEx |
| Brain_Hippocampus | RNA-seq | 11486 | GTEx |
| Brain_Hypothalamus | RNA-seq | 11783 | GTEx |
| Brain_Nucleus_accumbens_basal_ganglia | RNA-seq | 11564 | GTEx |
| Brain_Putamen_basal_ganglia.txt | RNA-seq | 11153 | GTEx |
| CMC_Brain_DLPFC | RNA-seq | 5419 | CMC |
| CMC_Brain_DLPFC_splicing | RNA-seq splicing | 7782 | CMC |

#### 3. Use TITANS to perform trio-based TWAS

The `TITANS.assoc.R` outputs the trio-based transcriptome-wide association test result **one gene in one tissue at a time**. `TITANS.assoc.R` takes the following inputs

```bash
Rscripts TITANS.assoc.R
  --tissue $TISSUE
  --chr $CHR
  --gene $GENE
  --fam $FAMPATH
  --pred $PRED
  --vcf $VCF
```
where the inputs are

| Flag | Description |
|-----|-------------|
| --tissue      | One of the above brain tissues. |
| --chr         | The chromosomal location of the gene. |
| --gene        | The name of the Gene. |                                                    
| --fam     | The path to the fam file; the fam file must follow the [PLINK format](https://www.cog-genomics.org/plink/1.9/formats#fam). The sample IDs in the fam file should also be consistent with the sample IDs in the vcf file. |
| --pred        | The path to the prediciton matrix. |
| --vcf         | The path to the phased and QCed vcf file. |

The final result has the following fields:

| Column | Description |
|-----|-------------|
| CHR | The chromosomal location of the gene |
| Nsnps | Number of SNPs in the GE prediction model |
| Nsnps.used | Number of SNPs used in building the association test |                                                    
| Gene | The name of the gene. |
| Beta | The estimated effect size. |
| SE | The estimated standard error of Beta. |
| Z | The Z test statistic for testing transmission disequilibrium. |
| P | The P-value for testing transmission disequilibrium. |

## Authors

See also the list of [contributors](##) who participated in this project.

## License

All rights reserved for Lu-Laboratory

## Acknowledgments
[bcftools](http://samtools.github.io/bcftools/bcftools.html)

[FUSION](http://gusevlab.org/projects/fusion/)

[UTMOST](https://github.com/Joker-Jerome/UTMOST)

[GTEx](https://www.gtexportal.org/home/)

[CMC](https://www.nimhgenetics.org/resources/commonmind)


