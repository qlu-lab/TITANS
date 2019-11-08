# TITANS

### Introduction

TITANS (TrIo-based Transcriptome-wide AssociatioN Study) is a statistical framework to conduct TWAS in proband-parent trios. The preprint is available at ???.


### Prerequisites

The software is developed and tested in Linux and Mac OS environments. The statistical computing software [R](https://www.r-project.org/) (>=3.5.1) and the following R packages are required for association tests:

* [data.table](https://cran.r-project.org/web/packages/data.table/index.html) (>=1.11.8)
* [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html) (>=0.8.3)
* [optparse](https://cran.r-project.org/web/packages/optparse/index.html) (>=1.6.0)
* [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html) (>=1.4.3)
* [survival](https://cran.r-project.org/web/packages/survival/index.html) (>=2.44-1.1)
* [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html) (>=1.2.1)

### Data preparation

`TITANS` is a statistical tool for analysing phased trio genotyped data. We recommand performing phasing and imputation using [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!). You can use [bcftools](http://samtools.github.io/bcftools/bcftools.html) and the following codes to perform post-imputation quality control

```bash
$ bcftools filter -i 'INFO/MAF>0.01 && INFO/R2>0.8' $INPUT -o $OUTPUT
```
where `$INPUT` is the [vcf.gz](https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/) file processed by Michigan Imputation Server. The `bcftools` will output a QCed vcf file `$OUTPUT`. 

**Note:** For vcf.gz or vcf files not from Michigan Imputation Server, they need fields `INFO/MAF` and `INFO/R2` to perform quality control.

### Getting started

#### 1. Clone the git repository

```bash
$ git clone https://github.com/qlu-lab/TITANS.git
cd TITANS
```

#### 2. Download and upzip the precomputed imputation model (99Mb before unzipping)

```bash
$ wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/TITANS/Software/GEmodel.zip
$ unzip GEmodel.zip
```

After upzipping, there will be directories containing gene expression imputation models from 12 brain tissues from [Genotype-Tissue Expression (GTEx)](https://www.gtexportal.org/home/) and [CommonMind consortium (CMC)](https://www.nimhgenetics.org/resources/commonmind).

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

For detailed instructions of gene expression imputation models, please visit the [wiki](https://github.com/qlu-lab/TITANS/wiki/1.-Installation).

#### 3. Tutorial

We provide a walkthrough using fake data in tutorial. First, create a dirctory for the output

```bash
$ mkdir Results ## under TITANS/
```

and use the toy data under `Example` to conduct a trio-based TWAS on a `FAKE_GENE` on chromosome 1 in `FAKE_TISSUE`, the gene expression imputation model for this fake gene is under `Example/weight.FAKE_TISSUE.FAKE_GENE.txt`, and the QCed vcf file and fam file are under `./Example/FAKEDATA.vcf` and `./Example/FAKEDATA.fam` respectively.

```bash
$ Rscript TITANS.assoc.R \
  --tissue FAKE_TISSUE \
  --chr 1 \
  --gene FAKE_GENE \
  --fam ./Example/FAKEDATA.fam \
  --pred ./Example/weight.FAKE_TISSUE.FAKE_GENE.txt \
  --vcf ./Example/FAKEDATA.vcf \
  --out ./Results/FAKE_TISSUE.FAKE_GENE.txt
```

`TITANS` will write out a result table `./Results/FAKE_TISSUE.FAKE_GENE.txt`

| Column | Value | Description |
|-----|-------|------|
| CHR |  1  | The chromosomal location of the gene |
| Nsnps | 10 | Number of SNPs in the gene expression imputation model |
| Nsnps.used | 10 | Number of SNPs used in building the association test |                                                    
| Gene | FAKE_GENE | The name of the gene. |
| Beta | -7.412979 | The estimated effect size. |
| SE | 7.313888 | The estimated standard error of Beta. |
| Z | -1.013548 | The Z test statistic for testing transmission disequilibrium. |
| P | 0.310798 | The P-value for testing transmission disequilibrium. |

**Interpretation:** This trio-based TWAS was based on `FAKE_GENE` on chromosome 1, with `10` SNPs in the gene expression imputation model. After post-imputation QC and checking for variant consistency between gene expression model and genotyping data, `10` SNPs remain in our analysis. The disease effect size of `FAKE_GENE` estimated by conditional logistic regression is `7.31388`, the direction of `Beta` indicates that the higher expression raises the disease risk. The P-value is `0.310798`, indicating that the gene is not significantly associated with the disease. In other words, `FAKE_GENE` is not considered associating with the disease in `FAKE_TISSUE.`

#### 4. Use TITANS to perform trio-based TWAS

The `TITANS.assoc.R` outputs the trio-based transcriptome-wide association test result **one gene in one tissue at a time**. `TITANS.assoc.R` takes the following inputs

```bash
$ Rscripts TITANS.assoc.R \
    --tissue $TISSUE \
    --chr $CHR \
    --gene $GENE \
    --fam $FAMPATH \
    --pred $PRED \
    --vcf $VCF \
    --out $OUTPUT \
```

where the inputs are

| Flag | Description |
|-----|------------------------------------------------------------------------|
| --tissue      | One of the above brain tissues. |
| --chr         | The chromosomal location of the gene. |
| --gene        | The name of the Gene. |                                                    
| --fam     | The path to the fam file; the fam file must follow the [PLINK format](https://www.cog-genomics.org/plink/1.9/formats#fam).<br>The sample IDs in the fam file should also be consistent with the sample IDs in the vcf file. |
| --pred        | The path to the prediciton matrix. |
| --vcf         | The path to the phased and QCed vcf file. |
| --out       | The path to the output file. |

The final result has the following fields:

| Column | Description |
|-----|-------------|
| CHR | The chromosomal location of the gene |
| Nsnps | Number of SNPs in the gene expression imputation model |
| Nsnps.used | Number of SNPs used in building the association test |                                                    
| Gene | The name of the gene. |
| Beta | The estimated effect size. |
| SE | The estimated standard error of Beta. |
| Z | The Z test statistic for testing transmission disequilibrium. |
| P | The P-value for testing transmission disequilibrium. |

#### 5. Organize the output (optional)

`Tables/gene.coordinate.txt` containing gene name, splicing ID (genes in CMC splicing will have its unique slicing ID, otherwise splicing ID is set as gene name), chromosomal location, and gene coordinate information based on hg19. Use this table and the following code to organize your TITANS results.

```bash
$ Rscripts TITANS.organize.R \
    --input $INPUT \
    --result $RESULT \
    --out $OUTPUT \
```

**Note:** If users provide wrong chromosomal and gene information, the result table will still be generated, but the gene coordinate will not be available.

## Authors

Qiongshi Lu (University of Wisconsin-Madison, Department of Biostatistics and Medical Informatics)

Kunling Huang (University of Wisconsin-Madison, Department of Statistics)

Yuchang Wu (University of Wisconsin-Madison, Department of Biostatistics and Medical Informatics)

Yupei Lin (University of Wisconsin-Madison, Department of Computer Sciences)

## License

All rights reserved for [Lu-Laboratory](https://qlu-lab.org/)

## Acknowledgments
The imputation models and part of the codes are modified from the precomputed gene expression imputation model [UTMOST](https://github.com/Joker-Jerome/UTMOST) and [FUSION](http://gusevlab.org/projects/fusion/), we thank the authors for sharing their gene expression imputation model and codes. 
