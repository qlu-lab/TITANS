# TITANS

## Introduction

TITANS (TrIo-based Transcriptome-wide AssociatioN Study) is a statistical framework to conduct TWAS in proband-parent trios. The preprint is available at ???.


### Prerequisites

TITANS is a statistical analysis based on phased trio genotyped data. We recommand performing phasing and imputation using Michigan Imputation Server. You can use [bcftools](http://samtools.github.io/bcftools/bcftools.html) and the following codes to perform post-imputation quality control

```
$ bcftools filter -i 'INFO/MAF>0.01 && INFO/R2>0.8' $vcfPath -Oz -o $output
```

The following R paskages are required for association tests:

* data.table

* dplyr

* optparse

* reshape2

* survival

* tidyverse


### Getting started

1. Clone the git repository

```
$ git clone https://github.com/qlu-lab/TITANS.git
cd TITANS
```

Also the precomputed prediction model is available at 



Downloading the following libraries and create a directory called trio-twas and change to this destination directory.

```
wget ..
mkdir trio-twas
cd trio-twas
```
Download xx.gz and decompress xx.gz under folder trio-twas
```
unzip xx.gz

```

Change to Scripts folder 
```
cd scripts
```


End with an example of getting some data out of the system or using it for a little demo

## Input Data

| Parameter                   | Directory Name | Description                                                                  |
|----------------------------|----------------|------------------------------------------------------------------------------|
| Tissue            | --tissue      | The Tissue |
| Chromosone         | --chr          | The chromosome of the gene        |
| Gene               | --gene        | Gene name                        |                                                    
| Family path           |--famPath | The path of the fam file, the ids in the fam file should fit the ids from the vcf file  |
| Prediction          |--pred | The path of the original prediciton matrix  |
| VCF           |--vcf | The path of the vcf file, retrieved by step 0|

Explain how to run the automated tests for this system

### Step 1

Conduct Quality Check by using [bcftools](http://samtools.github.io/bcftools/bcftools.html) 


```
bcftools filter -i 'INFO/MAF>0.01 && INFO/R2>0.8' $vcfPath -Oz -o $output
```

### Step 2
Run clogit association test

```
R 
3-assoc.clogit.R
```

## Step 3

Add additional notes about how to deploy this on a live system



## Authors

See also the list of [contributors](##) who participated in this project.

## License

All rights reserved for Lu-Laboratory

## Acknowledgments
[bcftools](http://samtools.github.io/bcftools/bcftools.html) \
[plink2R](https://github.com/gusevlab/fusion_twas/issues/13)

