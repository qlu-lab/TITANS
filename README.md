# ASD-TWAS

## Introduction

Our project gives...

### Prerequisites

The following R packages are required:

reshape2

survival

data.table

optparse

tidyverse

dplyr

```
Give examples
```

### Installing and step 0


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

Ask users to create folder results
```
mkdir results
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


