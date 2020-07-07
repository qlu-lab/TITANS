suppressMessages(library(reshape2))
suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(survival))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))


options(stringsAsFactors = F)
option_list = list(
  make_option("--tissue", action = "store", default = NA, type = "character", help = "The brain tissues"),
  make_option("--chr", action = "store", default = NA, type = "integer", help = "The chromosomal location of the gene"),
  make_option("--gene", action = "store", default = NA, type = "character", help = "The name of the Gene"),
  make_option("--fam", action = "store", default = NA, type = "character", help = "The path to the fam file"),
  make_option("--pred", action = "store", default = NA, type = "character", help = "The path to the prediciton matrix."),
  make_option("--vcf", action = "store", default = NA, type = "character",
              help="The path to the phased and QCed vcf file."),
  make_option("--out", action = "store", default = NA, type = "character",
              help="The path to the output file."),
  make_option("--qc", action = "store", default = T,
              help="Removing SNPs with imputation quality less than 0.8 or MAF less than 0.01"),
  make_option("--matching", action = "store", default = "3sib", type = "character",
              help="The pseudo sibling matching method, choose among 3sib, 1sib, and 2parent")
  
)

opt = parse_args(OptionParser(option_list=option_list))

# --- 0. I/O check
if(is.na(opt$tissue) | is.na(opt$chr) |  is.na(opt$gene) | is.na(opt$fam) | is.na(opt$pred) | is.na(opt$vcf) | is.na(opt$out)){
  cat("ERROR: Missing essential inputs\n")
  q("no")
}

if(!(dir.exists(dirname(opt$out)))){
  cat("ERROR: Output directory does not exist\n")
  q("no")
}

if(!(opt$matching == "3sib" | opt$matching == "1sib" | opt$matching == "2parent")){
  cat("ERROR: Matching method does not exist\n")
  q("no")
}

tissue = opt$tissue
chr = opt$chr
gene.name = opt$gene
output = opt$out
# --- 1. Performing strandflip

source("TITANS.strand_flip.R")

pred = read.table(opt$pred, header = T)
weight = strand.flip(tissue, chr, gene.name, pred) # weight matrix after strandflipping
if(nrow(weight) == 0){
  cat("WARNING: no snp left in the imputation model after strand flip.\n")
  Result = data.frame(chr, nrow(pred), nrow(weight), gene.name, opt$matching, t(rep(NA, 4))) %>% 
    `colnames<-`(c("CHR", "Nsnps", "Nsnps.used", "Gene", "Matching", "Beta", "Exp.beta", "SE", "Z", "P"))
  write.table(Result, file = output, row.names = F, col.names = T, append = F, quote = F, sep = "\t")
  q("no")
}

vcf = as.data.frame(fread(opt$vcf, header = T, sep = "\t"))
if(nrow(vcf) == 0){
  cat( "WARNING: No lines available in the vcf file\n")
  Result = data.frame(chr, nrow(pred), nrow(weight), gene.name, opt$matching, t(rep(NA, 4))) %>% 
    `colnames<-`(c("CHR", "Nsnps", "Nsnps.used", "Gene", "Matching", "Beta", "Exp.beta", "SE", "Z", "P"))
  write.table(Result, file = output, row.names = F, col.names = T, append = F, quote = F, sep = "\t")
  q("no")
}

fam = read.table(opt$fam, header = F) %>% 
  `colnames<-`(c("fam", "sample.id", "p1", "p2", "sex", "disease"))
sample.id = Reduce(intersect, list(colnames(vcf), fam[,2]))

if(length(sample.id) == 0){
  cat("ERROR: samples in vcf file and fam file do not match\n")
  q("no")
}
pos = c(1:9, which(colnames(vcf) %in% fam[,2]))

# perform post-imputation QC
if(opt$qc){
  # subtract snp info
  info = colsplit(vcf$INFO, ";", c("AF", "MAF", "R2", "ER2"))
  info$MAF = colsplit(info$MAF, "=", c("MAF", "value"))
  info$R2 = colsplit(info$R2, "=", c("R2", "value"))
  info.data = data.frame(info$MAF[,2], info$R2[,2])
  index = which((info.data$MAF < 0.01) | (info.data$R2 < 0.8))
  if(length(index) > 0){
    print("bad snps have been found")
    vcf = vcf[-index, ]
  }
  if(nrow(vcf) == 0){
    cat( "WARNING: No snp left after QC\n")
    Result = data.frame(chr, nrow(pred), nrow(weight), gene.name, opt$matching, t(rep(NA, 4))) %>% 
      `colnames<-`(c("CHR", "Nsnps", "Nsnps.used", "Gene", "Matching", "Beta", "Exp.beta", "SE", "Z", "P"))
    write.table(Result, file = output, row.names = F, col.names = T, append = F, quote = F, sep = "\t")
    q("no")
  }
}

# extract samples only in fam file
snp.intersect = intersect(vcf$POS, weight$POS)
if(length(snp.intersect) == 0){
  cat( "WARNING: No common snp in prediction model and vcf file\n")
  Result = data.frame(chr, nrow(pred), nrow(weight), gene.name, opt$matching, t(rep(NA, 4))) %>% 
    `colnames<-`(c("CHR", "Nsnps", "Nsnps.used", "Gene", "Matching", "Beta", "Exp.beta", "SE", "Z", "P"))
  write.table(Result, file = output, row.names = F, col.names = T, append = F, quote = F, sep = "\t")
  q("no")
}
w = weight
w = w %>% 
  filter(POS %in% snp.intersect)
vcf = vcf %>% 
  filter(POS %in% snp.intersect)


# --- 2. Calculating childs' imputed gene expression
child.ID = fam[as.character(fam[,3]) != "0",1:2]
child = cbind(vcf[,1:9], vcf[, c(colnames(vcf) %in% child.ID[,2])])
pos.child = which(colnames(vcf) %in% child.ID[,2])
child = vcf[, c(1:9, pos.child)]
child.g = child
weight.idx = match(child$POS, w$POS)
for(i in 10: ncol(child.g)){
  child.g[,i] = colsplit(child.g[, i], ":", c("GT","GQ","DP;HQ"))[1]
}
child = child.g
child[, 10: ncol(child.g)] = matrix(0, nrow = nrow(child.g), ncol = ncol(child.g)-9)
for(i in 10: ncol(child)){
  if(length(which(child.g[,i] == "0|0")) != 0){
    child[which(child.g[,i] == "0|0"), i] = rep(0, length(which(child.g[,i] == "0|0")))
  }
  
  if(length(which(child.g[,i] == "0|1")) != 0){
    child[which(child.g[,i] == "0|1"), i] = rep(1, length(which(child.g[,i] == "0|1")))
  }
  
  if(length(which(child.g[,i] == "1|0")) != 0){
    child[which(child.g[,i] == "1|0"), i] = rep(1, length(which(child.g[,i] == "1|0")))
  }
  
  if(length(which(child.g[,i] == "1|1")) != 0){
    child[which(child.g[,i] == "1|1"), i] = rep(2, length(which(child.g[,i] == "1|1")))
  }
}

ge.c = t(as.matrix(w$weight[weight.idx])) %*% as.matrix(child[,10:ncol(child)])
ge.c = c(ge.c)

idx.child = match(colnames(child)[10:ncol(child)], fam[,2])
child.fam.id = fam[idx.child, ] %>% 
  separate(sample.id, sep = "_", into = c("fam", "sample")) %>% 
  select(fam, sample)
ge.c = data.frame(child.fam.id, ge.c, rep(1, length(ge.c))) %>% 
  `colnames<-`(c("fam","sample","expr", "disease"))

# --- 3. Calculating pseudo siblings' imputed gene expressions
source("TITANS.matching.R")
parent_info = fam[as.character(fam[,3]) == "0", 1:2]
parent = cbind(vcf[,1:9], vcf[, c(colnames(vcf) %in% parent_info[,2])])
nctrls = (ncol(parent)-9) * 2
nparent = ncol(parent)-9
ctrl = data.frame()

# Substract parental strand information
for(i in 1:length(unique(parent_info[,1]))) {
  fam.ID = unique(parent_info[,1])[i]
  fam.par.ID = parent_info[as.character(parent_info[,1]) %in% fam.ID, ]
  fam.par = parent[, c(colnames(parent) %in% as.character(fam.par.ID[,2]))]
  #subtract genotype info
  tmp1 = colsplit(fam.par[, 1], '\\|', c("p1", "p2"))
  tmp1$p2 = colsplit(tmp1$p2, ':', c("p2", "info1", "info2"))
  tmp1 = cbind(tmp1$p1, tmp1$p2[,1])
  tmp2 = colsplit(fam.par[, 2], '\\|', c("m1", "m2"))
  tmp2$m2 = colsplit(tmp2$m2, ':', c("m2", "info1", "info2"))
  tmp2 = cbind(tmp2$m1, tmp2$m2[,1])
  if(i != 1){
    ctrl = cbind(ctrl, tmp1, tmp2)
  } else{
    ctrl = cbind(tmp1, tmp2)
  }
}
# Simulate the pseudo siblings and return the family data
GE.final = matching(opt$matching, ctrl)

# --- 4. Perform the association test
pred = read.table(opt$pred, header = T)
result = tryCatch({
  clogit(disease~expr+strata(fam), data = GE.final)
}, error = function(e) {
  return(NA)
}
)

if(length(result) == 1){
  Result = data.frame(chr, nrow(pred), nrow(weight), gene.name, opt$matching, t(rep(NA, 4))) %>% 
    `colnames<-`(c("CHR", "Nsnps", "Nsnps.used", "Gene", "Matching", "Beta", "Exp.beta", "SE", "Z", "P"))
  write.table(Result, file = output, row.names = F, col.names = T, append = F, quote = F, sep = "\t")
} else{
  Result = data.frame(chr, nrow(pred), nrow(weight), gene.name, opt$matching, summary(result)$coefficients) %>% 
    `colnames<-`(c("CHR", "Nsnps", "Nsnps.used", "Gene", "Matching", "Beta", "Exp.beta", "SE", "Z", "P")) %>% 
    select(-Exp.beta)
  write.table(Result, file = output, row.names = F, col.names = T, append = F, quote = F, sep = "\t")
  
}

