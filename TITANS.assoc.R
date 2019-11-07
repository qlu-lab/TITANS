suppressMessages(library("reshape2"))
suppressMessages(library("survival"))
suppressMessages(library("data.table"))
suppressMessages(library("optparse"))
suppressMessages(library("survival"))
suppressMessages(library("dplyr"))

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
              help="The path to the output file.")
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

tissue = opt$tissue
chr = opt$chr
gene.name = opt$gene
output = opt$out

# --- 1. Performing strandflip
pred = read.table(opt$pred, header = T)

#read in files
source("TITANS.strand_flip.R")
weight = strand.flip(tissue, chr, gene.name, pred) # weight matrix after strandflipping
if(nrow(weight) == 0){
  cat("WARNING: no snp left in the imputation model after strand flip.\n")
  Result = data.frame(chr, nrow(pred), NA, gene.name, t(rep(NA, 4))) %>% 
    `colnames<-`(c("CHR", "Nsnps", "Nsnps.used", "Gene", "Beta", "SE", "Z", "P"))
  write.table(Result, file = output, row.names = F, col.names = T, append = F, quote = F, sep = "\t")
  q("no")
}

vcf = as.data.frame(fread(opt$vcf, header = T, sep = "\t"))

if(nrow(vcf) == 0){
  cat( "WARNING: No lines available in the vcf file\n")
  Result = data.frame(chr, nrow(pred), NA, gene.name, t(rep(NA, 4))) %>% 
    `colnames<-`(c("CHR", "Nsnps", "Nsnps.used", "Gene", "Beta", "SE", "Z", "P"))
  write.table(Result, file = output, row.names = F, col.names = T, append = F, quote = F, sep = "\t")
  q("no")
}
fam = read.table(opt$fam, header = F)
sample.id = Reduce(intersect, list(colnames(vcf), fam[,2]))

if(length(sample.id) == 0){
  cat("ERROR: samples in vcf file and fam file do not match\n")
  q("no")
}
pos = c(1:9, which(colnames(vcf) %in% fam[,2]))

# perform post-imputation QC
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
  cat( "WARNING: No snp left after post-imputation QC\n")
  Result = data.frame(chr, nrow(pred), NA, gene.name, t(rep(NA, 4))) %>% 
    `colnames<-`(c("CHR", "Nsnps", "Nsnps.used", "Gene", "Beta", "SE", "Z", "P"))
  write.table(Result, file = output, row.names = F, col.names = T, append = F, quote = F, sep = "\t")
  q("no")
}
# extract samples only in fam file
snp.intersect = intersect(vcf$POS, weight$POS)
if(length(snp.intersect) == 0){
  cat( "WARNING: No common snp in weight model and vcf file\n")
  Result = data.frame(chr, nrow(pred), NA, gene.name, t(rep(NA, 4))) %>% 
    `colnames<-`(c("CHR", "Nsnps", "Nsnps.used", "Gene", "Beta", "SE", "Z", "P"))
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
tmp = child
for(i in 10: ncol(child.g)){
  child.g[,i] = colsplit(child.g[, i], ":", c("GT","GQ","DP;HQ"))[1]
}
child = child.g
# FIXME
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
# FIXME
idx.child = match(colnames(child)[10:ncol(child)], fam[,2])
child.fam.id = fam[idx.child, 1]
child.sample.id = fam[idx.child, 2]
ge.c = data.frame(child.fam.id, child.sample.id, ge.c, rep(1, length(ge.c)))

# --- 3. Calculating parents' imputed gene expression
parent_info = fam[as.character(fam[,3]) == "0", 1:2]
parent = cbind(vcf[,1:9], vcf[, c(colnames(vcf) %in% parent_info[,2])])
nctrls = (ncol(parent)-9) * 2
nparent = ncol(parent)-9
ctrl = data.frame()
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

# Pseudo-sibling simulation
ge.gamete = t(w$weight[weight.idx]) %*% ctrl
ge.sib = data.frame()
for(i in seq(from = 1, to = (ncol(ge.gamete)-3), by = 4)){
  if(i != 1){
    ge.sib = cbind(ge.sib,
                   ge.gamete[i] + ge.gamete[i+2], 
                   ge.gamete[i+1] + ge.gamete[i+2], 
                   ge.gamete[i] + ge.gamete[i+3],
                   ge.gamete[i+1] + ge.gamete[i+3])
  } else{
    ge.sib = cbind(ge.gamete[i] + ge.gamete[i+2], 
                   ge.gamete[i+1] + ge.gamete[i+2], 
                   ge.gamete[i] + ge.gamete[i+3],
                   ge.gamete[i+1] + ge.gamete[i+3])
  }
}

# --- 4. Create family data
nchild = ncol(child.g)-9
ge.sib = c(ge.sib)
names(ge.sib) = rep(unique(parent_info[,1]), each = 4)
ge.sib = data.frame(rep(unique(parent_info[,1]), each = 4), rep(unique(parent_info[,1]), each = 4), ge.sib, rep(0, length(ge.sib)))
colnames(ge.sib) = c("fam", "sample", "expr","disease")
colnames(ge.c) = colnames(ge.sib)
GE = rbind(ge.c, ge.sib)
GE = GE[order(GE[,1]), ]

frq.table = as.data.frame(table(GE[GE[,4] == "1",1]))
quad_list = frq.table[frq.table[,2] == 2, 1]
if(length(quad_list) > 0){
  GE.quad = GE[GE[,1] %in% quad_list, ]
  GE.trio = GE[!(GE[,1] %in% quad_list), ]
  duplicate = GE.quad #(1,1,0,0,0,0) for every family
  duplicate.1 = duplicate[-seq(from = 1, to = 6 * length(quad_list), by = 6), ] #(1,0,0,0)
  duplicate.2 = GE.quad[-seq(from = 2, to = (6 * length(quad_list) - 1), by = 6),]
  duplicate.1[, 1] = paste0(duplicate.1[,1], "_quad")
  GE = rbind(GE.trio, duplicate.1, duplicate.2)
}
GE = GE[order(GE[,1]), ]
GE.sample.id = GE[seq(from = 1, to = nrow(GE), by = 5), ]$sample
GE$sample = rep(GE.sample.id, each = 5)
# --- 5. Delete the pseudo sibling which has the same genotype as the proband in each family

remove = list()
for(i in seq(from = 1, to = nrow(GE)-4, by = 5)){
  index = which(abs(GE[(i+1):(i+4),]$expr - GE[i,]$expr) == min(abs(GE[(i+1):(i+4),]$expr-GE[i,]$expr)))[1] + i
  remove = append(remove, index)
}
remove = unlist(remove) #index that the sibling having the same genotype as child's
GE.final = GE[-remove, ]

# --- 6. Perform the association test
pred = read.table(opt$pred, header = T)
result = tryCatch({
  clogit(disease~expr+strata(fam), data = GE.final)
}, error = function(e) {
  return(NA)
}
)

if(is.na(result)){
  Result = data.frame(chr, nrow(pred), nrow(weight), gene.name, t(rep(NA, 4))) %>% 
    `colnames<-`(c("CHR", "Nsnps", "Nsnps.used", "Gene", "Beta", "Exp.beta", "SE", "Z", "P"))
  write.table(Result, file = output, row.names = F, col.names = T, append = F, quote = F, sep = "\t")
} else{
  Result = data.frame(chr, nrow(pred), nrow(weight), gene.name, summary(result)$coefficients) %>% 
    `colnames<-`(c("CHR", "Nsnps", "Nsnps.used", "Gene", "Beta", "Exp.beta", "SE", "Z", "P")) %>% 
    select(-Exp.beta)
  write.table(Result, file = output, row.names = F, col.names = T, append = F, quote = F, sep = "\t")
  
}

