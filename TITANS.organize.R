suppressMessages(library("reshape2"))
suppressMessages(library("survival"))
suppressMessages(library("data.table"))
suppressMessages(library("optparse"))
suppressMessages(library("survival"))
suppressMessages(library("dplyr"))

options(stringsAsFactors = F)
option_list = list(
  make_option("--result", action = "store", default = NA, type = "character",
              help="The path to the result file."),
  make_option("--infoTable", action = "store", default = NA, type = "character",
              help="The path to the gene info file."),
  make_option("--out", action = "store", default = NA, type = "character",
              help="The path to the output file.")
)

opt = parse_args(OptionParser(option_list=option_list))

if(is.na(opt$result) | is.na(opt$infoTable) | is.na(opt$out)){
  cat("ERROR: Missing essential inputs\n")
  q("no")
}

result = as.data.frame(fread(opt$result, header = T))
#FIXME
result = as.data.frame(fread("/Users/kunlinghuang/Box Sync/lu_lab/software_development/git_repo/TITANS/Results/FAKE_TISSUE.FAKE_GENE.txt", header = T))
result$Gene = "ABCB10"
result$CHR = 9
opt$infoTable = "/Users/kunlinghuang/Box Sync/lu_lab/software_development/git_repo/TITANS/Tables/gene.coordinate.txt"
info.table = as.data.frame(fread(opt$infoTable))

result = left_join(result, info.table, by = c("Gene" = "Splicing.ID", "CHR")) %>% 
  select(CHR, Gene.y, Gene, P0, Nsnps, Nsnps.used, Beta, SE, Z, P) %>% 
  `colnames<-`(c("CHR", "Gene", "Splicing.ID", "P0", "Nsnps", "Nsnps.used", "Beta", "SE", "Z", "P")) %>% 
  mutate(Splicing.ID = ifelse(!is.na(Gene) & Splicing.ID == Gene, "", Splicing.ID)) %>% 
  mutate(tmp = Gene) %>% 
  mutate(Gene = ifelse(is.na(Gene) & !grepl("clu", Splicing.ID), Splicing.ID, Gene)) %>% 
  mutate(Splicing.ID = ifelse(is.na(tmp) & !grepl("clu", Splicing.ID), tmp, Splicing.ID)) %>% 
  select(-tmp)

write.table(result, file = opt$out, col.names = T, row.names = F, quote = F, sep = "\t")


