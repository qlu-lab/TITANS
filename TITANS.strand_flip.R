strand.flip = function(tissue, chr, gene.name, pred){
  pred = read.table(opt$pred, header = T) %>% 
    unique %>% 
    `colnames<-`(c("rsid", "gene", "weight", "REF", "EFF", "CHR", "POS"))
  snp.dup.pred = pred[duplicated(pred$POS), ]$POS
  if(length(snp.dup.pred) > 0){
    pred = pred %>% 
      filter(!POS %in% snp.dup.pred)
  }
  vcf = as.data.frame(fread(opt$vcf, header = T, select = c(1:5), sep = "\t")) %>% 
    unique %>% 
    `colnames<-`(c("CHR", "POS", "ID", "REF", "ALT"))
  if(nrow(vcf) == 0){
    cat( "WARNING: No lines available in the vcf file, but still prints out weight file\n")
    #FIXME
    tmp.weight = pred %>% 
      filter(!(1:nrow(pred)))
    return(tmp.weight)
  }
  # remove duplicated snps
  snp.dup.vcf = vcf[duplicated(vcf$POS), ]$POS
  if(length(snp.dup.vcf) > 0){
    vcf = vcf %>% 
      filter(!POS %in% snp.dup.vcf)
  }
  
  snp.intersect = intersect(vcf$POS, pred$POS)
  vcf = vcf %>% 
    filter(POS %in% snp.intersect)
  pred = pred %>% 
    filter(POS %in% snp.intersect)
  # step 0: trim out vague strands
  vcf.cp = vcf
  vcf.mapped = vcf.cp %>% 
    left_join(pred, by = c("CHR", "POS")) %>% 
    `colnames<-`(c("CHR", "POS", "ID", "REF.vcf", "ALT.vcf", "rsid", "gene", "weight", "REF.pred", "ALT.pred")) %>% 
    mutate_at(vars(starts_with("REF")), as.character) %>% 
    mutate_at(vars(starts_with("ALT")), as.character) %>% 
    filter(REF.vcf == "A" | REF.vcf == "C" | REF.vcf == "G" | REF.vcf == "T") %>% 
    filter(ALT.vcf == "A" | ALT.vcf == "C" | ALT.vcf == "G" | ALT.vcf == "T") %>% 
    filter(!((REF.vcf == "A" & ALT.vcf == "T") | (REF.vcf == "T" & ALT.vcf == "A") | (REF.vcf == "C" & ALT.vcf == "G") | (REF.vcf == "G" & ALT.vcf == "C")))
  
  # step 1: delete same REF but different ALT or same ALT but dirrerent REF
  vcf.mapped = vcf.mapped %>% 
    filter(!(REF.vcf == REF.pred & ALT.vcf != ALT.pred)) %>% 
    filter(!(REF.vcf != REF.pred & ALT.vcf == ALT.pred)) %>% 
    filter(!(REF.vcf == ALT.pred & ALT.vcf != REF.pred)) %>% 
    filter(!(REF.vcf != ALT.pred & ALT.vcf == REF.pred))
  
  # step 2: strand flip
  # flip the reference allele in predction model
  ref = vcf.mapped$REF.pred
  ref.flip = ref
  ref.flip[ref == "A"] = "T"
  ref.flip[ref == "T"] = "A"
  ref.flip[ref == "C"] = "G"
  ref.flip[ref == "G"] = "C"
  
  # flip the alternative (effect) allele in predction model
  alt = vcf.mapped$ALT.pred
  alt.flip = alt
  alt.flip[alt == "A"] = "T"
  alt.flip[alt == "T"] = "A"
  alt.flip[alt == "C"] = "G"
  alt.flip[alt == "G"] = "C"
  
  vcf.mapped = vcf.mapped %>% 
    mutate(REF.flip.pred = ref.flip, ALT.flip.pred = alt.flip) %>% 
    filter(!(REF.vcf == REF.flip.pred & ALT.vcf != ALT.flip.pred)) %>% 
    filter(!(REF.vcf != REF.flip.pred & ALT.vcf == ALT.flip.pred)) %>% 
    filter(!(REF.vcf == ALT.flip.pred & ALT.vcf != REF.flip.pred)) %>% 
    filter(!(REF.vcf != ALT.flip.pred & ALT.vcf == REF.flip.pred)) %>% 
    mutate(flip = (REF.vcf == ALT.flip.pred & ALT.vcf == REF.flip.pred) | (REF.vcf == ALT.pred & ALT.vcf == REF.pred))
  
  if(sum(vcf.mapped$flip == T) > 0){
    vcf.mapped[vcf.mapped$flip == T, ]$weight = -1 * vcf.mapped[vcf.mapped$flip == T, ]$weight
  }
  
  vcf.mapped = vcf.mapped %>% 
    select(rsid, gene, weight, REF.vcf, ALT.vcf, CHR, POS)
  
  return(vcf.mapped)
}











