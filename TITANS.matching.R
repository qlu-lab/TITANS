matching.3sib = function(parental.strand){
  ge.gamete = t(w$weight[weight.idx]) %*% parental.strand
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
  nchild = ncol(child.g)-9
  ge.sib = c(ge.sib)
  names(ge.sib) = rep(unique(parent_info[,1]), each = 4)
  ge.sib = data.frame(rep(unique(parent_info[,1]), each = 4), ge.sib, rep(0, length(ge.sib))) %>% 
    `colnames<-`(c("fam","expr","disease"))
  ge.sib.wide = ge.sib %>% 
    mutate(sample.id = rep(c("s1", "s2", "s3", "s4"), nrow(ge.sib)/4), fam = as.character(fam)) %>% 
    select(fam, sample.id, expr) %>% 
    spread(., sample.id, expr)
  
  GE = left_join(ge.c, ge.sib.wide, by = "fam") %>% 
    select(fam, sample, expr, s1:s4) %>% 
    gather(disease, ge, expr:s4) %>% 
    arrange(sample) %>% 
    mutate(disease = ifelse(disease == "expr", 1, 0)) %>% 
    select(sample, ge, disease) %>% 
    `colnames<-`(c("fam", "expr", "disease"))
  
  remove = list()
  for(i in seq(from = 1, to = nrow(GE)-4, by = 5)){
    index = which(abs(GE[(i+1):(i+4),]$expr - GE[i,]$expr) == min(abs(GE[(i+1):(i+4),]$expr-GE[i,]$expr)))[1] + i
    remove = append(remove, index)
  }
  remove = unlist(remove) #index that the sibling having the same genotype as child's
  GE.final = GE[-remove, ]
  return(GE.final)
}

matching.1sib = function(parental.strand){
  ge.gamete = t(w$weight[weight.idx]) %*% parental.strand
  ge.sib = data.frame()
  for(i in seq(from = 1, to = (ncol(ge.gamete)-3), by = 4)){
    if(i != 1){
      ge.sib = cbind(ge.sib,
                     ge.gamete[i] + ge.gamete[i+1] + ge.gamete[i+2] + ge.gamete[i+3])
    } else{
      ge.sib = cbind(ge.gamete[i] + ge.gamete[i+1] + ge.gamete[i+2] + ge.gamete[i+3])
    }
  }
  nchild = ncol(child.g)-9
  ge.sib = c(ge.sib)
  names(ge.sib) = rep(unique(parent_info[,1]), each = 1)
  ge.sib = data.frame(rep(unique(parent_info[,1]), each = 1), ge.sib, rep(0, length(ge.sib))) %>% 
    `colnames<-`(c("fam","expr","disease"))
  ge.sib.wide = ge.sib %>% 
    mutate(sample.id = rep(c("s1"), nrow(ge.sib)), fam = as.character(fam)) %>% 
    select(fam, sample.id, expr) %>% 
    spread(., sample.id, expr)
  GE.final = left_join(ge.c, ge.sib.wide, by = "fam") %>% 
    mutate(s1 = s1 - expr) %>% 
    select(fam, sample, expr, s1) %>% 
    gather(disease, ge, expr:s1) %>% 
    arrange(sample) %>% 
    mutate(disease = ifelse(disease == "expr", 1, 0)) %>% 
    select(sample, ge, disease) %>% 
    `colnames<-`(c("fam", "expr", "disease"))
  
  return(GE.final)
}

matching.2parent = function(parental.strand){
  ge.gamete = t(w$weight[weight.idx]) %*% parental.strand
  ge.sib = data.frame() # match with 2 parents
  for(i in seq(from = 1, to = (ncol(ge.gamete)-3), by = 4)){
    if(i != 1){
      ge.sib = cbind(ge.sib,
                     ge.gamete[i] + ge.gamete[i+1],
                     ge.gamete[i+2] + ge.gamete[i+3])
    } else{
      ge.sib = cbind(ge.gamete[i] + ge.gamete[i+1],
                     ge.gamete[i+2] + ge.gamete[i+3])
    }
  }
  nchild = ncol(child.g)-9
  ge.sib = c(ge.sib)
  names(ge.sib) = rep(unique(parent_info[,1]), each = 2)
  ge.sib = data.frame(rep(unique(parent_info[,1]), each = 2), ge.sib, rep(0, length(ge.sib))) %>% 
    `colnames<-`(c("fam","expr","disease"))
  ge.sib.wide = ge.sib %>% 
    mutate(sample.id = rep(c("s1", "s2"), nrow(ge.sib)/2), fam = as.character(fam)) %>% 
    select(fam, sample.id, expr) %>% 
    spread(., sample.id, expr)
  GE.final = left_join(ge.c, ge.sib.wide, by = "fam") %>% 
    select(fam, sample, expr, s1:s2) %>% 
    gather(disease, ge, expr:s2) %>% 
    arrange(sample) %>% 
    mutate(disease = ifelse(disease == "expr", 1, 0)) %>% 
    select(sample, ge, disease) %>% 
    `colnames<-`(c("fam", "expr", "disease"))
  
  return(GE.final)
}


matching = function(matching, parental.strand){
  if(matching == "3sib"){
    GE = matching.3sib(parental.strand)
  } else if (matching == "1sib"){
    GE = matching.1sib(parental.strand)
  } else if (matching == "2parent"){
    GE = matching.2parent(parental.strand)
  }
  return(GE)
}