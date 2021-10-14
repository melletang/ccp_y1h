library(tidyverse)
library(readxl)
network <- readxl::read_excel("MSB-2021-10625-DatasetEV3-Network.xlsx")

pathtflist<- lapply(unique(network$Target_Pathway), function(x){
  unique(network$TF_AGI[which(network$Target_Pathway==x)])
})
names(pathtflist) <- unique(network$Target_Pathway)


pathtfsize<- lapply(unique(network$Target_Pathway), function(x){
  unique(network$TF_AGI[which(network$Target_Pathway==x)]) %>% length()
})
names(pathtfsize) <- unique(network$Target_Pathway)
pathtfsize <- do.call(rbind, pathtfsize)
names(pathtfsize) <- c("Pathway", "N")

#write.table(pathtfsize, "pathway_unique_tfsize.txt", sep = "\t",  quote = F)


pathtfcombn <- data.frame(t(combn(unique(network$Target_Pathway), 2)))

names(pathtfcombn) <- c("Var1", "Var2")

pathtfcombn$Var1 <- as.vector(pathtfcombn$Var1)
pathtfcombn$Var2 <- as.vector(pathtfcombn$Var2)

pathtfcombn$TF_intersect <- apply(pathtfcombn, 1, function(x){
  a <- which(names(pathtflist)==x[1])
  b <- which(names(pathtflist)==x[2])
  intersect(pathtflist[[a]], pathtflist[[b]]) %>% unlist() %>% length()
})


pathtfcombn <- pathtfcombn %>% filter(Var1 != Var2)

pathtfcombn$Var1Size <- apply(pathtfcombn, 1, function(x){
  pathtflist[which(names(pathtflist)==x[1])] %>% unlist() %>% length()
})
pathtfcombn$Var2Size <- apply(pathtfcombn, 1, function(x){
  pathtflist[which(names(pathtflist)==x[2])] %>% unlist() %>% length()
})

pathtfcombn$fisherpval <- apply(pathtfcombn, 1, function(x){
  a <- as.numeric(x[3]) #overlap
  b <- as.numeric(x[4]) - a
  c <- as.numeric(x[5]) - a
  d <- 2039 - (a + b + c)
  m <- matrix(c(a,b,c,d), nrow = 2)
  t <- fisher.test(m, alternative = "g")
  return(t[[1]])
})

pathtfcombn$odds_ratio <- apply(pathtfcombn, 1, function(x){
  a <- as.numeric(x[3]) #overlap
  b <- as.numeric(x[4]) - a
  c <- as.numeric(x[5]) - a
  d <- 2039 - (a + b + c)
  m <- matrix(c(a,b,c,d), nrow = 2)
  t <- fisher.test(m, alternative = "g")
  return(t[[3]])
})

pathtfcombn$adjpval <- p.adjust(pathtfcombn$fisherpval, method = "BH")
