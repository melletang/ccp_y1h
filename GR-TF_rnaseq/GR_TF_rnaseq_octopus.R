library(tidyverse)
library(octopus)
library(Biobase)

count_sum <- read.csv("ATOEMT_star_readcountsum_id.csv", row.names = 1)

count_sum<-count_sum[,colnames(count_sum)[order(colnames(count_sum))]]

count_sum <- count_sum[rowSums(count_sum) > 0,]
count_sum <- as.matrix(count_sum)

y = DGEList(count_sum)

count_sumannot <- data.frame(sample = colnames(count_sum))
count_sumannot$TF <- as.factor(str_extract(count_sumannot$sample, "GR_[A-Z0-9]+"))
count_sumannot$Allele <- as.factor(str_extract(count_sumannot$sample, "L[0-9]+"))
count_sumannot$Treatment <- as.factor(str_extract(count_sumannot$sample, "CONTROL|DEX"))

#####GR-CHA19######

GRCHA19 <- count_sum[,which(str_detect(colnames(count_sum), "GR_CHA19"))]
GRCHA19 <- as.matrix(GRCHA19[Biobase::rowMedians(as.matrix(GRCHA19))>0,])
GRCHA19 <- octopus.normalize(GRCHA19, dist_dir = "results/GRCHA19_results/")

GRCHA19_key <- data.frame(row.names = colnames(GRCHA19),
                          Allele = str_extract(colnames(GRCHA19), "L[0-9]+"),
                          Treatment = str_extract(colnames(GRCHA19), "CONTROL|DEX"),
                          Rep = str_extract(colnames(GRCHA19),"A|B"))

GRCHA19 <- t(GRCHA19)

GRCHA19results <- octopus.glm.nb(~Allele*Treatment, GRCHA19, GRCHA19_key, specs4lsmeans = ~Allele*Treatment)

GRCHA19anova <- data.frame(t(GRCHA19results$anova))
GRCHA19anova$`adj.p.Allele` <- p.adjust(GRCHA19anova$Allele_Pr..Chi., method = "BH")
GRCHA19anova$`adj.p.Treatment` <- p.adjust(GRCHA19anova$Treatment_Pr..Chi., method = "BH")
GRCHA19anova$`adj.p.Allele.Treatment` <- p.adjust(GRCHA19anova$Allele.Treatment_Pr..Chi., method = "BH")
write.csv(GRCHA19anova, "results/GRCHA19_results/GRCHA19_anova_padj.csv", row.names = T)
write.csv(data.frame(t(GRCHA19results$lsm)), "results/GRCHA19_results/GRCHA19_lsm.csv", row.names = T)
write.csv(data.frame(t(GRCHA19results$se)), "results/GRCHA19_results/GRCHA19_se.csv", row.names = T)


#####GR-ENAP1######


GRENAP1 <- count_sum[,which(str_detect(colnames(count_sum), "GR_ENAP1"))]
GRENAP1 <- as.matrix(GRENAP1[Biobase::rowMedians(as.matrix(GRENAP1))>0,])
GRENAP1 <- octopus.normalize(GRENAP1, dist_dir = "results/GRENAP1_results/")

GRENAP1_key <- data.frame(row.names = colnames(GRENAP1),
                          Allele = str_extract(colnames(GRENAP1), "L[0-9]+"),
                          Treatment = str_extract(colnames(GRENAP1), "CONTROL|DEX"),
                          Rep = str_extract(colnames(GRENAP1),"A|B"))

GRENAP1 <- t(GRENAP1)

GRENAP1results <- octopus.glm.nb(~Allele*Treatment, GRENAP1, GRENAP1_key, specs4lsmeans = ~Allele*Treatment)

GRENAP1anova <- data.frame(t(GRENAP1results$anova))
GRENAP1anova$`adj.p.Allele` <- p.adjust(GRENAP1anova$Allele_Pr..Chi., method = "BH")
GRENAP1anova$`adj.p.Treatment` <- p.adjust(GRENAP1anova$Treatment_Pr..Chi., method = "BH")
GRENAP1anova$`adj.p.Allele.Treatment` <- p.adjust(GRENAP1anova$Allele.Treatment_Pr..Chi., method = "BH")
write.csv(GRENAP1anova, "results/GRENAP1_results/GRENAP1_anova_padj.csv", row.names = T)
write.csv(data.frame(t(GRENAP1results$lsm)), "results/GRENAP1_results/GRENAP1_lsm.csv", row.names = T)
write.csv(data.frame(t(GRENAP1results$se)), "results/GRENAP1_results/GRENAP1_se.csv", row.names = T)


#####GR-LBD16######


GRLBD16 <- count_sum[,which(str_detect(colnames(count_sum), "GR_LBD16"))]
GRLBD16 <- as.matrix(GRLBD16[Biobase::rowMedians(as.matrix(GRLBD16))>0,])
GRLBD16 <- octopus.normalize(GRLBD16, dist_dir = "results/GRLBD16_results/")

GRLBD16_key <- data.frame(row.names = colnames(GRLBD16),
                          Allele = str_extract(colnames(GRLBD16), "L[0-9]+"),
                          Treatment = str_extract(colnames(GRLBD16), "CONTROL|DEX"),
                          Rep = str_extract(colnames(GRLBD16),"A|B"))

GRLBD16 <- t(GRLBD16)

GRLBD16results <- octopus.glm.nb(~Allele*Treatment, GRLBD16, GRLBD16_key, specs4lsmeans = ~Allele*Treatment)

GRLBD16anova <- data.frame(t(GRLBD16results$anova))
GRLBD16anova$`adj.p.Allele` <- p.adjust(GRLBD16anova$Allele_Pr..Chi., method = "BH")
GRLBD16anova$`adj.p.Treatment` <- p.adjust(GRLBD16anova$Treatment_Pr..Chi., method = "BH")
GRLBD16anova$`adj.p.Allele.Treatment` <- p.adjust(GRLBD16anova$Allele.Treatment_Pr..Chi., method = "BH")
write.csv(GRLBD16anova, "results/GRLBD16_results/GRLBD16_anova_padj.csv", row.names = T)
write.csv(data.frame(t(GRLBD16results$lsm)), "results/GRLBD16_results/GRLBD16_lsm.csv", row.names = T)
write.csv(data.frame(t(GRLBD16results$se)), "results/GRLBD16_results/GRLBD16_se.csv", row.names = T)


#####GR-WRI3######


GRWRI3 <- count_sum[,which(str_detect(colnames(count_sum), "GR_WRI3"))]
GRWRI3 <- as.matrix(GRWRI3[Biobase::rowMedians(as.matrix(GRWRI3))>0,])
GRWRI3 <- octopus.normalize(GRWRI3, dist_dir = "results/GRWRI3_results/")

GRWRI3_key <- data.frame(row.names = colnames(GRWRI3),
                          Allele = str_extract(colnames(GRWRI3), "L[0-9]+"),
                          Treatment = str_extract(colnames(GRWRI3), "CONTROL|DEX"),
                          Rep = str_extract(colnames(GRWRI3),"A|B"))

GRWRI3 <- t(GRWRI3)

GRWRI3results <- octopus.glm.nb(~Allele*Treatment, GRWRI3, GRWRI3_key, specs4lsmeans = ~Allele*Treatment)

GRWRI3anova <- data.frame(t(GRWRI3results$anova))
GRWRI3anova$`adj.p.Allele` <- p.adjust(GRWRI3anova$Allele_Pr..Chi., method = "BH")
GRWRI3anova$`adj.p.Treatment` <- p.adjust(GRWRI3anova$Treatment_Pr..Chi., method = "BH")
GRWRI3anova$`adj.p.Allele.Treatment` <- p.adjust(GRWRI3anova$Allele.Treatment_Pr..Chi., method = "BH")
write.csv(GRWRI3anova, "results/GRWRI3_results/GRWRI3_anova_padj.csv", row.names = T)
write.csv(data.frame(t(GRWRI3results$lsm)), "results/GRWRI3_results/GRWRI3_lsm.csv", row.names = T)
write.csv(data.frame(t(GRWRI3results$se)), "results/GRWRI3_results/GRWRI3_se.csv", row.names = T)


