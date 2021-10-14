library(tidyverse)
library(readxl)
network <- readxl::read_excel("MSB-2021-10625-DatasetEV3-Network.xlsx")

binding_summary <- network %>% select(Promoter_AGI, Target_Pathway) %>% unique() %>% group_by(Target_Pathway) %>% 
  tally() %>% rename(num_gene = n)
binding_summary <- left_join(binding_summary, 
                             network %>% select(TF_AGI, Target_Pathway) %>%
                               unique() %>% group_by(Target_Pathway) %>% tally() %>%
                               rename(num_tf = n))

binding_summary <- left_join(binding_summary, 
                             network %>% select(TF_AGI, Promoter_AGI, Target_Pathway) %>%
                               unique() %>% group_by(Target_Pathway) %>% tally() %>%
                               rename(num_int = n))
panel_b <- ggplot(binding_summary, aes(reorder(Target_Pathway,num_gene), num_gene)) + geom_bar(stat = "identity", fill = "black") + coord_flip() + theme_bw() +
  ylab("Number of genes") + xlab("Pathway") + theme(
    axis.text = element_text(color = "black", size = "10"),
    axis.title = element_text(color = "black", size = "10")
  )

panel_c <- ggplot(binding_summary, aes(reorder(Target_Pathway,num_gene), num_tf)) + geom_bar(stat = "identity", fill = "black") + coord_flip() + theme_bw() +
  ylab("Number of TFs") + xlab("Pathway") + theme(
    axis.text = element_text(color = "black", size = "10"),
    axis.title = element_text(color = "black", size = "10"),
    plot.margin = unit(c(0, 0.5, 0, 0), "cm")
  ) 

panel_d <- ggplot(binding_summary, aes(reorder(Target_Pathway,num_gene), num_int)) + geom_bar(stat = "identity", fill = "black") + 
  coord_flip() + theme_bw() + ylab("Number of interactions") + xlab("Pathway") + theme(
    axis.text = element_text(color = "black", size = "10"),
    axis.title = element_text(color = "black", size = "10")) 

num_path <- network %>% select(TF_AGI, Target_Pathway) %>% unique() %>% group_by(TF_AGI) %>% tally()

numpathbar <- num_path %>% group_by(n) %>% tally()

panel_e <- ggplot(numpathbar, aes(n, nn)) + geom_bar(stat = "identity", fill = "black")+ theme_bw() + ylab("Number of TFs") + xlab("Number of Pathways") + theme(
  axis.text = element_text(color = "black", size = "10"),
  axis.title = element_text(color = "black", size = "10")) + scale_x_continuous(breaks=seq(0,12,1))


