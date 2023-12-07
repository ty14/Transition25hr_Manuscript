library(annotables)
library(tidyverse)
grcm38 # mouse genes




asc <- read.csv("manuscript/brain/results_tables/asc_genes_mPFC.csv")
head(asc)

a_go <- asc %>% select(logFC = as_logFC,symbol, time )

a_go70 <- a_go %>% filter(time == 70)
a_go25 <- a_go%>% filter(time == 25)
source("functions/gettop10GO.R")

gettop10GO(a_go70, my_showCategory) %>% 
  mutate(comparison = "ASC vs. Stable", time = "70") -> top10go1

gettop10GO(a_go25, my_showCategory) %>% 
  mutate(comparison = "ASC vs. Stable", time = "25") -> top10go2
#descenders: 
des <- read.csv("manuscript/brain/results_tables/des_genes_mPFC.csv")
head(des)

d_go <- des %>% select(logFC = dd_logFC,symbol, time )


d_go70 <- d_go %>% filter(time == 70)
d_go25 <- d_go %>% filter(time == 25)


gettop10GO(d_go70, my_showCategory) %>% 
  mutate(comparison = "DES vs. Stable", time = "70") -> top10go3

gettop10GO(d_go25, my_showCategory) %>% 
  mutate(comparison = "DES vs. Stable", time = "25") -> top10go4


rbind(top10go1,top10go2, top10go3,top10go4) -> top10_GOterms

write.csv(top10_GOterms,"manuscript/brain/results_tables/topBP_GOterms_ASC_DES_stableXtime.csv", row.names = F)

# view
top10_GOterms$Description[1:5]
top10_GOterms$Description[6:10]
top10_GOterms$Description[11:15]
top10_GOterms$Description[16:20]
top10_GOterms$Description[21:25]

