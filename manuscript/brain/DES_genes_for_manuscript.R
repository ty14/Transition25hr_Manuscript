library(tidyverse)
#70min
# desvsDOM = 813 genes, 366 up & 447 down
# desvsSUB =  1637 genes, 842 up & 795 down

des <- read_csv("manuscript/brain/results_tables/des_genes_mPFC.csv")
colnames(des)

des70_up <-  des %>% filter(time == 70) %>% 
  arrange(-dd_logFC) %>% filter(dd_logFC > 0.2) %>% 
  select(symbol,time) 
#168 genes
unique(des70_up$symbol)

des70_down <-  des %>% filter(time == 70) %>% 
  arrange(dd_logFC) %>% filter(dd_logFC < 0.2) %>% 
  select(symbol,time)
#204 genes
unique(des70_down$symbol)

des70 <- des70_up %>% rbind(des70_down)

#25hr
# desvsDOM = 385 genes, 216 up & 169 down
# desvsSUB = 818 genes, 452 up & 366 down

des25_up <-  des %>% filter(time == 25) %>% 
  arrange(-dd_logFC) %>% filter(dd_logFC > 0.2) %>% 
  select(symbol,time)
#72 genes 
unique(des25_up$symbol)

des25_down <-  des %>% filter(time == 25) %>% 
  arrange(dd_logFC) %>% filter(dd_logFC < 0.2) %>% 
  select(symbol,time)
#58 genes 
unique(des25_down$symbol)

my_logFC_threshold = 0.2

limma_list<- readRDS("manuscript/brain/results_use/limma_PFC70min_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

dd <- limma_list$desdom
ds <- limma_list$dessub

x <- dd$symbol[dd$symbol %in% ds$symbol] 
xx <- unique(x) #372

dd_up <- add %>% filter(logFC > 0.2)
ds_up <- ds %>% filter(logFC > 0.2)

dd_down <- dd %>% filter(logFC < 0.2)
ds_down <- ds %>% filter(logFC < 0.2)

# upsetter Plot 
library(UpSetR)
library(workflowr)
library(ComplexUpset)
listInput <- list(dd_up= dd_up$symbol, ds_up=ds_up$symbol, dd_down = dd_down$symbol, 
                  ds_down = ds_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

mat <- matrix(c(168,0,0,204),ncol=2)
# Perform chi-squared test
chi_square <- chisq.test(mat)
chi_square


dd_down$symbol[dd_down$symbol %in% ds_down$symbol]

limma_list<- readRDS("manuscript/brain/results_use/limma_PFC70min_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

dd <- limma_list$desdom
ds <- limma_list$dessub

x <- dd$symbol[dd$symbol %in% ds$symbol] 
xx <- unique(x) #372

dd_up <- add %>% filter(logFC > 0.2)
ds_up <- ds %>% filter(logFC > 0.2)

dd_down <- dd %>% filter(logFC < 0.2)
ds_down <- ds %>% filter(logFC < 0.2)

# upsetter Plot 
library(UpSetR)
library(workflowr)
library(ComplexUpset)
listInput <- list(dd_up= dd_up$symbol, ds_up=ds_up$symbol, dd_down = dd_down$symbol, 
                  ds_down = ds_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

mat <- matrix(c(72,0,0,58),ncol=2)
# Perform chi-squared test
chi_square <- chisq.test(mat)
chi_square


# 25hr 
limma_list<- readRDS("manuscript/brain/results_use/limma_PFC25hr_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

dd25 <- limma_list$desdom
ds25 <- limma_list$dessub

x <- dd25$symbol[dd25$symbol %in% ds25$symbol] 
xx <- unique(x) #130


adom_up <- dd25 %>% filter(logFC > 0.2)
as_up <- ds25 %>% filter(logFC > 0.2)

adom_down <- dd25 %>% filter(logFC < 0.2)
as_down <- ds25 %>% filter(logFC < 0.2)


# upsetter Plot 
library(UpSetR)
library(workflowr)
library(ComplexUpset)
listInput <- list(adom_up= adom_up$symbol, as_up=as_up$symbol, adom_down = adom_down$symbol, 
                  as_down = as_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

mat <- matrix(c(72,0,0,58),ncol=2)
# Perform chi-squared test
chi_square <- chisq.test(mat)
chi_square

mat <- matrix(c(58,11,18,65),ncol=2)
# Perform chi-squared test
chi_square <- chisq.test(mat)
chi_square


dd %>% filter(symbol %in% c("Pyroxd2", "Etv6", "Pitpnm3", "Cxxc5", "Strip2"))
