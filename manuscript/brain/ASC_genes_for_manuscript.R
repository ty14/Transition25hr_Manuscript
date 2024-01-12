library(tidyverse)
#70min
# ASCvsDOM = 710 genes, 385 up & 325 down
# ASCvsSUB = 863 genes, 396 up & 467 down

asc <- read_csv("manuscript/brain/results_tables/asc_genes_mPFC.csv")
colnames(asc)

asc70_up <-  asc %>% filter(time == 70) %>% 
  arrange(-as_logFC) %>% filter(as_logFC > 0.2) %>% 
  select(symbol,time) 
#94 genes
unique(asc70_up$symbol)

asc70_down <-  asc %>% filter(time == 70) %>% 
  arrange(as_logFC) %>% filter(as_logFC < 0.2) %>% 
  select(symbol,time)
#92 genes
unique(asc70_down$symbol)

asc70 <- asc70_up %>% rbind(asc70_down)

#25hr
# ASCvsDOM = 493 genes, 252 up & 241 down
# ASCvsSUB = 630 genes, 331 up & 299 down


asc25_up <-  asc %>% filter(time == 25) %>% 
  arrange(-as_logFC) %>% filter(as_logFC > 0.2) %>% 
  select(symbol,time)
#63 genes 
unique(asc25_up$symbol)

asc25_down <-  asc %>% filter(time == 25) %>% 
  arrange(as_logFC) %>% filter(as_logFC < 0.2) %>% 
  select(symbol,time)
#54 genes 
unique(asc25_down$symbol)

my_logFC_threshold = 0.2

limma_list<- readRDS("manuscript/brain/results_use/limma_PFC70min_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

adom <- limma_list$ascdom
as <- limma_list$ascsub

x <- adom$symbol[adom$symbol %in% as$symbol] 
xx <- unique(x) #187

xx[129] #"Sash1" opposite direction. 

adom %>% filter(symbol =="Sash1")
as %>% filter(symbol =="Sash1")

adom_up <- adom %>% filter(logFC > 0.2)
as_up <- as %>% filter(logFC > 0.2)

adom_down <- adom %>% filter(logFC < 0.2)
as_down <- as %>% filter(logFC < 0.2)



# upsetter Plot 
library(UpSetR)
library(workflowr)
library(ComplexUpset)
listInput <- list(adom_up= adom_up$symbol, as_up=as_up$symbol, adom_down = adom_down$symbol, 
                  as_down = as_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)
mat <- matrix(c(94,0,1,92),ncol=2)
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

adom25 <- limma_list$ascdom
as25 <- limma_list$ascsub

x <- adom25$symbol[adom25$symbol %in% as25$symbol] 
xx <- unique(x) #117


adom_up <- adom25 %>% filter(logFC > 0.2)
as_up <- as25 %>% filter(logFC > 0.2)

adom_down <- adom25 %>% filter(logFC < 0.2)
as_down <- as25 %>% filter(logFC < 0.2)



# upsetter Plot 
library(UpSetR)
library(workflowr)
library(ComplexUpset)
listInput <- list(adom_up= adom_up$symbol, as_up=as_up$symbol, adom_down = adom_down$symbol, 
                  as_down = as_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

mat <- matrix(c(63,0,0,54),ncol=2)
# Perform chi-squared test
chi_square <- chisq.test(mat)
chi_square
