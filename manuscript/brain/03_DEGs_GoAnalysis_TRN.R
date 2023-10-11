# libraries 
library(limma)
library(edgeR)
library(Mus.musculus)
organism = 'org.Mm.eg.db'
library(organism, character.only = TRUE)
library(biomaRt)
library(AnnotationDbi)
library(pheatmap)
library(annotables)
library(clusterProfiler)
library(enrichplot)
library(organism, character.only = TRUE)
library(DOSE)
library(tidyverse)
grcm38 # mouse genes

#TRN GROUPS
my_logFC_threshold = 0.2

limma_list<- readRDS("manuscript/brain/results/limma_PFC70_TRN.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

y1a <- limma_list$tdom

y2a <- limma_list$tsub

y3a <- limma_list$ds


source("functions/gettop10GO.R")


gettop10GO(y1a, my_showCategory) %>% 
  mutate(comparison = "TRN-DOM") -> top10go1

gettop10GO(y2a, my_showCategory ) %>% 
  mutate(comparison = "TRN-SUB") -> top10go2

gettop10GO(y3a, my_showCategory ) %>% 
  mutate(comparison = "DOM-SUB") -> top10go3


rbind(top10go1,top10go2,top10go3) -> top10_GOterms


write.csv(top10_GOterms,"manuscript/brain/results_tables/topBP_GOterms_mPFC70min_TRN.csv", row.names = F)

#25 hour data TRN GROUPS

limma_list<- readRDS("manuscript/brain/results/limma_PFC_TRN.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

y1a <- limma_list$tdom

y2a <- limma_list$tsub

y3a <- limma_list$ds



gettop10GO(y1a, my_showCategory) %>% 
  mutate(comparison = "TRN-DOM") -> top10go1

gettop10GO(y2a, my_showCategory ) %>% 
  mutate(comparison = "TRN-SUB") -> top10go2

gettop10GO(y3a, my_showCategory ) %>% 
  mutate(comparison = "DOM-SUB") -> top10go3


rbind(top10go1,top10go2,top10go3) -> top10_GOterms

write.csv(top10_GOterms,"manuscript/brain/results_tables/topBP_GOterms_mPFC25hr_TRN.csv", row.names = F)
