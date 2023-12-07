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

#ALL GROUPS 
my_logFC_threshold = 0.2

limma_list<- readRDS("manuscript/brain/results_use/limma_PFC70min_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

y1a <- limma_list$ascdom

y2a <- limma_list$domsub

y3a <- limma_list$desdom

y4a <- limma_list$desasc

y5a <- limma_list$ascsub

y6a <- limma_list$dessub



source("functions/gettop10GO.R")

gettop10GO(y1a, my_showCategory) %>% 
  mutate(comparison = "ASC-DOM") -> top10go1

gettop10GO(y2a, my_showCategory ) %>% 
  mutate(comparison = "DOM-SUB") -> top10go2

gettop10GO(y3a, my_showCategory ) %>% 
  mutate(comparison = "DES-DOM") -> top10go3

gettop10GO(y4a, my_showCategory) %>% 
  mutate(comparison = "DES-ASC") -> top10go4

gettop10GO(y5a, my_showCategory ) %>% 
  mutate(comparison = "ASC-SUB") -> top10go5

gettop10GO(y6a, my_showCategory ) %>% 
  mutate(comparison = "DES-SUB") -> top10go6


rbind(top10go1,top10go2,top10go3,top10go4,top10go5,top10go6) -> top10_GOterms

write.csv(top10_GOterms,"manuscript/brain/results_tables/topBP_GOterms_mPFC70min_NORM.csv", row.names = F)

#25 hour data all groups 

limma_list<- readRDS("manuscript/brain/results_use/limma_PFC25hr_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

y1a <- limma_list$ascdom

y2a <- limma_list$domsub

y3a <- limma_list$desdom

y4a <- limma_list$desasc

y5a <- limma_list$ascsub

y6a <- limma_list$dessub

gettop10GO(y1a, my_showCategory) %>% 
  mutate(comparison = "ASC-DOM") -> top10go1

gettop10GO(y2a, my_showCategory ) %>% 
  mutate(comparison = "DOM-SUB") -> top10go2

gettop10GO(y3a, my_showCategory ) %>% 
  mutate(comparison = "DES-DOM") -> top10go3

gettop10GO(y4a, my_showCategory) %>% 
  mutate(comparison = "DES-ASC") -> top10go4

gettop10GO(y5a, my_showCategory ) %>% 
  mutate(comparison = "ASC-SUB") -> top10go5

gettop10GO(y6a, my_showCategory ) %>% 
  mutate(comparison = "DES-SUB") -> top10go6

rbind(top10go1,top10go2,top10go3,top10go4,top10go5,top10go6) -> top10_GOterms

write.csv(top10_GOterms,"manuscript/brain/results_tables/topBP_GOterms_mPFC25hr_NORM.csv", row.names = F)
