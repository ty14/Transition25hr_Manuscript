
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

#aggression limma models 
my_logFC_threshold = 0.2

#agg rec 70 min 
limma_list <- readRDS("manuscript/brain/results/limma_mPFC_CORT_AGGREC.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= 0.2)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

cort_r <- limma_list$cort

ar<- limma_list$aggrec

c_ar<- limma_list$cort_aggrec


des <- read_csv("manuscript/brain/results_tables/des_genes_mPFC.csv")

dar <- des  %>%  filter(time == 70) %>% filter(symbol %in% ar$symbol)
dar
#one gene Lmcd1


dcar <- des  %>%  filter(time == 70) %>% filter(symbol %in% c_ar$symbol)
dcar
#4 genes #Wfdc18, Cmah, Plin2, Cndp1


#agg_given 70 min 
limma_list <- readRDS("manuscript/brain/results/limma_mPFC_CORT_AGGgiven.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= 0.2)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

cort_g <- limma_list$cort

ag<- limma_list$agiven

c_ag<- limma_list$cort_agiven


asc <- read_csv("manuscript/brain/results_tables/asc_genes_mPFC.csv")

ag <- asc %>%  filter(time == 70) %>% filter(symbol %in% ag$symbol)
ag
#9 genes Gm10130, Med30 , Xcr1, Pop4, Afap1l1, Pclaf, Pus10, Olfr715, Proz


cag <- asc %>%  filter(time == 70) %>% filter(symbol %in% c_ag$symbol)
cag
#3 genes Tfap2a, Ddt, Kcns3


#agg rec 25hr 
limma_list <- readRDS("manuscript/brain/results/limma_mPFC25_CORT_AGGREC.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= 0.2)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

cort_r <- limma_list$cort

ar<- limma_list$aggrec

c_ar<- limma_list$cort_aggrec


des <- read_csv("manuscript/brain/results_tables/des_genes_mPFC.csv")

dar <- des  %>%  filter(time == 25) %>% filter(symbol %in% ar$symbol)
dar
#2 genes Cbr3, Vm2r31 


dcar <- des  %>%  filter(time == 25) %>% filter(symbol %in% c_ar$symbol)
dcar
#8 genes Fam183b, Rem2, Lrrc10b, Rxrg, Neurdo1, Sulf1, Vmn1r88, Inf2


#agg_given 25hr
limma_list <- readRDS("manuscript/brain/results/limma_mPFC25_CORT_AGGiven.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= 0.2)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

cort_g <- limma_list$cort

ag<- limma_list$aggiven

c_ag<- limma_list$cort_aggiven


asc <- read_csv("manuscript/brain/results_tables/asc_genes_mPFC.csv")

agx <- asc %>%  filter(time == 25) %>% filter(symbol %in% ag$symbol)
agx
#9 genes Gm10130, Med30 , Xcr1, Pop4, Afap1l1, Pclaf, Pus10, Olfr715, Proz


cag <- asc %>%  filter(time == 25) %>% filter(symbol %in% c_ag$symbol)
# one gene Bhlhe22
