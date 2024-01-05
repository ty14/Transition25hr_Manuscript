# libraries
library(limma)
library(edgeR)
library(Mus.musculus)
organism = 'org.Mm.eg.db'
organism1 = 'org.Hs.eg.db'
library(organism, character.only = TRUE)
library(biomaRt)
library(AnnotationDbi)
library(pheatmap)
library(annotables)
library(clusterProfiler)
library(enrichplot)
library(organism, character.only = TRUE)
library(organism1, character.only = TRUE)
library(DOSE)
library(EnhancedVolcano)
library(Orthology.eg.db)
library(tidyverse)
library(msigdbr)
msigdbr_collections() %>% as.data.frame()
#immune =  C7     IMMUNESIGDB
# biological pathways = C2 CP:WIKIPATHWAYS  
#hallmark = H
#reactome = immune, metabolic 


msigdbr(species = "Mus musculus", subcategory ="CP:REACTOME") -> reactome_sets
unique(reactome_sets$gs_name)


grcm38 # mouse genes
my_logFC_threshold = 0.2
limma_list<- readRDS("manuscript/brain/results_use/limma_PFC70min_NORM_RG.RDS") %>%
  map(~distinct(.)) %>%
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>%
  map(~filter(.,!is.na(entrez)))
y4a <- limma_list$desdom
y3a <- limma_list$dessub
y2a <- limma_list$ascdom
y1a <- limma_list$ascsub

# 
# #1 ATP metablic Energy 
# msigdbr(species = "Mus musculus", subcategory ="GO:BP") -> bp_sets
# # [2] "GOBP_ATP_METABOLIC_PROCESS"
# bp_sets %>%
#   filter(grepl("GOBP_ATP_METABOLIC_PROCESS",gs_name,ignore.case = T)) %>%
#   select(gene_symbol) %>%  as_tibble() %>%
#   .$gene_symbol -> enx
# 
# y4a  %>%
#   filter(symbol %in% enx) %>%
#   dplyr::select(symbol, logFC, P.Value) %>%
#   dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
#                              ifelse(logFC>my_logFC_threshold,"DES genes","DOM genes"))) %>%
#   dplyr::mutate(Sig = factor(Sig, levels = c("DES genes","DOM genes", "N.S."))) %>%
#   dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
#   unique() -> dd_atp
# 
# 
# y3a  %>%
#   filter(symbol %in% enx) %>%
#   dplyr::select(symbol, logFC, P.Value) %>%
#   dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
#                              ifelse(logFC>my_logFC_threshold,"DES genes","SUB genes"))) %>%
#   dplyr::mutate(Sig = factor(Sig, levels = c("DES genes","SUB genes", "N.S."))) %>%
#   dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
#   unique() -> ds_atp
# 
# y2a  %>%
#   filter(symbol %in% enx) %>%
#   dplyr::select(symbol, logFC, P.Value) %>%
#   dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
#                              ifelse(logFC>my_logFC_threshold,"ASC genes","DOM genes"))) %>%
#   dplyr::mutate(Sig = factor(Sig, levels = c("ASC genes","DOM genes", "N.S."))) %>%
#   dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
#   unique() -> ad_apt
# 
# y1a  %>%
#   filter(symbol %in% enx) %>%
#   dplyr::select(symbol, logFC, P.Value) %>%
#   dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
#                              ifelse(logFC>my_logFC_threshold,"ASC genes","SUB genes"))) %>%
#   dplyr::mutate(Sig = factor(Sig, levels = c("ASC genes","SUB genes", "N.S."))) %>%
#   dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
#   unique() -> as_apt
# 
# 
# 
# des1 <- dd_atp %>% filter(Sig == "DES genes")
# des2 <- ds_atp %>% filter(Sig == "DES genes")
# 
# des1$symbol %in% des2$symbol
# # Cox5a  Uqcc2 Tmsb4x
# des <- des1 %>% rbind(des2) 
# des 
# x <- unique(des$symbol)
# des %>% filter(symbol %in% x) %>% arrange(-logFC)
# 
# 
# 
# dom1 <- dd_atp %>% filter(Sig == "DOM genes")
# sub1 <- ds_atp %>% filter(Sig == "SUB genes")
# 
# 
# 
# 
# asc1 <- ad_apt %>% filter(Sig == "ASC genes") #11
# asc2 <- as_apt %>% filter(Sig == "ASC genes") #15
# 
# asc1$symbol %in% asc2$symbol
# # "Uqcc2",  "Eno1"    "Tefm"    "Atp5g3" , "Ndufb8" 
# asc <- asc1 %>% rbind(asc2) 
# asc 
# x <- unique(asc$symbol) #21 
# asc %>% filter(symbol %in% x) %>% arrange(-logFC)
# 
# 
# dom2 <- ad_apt %>% filter(Sig == "DOM genes") #4
# sub2 <- as_apt %>% filter(Sig == "SUB genes") #1
# 
# #trn dom
# ad_apt$symbol %in% dd_atp$symbol
# # "Uqcc2"  "Ndufaf1"
# #trn sub 
# as_apt$symbol %in% ds_atp$symbol
# # "Chchd10"  "Atp5g3"   "Cox6a1" "Ndufc2"   "Slc25a33" "Ndufa8"   "Ndufb8" "Uqcc2" 
# 
# 
#  # 2 Cytokine Signaling
# msigdbr(species = "Mus musculus", subcategory ="CP:REACTOME") -> reactome_sets
# reactome_sets %>%
# filter(grepl("REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",gs_name,ignore.case = T))  %>%
# select(gene_symbol) %>%  as_tibble() %>%
#  .$gene_symbol -> enx
# 
# 
# y4a  %>%
#   filter(symbol %in% enx) %>%
#   dplyr::select(symbol, logFC, P.Value) %>%
#   dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
#                              ifelse(logFC>my_logFC_threshold,"DES genes","DOM genes"))) %>%
#   dplyr::mutate(Sig = factor(Sig, levels = c("DES genes","DOM genes", "N.S."))) %>%
#   dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
#   unique() -> dd_c
# 
# 
# y3a  %>%
#   filter(symbol %in% enx) %>%
#   dplyr::select(symbol, logFC, P.Value) %>%
#   dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
#                              ifelse(logFC>my_logFC_threshold,"DES genes","SUB genes"))) %>%
#   dplyr::mutate(Sig = factor(Sig, levels = c("DES genes","SUB genes", "N.S."))) %>%
#   dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
#   unique() -> ds_c
# 
# y2a  %>%
#   filter(symbol %in% enx) %>%
#   dplyr::select(symbol, logFC, P.Value) %>%
#   dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
#                              ifelse(logFC>my_logFC_threshold,"ASC genes","DOM genes"))) %>%
#   dplyr::mutate(Sig = factor(Sig, levels = c("ASC genes","DOM genes", "N.S."))) %>%
#   dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
#   unique() -> ad_c
# 
# y1a  %>%
#   filter(symbol %in% enx) %>%
#   dplyr::select(symbol, logFC, P.Value) %>%
#   dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
#                              ifelse(logFC>my_logFC_threshold,"ASC genes","SUB genes"))) %>%
#   dplyr::mutate(Sig = factor(Sig, levels = c("ASC genes","SUB genes", "N.S."))) %>%
#   dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
#   unique() -> as_c
# 
# 
# 
# 
# des1 <- dd_c %>% filter(Sig == "DES genes") # 14 
# des2 <- ds_c %>% filter(Sig == "DES genes") # 29 
# 
# des1$symbol %in% des2$symbol
# # Ccl12"  "Vim"  "Ptgs2"  "Gab2" "Ltbr" "Ptk2b"
# des <- des1 %>% rbind(des2) 
# des 
# x <- unique(des$symbol)
# des %>% filter(symbol %in% x) %>% arrange(-logFC)
# 
# 
# 
# 
# dom1 <- dd_c%>% filter(Sig == "DOM genes") #7
# sub1 <- ds_c %>% filter(Sig == "SUB genes") #13
# dom1$symbol %in% sub1$symbol
# # "Rps6ka2"  "Flnb"    
# 
# asc1 <- ad_c %>% filter(Sig == "ASC genes") #13
# asc2 <- as_c%>% filter(Sig == "ASC genes") #9 
# 
# asc1$symbol %in% asc2$symbol
# # "Aaas" "Ube2m"  "Rnf7" "Arf1"   "Psmb7" 
# asc <- asc1 %>% rbind(asc2) 
# asc 
# x <- unique(asc$symbol) #22
# asc %>% filter(symbol %in% x) %>% arrange(-logFC)
# 
# 
# dom2 <- ad_c %>% filter(Sig == "DOM genes") #6
# sub2 <- as_c %>% filter(Sig == "SUB genes") #11
# 
# dom2$symbol %in% sub2$symbol
# "Birc5"
# 
# #trn dom
# ad_c$symbol %in% dd_c$symbol
# # "Aaas" "Seh1l" "Camk2d" "Ppp2r1b" "Ifnar2"  "Akt2"  
# #trn sub 
# as_c$symbol %in% ds_c$symbol
# # "Pik3r1"   "Flnb" "Rps6ka2""Birc5"



#3. Myelination Regulation 
my <- read_csv("manuscript/brain/gene_sets/myelin.csv")
my %>%  select(gene_symbol) %>%  as_tibble() %>%
  .$gene_symbol -> enx



y4a  %>%
  filter(symbol %in% enx) %>%
  dplyr::select(symbol, logFC, P.Value) %>%
  dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                             ifelse(logFC>my_logFC_threshold,"DES genes","DOM genes"))) %>%
  dplyr::mutate(Sig = factor(Sig, levels = c("DES genes","DOM genes", "N.S."))) %>%
  dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
  unique() -> dd_m


y3a  %>%
  filter(symbol %in% enx) %>%
  dplyr::select(symbol, logFC, P.Value) %>%
  dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                             ifelse(logFC>my_logFC_threshold,"DES genes","SUB genes"))) %>%
  dplyr::mutate(Sig = factor(Sig, levels = c("DES genes","SUB genes", "N.S."))) %>%
  dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
  unique() -> ds_m

y2a  %>%
  filter(symbol %in% enx) %>%
  dplyr::select(symbol, logFC, P.Value) %>%
  dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                             ifelse(logFC>my_logFC_threshold,"ASC genes","DOM genes"))) %>%
  dplyr::mutate(Sig = factor(Sig, levels = c("ASC genes","DOM genes", "N.S."))) %>%
  dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
  unique() -> ad_m

y1a  %>%
  filter(symbol %in% enx) %>%
  dplyr::select(symbol, logFC, P.Value) %>%
  dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                             ifelse(logFC>my_logFC_threshold,"ASC genes","SUB genes"))) %>%
  dplyr::mutate(Sig = factor(Sig, levels = c("ASC genes","SUB genes", "N.S."))) %>%
  dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
  unique() -> as_m




des1 <- dd_m %>% filter(Sig == "DES genes") # 1
des2 <- ds_m %>% filter(Sig == "DES genes") # 18 

des1$symbol %in% des2$symbol
des <- des1 %>% rbind(des2) 
des 
x <- unique(des$symbol)
des %>% filter(symbol %in% x) %>% arrange(-logFC)




dom1 <- dd_m%>% filter(Sig == "DOM genes") #1
sub1 <- ds_m %>% filter(Sig == "SUB genes") #2
dom1$symbol %in% sub1$symbol
    

asc1 <- ad_m %>% filter(Sig == "ASC genes") #1
asc2 <- as_m%>% filter(Sig == "ASC genes") #12

asc1$symbol %in% asc2$symbol
# "Aaas" "Ube2m"  "Rnf7" "Arf1"   "Psmb7" 
asc <- asc1 %>% rbind(asc2) 
asc 
x <- unique(asc$symbol) #22
asc %>% filter(symbol %in% x) %>% arrange(-logFC)


dom2 <- ad_m %>% filter(Sig == "DOM genes") #1
sub2 <- as_m %>% filter(Sig == "SUB genes") #2

dom2$symbol %in% sub2$symbol
"Birc5"

#trn dom
ad_m$symbol %in% dd_m$symbol
#Mtmr2
#trn sub 
 as_m$symbol %in% ds_m$symbol
# [1] "Sox10"   "Pikfyve" "Zfp488" "Ugt8a"    "Mobp"    "Plp1"    "Bcas1"  
# [11] "Cnp"     "Mog"     "Ctsc"    "Mal"   

# unique to ASC
# Itgb4",  "Opalin"  "Tspan2"

 
