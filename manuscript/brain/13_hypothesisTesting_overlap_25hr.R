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
limma_list<- readRDS("manuscript/brain/results_use/limma_PFC25hr_NORM_RG.RDS") %>%
  map(~distinct(.)) %>%
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>%
  map(~filter(.,!is.na(entrez)))
y4a <- limma_list$desdom
y3a <- limma_list$dessub
y2a <- limma_list$ascdom
y1a <- limma_list$ascsub


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
#   unique() -> dd_atp25
# 
# 
# y3a  %>%
#   filter(symbol %in% enx) %>%
#   dplyr::select(symbol, logFC, P.Value) %>%
#   dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
#                              ifelse(logFC>my_logFC_threshold,"DES genes","SUB genes"))) %>%
#   dplyr::mutate(Sig = factor(Sig, levels = c("DES genes","SUB genes", "N.S."))) %>%
#   dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
#   unique() -> ds_atp25
# 
# y2a  %>%
#   filter(symbol %in% enx) %>%
#   dplyr::select(symbol, logFC, P.Value) %>%
#   dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
#                              ifelse(logFC>my_logFC_threshold,"ASC genes","DOM genes"))) %>%
#   dplyr::mutate(Sig = factor(Sig, levels = c("ASC genes","DOM genes", "N.S."))) %>%
#   dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
#   unique() -> ad_apt25
# 
# y1a  %>%
#   filter(symbol %in% enx) %>%
#   dplyr::select(symbol, logFC, P.Value) %>%
#   dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
#                              ifelse(logFC>my_logFC_threshold,"ASC genes","SUB genes"))) %>%
#   dplyr::mutate(Sig = factor(Sig, levels = c("ASC genes","SUB genes", "N.S."))) %>%
#   dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
#   unique() -> as_apt25
# 
# 
# 
# des125 <- dd_atp25 %>% filter(Sig == "DES genes")
# des225 <- ds_atp25 %>% filter(Sig == "DES genes")
# 
# des125$symbol %in% des225$symbol
# # Cox5a  Uqcc2 Tmsb4x
# des25 <- des125 %>% rbind(des225) 
# des25 
# x25 <- unique(des25$symbol)
# des25 %>% filter(symbol %in% x25) %>% arrange(-logFC)
# des$symbol %in% des25$symbol
# 
# 
# dom125 <- dd_atp25 %>% filter(Sig == "DOM genes")
# sub125 <- ds_atp25 %>% filter(Sig == "SUB genes")
# 
# 
# asc125 <- ad_apt25 %>% filter(Sig == "ASC genes") #11
# asc225 <- as_apt25 %>% filter(Sig == "ASC genes") #15
# 
# asc125$symbol %in% asc225$symbol
# # "Uqcc2",  "Eno1"    "Tefm"    "Atp5g3" , "Ndufb8" 
# asc25 <- asc125 %>% rbind(asc225) 
# asc25 
# x25 <- unique(asc25$symbol) #21 
# asc25 %>% filter(symbol %in% x25) %>% arrange(-logFC)
# asc$symbol %in% asc25$symbol
# # "Sirt6" "Tefm"   
# dom225 <- ad_apt25 %>% filter(Sig == "DOM genes") #4
# sub225 <- as_apt25 %>% filter(Sig == "SUB genes") #1
# 
# 
# #trn dom
# ad_apt25$symbol %in% dd_atp25$symbol
# #no overlap
# #trn sub 
# as_apt25$symbol %in% ds_atp25$symbol
#  #no overlap 
# 
# # 2 Cytokine Signaling
# msigdbr(species = "Mus musculus", subcategory ="CP:REACTOME") -> reactome_sets
# reactome_sets %>%
#   filter(grepl("REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",gs_name,ignore.case = T))  %>%
#   select(gene_symbol) %>%  as_tibble() %>%
#   .$gene_symbol -> enx
# 
# 
# y4a  %>%
#   filter(symbol %in% enx) %>%
#   dplyr::select(symbol, logFC, P.Value) %>%
#   dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
#                              ifelse(logFC>my_logFC_threshold,"DES genes","DOM genes"))) %>%
#   dplyr::mutate(Sig = factor(Sig, levels = c("DES genes","DOM genes", "N.S."))) %>%
#   dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
#   unique() -> dd_c25
# 
# 
# y3a  %>%
#   filter(symbol %in% enx) %>%
#   dplyr::select(symbol, logFC, P.Value) %>%
#   dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
#                              ifelse(logFC>my_logFC_threshold,"DES genes","SUB genes"))) %>%
#   dplyr::mutate(Sig = factor(Sig, levels = c("DES genes","SUB genes", "N.S."))) %>%
#   dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
#   unique() -> ds_c25
# 
# y2a  %>%
#   filter(symbol %in% enx) %>%
#   dplyr::select(symbol, logFC, P.Value) %>%
#   dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
#                              ifelse(logFC>my_logFC_threshold,"ASC genes","DOM genes"))) %>%
#   dplyr::mutate(Sig = factor(Sig, levels = c("ASC genes","DOM genes", "N.S."))) %>%
#   dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
#   unique() -> ad_c25
# 
# y1a  %>%
#   filter(symbol %in% enx) %>%
#   dplyr::select(symbol, logFC, P.Value) %>%
#   dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
#                              ifelse(logFC>my_logFC_threshold,"ASC genes","SUB genes"))) %>%
#   dplyr::mutate(Sig = factor(Sig, levels = c("ASC genes","SUB genes", "N.S."))) %>%
#   dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
#   unique() -> as_c25
# 
# 
# 
# 
# des125 <- dd_c25 %>% filter(Sig == "DES genes") # 6
# des225 <- ds_c25 %>% filter(Sig == "DES genes") # 17 
# 
# des125$symbol %in% des225$symbol
# # Adam17"    "Nfkb1"    "Lif"   
# des25 <- des125 %>% rbind(des225) 
# des25 
# x25 <- unique(des25$symbol)
# des25 %>% filter(symbol %in% x25) %>% arrange(-logFC)
# 
# des$symbol %in% des25$symbol
# 
# 
# dom125 <- dd_c25%>% filter(Sig == "DOM genes") #4
# sub125 <- ds_c25 %>% filter(Sig == "SUB genes") #11
# 
# dom125$symbol %in% sub125$symbol
# # "Foxo3" 
# 
# asc125 <- ad_c25 %>% filter(Sig == "ASC genes") #9
# asc225 <- as_c25%>% filter(Sig == "ASC genes") #12 
# 
# asc125$symbol %in% asc225$symbol
# # "Nfkb1""Vrk3" "Il5ra"
# asc25 <- asc125 %>% rbind(asc225) 
# asc25 
# x25 <- unique(asc25$symbol) #21 
# asc25 %>% filter(symbol %in% x25) %>% arrange(-logFC)
# asc$symbol %in% asc25$symbol
# 
# dom225 <- ad_c25 %>% filter(Sig == "DOM genes") #2
# sub225 <- as_c25 %>% filter(Sig == "SUB genes") #7
# 
# dom225$symbol %in% sub225$symbol
# 
# #trn dom
# ad_c25$symbol %in% dd_c25$symbol
# # "Nfkb1" "Mapkapk2"
# #trn sub 
# as_c25$symbol %in% ds_c25$symbol
# # "Vrk3" "Rapgef1""Nfkb1""Kpna4"    "Trim12c "Il18r1" 
# 



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
  unique() -> dd_m25


y3a  %>%
  filter(symbol %in% enx) %>%
  dplyr::select(symbol, logFC, P.Value) %>%
  dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                             ifelse(logFC>my_logFC_threshold,"DES genes","SUB genes"))) %>%
  dplyr::mutate(Sig = factor(Sig, levels = c("DES genes","SUB genes", "N.S."))) %>%
  dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
  unique() -> ds_m25

y2a  %>%
  filter(symbol %in% enx) %>%
  dplyr::select(symbol, logFC, P.Value) %>%
  dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                             ifelse(logFC>my_logFC_threshold,"ASC genes","DOM genes"))) %>%
  dplyr::mutate(Sig = factor(Sig, levels = c("ASC genes","DOM genes", "N.S."))) %>%
  dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
  unique() -> ad_m25

y1a  %>%
  filter(symbol %in% enx) %>%
  dplyr::select(symbol, logFC, P.Value) %>%
  dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                             ifelse(logFC>my_logFC_threshold,"ASC genes","SUB genes"))) %>%
  dplyr::mutate(Sig = factor(Sig, levels = c("ASC genes","SUB genes", "N.S."))) %>%
  dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
  unique() -> as_m25




des125 <- dd_m25 %>% filter(Sig == "DES genes") # 1
des225 <- ds_m25 %>% filter(Sig == "DES genes") # 5 

des125$symbol %in% des225$symbol
# Rxrg
des25 <- des125 %>% rbind(des225) 
des25 
x25 <- unique(des25$symbol)
des25 %>% filter(symbol %in% x25) %>% arrange(-logFC)

des$symbol %in% des25$symbol
# "Ugt8a" "Mobp"

dom125 <- dd_m25%>% filter(Sig == "DOM genes") #1
sub125 <- ds_m25 %>% filter(Sig == "SUB genes") #0

asc125 <- ad_m25 %>% filter(Sig == "ASC genes") #0
asc225 <- as_m25%>% filter(Sig == "ASC genes") #12 
asc2$symbol %in% asc225$symbol
asc125$symbol %in% asc225$symbol
# "Nfkb1""Vrk3" "Il5ra"
asc25 <- asc125 %>% rbind(asc225) 
asc25 
x25 <- unique(asc25$symbol) #21 
asc25 %>% filter(symbol %in% x25) %>% arrange(-logFC)
asc$symbol %in% asc25$symbol
# Ugt8a, Mobp

dom225 <- ad_m25 %>% filter(Sig == "DOM genes") #0
sub225 <- as_m25 %>% filter(Sig == "SUB genes") #0
dom225$symbol %in% sub225$symbol

#trn dom
ad_m25$symbol %in% dd_m25$symbol

#trn sub 
as_m25$symbol %in% ds_m25$symbol
# "Cntn2" 
# "Mobp" 
# "Ugt8a"