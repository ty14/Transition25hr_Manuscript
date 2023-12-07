
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
library(RRHO2)
library(tidyverse)
grcm38 # mouse genes


#ALL GROUPS 
my_logFC_threshold = 0.2

limma_list<- readRDS("manuscript/brain/results/limma_PFC70min_ReorganizedGroup.RDS") %>% 
  map(~distinct(.)) %>% 
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

adom <- limma_list$ascdom

dd <- limma_list$desdom

as <- limma_list$ascsub

dsub <- limma_list$dessub


sta <- adom %>% rbind(as)

std <- dd %>% rbind(dsub)

# outlier removed 

limma_list<- readRDS("manuscript/brain/results/limma_PFC_ReorganizedGroups_outlierremoved.RDS") %>% 
  map(~distinct(.)) %>% 
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

adomx <- limma_list$ascdom

ddx <- limma_list$desdom

asx <- limma_list$ascsub

dsubx <- limma_list$dessub


stax <- adomx %>% rbind(asx)

stdx <- ddx %>% rbind(dsubx)

#ASC vs DOM both timepoints 

l <- list(stax,stdx)
names(l) <- c("adom_70", "dd_70")
lapply(l, head)


all <- l %>% map2_df(.,names(.), ~mutate(.x, group = .y)) 

head(all)
tail(all)

# adom <- adom[adom$symbol %in% adomx$symbol,]
# adomx <- adomx[adomx$symbol %in% adom$symbol,]

# adom <- adom %>% filter(symbol != "")
# adomx <- adomx%>% filter(symbol != "")

aDDDS <- stax %>% mutate(gene = row_number())
pDDDS <- stdx %>% mutate(gene = row_number())

aDDDS

pDDDS$gene <- paste0("Gene", pDDDS$gene)
head(pDDDS)
aDDDS$gene <- paste0("Gene", aDDDS$gene)
head(aDDDS)
#List 1
Gene = aDDDS$gene

list1_pvalue_1_200 <- aDDDS %>% filter(logFC > 0 & P.Value < 0.01) 
list1_pvalue_1_200 <-list1_pvalue_1_200$P.Value## up-regulated genes
list1_pvalue_201_400 <- aDDDS %>% filter(logFC < 0 & P.Value < 0.01)
list1_pvalue_201_400 <-list1_pvalue_201_400$P.Value ## down-regulated genes
an <- subset(aDDDS, !(P.Value %in% list1_pvalue_1_200))

an <- subset(an, !(P.Value %in% list1_pvalue_201_400))


list1_pvalue_401_2000 <- an$P.Value  ## non-changed genes
list1_DDE <- c(-log10(list1_pvalue_1_200), 
               -log10(list1_pvalue_201_400) * (-1), 
               -log10(list1_pvalue_401_2000) * sample(c(1,-1), length(list1_pvalue_401_2000), replace = TRUE))

lapply(list1_DDE, head)

gene_list1 <- data.frame(Genes=Gene,DDE = list1_DDE, stringsAsFactors = FALSE)

#List2

list2_pvalue_1_200 <- pDDDS %>% filter(logFC > 0 & P.Value < 0.01) 
list2_pvalue_1_200 <-list2_pvalue_1_200$P.Value## up-regulated genes
list2_pvalue_201_400 <- pDDDS %>% filter(logFC < 0 & P.Value < 0.01)
list2_pvalue_201_400 <-list2_pvalue_201_400$P.Value ## down-regulated genes
pn <- subset(pDDDS, !(P.Value %in% list2_pvalue_1_200))
pn <- subset(pn, !(P.Value %in% list2_pvalue_201_400))
list2_pvalue_401_2000 <- pn$P.Value  ## non-changed genes
list2_DDE <- c(-log10(list2_pvalue_1_200), 
               -log10(list2_pvalue_201_400) * (-1), 
               -log10(list2_pvalue_401_2000) * sample(c(1,-1), length(list2_pvalue_401_2000), replace = TRUE))


gene_list2 <- data.frame(Genes=Gene,DDE = list2_DDE, stringsAsFactors = FALSE)

library(RRHO2)
RRHO_obj <-  RRHO2_initialize(gene_list1, gene_list2, labels = c("ASC vs Stable", "DES vs Stable"), log10.ind=TRUE)


RRHO2_heatmap(RRHO_obj)

# 
# dd: down regulation in list1 and down regulation in list 2
# uu: up regulation in list1 and up regulation in list 2
# du: down regulation in list1 and up regulation in list 2
# ud: up regulation in list1 and down regulation in list 2
RRHO2_vennDiagram(RRHO_obj, type="du")
