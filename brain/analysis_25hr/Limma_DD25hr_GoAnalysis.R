# libraries 
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
organism = 'org.Mm.eg.db'
library(organism, character.only = TRUE)
library(biomaRt)
library(AnnotationDbi)
library(pheatmap)
library(annotables)
library(goseq)
library(clusterProfiler)
library(enrichplot)
library(organism, character.only = TRUE)
library(DOSE)
library(tidyverse)
grcm38 # mouse genes


library(tidyverse)

my_logFC_threshold = 0.2

limma_list<- readRDS("manuscript/brain/results/limma_MeA_25hrDD.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 


x<- limma_list$status 
x1 <- limma_list$cort




limma_list<- readRDS("manuscript/brain/results/limma_mPFC_25hrDD.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 


y<- limma_list$status 
y1 <- limma_list$cort


x %>% filter(grepl("Chr",symbol))
x %>% filter(symbol == "Nr3c2")
x %>% filter(symbol == "Htr5b")
x %>% filter(symbol == "Gal")
x %>% filter(symbol == "Npy")
x %>% filter(symbol == "Agrp")

x %>% filter(symbol == "Chat")
x %>% filter(symbol == "Pomc")



limma_cort_DEG <- limma_list$cort 


#upreg with status 
upstatus <- limma_list$status %>% 
  filter(logFC>0.2)

ups <- upstatus %>% arrange(-logFC) %>% head(25) %>% tibble()


downstatus <- x %>% 
  filter(logFC<0.2) 

downs <- downstatus %>% arrange(logFC) %>% head(25) %>% tibble()


##CORT
upcort <- limma_list$cort %>% 
  filter(logFC>0.2) 
upc <- upcort %>% arrange(-logFC) %>% head(25) %>% tibble()


downcort <- limma_list$cort %>% 
  filter(logFC<0.2) 

dc <- downcort %>% arrange(logFC) %>% head(25) %>% tibble()


my_ont = "BP"
my_showCategory = 100

gettop10GO <- function(limma_df,my_showCategory = 10){
  go_df_up <- limma_df %>% 
    filter(logFC>0) %>% 
    left_join(grcm38 %>% dplyr::select(symbol, entrez)) %>% 
    mutate(entrez = as.character(entrez)) %>% 
    dplyr::select(entrez, logFC) %>% 
    filter(!is.na(entrez)) %>% 
    arrange(desc(logFC))
  
  
  go_df_down <- limma_df %>% 
    filter(logFC<0) %>% 
    left_join(grcm38 %>% dplyr::select(symbol, entrez)) %>% 
    mutate(entrez = as.character(entrez)) %>% 
    dplyr::select(entrez, logFC) %>% 
    filter(!is.na(entrez)) %>% 
    arrange(desc(logFC))
  
  ggo_up <- enrichGO(gene = go_df_up$entrez %>% unique(),
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                     keyType = "ENTREZID",
                     ont = my_ont,
                     readable = T,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.2,
                     qvalueCutoff  = 0.20)
  
  ggo_down <- enrichGO(gene = go_df_down$entrez %>% unique() ,
                       OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = my_ont,
                       readable = T,
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.2, # change to 0.2 for Liver-status only
                       qvalueCutoff  = 0.20) 
  
  fortify(
    ggo_up,
    showCategory = my_showCategory,
    by = "GeneRatio",
    split = NULL,
    includeAll = TRUE
  ) %>% 
    arrange(desc(GeneRatio)) %>% 
    mutate(direction = "Up") -> temp1
  
  
  fortify(
    ggo_down,
    showCategory = my_showCategory,
    by = "GeneRatio",
    split = NULL,
    includeAll = TRUE
  ) %>% 
    arrange(desc(GeneRatio)) %>% 
    mutate(direction = "Down") -> temp2
  return(rbind(temp1,temp2))
  
}


gettop10GO(x, my_showCategory) %>% 
  mutate(comparison = "Social status") -> top10go1

gettop10GO(x1, my_showCategory ) %>% 
  mutate(comparison = "Post Corticosterone") -> top10go2


rbind(top10go1,top10go2) -> top10_GOterms

 write.csv(top10_GOterms,"manuscript/brain/results/RNAseqTable/topBP_GOterms_abovelogFC0.2_MeA_DD25hr.csv", row.names = F)


gettop10GO(y, my_showCategory) %>% 
  mutate(comparison = "Social status") -> top10go1

gettop10GO(y1, my_showCategory ) %>% 
  mutate(comparison = "Post Corticosterone") -> top10go2

rbind(top10go1,top10go2) -> top10_GOterms
write.csv(top10_GOterms,"manuscript/brain/results/RNAseqTable/topBP_GOterms_abovelogFC0.2_mPFC_DD25hr.csv", row.names = F)



top10_GOterms %>% arrange(qvalue) %>% select(Description, direction, comparison) %>%  head(20)
# view

top10_GOterms$Description[1:5]
top10_GOterms$Description[6:10]
top10_GOterms$Description[11:15]
top10_GOterms$Description[16:20]
top10_GOterms$Description[21:25]

limma_list<- readRDS("manuscript/brain/results/limma_mPFC_25hrDD.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 
x<- limma_list$status 

# x %>% filter(symbol == "Efhd2")
# x %>% filter(symbol == "Rnf13")
# x %>% filter(symbol == "Cap1")
# x %>% filter(symbol == "Gins4")
# x %>% filter(symbol == "Fos")
# x %>% filter(symbol == "Sdc4")
# x %>% filter(symbol == "Cdkl1")
# x %>% filter(symbol == "Tmem69")
# x %>% filter(symbol == "Tmem14c")
# x %>% filter(symbol == "Casp9")
# x %>% filter(symbol == "Chat")
# x %>% filter(symbol == "Dusp1")
# x %>% filter(symbol == "Hdac4")
# x %>% filter(symbol == "Ier2")
# x %>% filter(symbol == "Bdnf")
# x %>% filter(symbol == "Btg2")
# x %>% filter(symbol == "Nr4a1")
# x %>% filter(symbol == "Crh")
# x %>% filter(symbol == "Avp")
# x %>% filter(symbol == "Oxt")
# x %>% filter(symbol == "Cyp11a1")
# x %>% filter(symbol == "Star")
x %>% filter(symbol == "Pomc")
 x %>% filter(symbol == "ERK")
