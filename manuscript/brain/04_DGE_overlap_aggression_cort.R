
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

#upregulated 
cortr_up <- cort_r %>% filter(logFC >= 0.2) %>% arrange(-logFC)
ar_up <- ar %>% filter(logFC >= 0.2) %>% arrange(-logFC)
car_up <- c_ar %>% filter(logFC >= 0.2) %>% arrange(-logFC)
cortg_up <- cort_g %>% filter(logFC >= 0.2) %>% arrange(-logFC)
ag_up <- ag %>% filter(logFC >= 0.2) %>% arrange(-logFC)
cag_up <- c_ag %>% filter(logFC >= 0.2) %>% arrange(-logFC)


#downregulated 
cortr_down <- cort_r %>% filter(logFC <= 0.2) %>% arrange(logFC)
ar_down <- ar %>% filter(logFC <= 0.2) %>% arrange(logFC)
car_down <- c_ar %>% filter(logFC <= 0.2) %>% arrange(logFC)
cortg_down <- cort_g %>% filter(logFC <= 0.2) %>% arrange(logFC)
ag_down <- ag %>% filter(logFC <= 0.2) %>% arrange(logFC)
cag_down <- c_ag %>% filter(logFC <= 0.2) %>% arrange(logFC)



library(UpSetR)
library(workflowr)
library(ComplexUpset)

#aggression 
listInput <- list(ag_up= ag_up$symbol,ar_up_up=ar_up$symbol,ar_down = ar_down$symbol, 
                  ag_down = ag_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

#getting overlap genes
ag_up.ar_down <- ag_up$symbol[ag_up$symbol %in% ar_down$symbol] %>% as.data.frame()
ar_up.ag_down <- ag_down$symbol[ag_down$symbol %in% ar_up$symbol]%>% as.data.frame()



#agggression and cort interaction 
listInput <- list(cag_up= cortg_up$symbol,car_up_up=cortr_up$symbol,car_down = cortr_down$symbol, 
                  cag_down = cortg_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

#getting overlap genes
cag_up.car_down <- cag_up$symbol[cag_up$symbol %in% car_down$symbol] %>% as.data.frame()
car_up.cag_down <- cag_down$symbol[cag_down$symbol %in% car_up$symbol]%>% as.data.frame()

#just cort 
listInput <- list(cag_up= cag_up$symbol,car_up_up=car_up$symbol,car_down = car_down$symbol, 
                  cag_down = cag_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)


#getting overlap genes 
cort_down <- cortr_down$symbol[cortr_down$symbol %in% cortg_down$symbol] %>% as.data.frame()
cort_up <- cortr_up$symbol[cortr_up$symbol %in% cortg_up$symbol]%>% as.data.frame()


######## On looking at 25 hr data 
#rec at 25hr 
limma_list <- readRDS("manuscript/brain/results/limma_mPFC25_CORT_AGGREC.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= 0.2)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

cort_r <- limma_list$cort

ar<- limma_list$aggrec

c_ar<- limma_list$cort_aggrec


#given 25 hr 
limma_list <- readRDS("manuscript/brain/results/limma_mPFC25_CORT_AGGiven.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= 0.2)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

cort_g <- limma_list$cort

ag<- limma_list$aggiven

c_ag<- limma_list$cort_aggiven

#upregulated 
cortr_up <- cort_r %>% filter(logFC >= 0.2) %>% arrange(-logFC)
ar_up <- ar %>% filter(logFC >= 0.2) %>% arrange(-logFC)
car_up <- c_ar %>% filter(logFC >= 0.2) %>% arrange(-logFC)
cortg_up <- cort_g %>% filter(logFC >= 0.2) %>% arrange(-logFC)
ag_up <- ag %>% filter(logFC >= 0.2) %>% arrange(-logFC)
cag_up <- c_ag %>% filter(logFC >= 0.2) %>% arrange(-logFC)


#downregulated 
cortr_down <- cort_r %>% filter(logFC <= 0.2) %>% arrange(logFC)
ar_down <- ar %>% filter(logFC <= 0.2) %>% arrange(logFC)
car_down <- c_ar %>% filter(logFC <= 0.2) %>% arrange(logFC)
cortg_down <- cort_g %>% filter(logFC <= 0.2) %>% arrange(logFC)
ag_down <- ag %>% filter(logFC <= 0.2) %>% arrange(logFC)
cag_down <- c_ag %>% filter(logFC <= 0.2) %>% arrange(logFC)



#aggression
listInput <- list(ag_up= ag_up$symbol,ar_up_up=ar_up$symbol,ar_down = ar_down$symbol, 
                  ag_down = ag_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

#getting overlap genes
ag_up.ar_down25 <- ag_up$symbol[ag_up$symbol %in% ar_down$symbol] %>% as.data.frame()
ar_up.ag_down25 <- ag_down$symbol[ag_down$symbol %in% ar_up$symbol]%>% as.data.frame()



#aggresion and cort interaction 
listInput <- list(cag_up= cag_up$symbol,car_up_up=car_up$symbol,car_down = car_down$symbol, 
                  cag_down = cag_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

#getting overlap genes
cag_up.car_down25 <- cag_up$symbol[cag_up$symbol %in% car_down$symbol] %>% as.data.frame()
car_up.cag_down25 <- cag_down$symbol[cag_down$symbol %in% car_up$symbol]%>% as.data.frame()


#just cort 
listInput <- list(cag_up= cortg_up$symbol,car_up=cortr_up$symbol,car_down = cortr_down$symbol, 
                  cag_down = cortg_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

#opposite directions. use RROH analysis. 






