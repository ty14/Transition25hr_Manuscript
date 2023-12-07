
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


limma_list<- readRDS("manuscript/brain/results/limma_PFC70_TRN.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

dt <- limma_list$tdom 
ts <- limma_list$tsub

#upreg
dt_up <- dt %>% filter(logFC >= 0.2) %>% arrange(-logFC)
ds_up <- ts %>% filter(logFC >= 0.2)%>% arrange(-logFC)

#downreg
dt_down <- dt %>% filter(logFC <= -0.2)%>% arrange(logFC)
ds_down <- ts %>% filter(logFC <= -0.2)%>% arrange(logFC)



# start with TRN groups at 70 min 
# upsetter Plot 
library(UpSetR)
library(workflowr)
library(ComplexUpset)
listInput <- list(dt_up= dt_up$symbol, ds_up=ds_up$symbol, ds_down = ds_down$symbol, 
                  dt_down = dt_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

#getting overlap genes with for dom and trn groups 
dt_upx <- dt_up$symbol[dt_up$symbol %in% ds_up$symbol] %>% as.data.frame()#69
dt_downx <- dt_down$symbol[dt_down$symbol %in% ds_down$symbol]%>% as.data.frame()#72


#significant test
mat <- matrix(c(69,0,0,72),ncol=2)
# Perform chi-squared test
chi_square <- chisq.test(mat)
chi_square

## at 25 hours 

#ALL GROUPS 
my_logFC_threshold = 0.2

limma_list<- readRDS("manuscript/brain/results/limma_PFC_TRN_outlierRemoved.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

dt <- limma_list$tdom 
ts <- limma_list$tsub

#upreg
dt25_up <- dt %>% filter(logFC >= 0.2) %>% arrange(-logFC)
ds25_up <- ts %>% filter(logFC >= 0.2)%>% arrange(-logFC)

#downreg
dt25_down <- dt %>% filter(logFC <= -0.2)%>% arrange(logFC)
ds25_down <- ts %>% filter(logFC <= -0.2)%>% arrange(logFC)


library(UpSetR)
library(workflowr)
library(ComplexUpset)
listInput <- list(dt_up= dt25_up$symbol, ds_up=ds25_up$symbol, ds_down = ds25_down$symbol, 
                  dt_down = dt25_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

#getting overlap genes with for dom and trn groups 
dt25_upx <- ds25_up$symbol[ds25_up$symbol %in% dt25_up$symbol] %>% as.data.frame()#46
dt25_downx <- dt25_down$symbol[dt25_down$symbol %in% ds25_down$symbol]%>% as.data.frame()#47


#significant test
mat <- matrix(c(40,0,0,46),ncol=2)
# Perform chi-squared test
chi_square <- chisq.test(mat)
chi_square


listInput <- list(dt25_up= dt25_upx$., dt_up=dt_upx$., ds25_down = dt25_downx$., 
                  ds_down = dt_downx$.)

