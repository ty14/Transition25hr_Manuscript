#social transitions withcontrols. 

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

my_logFC_threshold = 0.2

#ALL GROUPS 
my_logFC_threshold = 0.2

limma_list<- readRDS("manuscript/brain/results/limma_PFC70min_ReorganizedGroup.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

adom <- limma_list$ascdom

doms <- limma_list$domsub

dd <- limma_list$desdom

da <- limma_list$desasc

as <- limma_list$ascsub

dsub <- limma_list$dessub


#upreg
adom_up <- adom %>% filter(logFC >= 0.2) %>% arrange(-logFC)
doms_up <- doms %>% filter(logFC >= 0.2)%>% arrange(-logFC)
dd_up <- dd %>% filter(logFC >= 0.2)%>% arrange(-logFC) 
da_up <- da %>% filter(logFC >= 0.2)%>% arrange(-logFC)
as_up <- as %>% filter(logFC >= 0.2)%>% arrange(-logFC)
dsub_up <- dsub %>% filter(logFC >= 0.2)%>% arrange(-logFC)

#downreg
adom_down <- adom %>% filter(logFC <= -0.2)%>% arrange(logFC)
doms_down <- doms %>% filter(logFC <= -0.2)%>% arrange(logFC)
dd_down <- dd %>% filter(logFC <= -0.2)%>% arrange(logFC)
da_down <- da %>% filter(logFC <= -0.2)%>% arrange(logFC)
as_down <- as %>% filter(logFC <= -0.2)%>% arrange(logFC)
dsub_down <- dsub %>% filter(logFC <= -0.2)%>% arrange(logFC)




# start with TRN groups with dom
# upsetter Plot 
library(UpSetR)
library(workflowr)
library(ComplexUpset)
listInput <- list(adom_up= adom_up$symbol, dd_up=dd_up$symbol, adom_down = adom_down$symbol, 
                  dd_down = dd_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

#getting overlap genes with for dom and trn groups 
dt_up <- adom_up$symbol[adom_up$symbol %in% dd_up$symbol] %>% as.data.frame()
dt_down %>% head(.,20)
dt_down <- adom_down$symbol[adom_down$symbol %in% dd_down$symbol]%>% as.data.frame()


#significant test
mat <- matrix(c(66,0,0,57),ncol=2)
# Perform chi-squared test
chi_square <- chisq.test(mat)


# subs
listInput <- list(dsub_up= dsub_up$symbol, as_up=as_up$symbol, dsub_down = dsub_down$symbol, 
                  as_down = as_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)


#getting overlap genes 
st_up <- as_up$symbol[as_up$symbol %in% dsub_up$symbol] %>% as.data.frame()
st_down <- as_down$symbol[as_down$symbol %in% dsub_down$symbol]%>% as.data.frame()

#significant test
mat <- matrix(c(138,0,0,134),ncol=2)
# Perform chi-squared test
chi_square <- chisq.test(mat)



# looking for overlap in transition groups between sub and dom comparision 
listInput <- list(st_up= st_up$., dt_up=dt_up$., st_down = st_down$.,dt_down = dt_down$.)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

#getting overlap genes 
t_up <- st_up$.[st_up$. %in% dt_up$.] %>% as.data.frame()
t_down <- st_down$.[st_down$. %in% dt_down$.]%>% as.data.frame()

#significant test
mat <- matrix(c(8,0,0,11),ncol=2)
# Perform chi-squared test
chi_square <- chisq.test(mat)
