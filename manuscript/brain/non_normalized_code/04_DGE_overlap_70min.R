
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

limma_list<- readRDS("manuscript/brain/results/limma_PFC70min_ReorganizedGroup.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

adom <- limma_list$ascdom #666

dd <- limma_list$desdom #763


as <- limma_list$ascsub #807

dsub <- limma_list$dessub #1580



xx <- adom %>% rbind(dd,as,dsub)

all <- unique(xx$entrez)

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


### Number of genes in each cut off 
#0.2
adom %>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #309
adom %>% filter(between(logFC, .2, .5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 276
adom %>% filter(between(logFC, .5, 1))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #36 
adom %>% filter(between(logFC, 1, 3))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #6

adom %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 322 
adom %>% filter(between(logFC, -0.5, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #248
adom %>% filter(between(logFC, -1, -0.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 62
adom %>% filter(between(logFC, -3, -1)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 12



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



colnames(dt_up)[1]<- "symbol"
dt_up <- dt_up %>% mutate(regulated = "up") %>% mutate(condition = "TRN vs DOM") %>% mutate(time = 70)

colnames(dt_down)[1]<- "symbol"
dt_down <- dt_down %>% mutate(regulated = "down") %>% mutate(condition = "TRN vs DOM") %>% mutate(time = 70)

dt <- dt_up %>%  rbind(dt_down)

df_dd <- dd %>% select(symbol, dd_logFC= logFC, dd_pv= P.Value) %>%
  mutate(time = 70) %>%  filter(symbol %in% dt$symbol)
df_adom <- adom %>% select(symbol, adom_logFC= logFC, adom_pv= P.Value) %>%
  mutate(time = 70) %>% filter(symbol %in% dt$symbol)

dd_dt <- df_dd %>% full_join(df_adom)

dd_dt70 <- dd_dt %>% full_join(dt)



colnames(st_up)[1]<- "symbol"
st_up <- st_up %>% mutate(regulated = "up") %>% mutate(condition = "TRN vs SUB") %>% mutate(time = 70)

colnames(st_down)[1]<- "symbol"
st_down <- st_down %>% mutate(regulated = "down") %>% mutate(condition = "TRN vs SUB") %>% mutate(time = 70)

st <- st_up %>%  rbind(st_down)

df_as <- as %>% select(symbol, as_logFC= logFC, as_pv= P.Value) %>%
  mutate(time = 70) %>%  filter(symbol %in% st$symbol)
df_ds<- dsub %>% select(symbol, ds_logFC= logFC, ds_pv= P.Value) %>%
  mutate(time = 70) %>% filter(symbol %in% st$symbol)

dd_st <- df_as %>% full_join(df_ds)

dd_st70 <- dd_st %>% full_join(st)
