
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
library(ComplexUpset)
library(UpSetR)
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

adom <- limma_list$ascdom

dd <- limma_list$desdom

as <- limma_list$ascsub

dsub <- limma_list$dessub

# outlier removed 

limma_list<- readRDS("manuscript/brain/results_use/limma_PFC25hr_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

adomx <- limma_list$ascdom

ddx <- limma_list$desdom

asx <- limma_list$ascsub

dsubx <- limma_list$dessub


#upreg 70min
adom_up <- adom %>% filter(logFC >= 0.2) %>% arrange(-logFC)
dd_up <- dd %>% filter(logFC >= 0.2)%>% arrange(-logFC) 
as_up <- as %>% filter(logFC >= 0.2)%>% arrange(-logFC)
dsub_up <- dsub %>% filter(logFC >= 0.2)%>% arrange(-logFC)


#upreg 25hr
adom_upx <- adomx %>% filter(logFC >= 0.2) %>% arrange(-logFC)
dd_upx <- ddx %>% filter(logFC >= 0.2)%>% arrange(-logFC) 
as_upx <- asx %>% filter(logFC >= 0.2)%>% arrange(-logFC)
dsub_upx <- dsubx %>% filter(logFC >= 0.2)%>% arrange(-logFC)


#downreg 70 min
adom_down <- adom %>% filter(logFC <= -0.2)%>% arrange(logFC)
dd_down <- dd %>% filter(logFC <= -0.2)%>% arrange(logFC)
as_down <- as %>% filter(logFC <= -0.2)%>% arrange(logFC)
dsub_down <- dsub %>% filter(logFC <= -0.2)%>% arrange(logFC)


#downreg 25hr 
adom_downx <- adomx %>% filter(logFC <= -0.2)%>% arrange(logFC)
dd_downx <- ddx %>% filter(logFC <= -0.2)%>% arrange(logFC)
as_downx <- asx %>% filter(logFC <= -0.2)%>% arrange(logFC)
dsub_downx <- dsubx %>% filter(logFC <= -0.2)%>% arrange(logFC)


# compare across time periods 
# ASC vs DOM 
listInput <- list(adom_up= adom_up$symbol, adom_up25=adom_upx$symbol, adom_down = adom_down$symbol, 
                  adom_down25 = adom_downx$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)


#getting overlap genes 
adt_up<- adom_up$symbol[adom_up$symbol %in% adom_upx$symbol] %>% as.data.frame()
adt_down <- adom_down$symbol[adom_down$symbol %in% adom_downx$symbol]%>% as.data.frame()

adtx_up<- adom_up$symbol[adom_up$symbol %in% adom_downx$symbol] %>% as.data.frame()
adtx_down <- adom_down$symbol[adom_down$symbol %in% adom_upx$symbol]%>% as.data.frame()

#significant test
mat <- matrix(c(9,8,4,4),ncol=2)
# Perform chi-squared test
chi_square <- chisq.test(mat)
chi_square

# DES vs DOM 
listInput <- list(dd_up= dd_up$symbol, dd_up25=dd_upx$symbol, dd_down = dd_down$symbol, 
                  dd_down25 = dd_downx$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)


#getting overlap genes 
ddt_up<- dd_up$symbol[dd_up$symbol %in% dd_upx$symbol] %>% as.data.frame()
ddt_down <- dd_down$symbol[dd_down$symbol %in% dd_downx$symbol]%>% as.data.frame()

ddtx_up<- dd_up$symbol[dd_up$symbol %in% dd_downx$symbol] %>% as.data.frame()
ddtx_down <- dd_down$symbol[dd_down$symbol %in% dd_upx$symbol]%>% as.data.frame()

#significant test
mat <- matrix(c(11,7,5,3),ncol=2)
# Perform chi-squared test
chi_square <- chisq.test(mat)
chi_square

# ASC vs SUB 
listInput <- list(as_up= as_up$symbol, as_up25=as_upx$symbol, as_down = as_down$symbol, 
                  as_down25 = as_downx$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)
#significant test
mat <- matrix(c(16,12,10,8),ncol=2)
# Perform chi-squared test
chi_square <- chisq.test(mat)
chi_square


#getting overlap genes 
ast_up<- as_up$symbol[as_up$symbol %in% as_upx$symbol] %>% as.data.frame()
ast_down <- as_down$symbol[as_down$symbol %in% as_downx$symbol]%>% as.data.frame()

astx_up<- as_up$symbol[as_up$symbol %in% as_downx$symbol] %>% as.data.frame()
astx_down <- as_down$symbol[as_down$symbol %in% as_upx$symbol]%>% as.data.frame()



# DES vs SUB 
listInput <- list(dsub_up= dsub_up$symbol, dsub_up25=dsub_upx$symbol, dsub_down = dsub_down$symbol, 
                  dsub_down25 = dsub_downx$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)


#getting overlap genes 
dst_up<- dsub_up$symbol[dsub_up$symbol %in% dsub_upx$symbol] %>% as.data.frame()
dst_down <- dsub_down$symbol[dsub_down$symbol %in% dsub_downx$symbol]%>% as.data.frame()

dstx_up<- dsub_up$symbol[dsub_up$symbol %in% dsub_downx$symbol] %>% as.data.frame()
dstx_down <- dsub_down$symbol[dsub_down$symbol %in% dsub_upx$symbol]%>% as.data.frame()

#significant test
mat <- matrix(c(63,58,18,11),ncol=2)
# Perform chi-squared test
chi_square <- chisq.test(mat)
chi_square

#now descenders in dom and sub at 70 min and 25 hr 
# 
listInput <- list(dsub_up= dsub_up$symbol, dsub_up25=dsub_upx$symbol, dsub_down = dsub_down$symbol, 
                  dsub_down25 = dsub_downx$symbol,dd_up= dd_up$symbol, dd_up25=dd_upx$symbol, dd_down = dd_down$symbol, 
                  dd_down25 = dd_downx$symbol)

UpSetR::upset(fromList(listInput), nsets = 8, order.by = "freq", keep.order = F)





#descender 70 min
listInput <- list(dsub_up= dsub_up$symbol, dsub_down = dsub_down$symbol, 
                  dd_up= dd_up$symbol,  dd_down = dd_down$symbol)

UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

# des 25 hr
listInput <- list( dsub_up25=dsub_upx$symbol, 
                  dsub_down25 = dsub_downx$symbol, dd_up25=dd_upx$symbol,
                  dd_down25 = dd_downx$symbol)

UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

d_up<- dsub_up$symbol[dsub_up$symbol %in% dd_up$symbol] %>% as.data.frame()
d_down <- dsub_down$symbol[dsub_down$symbol %in% dd_down$symbol]%>% as.data.frame()

d_upx<- dsub_upx$symbol[dsub_upx$symbol %in% dd_upx$symbol] %>% as.data.frame()
d_downx <- dsub_downx$symbol[dsub_downx$symbol %in% dd_downx$symbol]%>% as.data.frame()

d_up$.[d_up$. %in% d_upx$.]
d_down$.[d_down$. %in% d_downx$.]

#now ascenders in dom and sub at 70 min and 25 hr 
# 
listInput <- list(as_up= as_up$symbol, as_up25=as_upx$symbol, as_down = as_down$symbol, 
                  as_down25 = as_downx$symbol,adom_up= adom_up$symbol, adom_up25=dd_upx$symbol, 
                  adom_down = adom_down$symbol, adom_down25 = adom_downx$symbol)

UpSetR::upset(fromList(listInput), nsets = 8, order.by = "freq", keep.order = F)


#just 70min
listInput <- list(as_up= as_up$symbol,as_down = as_down$symbol, 
                  adom_up= adom_up$symbol, 
                  adom_down = adom_down$symbol)

UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)


#25 hr
listInput <- list( as_up25=as_upx$symbol,
                  as_down25 = as_downx$symbol, adom_up25=dd_upx$symbol, 
                   adom_down25 = adom_downx$symbol)

 UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

 #one gene down in sub and up in dom 70min 
 
 adom_up$symbol[adom_up$symbol %in% as_down$symbol]%>% as.data.frame()
 #one gene up in sub and down in dom 25hr
 as_upx$symbol[as_upx$symbol %in% adom_downx$symbol]%>% as.data.frame()
 
 a_up<- as_up$symbol[as_up$symbol %in% adom_up$symbol] %>% as.data.frame()
 a_down <-as_down$symbol[as_down$symbol %in% adom_down$symbol]%>% as.data.frame()
 
 a_upx<- as_upx$symbol[dsub_upx$symbol %in% adom_upx$symbol] %>% as.data.frame()
 a_downx <- as_downx$symbol[as_downx$symbol %in% adom_downx$symbol]%>% as.data.frame()
 
 a_up$.[a_up$. %in% a_upx$.]
 a_down$.[a_down$. %in% a_downx$.]
 

 
 colnames(a_up)[1]<- "symbol"
a_up <- a_up %>% mutate(regulated = "up") %>% mutate(condition = "ascenders") %>% mutate(time = 70)
 
colnames(a_down)[1]<- "symbol"
a_down <- a_down %>% mutate(regulated = "down") %>% mutate(condition = "ascenders") %>% mutate(time = 70)

colnames(a_upx)[1]<- "symbol"
a_upx <- a_upx %>% mutate(regulated = "up") %>% mutate(condition = "ascenders") %>% mutate(time = 25)

colnames(a_downx)[1]<- "symbol"
a_downx <- a_downx %>% mutate(regulated = "down") %>% mutate(condition = "ascenders") %>% mutate(time = 25)

df_a <- as %>% select(symbol, as_logFC= logFC, as_pv= P.Value) %>%  mutate(time = 70)
df_ad <- adom %>% select(symbol, adom_logFC= logFC, adom_pv= P.Value) %>%  mutate(time = 70)
df_ax <- asx %>% select(symbol, as_logFC= logFC, as_pv= P.Value) %>%  mutate(time = 25)
df_axd <- adomx %>% select(symbol, adom_logFC= logFC, adom_pv= P.Value) %>%  mutate(time = 25)

a70 <- a_up %>% rbind(a_down)
df_a <- as %>% select(symbol, as_logFC= logFC, as_pv= P.Value) %>%  mutate(time = 70) %>% filter(symbol %in% a70$symbol)
df_ad <- adom %>% select(symbol, adom_logFC= logFC, adom_pv= P.Value) %>%  mutate(time = 70) %>% filter(symbol %in% a70$symbol)

asc_genes70 <- df_a %>% full_join(df_ad)


a25 <- a_upx %>% rbind(a_downx)
df_ax <- asx %>% select(symbol, as_logFC= logFC, as_pv= P.Value) %>%  mutate(time = 25) %>% filter(symbol %in% a25$symbol)
df_axd <- adomx %>% select(symbol, adom_logFC= logFC, adom_pv= P.Value) %>%  mutate(time = 25) %>% filter(symbol %in% a25$symbol)

asc_genes25 <- df_ax %>% full_join(df_axd)


asc_genes <- asc_genes70 %>% rbind(asc_genes25) %>% unique(.)
write.csv(asc_genes, "manuscript/brain/results_tables/asc_genes_mPFC.csv")




colnames(d_up)[1]<- "symbol"
d_up <- d_up %>% mutate(regulated = "up") %>% mutate(condition = "descenders") %>% mutate(time = 70)

colnames(d_down)[1]<- "symbol"
d_down <- a_down %>% mutate(regulated = "down") %>% mutate(condition = "descenders") %>% mutate(time = 70)

colnames(d_upx)[1]<- "symbol"
d_upx <- d_upx %>% mutate(regulated = "up") %>% mutate(condition = "descenders") %>% mutate(time = 25)

colnames(d_downx)[1]<- "symbol"
d_downx <- d_downx %>% mutate(regulated = "down") %>% mutate(condition = "descenders") %>% mutate(time = 25)

df_ds <- dsub %>% select(symbol, dsub_logFC= logFC, dsub_pv= P.Value) %>%  mutate(time = 70)
df_dd <- dd %>% select(symbol, dd_logFC= logFC, dd_pv= P.Value) %>%  mutate(time = 70)
df_dsx <- dsubx %>% select(symbol, dsub_logFC= logFC, dsub_pv= P.Value) %>%  mutate(time = 25)
df_ddx <- ddx %>% select(symbol, dd_logFC= logFC, dd_pv= P.Value) %>%  mutate(time = 25)

d70 <- d_up %>% rbind(d_down)
df_ds <- dsub %>% select(symbol, dsub_logFC= logFC, dsub_pv= P.Value) %>%  mutate(time = 70) %>% filter(symbol %in% d70$symbol)
df_dd <- dd %>% select(symbol, dd_logFC= logFC, dd_pv= P.Value) %>%  mutate(time = 70) %>% filter(symbol %in% d70$symbol)

des_genes70 <- df_ds %>% full_join(df_dd)


d25 <- d_upx %>% rbind(d_downx)
df_dsx <- dsubx %>% select(symbol, dsub_logFC= logFC, dsub_pv= P.Value) %>%  mutate(time = 25)%>% filter(symbol %in% d25$symbol)
df_ddx <- ddx %>% select(symbol, dd_logFC= logFC, dd_pv= P.Value) %>%  mutate(time = 25) %>% filter(symbol %in% d25$symbol)

des_genes25 <- df_dsx %>% full_join(df_ddx)


des_genes <- des_genes70 %>% rbind(des_genes25) %>% unique(.)
write.csv(des_genes, "manuscript/brain/results_tables/des_genes_mPFC.csv")


