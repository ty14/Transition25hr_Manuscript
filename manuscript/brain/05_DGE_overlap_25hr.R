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


#ALL GROUPS 
my_logFC_threshold = 0.2

# outlier removed 

limma_list<- readRDS("manuscript/brain/results_use/limma_PFC25hr_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

adom <- limma_list$ascdom

dd <- limma_list$desdom

as <- limma_list$ascsub

dsub <- limma_list$dessub


x <- adom %>% rbind(dd, as, dsub)
xy <- unique(x$symbol)
#upreg
adom_up <- adom %>% filter(logFC >= 0.2) %>% arrange(-logFC)
dd_up <- dd %>% filter(logFC >= 0.2)%>% arrange(-logFC) 
as_up <- as %>% filter(logFC >= 0.2)%>% arrange(-logFC)
dsub_up <- dsub %>% filter(logFC >= 0.2)%>% arrange(-logFC)


#downreg
adom_down <- adom %>% filter(logFC <= -0.2)%>% arrange(logFC)
dd_down <- dd %>% filter(logFC <= -0.2)%>% arrange(logFC)
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
dt_up25 <- adom_up$symbol[adom_up$symbol %in% dd_up$symbol] %>% as.data.frame() %>% unique(.)
dt_down25 <- adom_down$symbol[adom_down$symbol %in% dd_down$symbol]%>% as.data.frame()


#significant test
mat <- matrix(c(29,0,0,32),ncol=2)
# Perform chi-squared test
chi_square <- chisq.test(mat)

# data:  mat
# X-squared = 57.056, df = 1, p-value = 4.235e-14

# subs
listInput <- list(dsub_up= dsub_up$symbol, as_up=as_up$symbol, dsub_down = dsub_down$symbol, 
                  as_down = as_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)


#getting overlap genes 
st_up25 <- as_up$symbol[as_up$symbol %in% dsub_up$symbol] %>% as.data.frame()
st_down25 <- as_down$symbol[as_down$symbol %in% dsub_down$symbol]%>% as.data.frame()

#significant test
mat <- matrix(c(105,0,0,73),ncol=2)
# Perform chi-squared test
chi_square <- chisq.test(mat)
chi_square
# data:  mat
# X-squared = 173.89, df = 1, p-value < 2.2e-16

# looking for overlap in transition groups between sub and dom comparision 
listInput <- list(st_up25= st_up25$., dt_up25=dt_up25$., st_down25 = st_down25$.,dt_down25 = dt_down25$.)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

#getting overlap genes 
t_up25 <- st_up25$.[st_up25$. %in% dt_up25$.] %>% as.data.frame()
t_down25 <- st_down25$.[st_down25$. %in% dt_down25$.]%>% as.data.frame()

#significant test
mat <- matrix(c(5,0,0,5),ncol=2)
# Perform chi-squared test
chi_square <- chisq.test(mat)
chi_square
# data:  mat
# X-squared = 6.4, df = 1, p-value = 0.01141

colnames(dt_up25)[1]<- "symbol"
dt_up <- dt_up25 %>% mutate(regulated = "up") %>% mutate(condition = "TRN vs DOM") %>% mutate(time = 25)

colnames(dt_down25)[1]<- "symbol"
dt_down <- dt_down25 %>% mutate(regulated = "down") %>% mutate(condition = "TRN vs DOM") %>% mutate(time =25 )

dt <- dt_up %>%  rbind(dt_down) %>% unique(.)

df_dd <- dd %>% select(symbol, dd_logFC= logFC, dd_pv= P.Value) %>%
  mutate(time = 25) %>%  filter(symbol %in% dt$symbol) %>% filter(symbol != "")

df_adom <- adom %>% select(symbol, adom_logFC= logFC, adom_pv= P.Value) %>%
  mutate(time = 25) %>% filter(symbol %in% dt$symbol)%>% filter(symbol != "")

dd_dt <- df_dd %>% full_join(df_adom)

dd_dt_all  <-dd_dt %>% rbind(dd_dt70) %>% unique(.)



colnames(st_up25)[1]<- "symbol"
st_up <- st_up25 %>% mutate(regulated = "up") %>% mutate(condition = "TRN vs SUB") %>% mutate(time = 25)

colnames(st_down25)[1]<- "symbol"
st_down <- st_down25 %>% mutate(regulated = "down") %>% mutate(condition = "TRN vs SUB") %>% mutate(time = 25)

st <- st_up %>%  rbind(st_down) %>% unique(.)

df_as <- as %>% select(symbol, as_logFC= logFC, as_pv= P.Value) %>%
  mutate(time = 25) %>%  filter(symbol %in% st$symbol)
df_ds<- dsub %>% select(symbol, ds_logFC= logFC, ds_pv= P.Value) %>%
  mutate(time = 25) %>% filter(symbol %in% st$symbol)

dd_st <- df_as %>% full_join(df_ds)

dd_st_all <- dd_st %>% full_join(dd_st70) %>% unique(.)



source("manuscript/brain/05_DGE_overlap_70min.R")
# dt_all <- dd_dt70 %>% rbind(dd_dt25) %>% unique(.)
write.csv(dd_dt_all,"manuscript/brain/results_tables/TRNgenesvsDOM_bothtimepoints.csv", row.names = F)

# st_all <- dd_st70 %>% rbind(dd_st25) %>% unique(.)
write.csv(dd_st_all,"manuscript/brain/results_tables/TRNgenesvsSUB_bothtimepoints.csv", row.names = F)


dd_up <- dd_dt70 %>% filter(regulated == "up")
dd_down <- dd_dt70 %>% filter(regulated == "down")

dd25_up <-dd_dt25 %>%filter(regulated == "up")
dd25_down<- dd_dt25 %>% filter(regulated == "down")


# dominant vs trn in both time points 
listInput <- list(dt_up= dd_up$symbol, dt_up25=dd25_up$symbol, dt_down = dd_down$symbol,dt_down25 = dd25_down$symbol)

UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)
dd_dt70[dd_dt70$symbol %in% dd_dt25$symbol,]


# subordinates vs trn in both time points 
st_up <- dd_st70 %>% filter(regulated == "up")
st_down <- dd_st70 %>% filter(regulated == "down")
st25_up <-dd_st25 %>%filter(regulated == "up")
st25_down<- dd_st25 %>% filter(regulated == "down")



listInput <- list(st_up= st_up$symbol, st_up25=st25_up$symbol, st_down = st_down$symbol,st_down25 = st25_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

st_up[st_up$symbol %in% st25_up$symbol,]
st_down[st_down$symbol %in% st25_down$symbol,]


