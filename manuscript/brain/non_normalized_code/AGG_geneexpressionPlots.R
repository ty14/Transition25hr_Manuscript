
library(annotables)
library(tidyverse)
grcm38 # mouse genes




asc <- read.csv("manuscript/brain/results_tables/asc_genes_mPFC.csv")
head(asc)

ax <- asc %>% filter(time == 70) %>% select(logFC = as_logFC,symbol)


ax_up <- ax %>% filter(logFC > 0.2)
ax_down <- ax %>% filter(logFC < 0.2)


limma_list <- readRDS("manuscript/brain/results/limma_mPFC_CORT_AGGgiven.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= 0.2)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 


cort <- limma_list$cort
c_up <- cort %>% filter(logFC > 0.2)
c_down <- cort %>% filter(logFC < 0.2)


agiven <- limma_list$agiven
ag_up <- agiven%>% filter(logFC > 0.2)
ag_down <- agiven %>% filter(logFC < 0.2)


intx <- limma_list$cort_agiven
in_up <- intx %>% filter(logFC > 0.2)
in_down <- intx %>% filter(logFC < 0.2)


library(UpSetR)
library(workflowr)
library(ComplexUpset)
listInput <- list(a_up= ax_up$symbol, ag_up=ag_up$symbol, ag_down = ag_down$symbol, 
                  ag_down = ax_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

ag_upx <- ax_up$symbol[ax_up$symbol %in% ag_up$symbol] %>% as.data.frame()#3
# Med30
# Pop4
# Afap1l1
ag_downx <- ax_down$symbol[ax_down$symbol %in% ag_down$symbol]%>% as.data.frame()#5
# 1 Gm10130
# 2   Pclaf
# 3   Pus10
# 4 Olfr715
# 5    Proz
#8 that overlap with agg 
#no overlap in cort 
# 3 overlap with interaction with agg and cort 
listInput <- list(a_up= ax_up$symbol, c_up=in_up$symbol, c_down = in_down$symbol, 
                  a_down = ax_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

a_in <- ax_down$symbol[ax_down$symbol %in% in_up$symbol] %>% as.data.frame()#2
#Tfap2a
# 2    Ddt
# 3  Kcns3

agenes <- c("Med30", "Pop4","Afap1l1","Gm10130","Pclaf", "Pus10", "Olfr715", "Proz", "Tfap2a", "Ddt", "Kcns3")


###############now descenders 
des <- read.csv("manuscript/brain/results_tables/des_genes_mPFC.csv")
head(des)

des <- des %>% filter(time == 70) %>% select(logFC = dd_logFC,symbol)


dx_up <- des %>% filter(logFC > 0.2)
dx_down <- des %>% filter(logFC < 0.2)


limma_list <- readRDS("manuscript/brain/results/limma_mPFC_CORT_AGGREC.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= 0.2)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 


cort <- limma_list$cort
c_up <- cort %>% filter(logFC > 0.2)
c_down <- cort %>% filter(logFC < 0.2)


agiven <- limma_list$aggrec
ag_up <- agiven%>% filter(logFC > 0.2)
ag_down <- agiven %>% filter(logFC < 0.2)


intx <- limma_list$cort_aggrec
in_up <- intx %>% filter(logFC > 0.2)
in_down <- intx %>% filter(logFC < 0.2)

listInput <- list(a_up= dx_up$symbol, ag_up=ag_up$symbol, ag_down = ag_down$symbol, 
                  ag_down = dx_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)
# one gene with opp fold change

ar_upx <- dx_up$symbol[dx_up$symbol %in% ag_up$symbol] %>% as.data.frame()#3
# Lmcd1

listInput <- list(a_up= dx_up$symbol, c_up=in_up$symbol, c_down = in_down$symbol, 
                  a_down = dx_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)
d_in_up <- dx_up$symbol[dx_up$symbol %in% in_up$symbol] %>% as.data.frame()#2
# Wfdc18
# Plin2
d_in_down <- dx_down$symbol[dx_down$symbol %in% in_down$symbol] %>% as.data.frame()#2
# 1  Cmah
# 2 Cndp1




###
agenes <- c("Med30", "Pop4","Afap1l1","Gm10130","Pclaf", "Pus10", "Olfr715", "Proz", "Tfap2a", "Ddt", "Kcns3")
dgenes <- c("Lmcd1", "Wfdc18", "Plin2", "Cmah", "Cndp1")

ex <- readRDS("manuscript/brain/results/limma_vdl_mPFC_CORT_AGGiven.RDS")
head(ex)

x <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)


x %>% 
  filter(symbol %in% agenes) -> xex



xex2 <- xex %>% pivot_longer(cols = 2:29, names_to = "ids")

p <- xex2 %>% full_join(id)
p <- p %>% filter(group %in% c('DOM',"ASC")) 
# p$group <- factor(p$group, levels = c("ASC", "DOM", "SUB"))
p$ids
p$SampleID <- p$ids


agg <- read_csv("manuscript/brain/results_tables/coldata_ALLXAGG.csv")

head(agg)

agg70<- agg %>% filter(time == 70)

x <- agg70 %>% full_join(p) %>% na.omit(.)
x$SampleID


#Tfap2a
# 2    Ddt
# 3  Kcns3

colnames(x)
xx <- x %>% filter(symbol != "Tfap2a") %>%filter(symbol != "Ddt")%>% filter(symbol != "Kcns3")
xx %>% ggplot(., aes(post.given1,value, color = group))+
               geom_point(alpha=.4, size =2)+
               geom_smooth(method = "lm", se =F)+
  ylab("Normalized Expression")+
  xlab("Rate of Aggression Given")+
  facet_wrap(factor(symbol,levels = c("Med30", "Pop4","Afap1l1","Gm10130","Pclaf", "Pus10", "Olfr715", "Proz")) ~ ., scales = 'free', ncol =8)

  
xx2 <- x %>% filter(symbol%in% c("Tfap2a", "Ddt", "Kcns3"))


# colnames(xx2)
# xx3 <- xx2 %>% select(symbol,value,post.given1, CORT) %>% pivot_longer(cols = 3:4, values_to = "bval")

# xx3 %>% ggplot(., aes(bval,value, color = name))+
  # geom_point(alpha=.4, size =2)+
  # geom_smooth(method = "lm", se =F)+
  # ylab("Normalized Expression")+
  # xlab("Rate of Aggression Given")+
  # facet_wrap(~symbol)
    
    
    
    dgenes <- c("Lmcd1", "Wfdc18", "Plin2", "Cmah", "Cndp1")
  
  ex <- readRDS("manuscript/brain/results/limma_vdl_mPFC_CORT_AGGREC.RDS")
  head(ex)
  
  x <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
    left_join(grcm38) %>%
    filter(!is.na(symbol)) 
  
  id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)
  
  
  x %>% 
    filter(symbol %in% dgenes) -> xex
  
  
  
  xex2 <- xex %>% pivot_longer(cols = 2:29, names_to = "ids")
  
  p <- xex2 %>% full_join(id)
  p <- p %>% filter(group %in% c("SUB","DES")) 
  # p$group <- factor(p$group, levels = c("ASC", "DOM", "SUB"))
  p$ids
  p$SampleID <- p$ids
  
  
  agg <- read_csv("manuscript/brain/results_tables/coldata_ALLXAGG.csv")
  
  head(agg)
  
  agg70<- agg %>% filter(time == 70)
  
  x <- agg70 %>% full_join(p) %>% na.omit(.)
  x$SampleID

 
  x %>% ggplot(., aes(post.received1,value, color = group))+
    geom_point(alpha=.4, size =2)+
    geom_smooth(method = "lm", se =F)+
    ylab("Normalized Expression")+
    xlab("Rate of Aggression Received")+
    facet_wrap(factor(symbol,levels = c("Lmcd1", "Wfdc18", "Plin2", "Cmah", "Cndp1")) ~ ., scales = 'free', ncol =8)
  
