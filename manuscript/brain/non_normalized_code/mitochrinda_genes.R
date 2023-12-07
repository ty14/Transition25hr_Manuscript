
mt <- read_csv("manuscript/brain/gene_sets/Mouse.MitoCarta3.csv")
head(mt)
mt %>%
  as_tibble() %>%.$Symbol -> mtx

ros <- read_csv("manuscript/brain/gene_sets/OXPHOS.csv")
head(ros)
ros %>%
  as_tibble() %>%.$Symbol -> rosx

#70 min 
my_logFC_threshold = 0.2

limma_list<- readRDS("manuscript/brain/results/limma_PFC70min_ReorganizedGroup.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>%
  map(~filter(.,!is.na(entrez))) 

dd <- limma_list$desdom

ds <- limma_list$dessub

ad <- limma_list$ascdom

as <- limma_list$ascsub
trans70 <- c("Kcns3", "Crybg2", "Stx19", "Cmah", "Cndp1", "Dhx16", "Zfp280d", 'Chka', 'Kmt2a', 'Sfmbt1', 'Plekha6',"Crtac1", "Itga9", "D630023F18Rik", "Tmem147", "Slc7a4", "Egfl7", "Kcnip3", "Uqcc2")
trans25 <- c("Exosc9", "Gstm7", "Zscan18", "Synm", "Fam117b","Esrrg", "Fam169b", "Cdc23", "Mief1", "Nfkb1" )

trans70[trans70 %in% rosx]
# [1] "Uqcc2"

asc <- read_csv("manuscript/brain/results_tables/asc_genes_mPFC.csv")
asc %>% filter(symbol %in% rosx) #7, 6 at 70 min 
des <- read_csv("manuscript/brain/results_tables/des_genes_mPFC.csv")
des %>% filter(symbol %in% rosx) #2 at 70 min 

dd %>% filter(symbol %in% rosx)#57
ds %>% filter(symbol %in% rosx)#101
ad%>% filter(symbol %in% rosx)#63
as %>% filter(symbol %in% rosx)#


dd25 %>% filter(symbol %in% rosx)#57
ds25 %>% filter(symbol %in% rosx)#101
ad25%>% filter(symbol %in% rosx)#63
as25 %>% filter(symbol %in% rosx)#

ddmt <- dd %>% filter(symbol %in% mtx)#57
dsmt <- ds %>% filter(symbol %in% mtx)#101
admt <- ad%>% filter(symbol %in% mtx)#63
asmt <-  as %>% filter(symbol %in% mtx)#66

ddmt_up <- ddmt %>% filter(logFC > 0.2)
dsmt_up <- dsmt %>% filter(logFC > 0.2)
admt_up <- admt %>% filter(logFC > 0.2)
asmt_up <-  asmt %>% filter(logFC > 0.2)

ddmt_down <- ddmt %>% filter(logFC < 0.2)
dsmt_down <- dsmt %>% filter(logFC < 0.2)
admt_down <- admt %>% filter(logFC < 0.2)
asmt_down <-  asmt %>% filter(logFC < 0.2)


library(UpSetR)
library(workflowr)
library(ComplexUpset)
listInput <- list(ddmt_up = ddmt_up$symbol, ddmt_down=ddmt_down$symbol,
                  dsmt_up  = dsmt_up$symbol, dsmt_down = dsmt_down$symbol,
                  admt_up  = admt_up$symbol, admt_down = admt_down$symbol,
                  asmt_up  = asmt_up$symbol, asmt_down = asmt_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 8, order.by = "freq", keep.order = F)




#14 genes over lap in asmt and admt up
asmt_up$symbol[asmt_up$symbol %in% admt_up$symbol]
# 14 genes over lap in asmt and ddmt up
asmt_up$symbol[asmt_up$symbol %in% dsmt_up$symbol]
#9 genes overlap in ddmt and dsmt up
ddmt_up$symbol[ddmt_up$symbol %in% dsmt_up$symbol]
#6 genes over lap in admt and ddmt up 
admt_up$symbol[admt_up$symbol %in% ddmt_up$symbol]
# 5 genes overlap in ddmt and dsmt down
ddmt_down$symbol[ddmt_down$symbol %in% dsmt_down$symbol]
#4 genes overlap in asmt and admt down 
admt_down$symbol[admt_down$symbol %in% asmt_down$symbol]



#25 hour
my_logFC_threshold = 0.2

limma_list<- readRDS("manuscript/brain/results/limma_PFC_ReorganizedGroups_outlierRemoved.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

dd25 <- limma_list$desdom

ds25 <- limma_list$dessub

ad25 <- limma_list$ascdom

as25 <- limma_list$ascsub


dd25mt<- dd25 %>% filter(symbol %in% mtx)#22
ds25mt <- ds25 %>% filter(symbol %in% mtx)#53
ad25mt <- ad25 %>% filter(symbol %in% mtx)#37
as25mt <- as25 %>% filter(symbol %in% mtx)#53



ddmt_up <- dd25mt %>% filter(logFC > 0.2)
dsmt_up <- ds25mt %>% filter(logFC > 0.2)
admt_up <- ad25mt %>% filter(logFC > 0.2)
asmt_up <-  as25mt %>% filter(logFC > 0.2)

ddmt_down <- dd25mt %>% filter(logFC < 0.2)
dsmt_down <- ds25mt %>% filter(logFC < 0.2)
admt_down <- ad25mt %>% filter(logFC < 0.2)
asmt_down <-  as25mt %>% filter(logFC < 0.2)


library(UpSetR)
library(workflowr)
library(ComplexUpset)
listInput <- list(ddmt_up = ddmt_up$symbol, ddmt_down=ddmt_down$symbol,
                  dsmt_up  = dsmt_up$symbol, dsmt_down = dsmt_down$symbol,
                  admt_up  = admt_up$symbol, admt_down = admt_down$symbol,
                  asmt_up  = asmt_up$symbol, asmt_down = asmt_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 8, order.by = "freq", keep.order = F)
