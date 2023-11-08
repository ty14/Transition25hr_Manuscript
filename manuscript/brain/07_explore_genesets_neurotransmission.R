library(annotables)
library(tidyverse)
grcm38 # mouse genes


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


endo <- c("Oxtr","Esr1", "Esrrg", "Esrra",'Nr3c4','Nr3c2', 'Nr3c1','Nr2c2',
          'Nr4a1','Nr4a3', 'Nr6a1','Avpi1',"Crhr1", 'Crhbp',"Lhcgr", "Avp", "Oxt", "Ar", 
          "Crh", "Cort", "Gnrh", "Pomc", "Mt1", 'Mt2', 'Mtr', 'Srd5a1', 'Srd5a2', 'Hsd17b3', 
          'Cyp17a1', 'Cyp19a1', 'Srd5a3','Hsd3b1',"Hsd3b2", "Trh", "Trhr", "Thrb", "Thra", 
          "Dio2", "Crym", "Mybpc1")

dd %>% filter(symbol %in% endo)#2
ds %>% filter(symbol %in% endo)#6
ad %>% filter(symbol %in% endo)#2
as %>% filter(symbol %in% endo)#3

dd25 %>% filter(symbol %in% endo)#1
ds25 %>% filter(symbol %in% endo)#3
ad25 %>% filter(symbol %in% endo)#2
as25 %>% filter(symbol %in% endo)#1


dop <- c("Th", "Dbh", "Ddc", "Comt", "Maoa", "Maob", "Slc6a3", "Dat1", "Snca", "Sncb", "Drd1", 
         "Drd2",'Drd5', "Drd3", "Ppp1r1b", "Nr4a2")

dd %>% filter(symbol %in% dop)#0
ds %>% filter(symbol %in% dop)#2
ad %>% filter(symbol %in% dop)#1
as %>% filter(symbol %in% dop)#1

dd25 %>% filter(symbol %in% dop)#0
ds25 %>% filter(symbol %in% dop)#2
ad25 %>% filter(symbol %in% dop)#0
as25 %>% filter(symbol %in% dop)#0

ser <- c("Tph2","Ddc", "Maoa", "Maob", "Slc18a2", 'Slc6a4', "Aldh1b1","Ald3a1", 'Htr1a',
         'Htr1b', 'Htr2a', 'Htr2c', 'Htr3a', "Htr4", "Htr5b", "Htr6", "Htr7" )

dd %>% filter(symbol %in% ser)#0
ds %>% filter(symbol %in% ser)#3
ad %>% filter(symbol %in% ser)#0
as %>% filter(symbol %in% ser)#2

dd25 %>% filter(symbol %in% ser)#0
ds25 %>% filter(symbol %in% ser)#0
ad25 %>% filter(symbol %in% ser)#1
as25 %>% filter(symbol %in% ser)#0

glu <- c("Grin1", "Grin2a","Grin2b","Grin2c","Grin2d","Gria1","Gria2", "Gria3", "Gria4",
         "Slc1a1", "Slc1a2", "Slc1a6", "slc17a7","slc17a6", "Grid1", "Grid2", "Gad1", "Gad2",
         "Camk2a","Camk2b","Shank3", "Syp" )

dd %>% filter(symbol %in% glu)#0
ds %>% filter(symbol %in% glu)#3
ad %>% filter(symbol %in% glu)#0
as %>% filter(symbol %in% glu)#2

dd25 %>% filter(symbol %in% glu)#0
ds25 %>% filter(symbol %in% glu)#0
ad25 %>% filter(symbol %in% glu)#1
as25 %>% filter(symbol %in% glu)#0


gaba <- c("Gad1", "Gad2", "Gad65", "Gad67", "Slc6a1", "Gat1", "Gabra1","Gabra2","Gabra3",
         "Gabrb1","Gabrb2", "Gabrb3", "Abat", "Slc32a1")

dd %>% filter(symbol %in% gaba)#0
ds %>% filter(symbol %in% gaba)#3
ad %>% filter(symbol %in% gaba)#0
as %>% filter(symbol %in% gaba)#2

dd25 %>% filter(symbol %in% gaba)#0
ds25 %>% filter(symbol %in% gaba)#0
ad25 %>% filter(symbol %in% gaba)#1
as25 %>% filter(symbol %in% gaba)#0




