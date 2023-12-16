library(annotables)
library(tidyverse)
grcm38 # mouse genes




des <- read_csv("manuscript/brain/results_tables/des_genes_mPFC.csv")
colnames(des)

des70_up <-  des %>% filter(time == 70) %>% 
  arrange(-dd_logFC) %>% head(.,10) %>% 
  select(symbol,time) 

des70_down <-  des %>% filter(time == 70) %>% 
  arrange(dd_logFC) %>% head(.,10) %>% 
  select(symbol,time )

des70 <- des70_up %>% rbind(des70_down)

des25_up <-  des %>% filter(time == 25) %>% 
  arrange(-dd_logFC) %>% head(.,11) %>% 
  select(symbol,time)

des25_down <-  des %>% filter(time == 25) %>% 
  arrange(dd_logFC) %>% head(.,10) %>% 
  select(symbol,time)

des.genes <- des70_up %>% rbind(des70_down, des25_up,des25_down)



# my_logFC_threshold = 0.2

limma_list<- readRDS("manuscript/brain/results_use/limma_PFC70min_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

adom <- limma_list$desdom
as <- limma_list$dessub


ad70.log <- adom %>% filter(symbol %in% des.genes$symbol) %>% select(symbol,`DES vs. DOM 70min` = logFC)
as70.log <- as %>% filter(symbol %in% des.genes$symbol) %>% select(symbol,`DES vs. SUB 70min` = logFC)

des70 <- ad70.log %>% full_join(as70.log)%>% unique(.) 

limma_list<- readRDS("manuscript/brain/results_use/limma_PFC25hr_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

adom25 <- limma_list$desdom
as25 <- limma_list$dessub


ad25.log <- adom25 %>% filter(symbol %in% des.genes$symbol) %>% select(symbol,`DES vs. DOM 25hr` = logFC)
as25.log <- as25 %>% filter(symbol %in% des.genes$symbol) %>% select(symbol,`DES vs. SUB 25hr` = logFC)

des25 <- ad25.log %>% full_join(as25.log) %>% unique(.)

des.hm <- des70 %>% full_join(des25) %>% unique(.) 

des.hm$dir <-ifelse(des.hm$symbol %in% des70_up$symbol, 'up70',"") 
des.hm$dir <-ifelse(des.hm$symbol %in% des70_down$symbol, 'down70',des.hm$dir)
des.hm$dir <-ifelse(des.hm$symbol %in% des25_down$symbol, 'down25',des.hm$dir)
des.hm$dir <-ifelse(des.hm$symbol %in% des25_up$symbol, 'up25',des.hm$dir)
des.hm$dir <- factor(des.hm$dir, levels = c("up70", "down70", "up25", "down25"))


des.hm <- des.hm %>% arrange(dir) %>% column_to_rownames(., var = "symbol") %>% select(-dir)
# table(des.hm$dir)
des.hmx <- as.matrix(des.hm)
library(ComplexHeatmap)
th <- t(des.hmx)
colnames(des.hmx) <- substr(colnames(des.hmx),1,12)

png("manuscript/brain/imgs/des_StableHeatmap.png",width=3,height=12,units="in",res=1200)

Heatmap(des.hmx, name = "LogFC",cluster_rows = FALSE,
        show_heatmap_legend = FALSE, rect_gp = gpar(col = "white", lwd = 2))

dev.off()

ahx <- des.hm %>% select(2,1,4,3) %>% as.matrix(.)

th <- t(ahx)
png("manuscript/brain/imgs/ASC_StableHeatmap2.png",width=18,height=3,units="in",res=1200)

Heatmap(th, name = "LogFC",cluster_columns = FALSE,cluster_rows = FALSE,
        show_heatmap_legend = F, column_names_rot = 60, rect_gp = gpar(col = "white", lwd = 2))

dev.off()


# rownames(th) <- substr(rownames(th),1,12)
png("manuscript/brain/imgs/des_StableHeatmap2.png",width=18,height=2.5,units="in",res=1200)

Heatmap(th, name = "LogFC",cluster_columns = FALSE,cluster_rows = FALSE,
        show_heatmap_legend = F, column_names_rot = 60, rect_gp = gpar(col = "white", lwd = 2))

dev.off()

dev.off()
