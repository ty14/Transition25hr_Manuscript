library(annotables)
library(tidyverse)
grcm38 # mouse genes




asc <- read_csv("manuscript/brain/results_tables/asc_genes_mPFC.csv")
colnames(asc)

asc70_up <-  asc %>% filter(time == 70) %>% 
  arrange(-as_logFC) %>% head(.,10) %>% 
  select(symbol,time) 

asc70_down <-  asc %>% filter(time == 70) %>% 
  arrange(as_logFC) %>% head(.,10) %>% 
  select(symbol,time )

asc70 <- asc70_up %>% rbind(asc70_down)

asc25_up <-  asc %>% filter(time == 25) %>% 
  arrange(-as_logFC) %>% head(.,10) %>% 
  select(symbol,time)

asc25_down <-  asc %>% filter(time == 25) %>% 
  arrange(as_logFC) %>% head(.,10) %>% 
  select(symbol,time)

asc.genes <- asc70_up %>% rbind(asc70_down, asc25_up,asc25_down)



# my_logFC_threshold = 0.2

limma_list<- readRDS("manuscript/brain/results_use/limma_PFC70min_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

adom <- limma_list$ascdom
as <- limma_list$ascsub


ad70.log <- adom %>% filter(symbol %in% asc.genes$symbol) %>% select(symbol,`ASC vs. DOM 70min` = logFC)
as70.log <- as %>% filter(symbol %in% asc.genes$symbol) %>% select(symbol,`ASC vs. SUB 70min` = logFC)

asc70 <- ad70.log %>% full_join(as70.log)%>% unique(.) 

limma_list<- readRDS("manuscript/brain/results_use/limma_PFC25hr_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

adom25 <- limma_list$ascdom
as25 <- limma_list$ascsub


ad25.log <- adom25 %>% filter(symbol %in% asc.genes$symbol) %>% select(symbol,`ASC vs. DOM 25hr` = logFC)
as25.log <- as25 %>% filter(symbol %in% asc.genes$symbol) %>% select(symbol,`ASC vs. SUB 25hr` = logFC)

asc25 <- ad25.log %>% full_join(as25.log) %>% unique(.)

asc.hm <- asc70 %>% full_join(asc25) %>% unique(.) 

asc.hm$dir <-ifelse(asc.hm$symbol %in% asc70_up$symbol, 'up70',"") 
asc.hm$dir <-ifelse(asc.hm$symbol %in% asc70_down$symbol, 'down70',asc.hm$dir)
asc.hm$dir <-ifelse(asc.hm$symbol %in% asc25_down$symbol, 'down25',asc.hm$dir)
asc.hm$dir <-ifelse(asc.hm$symbol %in% asc25_up$symbol, 'up25',asc.hm$dir)
asc.hm$dir <- factor(asc.hm$dir, levels = c("up70", "down70", "up25", "down25"))


asc.hm <- asc.hm %>% arrange(dir) %>% column_to_rownames(., var = "symbol") %>% select(-dir)

asc.hmx <- as.matrix(asc.hm)
library(ComplexHeatmap)
th <- t(asc.hmx)
colnames(asc.hmx) <- substr(colnames(asc.hmx),1,12)

png("manuscript/brain/imgs/ASC_StableHeatmap.png",width=3,height=12,units="in",res=1200)

Heatmap(asc.hmx, name = "LogFC",cluster_rows = FALSE,
        show_heatmap_legend = FALSE, rect_gp = gpar(col = "white", lwd = 2))

dev.off()


# rownames(th) <- substr(rownames(th),1,12)

ahx <- asc.hm %>% select(2,1,4,3) %>% as.matrix(.)

th <- t(ahx)
png("manuscript/brain/imgs/ASC_StableHeatmap2.png",width=18,height=2.95,units="in",res=1200)

Heatmap(th, name = "LogFC",cluster_columns = FALSE,cluster_rows = FALSE,
        show_heatmap_legend = F, column_names_rot = 60, rect_gp = gpar(col = "white", lwd = 2))

dev.off()




# Heatmap(asc.hmx, name = "LogFC",cluster_rows = FALSE,
# column_names_rot = 45, rect_gp = gpar(col = "white", lwd = 2))

# dev.off()
