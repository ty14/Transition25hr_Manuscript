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


ad70.log <- adom %>% filter(symbol %in% asc.genes$symbol) %>% select(symbol,`asc vs. dom 70min` = logFC)
as70.log <- as %>% filter(symbol %in% asc.genes$symbol) %>% select(symbol,`asc vs. sub 70min` = logFC)

asc70 <- ad70.log %>% full_join(as70.log)%>% unique(.) 

limma_list<- readRDS("manuscript/brain/results_use/limma_PFC25hr_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

adom25 <- limma_list$ascdom
as25 <- limma_list$ascsub


ad25.log <- adom25 %>% filter(symbol %in% asc.genes$symbol) %>% select(symbol,`asc vs. dom 25hr` = logFC)
as25.log <- as25 %>% filter(symbol %in% asc.genes$symbol) %>% select(symbol,`asc vs. sub 25hr` = logFC)

asc25 <- ad25.log %>% full_join(as25.log) %>% unique(.)

asc.hm <- asc70 %>% full_join(asc25) %>% unique(.) 

asc.hm$dir <-ifelse(asc.hm$symbol %in% asc70_up$symbol, 'up70',"") 
asc.hm$dir <-ifelse(asc.hm$symbol %in% asc70_down$symbol, 'down70',asc.hm$dir)
asc.hm$dir <-ifelse(asc.hm$symbol %in% asc25_down$symbol, 'down25',asc.hm$dir)
asc.hm$dir <-ifelse(asc.hm$symbol %in% asc25_up$symbol, 'up25',asc.hm$dir)
asc.hm$dir <- factor(asc.hm$dir, levels = c("up70", "down70", "up25", "down25"))


asc.hm <- asc.hm %>% arrange(dir) %>% column_to_rownames(., var = "symbol")

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

asc.hmx <- as.matrix(asc.hm)
heatmap(asc.hmx,
        Rowv=NA, Colv=NA, col=rev(brewer.pal(9,"RdBu")))

# colnames(asc.hmx) <- substr(colnames(asc.hmx), 1,11)
library(plotly)
library(heatmaply)
heatmaply(th, 
          dendrogram = "none",
          xlab = "", ylab = "", 
          main = "",
          scale = "column",
          margins = c(60,100,40,20),
          grid_color = "white",
          grid_width = 0.00001,
          titleX = FALSE,
          hide_colorbar = TRUE,
          branches_lwd = 0.1,
          # label_names = c("Country", "Feature:", "Value"),
          fontsize_row = 5, fontsize_col = 5,
          # labCol = colnames(mat),
          # labRow = rownames(mat),
          plot_method = c("ggplot"),
          heatmap_layers = theme(axis.line=element_blank(), axis.title = element_text(size =100)))

library(ComplexHeatmap)
th <- t(asc.hmx)
Heatmap(th, name = "LogFC", cluster_columns = FALSE,
        column_names_rot = 45)
