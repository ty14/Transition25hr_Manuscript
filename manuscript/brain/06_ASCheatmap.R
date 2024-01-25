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

asc25 <- asc25_up %>% rbind(asc25_down)

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


### heatmap significant testing

adom %>% filter(symbol %in% asc70$symbol)
as %>% filter(symbol %in% asc70$symbol)

adom25 %>% filter(symbol %in% asc70$symbol) %>% filter(P.Value < 0.06)
as25 %>% filter(symbol %in% asc70$symbol) %>% filter(P.Value < 0.06)

adom25 %>% filter(symbol %in% asc25$symbol)
as25 %>% filter(symbol %in% asc25$symbol)


adom %>% filter(symbol %in% asc25$symbol) %>% filter(P.Value < 0.06)
as %>% filter(symbol %in% asc25$symbol) %>% filter(P.Value < 0.06)

###heatmap 
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


#############
#logfc scatter plots

 my_logFC_threshold = 0.2
#70min data 
limma_list<- readRDS("manuscript/brain/results_use/limma_PFC70min_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

adom <- limma_list$ascdom
as <- limma_list$ascsub

# 25 hr data 
limma_list<- readRDS("manuscript/brain/results_use/limma_PFC25hr_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>%
  map(~filter(.,!is.na(entrez))) 

adom25 <- limma_list$ascdom
as25 <- limma_list$ascsub


adom70_up <- adom %>% filter(logFC > 0.2)
adom25_up <- adom25 %>% filter(logFC > 0.2)
adom70_up$symbol[adom70_up$symbol %in% adom25_up$symbol] #381


adom70_down <- adom %>% filter(logFC < 0.2)
adom25_down <- adom25 %>% filter(logFC < 0.2)
adom70_down$symbol %in% adom25_down$symbol

library(UpSetR)
library(workflowr)
library(ComplexUpset)
listInput <- list(adom70_up= adom70_up$symbol, adom25_up=adom25_up$symbol, adom70_down = adom70_down$symbol, 
                  adom25_down = adom25_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)




as70_up <- as %>% filter(logFC > 0.2)
as25_up <- as25 %>% filter(logFC > 0.2)
as70_up$symbol[as70_up$symbol %in% as25_up$symbol] #8


as70_down <- as %>% filter(logFC < 0.2)
as25_down <- as25 %>% filter(logFC < 0.2)
as70_down$symbol %in% as25_down$symbol

library(UpSetR)
library(workflowr)
library(ComplexUpset)
listInput <- list(as70_up= as70_up$symbol, as25_up=as25_up$symbol, as70_down = as70_down$symbol, 
                  as25_down = as25_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)


#############
#logfc scatter plots

my_logFC_threshold = 0.2
#70min data 
limma_list<- readRDS("manuscript/brain/results_use/limma_PFC70min_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

adomx<- limma_list$ascdom
asx <- limma_list$ascsub

# 25 hr data 
limma_list<- readRDS("manuscript/brain/results_use/limma_PFC25hr_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>%
  map(~filter(.,!is.na(entrez))) 

adom25x <- limma_list$ascdom
as25x <- limma_list$ascsub

#2901 genes at 70 min
ad70 <- adomx %>% select(ad70_lf = logFC, ad70_pv = P.Value, symbol) %>% filter(.,abs(ad70_lf) >= my_logFC_threshold)
#2421 genes at 25 hr
ad25 <- adom25x %>% select(ad25_lf = logFC,ad25_pv = P.Value, symbol)%>% filter(.,abs(ad25_lf) >= my_logFC_threshold)
ad.df <- ad70 %>% full_join(ad25) %>%  na.omit(.)
colnames(ad.df)
ad.df$Sig <- ifelse(ad.df$ad70_pv <= 0.05 & ad.df$ad25_pv <= 0.05, "SIG", "N.S.")
ad.df <- ad.df %>% unique(.)
df.c <- ad.df %>% filter(Sig == "SIG")
df.c.same <- df.c %>% filter(ad70_lf >=0.2 & ad25_lf >= 0.2)
df.c.samex <- df.c %>% filter(ad70_lf <=0.2 & ad25_lf <= 0.2)
df.c.samexx <- df.c.same %>% rbind(df.c.samex)

df.c.diff <- df.c %>% filter(ad70_lf >=0.2 & ad25_lf <= 0.2)
df.c.diffx <- df.c %>% filter(ad70_lf <=0.2 & ad25_lf >= 0.2)
df.c.diffxx <- df.c.diff %>% rbind(df.c.diffx)

ad.same <- df.c.samexx

ad.diff <- df.c.diffxx

df.g <- ad.df %>% filter(Sig == "N.S.")

ggplot(ad.df, aes(ad70_lf, ad25_lf, group = Sig))+
  geom_point(data = df.g, alpha = 0.1, shape = 21, size = 3, fill = "grey50") +
  geom_point(data = df.c.samexx, alpha = .4, shape = 21, size = 4, fill = "blue") + 
  geom_point(data = df.c.diffxx, alpha = .4, shape = 21, size = 4, fill = "red") + 
  scale_x_continuous(limits = c(-2,2))+
  geom_text_repel(data= df.c.samexx, aes(label = symbol),vjust=0.5,max.overlaps = Inf)+
  ylim(-2,2)+ theme_classic()+
  ylab("LogFC at 25 hr")+
  xlab("LogFC at 70 min")+
  ggtitle("ASC vs. DOM")+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15))

####### asc vs sub
#3131 genes at 70 min
ad70 <- asx %>% select(ad70_lf = logFC, ad70_pv = P.Value, symbol) %>% filter(.,abs(ad70_lf) >= my_logFC_threshold)
#2775 genes at 25 hr
ad25 <- as25x %>% select(ad25_lf = logFC,ad25_pv = P.Value, symbol)%>% filter(.,abs(ad25_lf) >= my_logFC_threshold)
ad.df <- ad70 %>% full_join(ad25) %>%  na.omit(.)
colnames(ad.df)
ad.df$Sig <- ifelse(ad.df$ad70_pv <= 0.05 & ad.df$ad25_pv <= 0.05, "SIG", "N.S.")
ad.df <- ad.df %>% unique(.)
df.c <- ad.df %>% filter(Sig == "SIG")
df.c.same <- df.c %>% filter(ad70_lf >=0.2 & ad25_lf >= 0.2)
df.c.samex <- df.c %>% filter(ad70_lf <=0.2 & ad25_lf <= 0.2)
df.c.samexx <- df.c.same %>% rbind(df.c.samex)

df.c.diff <- df.c %>% filter(ad70_lf >=0.2 & ad25_lf <= 0.2)
df.c.diffx <- df.c %>% filter(ad70_lf <=0.2 & ad25_lf >= 0.2)
df.c.diffxx <- df.c.diff %>% rbind(df.c.diffx)

x <- df.c.samex %>% arrange(ad70_lf) %>% head(.)
y <- df.c.samexx %>% arrange(-ad70_lf) %>% head(.)
for_label <- x %>% rbind(y)
df.g <- ad.df %>% filter(Sig == "N.S.")

ggplot(ad.df, aes(ad70_lf, ad25_lf, group = Sig))+
  geom_point(data = df.g, alpha = 0.1, shape = 21, size = 3, fill = "grey50") +
  geom_point(data = df.c.samexx, alpha = .5, shape = 21, size = 4, fill = "blue") + 
  geom_point(data = df.c.diffxx, alpha = .5, shape = 21, size = 4, fill = "red") + 
  scale_x_continuous(limits = c(-2,2))+
  geom_text_repel(data= for_label, aes(label = symbol), vjust =0.6)+
  ylim(-2,2)+ theme_classic()+
  ylab("LogFC at 25 hr")+
  xlab("LogFC at 70 min")+
  ggtitle("ASC vs. SUB")+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15))




ad.same$symbol[ad.same$symbol %in% df.c.samexx$symbol]
ad.diff$symbol[ad.diff$symbol %in% df.c.diffxx$symbol]
