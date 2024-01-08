library(annotables)
library(tidyverse)
grcm38 # mouse genes




des <- read_csv("manuscript/brain/results_tables/des_genes_mPFC.csv")
colnames(des)

des70_up <-  des %>% filter(time == 70) %>% 
  arrange(-dsub_logFC) %>% head(.,10) %>% 
  select(symbol,time) 

des70_down <-  des %>% filter(time == 70) %>% 
  arrange(dsub_logFC) %>% head(.,10) %>% 
  select(symbol,time )

des70 <- des70_up %>% rbind(des70_down)

des25_up <-  des %>% filter(time == 25) %>% 
  arrange(-dsub_logFC) %>% head(.,10) %>% 
  select(symbol,time)

des25_down <-  des %>% filter(time == 25) %>% 
  arrange(dsub_logFC) %>% head(.,10) %>% 
  select(symbol,time)

des25 <- des25_up %>% rbind(des25_down)

des.genes <- des70_up %>% rbind(des70_down, des25_up,des25_down)

x <- des %>% filter(time == 70)  
y <- des %>% filter(time == 25) 

x$symbol[x$symbol %in% y$symbol]
# my_logFC_threshold = 0.2

limma_list<- readRDS("manuscript/brain/results_use/limma_PFC70min_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

dd <- limma_list$desdom
as <- limma_list$dessub




limma_list<- readRDS("manuscript/brain/results_use/limma_PFC25hr_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

dd25 <- limma_list$desdom
as25 <- limma_list$dessub


### heatmap significant testing

dd %>% filter(symbol %in% des70$symbol)
as %>% filter(symbol %in% des70$symbol)

dd25 %>% filter(symbol %in% des70$symbol) %>% filter(P.Value < 0.05)
as25 %>% filter(symbol %in% des70$symbol) %>% filter(P.Value < 0.05)

dd25 %>% filter(symbol %in% des25$symbol)
as25 %>% filter(symbol %in% des25$symbol)


dd %>% filter(symbol %in% des25$symbol) %>% filter(P.Value < 0.06)
as %>% filter(symbol %in% des25$symbol) %>% filter(P.Value < 0.06)

###heatmap 

ad70.log <- dd %>% filter(symbol %in% des.genes$symbol) %>% select(symbol,`DES vs. DOM 70min` = logFC)
as70.log <- as %>% filter(symbol %in% des.genes$symbol) %>% select(symbol,`DES vs. SUB 70min` = logFC)

des70 <- ad70.log %>% full_join(as70.log)%>% unique(.) 


ad25.log <- dd25 %>% filter(symbol %in% des.genes$symbol) %>% select(symbol,`DES vs. DOM 25hr` = logFC)
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



# th <- t(ahx)
# png("manuscript/brain/imgs/ASC_StableHeatmap2.png",width=18,height=3,units="in",res=1200)
# 
# Heatmap(th, name = "LogFC",cluster_columns = FALSE,cluster_rows = FALSE,
#         show_heatmap_legend = F, column_names_rot = 60, rect_gp = gpar(col = "white", lwd = 2))
# 
# dev.off()
# 

# rownames(th) <- substr(rownames(th),1,12)
png("manuscript/brain/imgs/des_StableHeatmap2.png",width=18,height=2.75,units="in",res=1200)

Heatmap(th, name = "LogFC",cluster_columns = FALSE,cluster_rows = FALSE,
        show_heatmap_legend = F, column_names_rot = 60, rect_gp = gpar(col = "white", lwd = 2))

dev.off()

dev.off()



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

dd <- limma_list$desdom
ds <- limma_list$dessub

# 25 hr data 
limma_list<- readRDS("manuscript/brain/results_use/limma_PFC25hr_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>%
  map(~filter(.,!is.na(entrez))) 

dd25 <- limma_list$desdom
ds25 <- limma_list$dessub


dd70_up <- dd %>% filter(logFC > 0.2)
dd25_up <- dd25 %>% filter(logFC > 0.2)
dd70_up$symbol[dd70_up$symbol %in% dd25_up$symbol] #381


dd70_down <- dd %>% filter(logFC < 0.2)
dd25_down <- dd25 %>% filter(logFC < 0.2)
dd70_down$symbol %in% dd25_down$symbol

library(UpSetR)
library(workflowr)
library(ComplexUpset)
listInput <- list(dd70_up= dd70_up$symbol, dd25_up=dd25_up$symbol, dd70_down = dd70_down$symbol, 
                  dd25_down = dd25_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)




ds70_up <- ds %>% filter(logFC > 0.2)
ds25_up <- ds25 %>% filter(logFC > 0.2)
ds70_up$symbol[ds70_up$symbol %in% ds25_up$symbol] #8


ds70_down <- ds %>% filter(logFC < 0.2)
ds25_down <- ds25 %>% filter(logFC < 0.2)
ds70_down$symbol %in% ds25_down$symbol

library(UpSetR)
library(workflowr)
library(ComplexUpset)
listInput <- list(ds70_up= ds70_up$symbol, ds25_up=ds25_up$symbol, ds70_down = ds70_down$symbol, 
                  ds25_down = ds25_down$symbol)


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

ddx<- limma_list$desdom
dsx <- limma_list$dessub

# 25 hr data 
limma_list<- readRDS("manuscript/brain/results_use/limma_PFC25hr_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>%
  map(~filter(.,!is.na(entrez))) 

dd25x <- limma_list$desdom
ds25x <- limma_list$dessub

#2901 genes at 70 min
ad70 <- ddx %>% select(ad70_lf = logFC, ad70_pv = P.Value, symbol) %>% filter(.,abs(ad70_lf) >= my_logFC_threshold)
#2421 genes at 25 hr
ad25 <- dd25x %>% select(ad25_lf = logFC,ad25_pv = P.Value, symbol)%>% filter(.,abs(ad25_lf) >= my_logFC_threshold)
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


x <- df.c.samex %>% arrange(ad70_lf) %>% head(.,4)
y <- df.c.samexx %>% arrange(-ad70_lf) %>% head(.,4)
for_label <- x %>% rbind(y)
df.g <- ad.df %>% filter(Sig == "N.S.")



ggplot(ad.df, aes(ad70_lf, ad25_lf, group = Sig))+
  geom_point(data = df.g, alpha = 0.1, shape = 21, size = 3, fill = "grey50") +
  geom_point(data = df.c.samexx, alpha = .4, shape = 21, size = 4, fill = "blue") + 
  geom_point(data = df.c.diffxx, alpha = .4, shape = 21, size = 4, fill = "red") + 
  scale_x_continuous(limits = c(-2,2))+
  geom_text_repel(data= for_label, aes(label = symbol),vjust=.25,max.overlaps = Inf)+
  ylim(-2,2)+ theme_classic()+
  ylab("LogFC at 25 hr")+
  xlab("LogFC at 70 min")+
  ggtitle("DES vs. DOM")+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15))


#2901 genes at 70 min
ad70 <- dsx %>% select(ad70_lf = logFC, ad70_pv = P.Value, symbol) %>% filter(.,abs(ad70_lf) >= my_logFC_threshold)
#2421 genes at 25 hr
ad25 <- ds25x %>% select(ad25_lf = logFC,ad25_pv = P.Value, symbol)%>% filter(.,abs(ad25_lf) >= my_logFC_threshold)
ad.df <- ad70 %>% full_join(ad25) %>%  na.omit(.)
colnames(ad.df)
ad.df$Sig <- ifelse(ad.df$ad70_pv <= 0.05 & ad.df$ad25_pv <= 0.05, "SIG", "N.S.")
ad.df <- ad.df %>% unique(.)
ad.df$ad70_lf <- ifelse(ad.df$ad70_lf > 2.2, 2, ad.df$ad70_lf)
df.c <- ad.df %>% filter(Sig == "SIG")
df.c.same <- df.c %>% filter(ad70_lf >=0.2 & ad25_lf >= 0.2)
df.c.samex <- df.c %>% filter(ad70_lf <=0.2 & ad25_lf <= 0.2)
df.c.samexx <- df.c.same %>% rbind(df.c.samex)

df.c.diff <- df.c %>% filter(ad70_lf >=0.2 & ad25_lf <= 0.2)
df.c.diffx <- df.c %>% filter(ad70_lf <=0.2 & ad25_lf >= 0.2)
df.c.diffxx <- df.c.diff %>% rbind(df.c.diffx)


x <- df.c.samex %>% arrange(ad70_lf) %>% head(.,5)
y <- df.c.samexx %>% arrange(-ad70_lf) %>% head(.,5)
for_label <- x %>% rbind(y)
df.g <- ad.df %>% filter(Sig == "N.S.")

 

ggplot(ad.df, aes(ad70_lf, ad25_lf, group = Sig))+
  geom_point(data = df.g, alpha = 0.1, shape = 21, size = 3, fill = "grey50") +
  geom_point(data = df.c.samexx, alpha = .4, shape = 21, size = 4, fill = "blue") + 
  geom_point(data = df.c.diffxx, alpha = .4, shape = 21, size = 4, fill = "red") + 
  scale_x_continuous(limits = c(-2,2))+
  geom_text_repel(data= for_label, aes(label = symbol),vjust=.5,max.overlaps = Inf)+
  ylim(-2,2)+ theme_classic()+
  ylab("LogFC at 25 hr")+
  xlab("LogFC at 70 min")+
  ggtitle("DES vs. SUB")+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15))

