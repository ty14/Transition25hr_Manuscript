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
library(EnhancedVolcano)
library(tidyverse)
grcm38 # mouse genes


my_logFC_threshold = 0.2

limma_list<- readRDS("manuscript/brain/results_use/limma_PFC70min_Norm_RG.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~unique(.)) %>% 
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(.,symbol, entrez))) %>%
  map(~unique(.,)) %>% 
  map(~filter(.,!is.na(entrez)))  


y1a <- limma_list$ascdom
y2a <- limma_list$desdom
y3a <- limma_list$ascsub
y4a <- limma_list$dessub

###volano plots
head(y1a)
dc <- y1a %>% mutate(contrast = "ASC vs. DOM") %>% mutate(log10 = -log10(P.Value))

dc$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
dc$diffexpressed[dc$logFC > 0.2 & dc$P.Value < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
dc$diffexpressed[dc$logFC < -0.2 & dc$P.Value < 0.05] <- "DOWN"


dcx <- dc %>% filter(.,logFC >= 0.8)%>% filter(P.Value < 0.05)
dcxx <- dc %>% filter(.,logFC <= -1.25) %>% filter(P.Value < 0.05)
dc$log10 <- ifelse(dc$log10 == Inf, 4,dc$log10)
# dcx <- dc %>% filter(.,abs(logFC) >= 1.5) %>% filter(P.Value < 0.05) 

vp_ad <- ggplot(data = dc, 
                aes(x = logFC, 
                    y = log10, 
                    colour=diffexpressed)) +
  geom_point(alpha=0.25, size=3.5) +
  scale_color_manual(values=c("orange", "grey","purple4"))+
  xlim(c(-2.5, 2.5)) +
  ylim(0,4)+
  # geom_vline(xintercept=c(-0.2,0.2),lty=4,col="black",lwd=0.8) +
  geom_vline(xintercept=c(-.5,.5),lty=4,col="black",lwd=0.8) +
  geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
  geom_text_repel(data = dcx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", vjust = 0.1,  hjust = -1)+
  geom_text_repel(data = dcxx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", hjust = 1, vjust = -.65)+
  labs(x="log2 Fold Change",
       y=bquote(~-Log[10]~italic(eFDR)))  +
  theme_bw() +
  annotate(geom="text", x=2, y=.5, label=paste0(" ", "\n", "in ASC"),
           color="black", size = 5)+
  annotate(geom="text", x=-2, y=.5, label=paste0(" ", "\n", "in DOM"),
           color="black",, size = 5)+
  scale_x_continuous(limits = c(-2.5,2.5),breaks = c(-2.5,-1.5,-0.5,0,0.5,1.5,2.5))+
  theme(axis.text.x = element_text(vjust = 1,size = 20),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),   
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        text = element_text(size = 20)
        
  )


vp_ad


 ggsave("manuscript/brain/results_figures/vplot_ad_70min.png",vp_ad,height =5.5, width =6.5, dpi=600)

### Number of genes in each cut off 
#0.2
dc %>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #382
dc %>% filter(between(logFC, .2, .5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 357
dc %>% filter(between(logFC, .5, 1))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #24
dc %>% filter(between(logFC, 1, 5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #1

dc %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 325
dc %>% filter(between(logFC, -0.5, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #260
dc %>% filter(between(logFC, -1, -0.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 62
dc %>% filter(between(logFC, -3, -1)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 3

#top 10 genes up in asc
dc %>% filter(P.Value < 0.05) %>% arrange(-logFC) %>% head(., 10)
#top 10 genes down in down
dc %>% filter(P.Value < 0.05) %>%  arrange(logFC) %>% head(., 10)
################
head(y2a)

dc <- y2a %>% mutate(contrast = "DES vs. DOM") %>% mutate(log10 = -log10(P.Value))

dc$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
dc$diffexpressed[dc$logFC > 0.2 & dc$P.Value < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
dc$diffexpressed[dc$logFC < -0.2 & dc$P.Value < 0.05] <- "DOWN"


dcx <- dc %>% filter(.,logFC >= 1.6) %>% filter(P.Value < 0.05)
dcxx <- dc %>% filter(.,logFC <= -1.2) %>% filter(P.Value < 0.05)
dc$log10 <- ifelse(dc$log10 == Inf, 4,dc$log10)
# dcx <- dc %>% filter(.,abs(logFC) >= 1.5) %>% filter(P.Value < 0.05) 

vp_dd <- ggplot(data = dc, 
                  aes(x = logFC, 
                      y = log10, 
                      colour=diffexpressed)) +
  geom_point(alpha=0.25, size=3.5) +
  scale_color_manual(values=c("orange", "grey","purple4"))+
  xlim(c(-2.5, 2.5)) +
  ylim(0,4)+
  # geom_vline(xintercept=c(-0.2,0.2),lty=4,col="black",lwd=0.8) +
  geom_vline(xintercept=c(-.5,.5),lty=4,col="black",lwd=0.8) +
  geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
  geom_text_repel(data = dcx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", vjust = 0.65,  hjust = -0.65)+
  geom_text_repel(data = dcxx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", hjust = 1.5, vjust = -.65)+
  labs(x="log2 Fold Change",
       y=bquote(~-Log[10]~italic(eFDR)))  +
  theme_bw() +
  annotate(geom="text", x=2, y=.5, label=paste0(" ", "\n", "in DES"),
           color="black", size = 5)+
  annotate(geom="text", x=-2, y=.5, label=paste0(" ", "\n", "in DOM"),
           color="black",, size = 5)+
  scale_x_continuous(limits = c(-2.5,2.5),breaks = c(-2.5,-1.5,-0.5,0,0.5,1.5,2.5))+
  theme(axis.text.x = element_text(vjust = 1,size = 20),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),   
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        text = element_text(size = 20)
        
  )

vp_dd

ggsave("manuscript/brain/results_figures/vplot_dd_70min.png",vp_dd,height =5.5, width =6.5, dpi=600)

### Number of genes in each cut off 
#0.2
dc %>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #366
dc %>% filter(between(logFC, .2, .5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 268
dc %>% filter(between(logFC, .5, 1.))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #76
dc %>% filter(between(logFC, 1., 3))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #22

dc %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 447
dc %>% filter(between(logFC, -0.5, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #375
dc %>% filter(between(logFC, -1, -0.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 63
dc %>% filter(between(logFC, -3, -1)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #9


#top 10 genes up in des
dc %>% filter(P.Value < 0.05) %>% arrange(-logFC) %>% head(., 10)
#top 10 genes down in des
dc %>% filter(P.Value < 0.05) %>%  arrange(logFC) %>% head(., 10)
################


#################
head(y3a)

dc <- y3a %>% mutate(contrast = "ASC vs. SUB") %>% mutate(log10 = -log10(P.Value))

dc$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
dc$diffexpressed[dc$logFC > 0.2 & dc$P.Value < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
dc$diffexpressed[dc$logFC < -0.2 & dc$P.Value < 0.05] <- "DOWN"


dcx <- dc %>% filter(.,logFC >= 1.25) %>% filter(P.Value < 0.05)
dcxx <- dc %>% filter(.,logFC <= -1) %>% filter(P.Value < 0.05)
dc$log10 <- ifelse(dc$log10 == Inf, 4,dc$log10)
# dcx <- dc %>% filter(.,abs(logFC) >= 1.5) %>% filter(P.Value < 0.05) 

vp_as <- ggplot(data = dc, 
                aes(x = logFC, 
                    y = log10, 
                    colour=diffexpressed)) +
  geom_point(alpha=0.25, size=3.5) +
  scale_color_manual(values=c("orange", "grey","purple4"))+
  xlim(c(-2.5, 2.5)) +
  ylim(0,4)+
  # geom_vline(xintercept=c(-0.2,0.2),lty=4,col="black",lwd=0.8) +
  geom_vline(xintercept=c(-.5,.5),lty=4,col="black",lwd=0.8) +
  geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
  geom_text_repel(data = dcx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", vjust = 0.65,  hjust = -1)+
  geom_text_repel(data = dcxx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", hjust = 1.5, vjust = -1)+
  labs(x="log2 Fold Change",
       y=bquote(~-Log[10]~italic(eFDR)))  +
  theme_bw() +
  annotate(geom="text", x=2, y=.5, label=paste0(" ", "\n", "in ASC"),
           color="black", size = 5)+
  annotate(geom="text", x=-2, y=.5, label=paste0(" ", "\n", "in SUB"),
           color="black",, size = 5)+
  scale_x_continuous(limits = c(-2.5,2.5),breaks = c(-2.5,-1.5,-0.5,0,0.5,1.5,2.5))+
  theme(axis.text.x = element_text(vjust = 1,size = 20),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),   
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        text = element_text(size = 20)
        
  )


vp_as

ggsave("manuscript/brain/results_figures/vplot_as_70min.png",vp_as,height =5.5, width =6.5, dpi=600)

### Number of genes in each cut off 
#0.2
dc %>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #396
dc %>% filter(between(logFC, .2, .5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 327
dc %>% filter(between(logFC, .5, 1))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #61
dc %>% filter(between(logFC, 1, 3.5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #8

dc %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 467
dc %>% filter(between(logFC, -0.5, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #353
dc %>% filter(between(logFC, -1, -0.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 110
dc %>% filter(between(logFC, -3.5, -1.)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 4

#top 10 genes up in dom
dc %>% filter(P.Value < 0.05) %>% arrange(-logFC) %>% head(., 10)
#top 10 genes down in dom 
dc %>% filter(P.Value < 0.05) %>%  arrange(logFC) %>% head(., 10)


########
head(y4a)

dc <- y4a %>% mutate(contrast = "DES vs. SUB") %>% mutate(log10 = -log10(P.Value))

dc$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
dc$diffexpressed[dc$logFC > 0.2 & dc$P.Value < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
dc$diffexpressed[dc$logFC < -0.2 & dc$P.Value < 0.05] <- "DOWN"


dcx <- dc %>% filter(.,logFC >= 1.99) %>% filter(P.Value < 0.05)
dcxx <- dc %>% filter(.,logFC <= -1.5) %>% filter(P.Value < 0.05)
dc$log10 <- ifelse(dc$log10 == Inf, 4,dc$log10)
dc$logFC <- ifelse(dc$logFC >2.5, 2.495, dc$logFC)
# dcx <- dc %>% filter(.,abs(logFC) >= 1.5) %>% filter(P.Value < 0.05) 

vp_ds <-ggplot(data = dc, 
               aes(x = logFC, 
                   y = log10, 
                   colour=diffexpressed)) +
  geom_point(alpha=0.25, size=3.5) +
  scale_color_manual(values=c("orange", "grey","purple4"))+
  xlim(c(-2.5, 2.5)) +
  ylim(0,4)+
  # geom_vline(xintercept=c(-0.2,0.2),lty=4,col="black",lwd=0.8) +
  geom_vline(xintercept=c(-.5,.5),lty=4,col="black",lwd=0.8) +
  geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
  geom_text_repel(data = dcx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", vjust = 0.65,  hjust = -1)+
  geom_text_repel(data = dcxx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", hjust = 1.25, vjust = -1)+
  labs(x="log2 Fold Change",
       y=bquote(~-Log[10]~italic(eFDR)))  +
  theme_bw() +
  annotate(geom="text", x=2, y=.5, label=paste0(" ", "\n", "in DES"),
           color="black", size = 5)+
  annotate(geom="text", x=-2, y=.5, label=paste0(" ", "\n", "in SUB"),
           color="black",, size = 5)+
  scale_x_continuous(limits = c(-2.5,2.5),breaks = c(-2.5,-1.5,-0.5,0,0.5,1.5,2.5))+
  theme(axis.text.x = element_text(vjust = 1,size = 20),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),   
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        text = element_text(size = 20)
        
  )



vp_ds

ggsave("manuscript/brain/results_figures/vplot_ds70min.png",vp_ds,height =5.5, width =6.5, dpi=600)

### Number of genes in each cut off 
#0.2
dc %>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #842
dc %>% filter(between(logFC, .2, .5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 533
dc %>% filter(between(logFC, .5, 1))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #241
dc %>% filter(between(logFC, 1, 3))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #68

dc %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 795
dc %>% filter(between(logFC, -0.5, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #603
dc %>% filter(between(logFC, -1, -0.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 159
dc %>% filter(between(logFC, -3, -1.)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 33

#top 10 genes up in dom
dc %>% filter(P.Value < 0.05) %>% arrange(-logFC) %>% head(., 10)
#top 10 genes down in dom 
dc %>% filter(P.Value < 0.05) %>%  arrange(logFC) %>% head(., 10)
