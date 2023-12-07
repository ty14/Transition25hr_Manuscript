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

limma_list<- readRDS("manuscript/brain/results_use/limma_PFC25hr_Norm_RG.RDS") %>% 
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

dc$logFC <- ifelse(dc$logFC> 2.5, 2.5, dc$logFC)
dcx <- dc %>% filter(.,logFC >= 1.1)%>% filter(P.Value < 0.05)
dcxx <- dc %>% filter(.,logFC <= -0.75) %>% filter(P.Value < 0.05)
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
  geom_text_repel(data = dcx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", vjust = 0.1,  hjust = -1.2)+
  geom_text_repel(data = dcxx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", hjust = 0., vjust = 0)+
  labs(x="log2 Fold Change",
       y=bquote(~-Log[10]~italic(eFDR)))  +
  theme_bw() +
  annotate(geom="text", x=2.15, y=.5, label=paste0(" ", "\n", "ASC vs. DOM"),
           color="black", size = 5)+
  annotate(geom="text", x=-2.15, y=.5, label=paste0(" ", "\n", "ASC vs. DOM"),
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


ggsave("manuscript/brain/results_figures/vplot_ad_25hr.png",vp_ad,height =5.5, width =6.5, dpi=600)

### Number of genes in each cut off 
#0.2
dc %>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #252
dc %>% filter(between(logFC, .2, .5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 175
dc %>% filter(between(logFC, .5, 1))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #70
dc %>% filter(between(logFC, 1, 5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #7

dc %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 241
dc %>% filter(between(logFC, -0.5, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #219
dc %>% filter(between(logFC, -1, -0.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 22
dc %>% filter(between(logFC, -3, -1)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 0

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


dcx <- dc %>% filter(.,logFC >= 1.) %>% filter(P.Value < 0.05)
dcxx <- dc %>% filter(.,logFC <= -0.8) %>% filter(P.Value < 0.05)
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
  geom_text_repel(data = dcx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", vjust = 1,  hjust = -1)+
  geom_text_repel(data = dcxx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", hjust = 1.11, vjust = -.5)+
  labs(x="log2 Fold Change",
       y=bquote(~-Log[10]~italic(eFDR)))  +
  theme_bw() +
  annotate(geom="text", x=2.15, y=.5, label=paste0(" ", "\n", "DES vs. DOM"),
           color="black", size = 5)+
  annotate(geom="text", x=-2.15, y=.5, label=paste0(" ", "\n", "DES vs. DOM"),
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

ggsave("manuscript/brain/results_figures/vplot_dd_25hr.png",vp_dd,height =5.5, width =6.5, dpi=600)

### Number of genes in each cut off 
#0.2
dc %>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #216
dc %>% filter(between(logFC, .2, .5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 148
dc %>% filter(between(logFC, .5, 1.))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #63
dc %>% filter(between(logFC, 1., 3))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #5

dc %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 169
dc %>% filter(between(logFC, -0.5, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #137
dc %>% filter(between(logFC, -1, -0.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 30
dc %>% filter(between(logFC, -3, -1)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #2


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


dcx <- dc %>% filter(.,logFC >= 1) %>% filter(P.Value < 0.05)
dcxx <- dc %>% filter(.,logFC <= -0.85) %>% filter(P.Value < 0.05)
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
  geom_text_repel(data = dcx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", vjust = 0.8,  hjust = -0.2)+
  geom_text_repel(data = dcxx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", hjust = 0, vjust = 0.1)+
  labs(x="log2 Fold Change",
       y=bquote(~-Log[10]~italic(eFDR)))  +
  theme_bw() +
  annotate(geom="text", x=2.15, y=.5, label=paste0(" ", "\n", "ASC vs. SUB"),
           color="black", size = 5)+
  annotate(geom="text", x=-2.15, y=.5, label=paste0(" ", "\n", "ASC vs. SUB"),
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

ggsave("manuscript/brain/results_figures/vplot_as_25hr.png",vp_as,height =5.5, width =6.5, dpi=600)

### Number of genes in each cut off 
#0.2
dc %>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #331
dc %>% filter(between(logFC, .2, .5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 267
dc %>% filter(between(logFC, .5, 1))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #59
dc %>% filter(between(logFC, 1, 3.5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #5

dc %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 299
dc %>% filter(between(logFC, -0.5, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #245
dc %>% filter(between(logFC, -1, -0.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 53
dc %>% filter(between(logFC, -3.5, -1.)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 1

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

dc$logFC <- ifelse(dc$logFC >2.5, 2.495, dc$logFC)
dcx <- dc %>% filter(.,logFC >= 1.6) %>% filter(P.Value < 0.05)
dcxx <- dc %>% filter(.,logFC <= -1.4) %>% filter(P.Value < 0.05)
dc$log10 <- ifelse(dc$log10 == Inf, 4,dc$log10)

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
  geom_text_repel(data = dcx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", vjust = 0.2,  hjust = -1)+
  geom_text_repel(data = dcxx, aes(x = logFC, y = -log10(P.Value),label = symbol), color = "black", hjust = 0.25, vjust = -1)+
  labs(x="log2 Fold Change",
       y=bquote(~-Log[10]~italic(eFDR)))  +
  theme_bw() +
  annotate(geom="text", x=2.15, y=.5, label=paste0(" ", "\n", "DES vs. SUB"),
           color="black", size = 5)+
  annotate(geom="text", x=-2.15, y=.5, label=paste0(" ", "\n", "DES vs. SUB"),
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

ggsave("manuscript/brain/results_figures/vplot_ds25hr.png",vp_ds,height =5.5, width =6.5, dpi=600)

### Number of genes in each cut off 
#0.2
dc %>% filter(.,logFC >= .2)%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #452
dc %>% filter(between(logFC, .2, .5))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) # 354
dc %>% filter(between(logFC, .5, 1))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #81
dc %>% filter(between(logFC, 1, 3))%>% filter(P.Value < 0.05) %>% summarise(n = nrow(.)) #17

dc %>% filter(.,logFC <= -.2) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 366
dc %>% filter(between(logFC, -0.5, -0.2)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) #262
dc %>% filter(between(logFC, -1, -0.5)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.))# 88
dc %>% filter(between(logFC, -3, -1.)) %>% filter(P.Value < 0.05)%>% summarise(n = nrow(.)) # 16

#top 10 genes up in dom
dc %>% filter(P.Value < 0.05) %>% arrange(-logFC) %>% head(., 10)
#top 10 genes down in dom 
dc %>% filter(P.Value < 0.05) %>%  arrange(logFC) %>% head(., 10)
