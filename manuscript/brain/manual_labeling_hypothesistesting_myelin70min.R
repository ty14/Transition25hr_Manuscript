# libraries
library(limma)
library(edgeR)
library(Mus.musculus)
organism = 'org.Mm.eg.db'
organism1 = 'org.Hs.eg.db'
library(organism, character.only = TRUE)
library(biomaRt)
library(AnnotationDbi)
library(pheatmap)
library(annotables)
library(clusterProfiler)
library(enrichplot)
library(organism, character.only = TRUE)
library(organism1, character.only = TRUE)
library(DOSE)
library(EnhancedVolcano)
library(Orthology.eg.db)
library(tidyverse)
library(msigdbr)
msigdbr_collections() %>% as.data.frame()
#immune =  C7     IMMUNESIGDB
# biological pathways = C2 CP:WIKIPATHWAYS  
#hallmark = H
#reactome = immune, metabolic 
msigdbr(species = "Mus musculus", subcategory ="CP:REACTOME") -> reactome_sets
unique(reactome_sets$gs_name)


grcm38 # mouse genes
my_logFC_threshold = 0.2
limma_list<- readRDS("manuscript/brain/results_use/limma_PFC70min_NORM_RG.RDS") %>%
  map(~distinct(.)) %>%
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>%
  map(~filter(.,!is.na(entrez)))
y4a <- limma_list$desdom
y3a <- limma_list$dessub
y2a <- limma_list$ascdom
y1a <- limma_list$ascsub

# #3. Myelination Regulation 
my <- read_csv("manuscript/brain/gene_sets/myelin.csv")
my %>%  select(gene_symbol) %>%  as_tibble() %>%
 .$gene_symbol -> enx

###############################
head(y4a)

y4a  %>%
  filter(symbol %in% enx) %>%
  dplyr::select(symbol, logFC, P.Value) %>%
  dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                             ifelse(logFC>my_logFC_threshold,"DES genes","DOM genes"))) %>%
  dplyr::mutate(Sig = factor(Sig, levels = c("DES genes","DOM genes", "N.S."))) %>%
  dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
  unique() -> df

keyvals <- ifelse(
  df$logFC < -my_logFC_threshold & df$P.Value<0.05, 'purple4',
  ifelse(df$logFC > my_logFC_threshold & df$P.Value<0.05, 'red',
         'grey'))

table(keyvals) %>%
  as.data.frame() %>%
  mutate(prop = glue::glue("{round(Freq/sum(Freq),3)*100}")) %>%
  mutate(text = glue::glue("{Freq} genes ({prop}%)")) -> my_texts
if(is.na(my_texts[2,4])){my_texts[2,4] <- "0 genes"}
if(is.na(my_texts[3,4])){my_texts[3,4] <- "0 genes"}
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'red'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == 'purple4'] <- 'low'


df %>%
  filter(P.Value <0.05) %>%
  filter(abs(logFC) > my_logFC_threshold) %>%
  arrange(logFC) %>%
  mutate(effects = abs(sin(logFC)*-log10(P.Value))) %>%
  # top_n(10) %>%
  top_frac(.,1) %>%
  # filter(P.Value <0.01 | abs(logFC) >0.5) %>%
  .$symbol -> for_label

EnhancedVolcano(df,
                lab = df$symbol,
                selectLab= for_label,
                x = 'logFC',
                y = 'P.Value',
                title = "",
                pCutoff = 0.05,
                FCcutoff = 0.0,
                cutoffLineType = 'blank',
                drawConnectors = TRUE, # just for a few plots
                widthConnectors = 0.05, # just for a few plots
                vline = c(-0.15, 0.15),
                vlineCol = c('grey90'),
                vlineType = c( 'dashed'),
                vlineWidth = c(0.3),
                hline = c(0.05),
                hlineCol = c('grey90'),
                hlineType = c( 'dashed'),
                hlineWidth = c(0.5),
                colCustom = keyvals,
                shape = 21,
                pointSize = 2,
                labCol = "black",
                labSize = 2.1)+
  annotate("text", x = 0.8, y = 4.4,color = "red" , size = 2,
           label =glue::glue(" DES", "\n", "{my_texts[3,4]}"))+
  annotate("text", x = -0.8, y = 4.4,color = "purple4" , size = 2,
           label = glue::glue(" DOM","\n", "{my_texts[2,4]}"))+
  scale_x_continuous(limits = c(-1.5,1.5),breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5))+
  scale_y_continuous(limits = c(-0.1,4.8),breaks = c(0,1,2,3,4),expand=expansion(mult=c(0.0,0.0)))+
  theme_bw(base_size = 7)+
  labs(color = "",
       caption = paste0('total = ', nrow(df), ' genes'),
       y = bquote(~-Log[10]~italic(eFDR)))+
  theme(legend.position = "none",
        plot.title = element_text(size = 7),
        plot.subtitle = element_blank()) -> a
print(a)
invisible(dev.off())
###

head(y3a)
y3a  %>%
  filter(symbol %in% enx) %>%
  dplyr::select(symbol, logFC, P.Value) %>%
  dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                             ifelse(logFC>my_logFC_threshold,"DES genes","SUB genes"))) %>%
  dplyr::mutate(Sig = factor(Sig, levels = c("DES genes","SUB genes", "N.S."))) %>%
  dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
  unique() -> df
keyvals <- ifelse(
  df$logFC < -my_logFC_threshold & df$P.Value<0.05, 'purple4',
  ifelse(df$logFC > my_logFC_threshold & df$P.Value<0.05, 'red',
         'grey'))
table(keyvals) %>%
  as.data.frame() %>%
  mutate(prop = glue::glue("{round(Freq/sum(Freq),3)*100}")) %>%
  mutate(text = glue::glue("{Freq} genes ({prop}%)")) -> my_texts
if(is.na(my_texts[2,4])){my_texts[2,4] <- "0 genes"}
if(is.na(my_texts[3,4])){my_texts[3,4] <- "0 genes"}
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'red'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == 'purple4'] <- 'low'
df %>%
  filter(P.Value <0.05) %>%
  filter(abs(logFC) > my_logFC_threshold) %>%
  arrange(logFC) %>%
  mutate(effects = abs(sin(logFC)*-log10(P.Value))) %>%
  # top_n(10) %>%
  top_frac(.,.35) %>%
  # filter(P.Value <0.01 | abs(logFC) >0.5) %>%
  .$symbol -> for_label

for_label <-  c( "Bcas1", "Myrf",  "Mog",   "Plp1",  "Ctsc",  "Ugt8a", "Mobp", "Mal", "Cnp", "Sox10",
                             "Pikfyve", "Zfp488", "Myrf") 

EnhancedVolcano(df,
                lab = df$symbol,
                selectLab= for_label,
                x = 'logFC',
                y = 'P.Value',
                title = "",
                pCutoff = 0.05,
                FCcutoff = 0.0,
                cutoffLineType = 'blank',
                drawConnectors = TRUE, # just for a few plots
                widthConnectors = 0.04, # just for a few plots
                vline = c(-0.2, 0.2),
                vlineCol = c('grey90'),
                vlineType = c( 'dashed'),
                vlineWidth = c(0.3),
                hline = c(0.03),
                hlineCol = c('grey90'),
                hlineType = c( 'dashed'),
                hlineWidth = c(0.5),
                colCustom = keyvals,
                shape = 21,
                pointSize = 2,
                labCol = "black",
                labSize = 2.1)+
  annotate("text", x = 0.8, y = 4.4,color = "red" , size = 2,
           label =glue::glue(" DES", "\n", "{my_texts[3,4]}"))+
  annotate("text", x = -0.8, y = 4.4,color = "purple4" , size = 2,
           label = glue::glue(" SUB","\n", "{my_texts[2,4]}"))+
  scale_x_continuous(limits = c(-1.5,1.5),breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5))+
  scale_y_continuous(limits = c(-0.1,4.8),breaks = c(0,1,2,3,4),expand=expansion(mult=c(0.0,0.0)))+
  theme_bw(base_size = 7)+
  labs(color = "",
       caption = paste0('total = ', nrow(df), ' genes'),
       y = bquote(~-Log[10]~italic(eFDR)))+
  theme(legend.position = "none",
        plot.title = element_text(size = 7),
        plot.subtitle = element_blank()) -> b
print(b)
invisible(dev.off())
#####
head(y2a)
y2a  %>%
  filter(symbol %in% enx) %>%
  dplyr::select(symbol, logFC, P.Value) %>%
  dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                             ifelse(logFC>my_logFC_threshold,"ASC genes","DOM genes"))) %>%
  dplyr::mutate(Sig = factor(Sig, levels = c("ASC genes","DOM genes", "N.S."))) %>%
  dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
  unique() -> df
keyvals <- ifelse(
  df$logFC < -my_logFC_threshold & df$P.Value<0.05, 'purple4',
  ifelse(df$logFC > my_logFC_threshold & df$P.Value<0.05, 'red',
         'grey'))
table(keyvals) %>%
  as.data.frame() %>%
  mutate(prop = glue::glue("{round(Freq/sum(Freq),3)*100}")) %>%
  mutate(text = glue::glue("{Freq} genes ({prop}%)")) -> my_texts
if(is.na(my_texts[2,4])){my_texts[2,4] <- "0 genes"}
if(is.na(my_texts[3,4])){my_texts[3,4] <- "0 genes"}
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'red'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == 'purple4'] <- 'low'


df %>%
  filter(P.Value <0.05) %>%
  filter(abs(logFC) > my_logFC_threshold) %>%
  arrange(logFC) %>%
  mutate(effects = abs(sin(logFC)*-log10(P.Value))) %>%
  # top_n(10) %>%
  top_frac(.,1) %>%
  # filter(P.Value <0.01 | abs(logFC) >0.5) %>%
  .$symbol -> for_label


EnhancedVolcano(df,
                lab = df$symbol,
                selectLab= for_label,
                x = 'logFC',
                y = 'P.Value',
                title = "",
                pCutoff = 0.05,
                FCcutoff = 0.0,
                cutoffLineType = 'blank',
                drawConnectors = TRUE, # just for a few plots
                widthConnectors = 0.06, # just for a few plots
                vline = c(-0.2, 0.2),
                vlineCol = c('grey90'),
                vlineType = c( 'dashed'),
                vlineWidth = c(0.3),
                hline = c(0.05),
                hlineCol = c('grey90'),
                hlineType = c( 'dashed'),
                hlineWidth = c(0.5),
                colCustom = keyvals,
                shape = 21,
                pointSize = 2,
                labCol = "black",
                labSize = 2.1)+
  annotate("text", x = 0.8, y = 4.4,color = "red" , size = 2,
           label =glue::glue(" ASC", "\n", "{my_texts[3,4]}"))+
  annotate("text", x = -0.8, y = 4.4,color = "purple4" , size = 2,
           label = glue::glue(" DOM","\n", "{my_texts[2,4]}"))+
  scale_x_continuous(limits = c(-1.5,1.5),breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5))+
  scale_y_continuous(limits = c(-0.1,4.8),breaks = c(0,1,2,3,4),expand=expansion(mult=c(0.0005,0.0)))+
  theme_bw(base_size = 7)+
  labs(color = "",
       caption = paste0('total = ', nrow(df), ' genes'),
       y = bquote(~-Log[10]~italic(eFDR)))+
  theme(legend.position = "none",
        plot.title = element_text(size = 7),
        plot.subtitle = element_blank()) -> c
print(c)
invisible(dev.off())
#
head(y1a)
y1a  %>%
  filter(symbol %in% enx) %>%
  dplyr::select(symbol, logFC, P.Value) %>%
  dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                             ifelse(logFC>my_logFC_threshold,"ASC genes","SUB genes"))) %>%
  dplyr::mutate(Sig = factor(Sig, levels = c("ASC genes","SUB genes", "N.S."))) %>%
  dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>%
  unique() -> df
keyvals <- ifelse(
  df$logFC < -my_logFC_threshold & df$P.Value<0.05, 'purple4',
  ifelse(df$logFC > my_logFC_threshold & df$P.Value<0.05, 'red',
         'grey'))
table(keyvals) %>%
  as.data.frame() %>%
  mutate(prop = glue::glue("{round(Freq/sum(Freq),3)*100}")) %>%
  mutate(text = glue::glue("{Freq} genes ({prop}%)")) -> my_texts
if(is.na(my_texts[2,4])){my_texts[2,4] <- "0 genes"}
if(is.na(my_texts[3,4])){my_texts[3,4] <- "0 genes"}
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'red'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == 'purple4'] <- 'low'
df %>%
  filter(P.Value <0.05) %>%
  filter(abs(logFC) > my_logFC_threshold) %>%
  arrange(logFC) %>%
  mutate(effects = abs(sin(logFC)*-log10(P.Value))) %>%
  # top_n(10) %>%
  top_frac(.,.5) %>%
  # filter(P.Value <0.01 | abs(logFC) >0.5) %>%
  .$symbol -> for_label
df %>% filter(Sig != "N.S.")
for_label <-  c( "Bcas1", "Myrf",  "Mog",   "Plp1",  "Ctsc",  "Ugt8a", "Mobp", "Mal", "Cnp", "Sox10",
                 "Pikfyve", "Zfp488", "Dicer1","Itgb4") 


EnhancedVolcano(df,
                lab = df$symbol,
                selectLab= for_label,
                x = 'logFC',
                y = 'P.Value',
                title = "",
                pCutoff = 0.05,
                FCcutoff = 0.0,
                cutoffLineType = 'blank',
                drawConnectors = TRUE, # just for a few plots
                widthConnectors = 0.03, # just for a few plots
                vline = c(-0.2, 0.2),
                vlineCol = c('grey90'),
                vlineType = c( 'dashed'),
                vlineWidth = c(0.3),
                hline = c(0.05),
                hlineCol = c('grey90'),
                hlineType = c( 'dashed'),
                hlineWidth = c(0.5),
                colCustom = keyvals,
                shape = 21,
                pointSize = 2,
                labCol = "black",
                labSize = 2.1)+
  annotate("text", x = 0.8, y = 4.4,color = "red" , size = 2,
           label =glue::glue(" ASC", "\n", "{my_texts[3,4]}"))+
  annotate("text", x = -0.8, y = 4.4,color = "purple4" , size = 2,
           label = glue::glue(" SUB","\n", "{my_texts[2,4]}"))+
  scale_x_continuous(limits = c(-1.5,1.5),breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5))+
  scale_y_continuous(limits = c(-0.1,4.8),breaks = c(0,1,2,3,4),expand=expansion(mult=c(0.0005,0.0)))+
  theme_bw(base_size = 7)+
  labs(color = "",
       caption = paste0('total = ', nrow(df), ' genes'),
       y = bquote(~-Log[10]~italic(eFDR)))+
  theme(legend.position = "none",
        plot.title = element_text(size = 7),
        plot.subtitle = element_blank()) -> d
print(d)
invisible(dev.off())



# #3. Myelin Regulation 
# top<- grid::textGrob("Myelin Regulation at 70 min", gp = grid::gpar(fontsize = 10))
endo_plot <- gridExtra::grid.arrange(a,c,b,d, ncol =2)
ggsave("manuscript/brain/imgs/MyelinRegulation_70min_use_use2.png",endo_plot,height =5, width =5, dpi=1000)

