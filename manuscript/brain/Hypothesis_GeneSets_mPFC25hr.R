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
grcm38 # mouse genes
my_logFC_threshold = 0.1
limma_list<- readRDS("manuscript/brain/results/limma_PFC_ReorganizedGroups_OutlierRemoved.RDS") %>%
  map(~distinct(.)) %>%
  # map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  # map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>%
  map(~filter(.,!is.na(entrez)))

names(limma_list)


y4a <- limma_list$desdom
y3a <- limma_list$dessub
y2a <- limma_list$ascdom
y1a <- limma_list$ascsub
##### neuroendo gene list
# en <- readRDS("manuscript/brain/gene_sets/neuroendo_mouse_geneset.RDS")
# head(en)
# # 1. Neuroendorcine gene set = Donkelaar 2020
# en %>%
#   as_tibble() %>%
#   .$symbol -> enx
# enx_names <- c("Neuroendocrine gene set")
library(msigdbr)
msigdbr(species = "Mus musculus", subcategory ="CP:REACTOME") -> reactome_sets
msigdbr(species = "Mus musculus", category ="H") -> hallmark_sets
msigdbr(species = "Mus musculus", subcategory ='CP:WIKIPATHWAYS' ) -> wp_sets
######2. Energy
# REACTOME_RESPIRATORY_ELECTRON_TRANSPORT
# REACTOME_THE_CITRIC_ACID_TCA_CYCLE_AND_RESPIRATORY_ELECTRON_TRANSPORT
# reactome_sets %>%
#   filter(grepl("REACTOME_RESPIRATORY_ELECTRON_TRANSPORT" ,gs_name,ignore.case = T)) %>%
#   select(gene_symbol) %>%  as_tibble() %>%
#   .$gene_symbol -> enx
# enx_names <- c("Respiratory Electron Transport Chain Gene Set")
#
# [2] "GOBP_ATP_METABOLIC_PROCESS"
# bp_sets %>%
#   filter(grepl("GOBP_ATP_METABOLIC_PROCESS",gs_name,ignore.case = T)) %>%
#   select(gene_symbol) %>%  as_tibble() %>%
#             .$gene_symbol -> enx
######3. metabolic
# reactome_sets %>%
#   filter(grepl("REACTOME_METABOLISM_OF_RNA",gs_name,ignore.case = T)) %>%
#   select(gene_symbol) %>%  as_tibble() %>%
#    .$gene_symbol -> enx
# reactome_sets %>%
# filter(grepl("REACTOME_FATTY_ACID_METABOLISM",gs_name,ignore.case = T)) %>%
# select(gene_symbol) %>%  as_tibble() %>%
#   .$gene_symbol -> enx
# reactome_sets %>%
# filter(grepl("REACTOME_GLYCOLYSIS",gs_name,ignore.case = T)) %>%
#   select(gene_symbol) %>%  as_tibble() %>%
#   .$gene_symbol -> enx
#
######### 4. inflammation
# reactome_sets %>%
# filter(grepl("REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",gs_name,ignore.case = T))  %>%
# select(gene_symbol) %>%  as_tibble() %>%
#  .$gene_symbol -> enx
#
# hallmark_sets %>%
#  filter(grepl("HALLMARK_APOPTOSIS",gs_name,ignore.case = T)) %>%
# select(gene_symbol) %>%  as_tibble() %>%
#   .$gene_symbol -> enx
# wp_sets %>%
#    filter(grepl("WP_NEUROINFLAMMATION",gs_name,ignore.case = T)) %>%
#  select(gene_symbol) %>%  as_tibble() %>%
#    .$gene_symbol -> enx
# reactome_sets %>%
# filter(grepl("REACTOME_ADAPTIVE_IMMUNE_SYSTEM",gs_name,ignore.case = T)) %>%
  # select(gene_symbol) %>%  as_tibble() %>%
  # .$gene_symbol -> enx
# reactome_sets %>%
#   filter(grepl("REACTOME_INNATE_IMMUNE_SYSTEM",gs_name,ignore.case = T)) %>%
#   select(gene_symbol) %>%  as_tibble() %>%
#     .$gene_symbol -> enx
########5.  hormone pathways
# "WP_GLUCOCORTICOID_RECEPTOR_PATHWAY"
# "WP_WNT_SIGNALING"
# "WP_ANDROGEN_RECEPTOR_SIGNALING_PATHWAY"
# "WP_CORTICOTROPINRELEASING_HORMONE_SIGNALING_PATHWAY"
# wp_sets %>%
  # filter(grepl("WP_CORTICOTROPINRELEASING_HORMONE_SIGNALING_PATHWAY" ,gs_name,ignore.case = T))%>%
  # select(gene_symbol) %>%  as_tibble() %>%
  # .$gene_symbol -> enx
################## neurotransmission
# WP_NEUROINFLAMMATION_AND_GLUTAMATERGIC_SIGNALING
# wp_sets %>%
#   filter(grepl("WP_NEUROINFLAMMATION_AND_GLUTAMATERGIC_SIGNALING" ,gs_name,ignore.case = T))%>%
#   select(gene_symbol) %>%  as_tibble() %>%
#   .$gene_symbol -> enx
# GOBP_RESPONSE_TO_DOPAMINE
# GOBP_SYNAPTIC_TRANSMISSION_GABAERGIC
# "GOBP_CATECHOLAMINE_SECRETION
# GOBP_RESPONSE_TO_CATECHOLAMINE
# GOBP_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC
bp_sets %>%
  filter(grepl("GOBP_CATECHOLAMINE_SECRETION",gs_name,ignore.case = T)) %>%
  select(gene_symbol) %>%  as_tibble() %>%
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
  df$logFC < -my_logFC_threshold & df$P.Value<0.05, 'orange',
  ifelse(df$logFC > my_logFC_threshold & df$P.Value<0.05, 'purple4',
         'grey'))

table(keyvals) %>%
  as.data.frame() %>%
  mutate(prop = glue::glue("{round(Freq/sum(Freq),3)*100}")) %>%
  mutate(text = glue::glue("{Freq} genes ({prop}%)")) -> my_texts

if(is.na(my_texts[2,4])){my_texts[2,4] <- "0 genes"}
if(is.na(my_texts[3,4])){my_texts[3,4] <- "0 genes"}

keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'purple4'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == 'orange'] <- 'low'

df %>%
  filter(P.Value <0.05) %>%
  filter(abs(logFC) > my_logFC_threshold) %>%
  arrange(logFC) %>%
  mutate(effects = abs(sin(logFC)*-log10(P.Value))) %>%
  top_n(20) %>%
  # top_frac(.,1) %>%
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
                # drawConnectors = TRUE, # just for a few plots
                # widthConnectors = 0.1, # just for a few plots
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
                pointSize = 2.5,
                labCol = "black",
                labSize = 2.1)+
  annotate("text", x = 0.8, y = 4.4,color = "purple4" , size = 3,
           label =glue::glue("Relatively Upregulated in DES", "\n", "{my_texts[3,4]}"))+
  annotate("text", x = -0.8, y = 4.4,color = "orange" , size = 3,
           label = glue::glue("Relatively Upregulated in DOM","\n", "{my_texts[2,4]}"))+
  scale_x_continuous(limits = c(-1.5,1.5),breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5))+
  scale_y_continuous(limits = c(-0.1,4.8),breaks = c(0,1,2,3,4),expand=expansion(mult=c(0.0005,0.0)))+
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
  df$logFC < -my_logFC_threshold & df$P.Value<0.05, 'orange',
  ifelse(df$logFC > my_logFC_threshold & df$P.Value<0.05, 'purple4',
         'grey'))

table(keyvals) %>%
  as.data.frame() %>%
  mutate(prop = glue::glue("{round(Freq/sum(Freq),3)*100}")) %>%
  mutate(text = glue::glue("{Freq} genes ({prop}%)")) -> my_texts

if(is.na(my_texts[2,4])){my_texts[2,4] <- "0 genes"}
if(is.na(my_texts[3,4])){my_texts[3,4] <- "0 genes"}

keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'purple4'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == 'orange'] <- 'low'
df %>%
  filter(P.Value <0.05) %>%
  filter(abs(logFC) > my_logFC_threshold) %>%
  arrange(logFC) %>%
  mutate(effects = abs(sin(logFC)*-log10(P.Value))) %>%
  top_n(20) %>%
  # top_frac(.,1) %>%
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
                # drawConnectors = TRUE, # just for a few plots
                # widthConnectors = 0.1, # just for a few plots
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
                pointSize = 2.5,
                labCol = "black",
                labSize = 2.1)+
  annotate("text", x = 0.8, y = 4.4,color = "purple4" , size = 3,
           label =glue::glue("Relatively Upregulated in DES", "\n", "{my_texts[3,4]}"))+
  annotate("text", x = -0.8, y = 4.4,color = "orange" , size = 3,
           label = glue::glue("Relatively Upregulated in SUB","\n", "{my_texts[2,4]}"))+
  scale_x_continuous(limits = c(-1.5,1.5),breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5))+
  scale_y_continuous(limits = c(-0.1,4.8),breaks = c(0,1,2,3,4),expand=expansion(mult=c(0.0005,0.0)))+
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
  df$logFC < -my_logFC_threshold & df$P.Value<0.05, 'orange',
  ifelse(df$logFC > my_logFC_threshold & df$P.Value<0.05, 'purple4',
         'grey'))
table(keyvals) %>%
  as.data.frame() %>%
  mutate(prop = glue::glue("{round(Freq/sum(Freq),3)*100}")) %>%
  mutate(text = glue::glue("{Freq} genes ({prop}%)")) -> my_texts
if(is.na(my_texts[2,4])){my_texts[2,4] <- "0 genes"}
if(is.na(my_texts[3,4])){my_texts[3,4] <- "0 genes"}
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'purple4'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == 'orange'] <- 'low'
df %>%
  filter(P.Value <0.05) %>%
  filter(abs(logFC) > my_logFC_threshold) %>%
  arrange(logFC) %>%
  mutate(effects = abs(sin(logFC)*-log10(P.Value))) %>%
  top_n(20) %>%
  # top_frac(.,1) %>%
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
                # drawConnectors = TRUE, # just for a few plots
                # widthConnectors = 0.1, # just for a few plots
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
                pointSize = 2.5,
                labCol = "black",
                labSize = 2.1)+
  annotate("text", x = 0.8, y = 4.4,color = "purple4" , size = 3,
           label =glue::glue("Relatively Upregulated in ASC", "\n", "{my_texts[3,4]}"))+
  annotate("text", x = -0.8, y = 4.4,color = "orange" , size = 3,
           label = glue::glue("Relatively Upregulated in DOM","\n", "{my_texts[2,4]}"))+
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
  df$logFC < -my_logFC_threshold & df$P.Value<0.05, 'orange',
  ifelse(df$logFC > my_logFC_threshold & df$P.Value<0.05, 'purple4',
         'grey'))
table(keyvals) %>%
  as.data.frame() %>%
  mutate(prop = glue::glue("{round(Freq/sum(Freq),3)*100}")) %>%
  mutate(text = glue::glue("{Freq} genes ({prop}%)")) -> my_texts
if(is.na(my_texts[2,4])){my_texts[2,4] <- "0 genes"}
if(is.na(my_texts[3,4])){my_texts[3,4] <- "0 genes"}
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'purple4'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == 'orange'] <- 'low'
df %>%
  filter(P.Value <0.05) %>%
  filter(abs(logFC) > my_logFC_threshold) %>%
  arrange(logFC) %>%
  mutate(effects = abs(sin(logFC)*-log10(P.Value))) %>%
  top_n(20) %>%
  # top_frac(.,1) %>%
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
                # drawConnectors = TRUE, # just for a few plots
                # widthConnectors = 0.1, # just for a few plots
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
                pointSize = 2.5,
                labCol = "black",
                labSize = 2.1)+
  annotate("text", x = 0.8, y = 4.4,color = "purple4" , size = 3,
           label =glue::glue("Relatively Upregulated in ASC", "\n", "{my_texts[3,4]}"))+
  annotate("text", x = -0.8, y = 4.4,color = "orange" , size = 3,
           label = glue::glue("Relatively Upregulated in SUB","\n", "{my_texts[2,4]}"))+
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
##endocrine all
# 1.Neuroendocrine Gene Set
# top<- grid::textGrob("Neuroendocrine Gene Set", gp = grid::gpar(fontsize = 20))
# endo_plot <- gridExtra::grid.arrange(a,b,c,d, ncol =2, top = top)
# ggsave("manuscript/brain/imgs/endocrine_25hr.png",endo_plot,height =10, width =10, dpi=600)
#
#
# 2.REACTOME_RESPIRATORY_ELECTRON_TRANSPORT
# REACTOME_THE_CITRIC_ACID_TCA_CYCLE_AND_RESPIRATORY_ELECTRON_TRANSPORT
 # top<- grid::textGrob("ATP Metabolism", gp = grid::gpar(fontsize = 20))
 # endo_plot <- gridExtra::grid.arrange(a,b,c,d, ncol =2, top = top)
 # ggsave("manuscript/brain/imgs/Energy_25hr.png",endo_plot,height =10, width =10, dpi=600)
# ggsave("manuscript/brain/imgs/CAC_25hr.png",endo_plot,height =10, width =10, dpi=600)
# ggsave("manuscript/brain/imgs/ATP_metabolic_25hr.png",endo_plot,height =10, width =10, dpi=600)


#3 Metabolism
#RNA metabolism
# top<- grid::textGrob("RNA Metabolism Gene Set", gp = grid::gpar(fontsize = 20))
#  endo_plot <- gridExtra::grid.arrange(a,b,c,d, ncol =2, top = top)
#  ggsave("manuscript/brain/imgs/RNAMetabolism_25hr.png",endo_plot,height =10, width =10, dpi=600)
#fatty acid metabolism
# top<- grid::textGrob("Fatty Acid Metabolism Gene Set", gp = grid::gpar(fontsize = 20))
# endo_plot <- gridExtra::grid.arrange(a,b,c,d, ncol =2, top = top)
# ggsave("manuscript/brain/imgs/FattyAcidMetabolism_25hr.png",endo_plot,height =10, width =10, dpi=600)
#
#glycolysis
# top<- grid::textGrob("Glycolysis Gene Set", gp = grid::gpar(fontsize = 20))
# endo_plot <- gridExtra::grid.arrange(a,b,c,d, ncol =2, top = top)
# ggsave("manuscript/brain/imgs/GlycolysisMetabolism_25hr.png",endo_plot,height =10, width =10, dpi=600)


#4. inflammation
#cytokines
# top<- grid::textGrob("Cytokine Signaling Gene Set", gp = grid::gpar(fontsize = 20))
# endo_plot <- gridExtra::grid.arrange(a,b,c,d, ncol =2, top = top)
# ggsave("manuscript/brain/imgs/CytokineSignaling_25hr.png",endo_plot,height =10, width =10, dpi=600)

#inflammation
 # top<- grid::textGrob("Neuro Inflammation Gene Set", gp = grid::gpar(fontsize = 20))
 # endo_plot <- gridExtra::grid.arrange(a,b,c,d, ncol =2, top = top)
 # ggsave("manuscript/brain/imgs/inflammation_25hr.png",endo_plot,height =10, width =10, dpi=600)

# APOPTOSIS
 # top<- grid::textGrob("Apoptosis Gene Set", gp = grid::gpar(fontsize = 20))
 # endo_plot <- gridExtra::grid.arrange(a,b,c,d, ncol =2, top = top)
 # ggsave("manuscript/brain/imgs/apoptosis_25hr.png",endo_plot,height =10, width =10, dpi=600)

#adaptive immune system
# top<- grid::textGrob("Adaptive Immune Gene Set", gp = grid::gpar(fontsize = 20))
# endo_plot <- gridExtra::grid.arrange(a,b,c,d, ncol =2, top = top)
# ggsave("manuscript/brain/imgs/adapt_immune_25hr.png",endo_plot,height =10, width =10, dpi=600)

#innate 
# top<- grid::textGrob("Innate Immune Gene Set", gp = grid::gpar(fontsize = 20))
# endo_plot <- gridExtra::grid.arrange(a,b,c,d, ncol =2, top = top)
# ggsave("manuscript/brain/imgs/Innate_immune_25hr.png",endo_plot,height =10, width =10, dpi=600)


#5. hormones signaling pathways
# # "WP_GLUCOCORTICOID_RECEPTOR_PATHWAY"
# top<- grid::textGrob("Glucocorticoid Receptor Signaling Gene Set", gp = grid::gpar(fontsize = 20))
# endo_plot <- gridExtra::grid.arrange(a,b,c,d, ncol =2, top = top)
# ggsave("manuscript/brain/imgs/GluR_25hr.png",endo_plot,height =10, width =10, dpi=600)

# CORT
 # top<- grid::textGrob("CRHR Signaling Gene Set", gp = grid::gpar(fontsize = 20))
 # endo_plot <- gridExtra::grid.arrange(a,b,c,d, ncol =2, top = top)
 # ggsave("manuscript/brain/imgs/Crhr_25hr.png",endo_plot,height =10, width =10, dpi=600)

#wnt
# top<- grid::textGrob("Wnt Signaling Gene Set", gp = grid::gpar(fontsize = 20))
# endo_plot <- gridExtra::grid.arrange(a,b,c,d, ncol =2, top = top)
# ggsave("manuscript/brain/imgs/Wnt_25hr.png",endo_plot,height =10, width =10, dpi=600)

# "WP_ANDROGEN_RECEPTOR_SIGNALING_PATHWAY"
 # top<- grid::textGrob("AR Signaling Gene Set", gp = grid::gpar(fontsize = 20))
 # endo_plot <- gridExtra::grid.arrange(a,b,c,d, ncol =2, top = top)
 # ggsave("manuscript/brain/imgs/AR_25hr.png",endo_plot,height =10, width =10, dpi=600)


# 6. Neurotransmission: Glu, gaba, ser, dopamine, catecholamines,
# WP_NEUROINFLAMMATION_AND_GLUTAMATERGIC_SIGNALING
# top<- grid::textGrob("Neuroinflammation + GLU Signaling Gene Set", gp = grid::gpar(fontsize = 20))
# endo_plot <- gridExtra::grid.arrange(a,b,c,d, ncol =2, top = top)
# ggsave("manuscript/brain/imgs/inflamm_Glu_25hr.png",endo_plot,height =10, width =10, dpi=600)

# GOBP_RESPONSE_TO_DOPAMINE
# top<- grid::textGrob("Response to Dopamine Gene Set", gp = grid::gpar(fontsize = 20))
# endo_plot <- gridExtra::grid.arrange(a,b,c,d, ncol =2, top = top)
# ggsave("manuscript/brain/imgs/responseDopamine25hr.png",endo_plot,height =10, width =10, dpi=600)

# # GOBP_SYNAPTIC_TRANSMISSION_GABAERGIC
# top<- grid::textGrob("synaptic Transmission GABA Gene Set", gp = grid::gpar(fontsize = 20))
# endo_plot <- gridExtra::grid.arrange(a,b,c,d, ncol =2, top = top)
# ggsave("manuscript/brain/imgs/ST_GABA25hr.png",endo_plot,height =10, width =10, dpi =600)

# GOBP_SYNAPTIC_TRANSMISSION_GLU
# # top<- grid::textGrob("synaptic Transmission GLU Gene Set", gp = grid::gpar(fontsize = 20))
# endo_plot <- gridExtra::grid.arrange(a,b,c,d, ncol =2, top = top)
# ggsave("manuscript/brain/imgs/ST_GLU25hr.png",endo_plot,height =10, width =10, dpi =600)

# GOBP_RESPONSE_TO_CATECHOLAMINE
# top<- grid::textGrob("Response to Catecholamines Gene Set", gp = grid::gpar(fontsize = 20))
# endo_plot <- gridExtra::grid.arrange(a,b,c,d, ncol =2, top = top)
# ggsave("manuscript/brain/imgs/ResponsetoCat0min.png",endo_plot,height =10, width =10, dpi =600)

# "GOBP_CATECHOLAMINE_SECRETION
top<- grid::textGrob("Secretion of Catecholamines", gp = grid::gpar(fontsize = 20))
endo_plot <- gridExtra::grid.arrange(a,b,c,d, ncol =2, top = top)
ggsave("manuscript/brain/imgs/SecretionofCat25hr.png",endo_plot,height =10, width =10, dpi =600)
