library(annotables)
library(tidyverse)
grcm38 # mouse genes



chol <- c('Chrm2', 'Chrm4','Gnai1','Gnai2','Chrm1','Chrm3','Chrm5','Gna11', 'Gna14','Gnaq','Ache',
          'Chat','Grk2', 'Grk5', 'Rgs2', 'Rgs4', 'Rgs6', 'Slc18a3', 'Slc5a7', 'Nat1', 'Lhx8', 'Slc10a4', 
          'Gbx1','Chrna2', 'Chrna3', 'Chrna6', 'Chrna7', 'Chrnb4', 'Chrnb3', 'Agrn', 'Chrna1', 'Chrna10', 
          'Chrna4', 'Chrna5', 'Chrna9', 'Chrnb1', 'Chrnb2', 'Chrnd', 'Chrne', 'Chrng', 'Dok7', 'Lrp4', 'Musk',
          'Rapsn')


my <- c('Cnp', 'Mobp', 'Mbp', 'Myrf', 'Mal', 'Bcas1', 'Mog', 'Mag', 'Lpar1',  'Plp1', 'Tspan2', 'Cntn2', 'Opalin', 'Arhgef10')

# 70 min 
ex <- readRDS("manuscript/brain/results_use/limma_vdl_PFC70min_NormRG.RDS")
head(ex)

x <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

x %>% 
  filter(symbol %in% chol) -> xex
colnames(xex)

xex2 <- xex %>% pivot_longer(cols = 2:29, names_to = "ids")

p70 <- xex2 %>% full_join(id)


p70 <- p70 %>% select(symbol, ids, value, group) %>% mutate(time = "70 min")

#25 

ex <- readRDS("manuscript/brain/results_use/limma_vdl_PFC25hr_NormRG.RDS")
head(ex)

x <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

x %>% 
  filter(symbol %in% chol) -> xex

colnames(xex)

xex2 <- xex %>% pivot_longer(cols = 2:24, names_to = "ids")

p25 <- xex2 %>% full_join(id)


p25 <- p25 %>% select(symbol, ids, value, group) %>% mutate(time = "25 hr")


p <- p70 %>% full_join(p25)


p$group <- substr(p$group, 1,3)

p$group <- factor(p$group, levels = c("ASC","DES", "DOM", "SUB"))
p$time <- factor(p$time, levels = c("70 min", "25 hr" ))


source('functions/geom_boxjitter.R')
library(viridis)

p1 <- p %>% filter(time == "70 min") %>% 
  ggplot(., aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA, alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE, position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+
  facet_wrap(~symbol, scales= "free_y", ncol = 6)+
  ylab("Normalized Expression") +
  xlab("")+
  theme_bw()+
  theme(legend.position = "none", text =element_text(size = 15))+
  scale_y_continuous(expand = c(0, 0, 0.2, 0))
p1

ggsave("manuscript/brain/imgs/ache_genes_boxplots70min.png", height =7, width =12, dpi = 300)


p1 <- p  %>% filter(time == "25 hr") %>%  
  ggplot(., aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA, alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE, position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+
  facet_wrap(~symbol, scales= "free_y", ncol = 6)+
  ylab("Normalized Expression") +
  xlab("")+
  theme_bw()+
  theme(legend.position = "none", text =element_text(size = 15))+
  scale_y_continuous(expand = c(0, 0, 0.2, 0))
p1
ggsave("manuscript/brain/imgs/ache_genes_boxplots25hr.png", height =7, width =12, dpi = 300)


#myelin 
p1 <- ggplot(p, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA, alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE, position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+
  facet_grid(symbol~time, scales= "free")+
  ylab("Normalized Expression") +
  xlab("")+
  theme_bw()+
  theme(legend.position = "none", text =element_text(size = 15))+
  scale_y_continuous(expand = c(0, 0, 0.2, 0))
p1
# ggsave("manuscript/brain/imgs/myl_genes_boxplots.png", height =20, width = 4.5, dpi = 300)