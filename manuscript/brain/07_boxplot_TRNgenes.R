#boxplots for transition genes 


# 70 min 
ex <- readRDS("manuscript/brain/results_use/limma_vdl_PFC70min_NormRG.RDS")
head(ex)

x <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

table(id$group)

trans <- c("Kcns3", "Crybg2", "Stx19", "Cmah", "Cndp1", "Dhx16", "Zfp280d", 'Chka', 'Kmt2a', 'Plekha6',"Crtac1", "Itga9", "D630023F18Rik", "Tmem147", "Slc7a4", "Egfl7", "Kcnip3", "Uqcc2")
x %>% 
  filter(symbol %in% trans) -> xex
colnames(xex)

xex2 <- xex %>% pivot_longer(cols = 2:29, names_to = "ids")

p <- xex2 %>% full_join(id)

p$group <- substr(p$group, 1,3)

p$group <- factor(p$group, levels = c("ASC","DES", "DOM", "SUB"))
# p <- p %>% na.omit(.)

source('functions/geom_boxjitter.R')
library(viridis)

p1 <- ggplot(p, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA, alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE, position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+
  facet_wrap(factor(symbol,levels = c("Kcns3", "Crybg2", "Stx19", "Cmah", "Cndp1", "Dhx16", "Zfp280d", 'Chka', 'Kmt2a','Plekha6',
                                      "Crtac1", "Itga9", "D630023F18Rik", "Tmem147", "Slc7a4", "Egfl7", "Kcnip3", "Uqcc2")) ~ ., scales = 'free', ncol =11)+
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("")+
  theme_bw()+
  theme(legend.position = "none", text =element_text(size = 15))
p1

ggsave("manuscript/brain/imgs/TRN70_genes_boxplots.png", width =22 , height = 4.5, dpi = 300)

p70_down <- p %>% filter(symbol %in% c("Kcns3", "Crybg2", "Stx19", "Cmah", "Cndp1", "Dhx16", "Zfp280d", 'Chka', 'Kmt2a', 'Plekha6'))
p70_up <- p %>% filter(symbol %in% c("Crtac1", "Itga9", "D630023F18Rik", "Tmem147", "Slc7a4", "Egfl7", "Kcnip3", "Uqcc2"))


table(p70_up$group)

#up 
p1 <- ggplot(p70_up, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA, alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE, position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+
  facet_wrap(~symbol, scales = 'free', ncol =4)+
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("")+
  theme_bw()+
  theme(legend.position = "none", text =element_text(size = 15))
p1
print(p1)
ggsave("manuscript/brain/imgs/TRN70_UPreg_genes_boxplots.png",p1, width =9 , height = 4.5, dpi = 300)


#down
p1 <- ggplot(p70_down, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA, alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE, position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+
  facet_wrap(~symbol, scales = 'free', ncol =5)+
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("")+
  theme_bw()+
  theme(legend.position = "none", text =element_text(size = 15))
p1

ggsave("manuscript/brain/imgs/TRN70_DOWNreg_genes_boxplots.png", width =10 , height = 4.25, dpi = 300)


#25hr
ex <- readRDS("manuscript/brain/results_use/limma_vdl_PFC25hr_NormRG.RDS")
head(ex)

x <- ex$E %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) 

id <- ex$targets%>% as.data.frame() %>% rownames_to_column('ids') %>%  dplyr::select(ids, group)

trans <- c("Exosc9", "Gstm7", "Zscan18", "Synm", "Fam117b","Esrrg", "Fam169b", "Cdc23", "Mief1", "Nfkb1" )
x %>% 
  filter(symbol %in% trans) -> xex

colnames(xex)

xex2 <- xex %>% pivot_longer(cols = 2:24, names_to = "ids")

p <- xex2 %>% full_join(id)
p$group <- substr(p$group, 1,3)
p$group <- factor(p$group, levels = c("ASC","DES", "DOM", "SUB"))
# p <- p %>% na.omit(.)

table(p$group)
source('functions/geom_boxjitter.R')
library(viridis)


p1 <- ggplot(p, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA, alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE, position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+
  facet_wrap(factor(symbol,levels = c("Esrrg", "Fam169b", "Cdc23", "Mief1", "Nfkb1","Exosc9", "Gstm7", "Zscan18", "Synm", "Fam117b")) ~ ., scales = 'free', ncol =5)+
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("")+
  theme_bw()+
  theme(legend.position = "none", text =element_text(size = 15))
p1

ggsave("manuscript/brain/imgs/TRN25_genes_boxplots.png", width =10 , height = 4.25, dpi = 300)


p25_up <- p %>% filter(symbol %in% c("Esrrg", "Fam169b", "Cdc23", "Mief1", "Nfkb1"))
p25_down<- p %>%  filter(symbol %in% c("Exosc9", "Gstm7", "Zscan18", "Synm", "Fam117b"))


p1 <- ggplot(p25_up, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA, alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE, position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+
  facet_wrap(~symbol, scales = 'free', ncol =5)+
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("")+
  theme_bw()+
  theme(legend.position = "none", text =element_text(size = 15))
p1

ggsave("manuscript/brain/imgs/TRN25_UPreg_genes_boxplots.png", width =10 , height = 2.5, dpi = 300)


p1 <- ggplot(p25_down, aes(group,value, color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA, alpha = 0.4, jitter.size = 2,, size = .8,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE, position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+
  facet_wrap(~symbol, scales = 'free', ncol =5)+
  scale_y_continuous(expand = c(0, 1))+
  ylab("Normalized Expression") +
  xlab("")+
  theme_bw()+
  theme(legend.position = "none", text =element_text(size = 15))
p1

ggsave("manuscript/brain/imgs/TRN25_DOWNreg_genes_boxplots.png", width =10 , height = 2.5, dpi = 300)
