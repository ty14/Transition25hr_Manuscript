library(WGCNA)
library(tidyverse)
#########module number 

datExpr <- readRDS("manuscript/brain/results_use/WGCNA/WGCNA_datExpr25NORM_RG.RDS") 
net <- readRDS("manuscript/brain/results_use/WGCNA/WGCNA_net_mPFC25Power4.RDS")

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

moduleNumber = length(unique(moduleColors))

modNames = substring(names(MEs), 3)


moduleColors %>% 
  table() %>% 
  as.data.frame() %>% arrange(Freq)   -> modnum 
colnames(modnum) <- c("module","count")
modnum <- modnum %>% filter(module !="grey")

modnum$module <- factor(modnum$module,levels=unique(modnum$module[order(modnum$count)]))
mycol <- modnum %>% 
  .$module %>% as.character() 
# for some reason this has to be character again to be ordered properly in the figure...!! 


modnum %>%
  filter(module != "grey") %>% 
  ggplot(aes(y = module, x = count, fill = module))+
  geom_bar(stat = 'identity', color = NA)+
  geom_text(aes(label = count),hjust = -0.2, size = 4)+
  xlim(0, max(modnum$count)+75)+
  scale_fill_manual(values = mycol)+
  theme_minimal(base_size = 7)+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(margin = margin(r=-6)),
        panel.grid = element_blank())+
  labs(x = "",
       y = "",title = "Module Size")+
  theme_minimal(base_size = 15)+
  theme(legend.position = "none", 
        legend.title = element_text(size =rel(1.2)),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size =rel(1.5)),
        legend.key.size = unit(.9, 'cm'),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_blank(),
        # panel.background = element_blank(),
        # panel.grid.minor = element_line(color = 'grey50'),
        panel.grid = element_blank(),
        axis.text.y = element_text(hjust=1,vjust=0.25,size = rel(1.25))) -> temp_p
temp_p


########boxplots
coldata <- read.csv("manuscript/brain/non_normalized_code/results_tables/coldata_ALLXAGG.csv")
# get rid of controls
coldata <- coldata  %>% filter(time != 70)


MEs <- readRDS("manuscript/brain/results_use/WGCNA/WGCNA_MEs_mPFC25Power4.RDS")

ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleID") %>%
  pivot_longer(cols = 2:17, names_to = "Module") %>% 
  full_join(coldata) %>% filter(Module != "MEgrey") %>% 
  mutate(status = condition1)

head(ME_df)

ME_df$Module <- gsub("ME", "", ME_df$Module)


colnames(ME_df)
ME_df$condition2 <- ifelse(ME_df$condition1 == "SUB", "SD",ME_df$condition1)
ME_df$condition2 <- ifelse(ME_df$condition2 == "DOM", "SD",ME_df$condition2)

#boxplot of eigengenes and status 
source("functions/geom_boxjitter.R")
source("functions/newggtheme.R")

# purple midnightblue
#salmon, pink is 0.6

p <-  ME_df %>% filter(Module == "purple") %>% 
  ggplot(aes(condition2, value, fill = condition2, color = condition2))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  # facet_wrap(~Module, ncol =5)+
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+         
  labs(x = "Social condition",
       y = "Module eigengene", title = "Purple: 368 genes")+newggtheme_with_legends+ 
  theme(legend.position = "none", text =element_text(size =10))


ggsave(filename = "manuscript/brain/imgs/WGCNA25_Purple.png",
       p,
       height = 4, width = 4, dpi = 600)


r <- ME_df %>% filter(Module == "midnightblue") %>% 
  ggplot(aes(condition2, value, fill = condition2, color = condition2))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  # facet_wrap(~Module, ncol =5)+
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+         
  labs(x = "Social condition",
       y = "Module eigengene", title = "Midnightblue: 129 genes")+newggtheme_with_legends+ 
  theme(legend.position = "none", text =element_text(size =10))

ggsave(filename = "manuscript/brain/imgs/WGCNA25_midblue.png",
       r,
       height = 4, width = 4, dpi = 600)

g <- ME_df %>% filter(Module == "pink") %>% 
  ggplot(aes(condition2, value, fill = condition2, color = condition2))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  # facet_wrap(~Module, ncol =5)+
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+         
  labs(x = "Social condition",
       y = "Module eigengene", title = "Pink: 433 genes")+newggtheme_with_legends+ 
  theme(legend.position = "none", text =element_text(size =10))

ggsave(filename = "manuscript/brain/imgs/WGCNA25_pink.png",
       g,
       height = 4, width = 4, dpi = 600)


s <- ME_df %>% filter(Module == "salmon") %>% 
  ggplot(aes(condition2, value, fill = condition2, color = condition2))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  # facet_wrap(~Module, ncol =5)+
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+         
  labs(x = "Social condition",
       y = "Module eigengene", title = "Salmon: 278 genes")+newggtheme_with_legends+ 
  theme(legend.position = "none", text =element_text(size =10))

ggsave(filename = "manuscript/brain/imgs/WGCNA25_salmon.png",
       s,
       height = 4, width = 4, dpi = 600)

#estimates
d<- readRDS("manuscript/brain/results_use/WGCNA/WGCNA_lm_result_mPFC25.RDS")

head(d)
colnames(d)
d$`Pr(>|t|)` <- as.numeric(d$`Pr(>|t|)`)
d %>% 
  ggplot(aes(y = module, x =`Pr(>|t|)`)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")+
  geom_point(size = 3.5)+
  facet_wrap(~key)




moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

moduleNumber = length(unique(moduleColors))

modNames = substring(names(MEs), 3)
heatmap_df <- d %>% 
  mutate(my_alpha = ifelse(`Pr(>|t|)` < 0.05, 1, 0)) %>% 
  mutate(my_alpha2 = Estimate )

moduleColors %>% 
  table() %>% 
  as.data.frame() %>% arrange(Freq)  -> modnum 
colnames(modnum) <- c("module","count")

heatmap_df %>% as_tibble() %>% 
  left_join(modnum) -> heatmap_dfx


str(heatmap_df)

heatmap_df$Estimate <- as.numeric(heatmap_df$Estimate)
heatmap_df$`Std. Error` <- as.numeric(heatmap_df$`Std. Error`)
heatmap_df$module <- as.factor(heatmap_df$module)



color_above <- heatmap_df$`Pr(>|t|)` < 0.075
color_below <- heatmap_df$`Pr(>|t|)` > 0.075


heatmap_df$color <- ifelse(heatmap_df$`Pr(>|t|)`< color_above,heatmap_df$`Pr(>|t|)`,.075)
heatmap_df$key <- factor(heatmap_df$key, levels = c("ASC-DES", "ASC-Same", "DES-Same"))
heatmap_df$key2 <- ifelse(heatmap_df$key == "ASC-DES", "ASC vs. DES", heatmap_df$key)
heatmap_df$key2 <- ifelse(heatmap_df$key == "ASC-Same", "ASC vs. SD", heatmap_df$key2)
heatmap_df$key2 <- ifelse(heatmap_df$key == "DES-Same", "DES vs. SD", heatmap_df$key2)
heatmap_df$key2 <- factor(heatmap_df$key2, levels = c("ASC vs. DES", "ASC vs. SD", "DES vs. SD"))



str(heatmap_df)
p25_dot <- heatmap_df %>% 
  ggplot(aes(y = module, x = Estimate, color =color)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")+
  geom_point(size = 3)+
  facet_wrap(~key2)+
  labs(y="", color = "p-value")+
  xlim(-2,2)+
  scale_color_continuous(low = "red", high = "lightgray", breaks = c(0.02,0.04, 0.06, 0.08))+
  geom_errorbar(aes(xmin = Estimate-`Std. Error`, xmax = Estimate+`Std. Error`),width = 0.2)+
  scale_y_discrete(limits = rev(levels(heatmap_df$module)))+
  theme_bw()+
  theme(text = element_text(size = 15))

p25_dot 

ggsave("manuscript/brain/results_figures/WGCNA_pfc25_dotplot.png", p25_dot, width = 7, height = 6)





hub <- read.csv("manuscript/brain/results_tables/WGCNA/WCGNA_hubgene_list_mPFC25Power4.csv")
hub %>% filter(moduleName == "purple", moduleMembership > 0.8) %>% arrange(-moduleMembership) -> hp
hp

asx %>% filter(symbol %in% hp$symbol)
adomx %>% filter(symbol %in% hp$symbol)


d70 <- readRDS("manuscript/brain/results_use/KIN/kIN_dataframe_mPFC_power5.RDS")
colnames(d70)


d25 <- readRDS("manuscript/brain/results_use/KIN/kIN_dataframe_mPFC25_power4.RDS")


