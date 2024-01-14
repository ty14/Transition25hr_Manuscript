library(WGCNA)
library(tidyverse)
#########module number 

datExpr <- readRDS("manuscript/brain/results_use/WGCNA/WGCNA_datExpr70NORM_RG.RDS") 
net <- readRDS("manuscript/brain/results_use/WGCNA/WGCNA_net_mPFC70Power5.RDS")

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

moduleNumber = length(unique(moduleColors))

modNames = substring(names(MEs), 3)


table(module)

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
coldata <- coldata  %>% filter(time != 25)


MEs <- readRDS("manuscript/brain/results_use/WGCNA/WGCNA_MEs_mPFC70Power5.RDS")

ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleID") %>%
  pivot_longer(cols = 2:26, names_to = "Module") %>% 
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

#des pink(459) , red(648), grey60(262) 

p <-  ME_df %>% filter(Module == "pink") %>% 
  ggplot(aes(condition2, value, fill = condition2, color = condition2))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  # facet_wrap(~Module, ncol =5)+
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+         
  labs(x = "Social condition",
       y = "Module eigengene", title = "Pink: 459 genes")+newggtheme_with_legends+ 
  theme(legend.position = "none", text =element_text(size =15))+
  scale_y_continuous(expand = c(0, 0.1, 0.2, 0))

 p
  ggsave(filename = "manuscript/brain/imgs/WGCNA70_Pink.png",
         p,
         height = 4, width = 4, dpi = 600)
  

  r <- ME_df %>% filter(Module == "red") %>% 
  ggplot(aes(condition2, value, fill = condition2, color = condition2))+
    geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                   alpha = 0.3,
                   jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                   position = position_dodge(0.85)) +
    # facet_wrap(~Module, ncol =5)+
    scale_color_manual(values = viridis::viridis(4))+
    scale_fill_manual(values = viridis::viridis(4))+         
    labs(x = "Social condition",
         y = "Module eigengene", title = "Red: 648 genes")+newggtheme_with_legends+ 
    theme(legend.position = "none", text =element_text(size =15))+
    scale_y_continuous(expand = c(0, 0.1, 0.2, 0))
  
  
  ggsave(filename = "manuscript/brain/imgs/WGCNA70_Red.png",
         r,
         height = 4, width = 4, dpi = 600)

  g <- ME_df %>% filter(Module == "grey60") %>% 
  ggplot(aes(condition2, value, fill = condition2, color = condition2))+
    geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                   alpha = 0.3,
                   jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                   position = position_dodge(0.85)) +
    # facet_wrap(~Module, ncol =5)+
    scale_color_manual(values = viridis::viridis(4))+
    scale_fill_manual(values = viridis::viridis(4))+         
    labs(x = "Social condition",
         y = "Module eigengene", title = "Grey60: 262 genes")+newggtheme_with_legends+ 
    theme(legend.position = "none", text =element_text(size =15))+
    scale_y_continuous(expand = c(0, 0.1, 0.2, 0))
  
  ggsave(filename = "manuscript/brain/imgs/WGCNA70_grey60.png",
         g,
         height = 4, width = 4, dpi = 600)
  
  
#asc blue (1109), turquoise (1268), yellow (878)
  
  b <- ME_df %>% filter(Module == "blue") %>% 
  ggplot(aes(condition2, value, fill = condition2, color = condition2))+
    geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                   alpha = 0.3,
                   jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                   position = position_dodge(0.85)) +
    # facet_wrap(~Module, ncol =5)+
    scale_color_manual(values = viridis::viridis(4))+
    scale_fill_manual(values = viridis::viridis(4))+         
    labs(x = "Social condition",
         y = "Module eigengene", title = "Blue: 1109 genes")+newggtheme_with_legends+ 
    theme(legend.position = "none", text =element_text(size =15))+
    scale_y_continuous(expand = c(0, 0.1, 0.2, 0))
  ggsave(filename = "manuscript/brain/imgs/WGCNA70_Blue.png",
         b,
         height = 4, width = 4, dpi = 600)
  
  t <- ME_df %>% filter(Module == "turquoise") %>% 
  ggplot(aes(condition2, value, fill = condition2, color = condition2))+
    geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                   alpha = 0.3,
                   jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                   position = position_dodge(0.85)) +
    # facet_wrap(~Module, ncol =5)+
    scale_color_manual(values = viridis::viridis(4))+
    scale_fill_manual(values = viridis::viridis(4))+         
    labs(x = "Social condition",
         y = "Module eigengene", title = "Turquoise: 1268 genes")+newggtheme_with_legends+ 
    theme(legend.position = "none", text =element_text(size =15))+
    scale_y_continuous(expand = c(0, 0.1, 0.2, 0))
  ggsave(filename = "manuscript/brain/imgs/WGCNA70_Turquoise.png",
         t,
         height = 4, width = 4, dpi = 600)
  
  y <-  ME_df %>% filter(Module == "yellow") %>% 
  ggplot(aes(condition2, value, fill = condition2, color = condition2))+
    geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                   alpha = 0.3,
                   jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                   position = position_dodge(0.85)) +
    # facet_wrap(~Module, ncol =5)+
    scale_color_manual(values = viridis::viridis(4))+
    scale_fill_manual(values = viridis::viridis(4))+         
    labs(x = "Social condition",
         y = "Module eigengene", title = 'Yellow: 878 genes')+newggtheme_with_legends+ 
    theme(legend.position = "none", text =element_text(size =15))+
    scale_y_continuous(expand = c(0, 0.1, 0.2, 0))
  
  ggsave(filename = "manuscript/brain/imgs/WGCNA70_Yellow.png",
         y,
         height = 4, width = 4, dpi = 600)
  
  
  #estimates
  d<- readRDS("manuscript/brain/results_use/WGCNA/WGCNA_lm_result_mPFC70.RDS")
  d
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
  
  
  
  color_above <- heatmap_df$`Pr(>|t|)` < 0.072
  color_below <- heatmap_df$`Pr(>|t|)` > 0.072
  
  
  heatmap_df$color <- ifelse(heatmap_df$`Pr(>|t|)`< color_above,heatmap_df$`Pr(>|t|)`,.072)
  heatmap_df$key <- factor(heatmap_df$key, levels = c("ASC-DES", "ASC-Same", "DES-Same"))
  heatmap_df$key2 <- ifelse(heatmap_df$key == "ASC-DES", "ASC vs. DES", heatmap_df$key)
  heatmap_df$key2 <- ifelse(heatmap_df$key == "ASC-Same", "ASC vs. SD", heatmap_df$key2)
  heatmap_df$key2 <- ifelse(heatmap_df$key == "DES-Same", "DES vs. SD", heatmap_df$key2)
  heatmap_df$key2 <- factor(heatmap_df$key2, levels = c("ASC vs. DES", "ASC vs. SD", "DES vs. SD"))
 
  
  
  str(heatmap_df)
  p70_dot <- heatmap_df %>% 
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
  
  p70_dot 
  
  ggsave("manuscript/brain/results_figures/WGCNA_pfc70_dotplot.png", p70_dot, width = 7, height = 6)
  
   
  
  
  ###hub genes
  
  hub <- read.csv("manuscript/brain/results_tables/WGCNA/WCGNA_hubgene_list_mPFC70Power5.csv")
colnames(hub)  
hub$moduleName
library(tidyverse)
#red - 17 hub genes 
hub %>% filter(moduleName == "red", moduleMembership > 0.85) %>% arrange(-moduleMembership) -> hr
# ENSMUSG00000040536        0.9360687 -0.3114477   Necab1   4
# 2  ENSMUSG00000032890        0.9224119 -0.3230926    Rims3   4
# 3  ENSMUSG00000069171        0.9145851 -0.4130442    Nr2f1  13
# 4  ENSMUSG00000020173        0.9138274 -0.4277380     Cobl  11
# 5  ENSMUSG00000030067        0.9122423 -0.3222683    Foxp1   6
# 6  ENSMUSG00000021217        0.9082961 -0.3607931    Tshz3   7

dd %>% filter(symbol %in% hr$symbol)
ds %>% filter(symbol %in% hr$symbol)
hr

#grey60 - 7 hub genes
hub %>% filter(moduleName == "grey60", moduleMembership > 0.85) %>% arrange(-moduleMembership) ->hg
# 1 ENSMUSG00000020932        0.9625517 0.3939682   Gfap  11
# 2 ENSMUSG00000028583        0.9005563 0.2436228   Pdpn   4
# 3 ENSMUSG00000026728        0.8860994 0.2300441    Vim   2
# 4 ENSMUSG00000031762        0.8660734 0.3188989    Mt2   8

dd %>% filter(symbol %in% hg$symbol)
ds %>% filter(symbol %in% hg$symbol)

#pink - 12 hub genes
hub %>% filter(moduleName == "pink", moduleMembership > 0.85) %>% arrange(-moduleMembership) -> hp
# 1  ENSMUSG00000033717        0.9422263 0.4250606  Adra2a  19
# 2  ENSMUSG00000031557        0.9344391 0.4850286 Plekha2   8
# 3  ENSMUSG00000044071        0.9135104 0.3010287   Tafa2  10
# 4  ENSMUSG00000047714        0.9120924 0.2117752  Ppp1r2  16
# 5  ENSMUSG00000035168        0.9053383 0.4349292   Tanc1   2
# 6  ENSMUSG00000049001        0.9029729 0.5436078    Ndnf   6

dd %>% filter(symbol %in% hp$symbol)
ds %>% filter(symbol %in% hp$symbol)
#blue - 4 hub genes
hb <- hub %>% filter(moduleName == "blue", moduleMembership > 0.85) %>% arrange(-moduleMembership)
# ensgene moduleMembership   GS.trait symbol chr                              description moduleName
# 1 ENSMUSG00000043467        0.8716860 -0.2503836 Zbtb37   1 zinc finger and BTB domain containing 37       blue
# 2 ENSMUSG00000098905        0.8640126 -0.2303895 Zfp953  13                  zinc finger protein 953       blue
# 3 ENSMUSG00000039153        0.8558831 -0.4232184  Runx2  17      runt related transcription factor 2       blue
# 4 ENSMUSG00000028427        0.8525530 -0.3151227   Aqp7 

adom %>% filter(symbol %in% hb$symbol)
as %>% filter(symbol %in% hb$symbol)

#yellow - 9 hub genes
hy <- hub %>% filter(moduleName == "yellow", moduleMembership > 0.85) %>% arrange(-moduleMembership)
# 1 ENSMUSG00000113004        0.9442418 -0.2598936        Gm47951  13        predicted gene, 47951
# 2 ENSMUSG00000112122        0.9193776 -0.2193416        Gm47628  10        predicted gene, 47628
# 3 ENSMUSG00000090619        0.9104870 -0.2269341        Vmn2r60   7   vomeronasal 2, receptor 60
# 4 ENSMUSG00000069892        0.9044422 -0.2544997 9930111J21Rik2  11 RIKEN cDNA 9930111J21 gene 2
# 5 ENSMUSG00000117148        0.9034256 -0.2106667       Vmn1r229  17   vomeronasal 1 receptor 229

adom %>% filter(symbol %in% hy$symbol)
as %>% filter(symbol %in% hy$symbol)

dd%>% filter(symbol %in% hy$symbol)
ds %>% filter(symbol %in% hy$symbol)

#turquoise- 2 hub genes
ht <- hub %>% filter(moduleName == "turquoise", moduleMembership > 0.85) %>% arrange(-moduleMembership)
# 1 ENSMUSG00000020349        0.8924074 -0.2023527 Ppp2ca  11
# 2 ENSMUSG00000038690        0.8783653  0.2137025 Atp5j2   5


adom %>% filter(symbol %in% ht$symbol)
as %>% filter(symbol %in% ht$symbol)

dd%>% filter(symbol %in% ht$symbol)
ds %>% filter(symbol %in% ht$symbol)
