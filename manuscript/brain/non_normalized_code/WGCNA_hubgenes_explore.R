
#########module number 

datExpr <- readRDS("manuscript/brain/results/WGCNA_datExpr25.RDS") 
net <- readRDS("manuscript/brain/results/WGCNA_net_mPFC25_Power4_minmod75.RDS")

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

moduleNumber = length(unique(moduleColors))

modNames = substring(names(MEs), 3)


moduleColors %>% 
  table() %>% 
  as.data.frame() %>% arrange(Freq)  -> modnum 
colnames(modnum) <- c("module","count")


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


##hub genes: 
hub <- read_csv("manuscript/brain/results_tables/WCGNA_hubgene_list_mPFC25_Power4.csv")
head(hub)

source("manuscript/brain/03_DEGs_GoAnalysis_ASCvsstable.R")
asc <- read.csv("manuscript/brain/results_tables/asc_genes_mPFC.csv")
des <- read.csv("manuscript/brain/results_tables/des_genes_mPFC.csv")

b <- hub %>% filter(moduleName == "blue")
g <- hub %>% filter(moduleName == "green")

g60 <- hub %>% filter(moduleName == "grey60") %>% arrange(-moduleMembership)
p <- hub %>% filter(moduleName == "salmon") %>% arrange(-moduleMembership)

