# libraries 
library(limma)
library(edgeR)
library(Mus.musculus)
organism = 'org.Mm.eg.db'
library(organism, character.only = TRUE)
library(biomaRt)
library(AnnotationDbi)
library(annotables)
grcm38 # mouse genes
library(WGCNA)
library(tidyverse)


#Getting metadata ready 
coldata <- read_csv("brain/sample25hr_table.csv")
head(coldata)
str(coldata)

#group thing 
coldata$groupEX <- coldata$group
# Normalizing cort data
# df <- transform(df, N = (N - min(N)) / (max(N) - min(N))

coldata <-coldata %>% mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))
coldata <- coldata %>%  dplyr::select(-group, -period)
# coldata$post_Ncort[is.na(coldata$post_Ncort)] = mean(post_Ncort)

#getting condition1
table(coldata$condition)
# CDOM, RDOM to Descenders (DOM to SUB)(4->1)
# CSUB, SUB to Ascenders (Sub to DOM)  (1->4)

coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "DOM", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "SUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "DES", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "ASC", coldata$condition1)

coldata$condition1  <- coldata$condition1  %>% replace_na('SUB')
coldata$condition1  <- paste(coldata$condition1, coldata$time)

coldata <- coldata %>% 
  filter(region == "P")

coldata$Sampleid <- substr(coldata$SampleNames,6,12)

coldata25 <- coldata %>% select(condition1, Sampleid) 

#remove outlier
#b1.2.1 outlier DES
pfc_data25 <- coldata25 %>% filter(Sampleid != 'b1.2.1.')

#now 70 min data 
coldata <- read_csv("brain/sample70min_table.csv")
head(coldata)
str(coldata)

#group thing 
coldata$groupEX <- coldata$group
# Normalizing cort data
# df <- transform(df, N = (N - min(N)) / (max(N) - min(N))

coldata <-coldata %>% mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))
coldata <- coldata %>%  dplyr::select(-group, -period)
# coldata$post_Ncort[is.na(coldata$post_Ncort)] = mean(post_Ncort)

#getting condition1
table(coldata$condition)
# CDOM, RDOM to Descenders (DOM to SUB)(4->1)
# CSUB, SUB to Ascenders (Sub to DOM)  (1->4)

coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "DOM", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "SUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "DES", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "ASC", coldata$condition1)

coldata$condition1  <- coldata$condition1  %>% replace_na('SUB')
coldata$condition1 <- paste(coldata$condition1, coldata$time)
pfc_data1<- coldata %>% 
  filter(region != "AMY") %>% 
  filter(Postrank != 3) %>% 
  filter(condition1 != "control 1 hr") %>% 
  filter(condition1 != "ascenders 1 hr")

table(pfc_data1$condition1)
pfc_data1$Sampleid <- substr(pfc_data1$SampleName, 7,13)
pfc_data70 <- pfc_data1 %>% select(condition1, Sampleid)


pfc_data <- pfc_data70 %>% rbind(pfc_data25) 
pfc_data <- pfc_data %>% column_to_rownames(., var = "Sampleid")

### 70 min data 
# Expression values
dlNorm70 <-  read.csv("brain/PFC_counts.csv", row.names = 1)
#remove zeros
dlNorm70 <- dlNorm70[apply(dlNorm70[], 1, function(x) !all(x==0)),]
#trim sample ids
colnames(dlNorm70)[c(1:67)] <- substr(colnames(dlNorm70)[c(1:67)], 7, 13)


# Bring in count data for mPFC 25 hr 
dlNorm25 <-  read.csv("brain/PFC_25hcounts.csv", row.names = 1)
#remove zeros
dlNorm25 <- dlNorm25[apply(dlNorm25[], 1, function(x) !all(x==0)),]
#trim sample ids
colnames(dlNorm25)[c(1:24)] <- substr(colnames(dlNorm25)[c(1:24)], 6, 12)
#remove outlier B1.1.2.
dlNorm25 <- dlNorm25[,c(1:19,21:24)]


#combinf counts to normalize together 
dl70 <- dlNorm70 %>% rownames_to_column(., var = "gene")
dl25 <- dlNorm25 %>% rownames_to_column(., var = "gene")

dlNorm <- dl70 %>% full_join(dl25)
dlNorm <- dlNorm %>% column_to_rownames(., var= "gene")


colnames(dlNorm)
rownames(pfc_data)

#check before normalizing 
dlNorm<- dlNorm[ ,rownames(pfc_data)]
all(rownames(pfc_data) == colnames(dlNorm))

#normalize and filter with all groups 
# row.names(dlNorm) <- nrows
dlNorm <- dlNorm[!is.na(rowSums(dlNorm)),]

d = apply(dlNorm, 2, as.numeric)
dim(d)

d0= DGEList(d, group = pfc_data$condition1)
dim(d0)
rownames(d0) <- rownames(dlNorm)
d0 <- calcNormFactors(d0)

cutoff <- 10
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
# [1] 13367    51
# Now take out groups that you want
#70min first 
dge.dl$samples$group

dge.dl_dom <- dge.dl[, dge.dl$samples$group %in% c("ASC 1 hr", "DES 1 hr", "DOM 1 hr", "SUB 1 hr")]
dge.dl_dom$samples$group <- droplevels(dge.dl_dom$samples$group)
dge.dl_dom$samples$group
dge.dl<- dge.dl_dom
dge.dl$samples$group

pfc_data %>% 
  filter( condition1 %in% c("ASC 1 hr", "DES 1 hr", "DOM 1 hr", "SUB 1 hr") ) -> var_info  

dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check

var_info$condition1 %>% 
  factor(.,levels = c("DOM 1 hr","ASC 1 hr","DES 1 hr", "SUB 1 hr")) -> group.dl

group.dl<- gsub(" 1 hr", "", group.dl)

design.dl <- model.matrix(~ 0 + group.dl)
colnames(design.dl) -> mycolnames

v.dl = voom(dge.dl, design.dl, plot = F)

# Many functions expect the matrix to be transposed
datExpr <- t(v.dl$E)
## check rows/cols
nrow(datExpr)
ncol(datExpr)
rownames(datExpr)

# getting trait data
df <- read_csv("manuscript/brain/non_normalized_code/results_tables/coldata_ALLXAGG.csv")
colnames(df)
pfc <- df %>% filter(time == 70) %>% 
  select(batch,Postds, Preds, wt_d8, CORT, post.given1, post.received1, condition1, SampleID)
pfc <- pfc[c(1:25,27,28,29),]

ifelse(pfc$condition1 == "DOM", 3, pfc$condition1) -> pfc$condition2
ifelse(pfc$condition2 == "DES", 2, pfc$condition2) -> pfc$condition2
ifelse(pfc$condition2 == "ASC", 1, pfc$condition2) -> pfc$condition2
ifelse(pfc$condition2 == "SUB", 0, pfc$condition2) -> pfc$condition2


pfc1 <- pfc %>% 
  dplyr::select(-condition1, -SampleID)

#make everything numeric 
pfc1[1:8] <- lapply(pfc1[1:8], as.numeric)
table(pfc1$condition2)

sampleTree = hclust(dist(datExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)

par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 60, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 60, minSize = 10)
table(clust)
#outliers: B10.1.2.: ASC and B12.4.3: DES
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExprx = datExpr[keepSamples, ]
nGenes = ncol(datExprx)
nSamples = nrow(datExprx)


# Re-cluster samples with out any samples filtered 
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(pfc1, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(pfc1),
                    main = "mPFC:TRN dendrogram and trait heatmap")

collectGarbage()
#not going to filter because ascenders rreally are not grouping anyways?
saveRDS(datExpr,"manuscript/brain/results_use/WGCNA/WGCNA_datExpr70NORM_RG.RDS")



# Run WGCNA for each dom, descender, controls 
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, )
# Plot the results:
print(sft)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;


# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main ="mPFC 70 min: Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");


# this line corresponds to using an R^2 cut-off of h

abline(h=0.85,col="red")  #changed to 0.8 (originally 0.9 in the tutorial)
abline(h=0.90,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = "mPFC 70 min:Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
print(sft)



# TOMType = "signed", networkType = "signed",
# TOMType = "unsigned", networkType = "unsigned",
# TOMType = c("signed", networkType = "signed hybrid")
cor <- WGCNA::cor
WGCNA_get_net <- function(my_tissue = "mPFC70min",
                          my_power =  5,
                          my_TOMType ="signed", 
                          my_networkType = "signed hybrid"){
  
  x <- readRDS("manuscript/brain/results_use/WGCNA/WGCNA_datExpr70NORM_RG.RDS") 
  set.seed(312)
  net = blockwiseModules(x, 
                         power = my_power,
                         TOMType = my_TOMType, 
                         networkType = my_networkType,
                         minModuleSize = 75,
                         reassignThreshold = 0, 
                         mergeCutHeight = 0.25,
                         numericLabels = TRUE, 
                         pamRespectsDendro = FALSE,
                         saveTOMs = FALSE,
                         verbose = 3)
  
  
  saveRDS(net, "manuscript/brain/results_use/WGCNA/WGCNA_net_mPFC70Power5.RDS")
  
}

WGCNA_get_net("mPFC70min", 5, "signed", "signed hybrid")
# ========================================================================================

net <- readRDS(glue::glue("manuscript/brain/results_use/WGCNA/WGCNA_net_mPFC70Power5.RDS"))

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath

dev.off() # make sure you do this before AND after 
png(file = "manuscript/brain/imgs/cluster_dendo_mPFC70_Power5_minmod75.png",
    width=600, height=350)
plotDendroAndColors(net$dendrograms[[2]], mergedColors[net$blockGenes[[2]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "mPFC 70min:Power5")

dev.off()

MEs = net$MEs



##################################
set.seed(312)
datExpr <- readRDS("manuscript/brain/results_use/WGCNA/WGCNA_datExpr70NORM_RG.RDS") 
net <- readRDS("manuscript/brain/results_use/WGCNA/WGCNA_net_mPFC70Power5.RDS")
dim(datExpr)


#PART 1
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
moduleNumber = length(unique(moduleColors))
rownames(datExpr)


# getting trait data
pfc1
# p1 <- read_csv("manuscript/brain/results_tables/TraitsforWGCNA_mPFC70.csv")
# Re-cluster samples with out any samples filtered 
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(pfc1, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(pfc1),
                    main = "")

datTraits <-pfc1





#PART 2
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

saveRDS(MEs,"manuscript/brain/results_use/WGCNA/WGCNA_MEs_mPFC70Power5.RDS")


#PART3
# sizeGrWindow(10,6)

dev.off() # make sure you do this before AND after

png(file = "manuscript/brain/imgs/module_trait_mPFC70Power5use.png",
    width=2000, height=2600, res = 300)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = "mPFC 70min: Module-trait relationships")
#
## remember - they are correlation, not regression

dev.off() 
#PART 4 =====================================================================================
# Define variable David's score containing the David's score column of datTrait
status = as.data.frame(datTraits$condition2);
names(status) = "trait"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, status, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(status), sep="");
names(GSPvalue) = paste("p.GS.", names(status), sep="");


geneModuleMembership %>% 
  rownames_to_column("ensgene") %>% 
  left_join(geneTraitSignificance %>% rownames_to_column("ensgene")) -> gene_MM_TS_all



# to make holistic freaking dataframe ==================================================
gene_MM_TS_all$module <- moduleColors


get_mm <- function(x){
  x$moduleMembership <- x[colnames(x) == paste("MM",x$module,sep = "")] %>%
    unlist %>%
    as.numeric
  xx <- x %>%
    dplyr::select(ensgene,module,moduleMembership,GS.trait)
  return(x)
}


wgcna_whole <- get_mm(gene_MM_TS_all[1,])

for(i in 2:nrow(gene_MM_TS_all)){
  wgcna_whole <- rbind(wgcna_whole,get_mm(gene_MM_TS_all[i,]))
}

wgcna_whole %>%
  rename(GS.status = GS.trait) -> wgcna_whole


wgcna_whole %>%
  left_join(geneTraitSignificance %>% rownames_to_column("ensgene")) -> wgcna_all

saveRDS(wgcna_all, "manuscript/brain/results_use/WGCNA/WGCNA_WGCNA_MM_GS_mPFC70Power5.RDS")


##boxplots

ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleID") %>%
  pivot_longer(cols = 2:23, names_to = "Module") %>% 
  full_join(pfc) %>% filter(Module != "MEgrey") %>% 
  mutate(status = condition2)

head(ME_df)

ME_df$Module <- gsub("ME", "", ME_df$Module)
source("functions/geom_boxjitter.R")


p1 <- ME_df %>% 
  ggplot(aes(condition1, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+         
  labs(x = "Social condition",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol = 7) + theme_classic()+
  theme(legend.position = "none", text =element_text(size =10))



ggsave("manuscript/brain/imgs/boxplots_WGCNAPFC70NORM_RG.png", p1, width= 15, height = 5, dpi = 150)



ME_df$conditionx <- ifelse(ME_df$condition1 == "SUB", "Same", ME_df$condition1)
ME_df$conditionx <- ifelse(ME_df$conditionx == "DOM", "Same", ME_df$conditionx)


p2 <- ME_df %>% 
  ggplot(aes(conditionx, value, fill = conditionx, color = conditionx))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+         
  labs(x = "Social condition",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol = 7) + theme_classic()+
  theme(legend.position = "none", text =element_text(size =10))


p2
ggsave("manuscript/brain/imgs/boxplots_WGCNAPFC70NORM_RG.png", p1, width= 15, height = 5, dpi = 150)

#PART 5 =====================================================================================

# Module of interest via significance 
my_trait = "Status"
module = "green"
module = "black"
module = "greenyellow"
module = "tan"
module = 'blue'
module = 'turquoise'
module = 'yellow'
module = 'grey60'
module = 'royalblue'
module = 'brown'
module = 'cyan'
module = 'midnightblue'
module = 'purple'
module = 'pink'
module = 'lightcyan'
module = 'lightyellow'

module_list = c("green", "black", "greenyellow", "tan", "blue", "turquoise", "yellow", "grey60", 
                "royalblue", "brown", "cyan", "midnightblue", "purple", "pink", "lightcyan", "lightyellow")
hub_gene_list = vector('list', length = length(module_list))
names(hub_gene_list) <- module_list


column = match(module, modNames);
moduleGenes = moduleColors==module;
colnames(datExpr)[moduleColors==module] -> module_gene
# 
# 
# sizeGrWindow(7, 7);
# par(mfrow = c(1,1));
# verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                    abs(geneTraitSignificance[moduleGenes, 1]),
#                    xlab = paste("Module Membership in", module, "module"),
#                    ylab = glue::glue("Gene significance for {my_trait}"),
#                    main = paste("Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
# datExpr[moduleGenes,]

my_MM_threshold = 0.8
my_GS_threshold = 0.2


for (module in module_list){
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  
  colnames(datExpr)[moduleColors==module] -> module_gene
  
  
  grcm38 %>% 
    dplyr::select(ensgene, symbol, chr, description) %>% 
    filter(ensgene %in% module_gene) -> module_gene_info
  
  
  
  geneModuleMembership %>% 
    rownames_to_column("ensgene") %>% 
    left_join(geneTraitSignificance %>% rownames_to_column("ensgene")) -> gene_MM_TS
  
  
  
  gene_MM_TS %>% 
    filter(ensgene %in% module_gene) %>% 
    dplyr::select(ensgene, glue::glue("MM{module}"), GS.trait) %>% 
    filter(abs(GS.trait) >= my_GS_threshold) -> x
  
  x[x[,glue::glue("MM{module}")]>my_MM_threshold,] -> hub_genes
  
  hub_genes %>% 
    left_join(module_gene_info) %>% 
    mutate(moduleName = glue::glue("{module}")) %>% 
    rename(moduleMembership = glue::glue("MM{module}")) -> hub_genes
  hub_gene_list[[module]] <- hub_genes
  
}



hub_gene_list %>% 
  do.call(rbind,.)%>% 
  unique() -> hubgenes_df

hubgenes_df %>% 
  dplyr::arrange(desc(GS.trait)) %>% 
  dplyr::select(-ensgene, -chr) %>% 
  head(15)

write.csv(hubgenes_df, "manuscript/brain/results_tables/WGCNA/WCGNA_hubgene_list_mPFC70Power5.csv", row.names = F)


ADJ=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ, moduleColors) 

Alldegrees1 %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38 %>% dplyr::select(ensgene,symbol)) -> kIN_df

xx <- cbind(ensgene = colnames(datExpr), module = moduleColors) %>% 
  as.tibble()


kIN_df %>% 
  left_join(xx) -> kIN_df

saveRDS(kIN_df,"manuscript/brain/results_use/kIN/kIN_dataframe_mPFC_power5.RDS")


set.seed(312)
library(clusterProfiler)
library(enrichplot)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

moduleNumber = length(unique(moduleColors))

modNames = substring(names(MEs), 3)

gettop10GO_WGCNA <- function(module,my_showCategory = 10){
  
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  colnames(datExpr)[moduleColors==module] -> module_gene
  
  
  grcm38 %>% 
    filter(ensgene %in% module_gene) %>% 
    filter(!is.na(entrez)) %>% 
    dplyr::select(entrez) -> go_df_wgcna
  
  
  ggo <- enrichGO(gene = go_df_wgcna$entrez %>% unique(),
                  OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                  keyType = "ENTREZID",
                  ont = 'BP',
                  readable = T,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.50)
  
  
  fortify(
    ggo,
    showCategory = my_showCategory,
    by = "Count",
    split = NULL,
    includeAll = TRUE
  ) %>% 
    dplyr::arrange(desc(GeneRatio)) %>% 
    mutate(module = module) -> temp1
  
  return(rbind(temp1))
  
}

my_ont = "BP"
my_showCategory = 100


module_list %>% unique() -> allcolors

WGCNA_GOs <- vector('list', length(allcolors))

for(i in 1:length(allcolors)){
  gettop10GO_WGCNA(allcolors[i],my_showCategory) -> WGCNA_GOs[[i]]
}

WGCNA_GOs %>% 
  do.call(rbind,.) -> wgcna_all_gos

write.csv(wgcna_all_gos, 
          "manuscript/brain/results_tables/WGCNA/WGCNA_all_gos_catogeryBP_mPFC70_.csv",
          row.names = F)


######linear model 
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes

pfc
pfc$conditionx <- ifelse(pfc$condition1 == "SUB", "Same", pfc$condition1)
pfc$conditionx <- ifelse(pfc$conditionx == "DOM", "Same", pfc$conditionx)
pfc$conditionx <- factor(pfc$conditionx, levels = c("ASC", "DES", "Same"))


orderMEs(MEs0) %>% 
  rownames_to_column("SampleID") %>%
  left_join(pfc) %>% na.omit(.) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  mutate(condition1 = factor(conditionx, levels = c("DES", "ASC", 'Same'))) %>%
  relocate(conditionx,condition1, batch,Postds, Preds, post.given1, post.received1) %>% 
  dplyr::select(-SampleID, -condition2, -CORT, -wt_d8)-> ME_df

lm_result_list <- list()

library(lme4)
library(lmerTest)

lmer(MEbrown ~ conditionx +(1|batch) , data = ME_df) -> mod1
lmer(MEbrown ~ condition1 +(1|batch) , data = ME_df) -> mod2
summary(mod1)
summary(mod2)

for(x in 1:length(MEs0)){
  k = x + 7
  ME_df[,c(1:7,k)] -> df
  md <- gsub("ME","",colnames(df)[8])
  colnames(df)[8] <- "module"
  lmer(module ~ conditionx +(1|batch) , data = df) -> mod1
  lmer(module ~ condition1 +(1|batch) , data = df) -> mod2
  summary(mod1)
  summary(mod2)
  df
  rbind(summary(lm(df[,8] ~ df[,1]))$coefficients[2,],
        summary(lm(df[,8] ~ df[,1]))$coefficients[3,],
        summary(lm(df[,8] ~ df[,2]))$coefficients[3,])%>%
    as.data.frame() %>%
    cbind(key = c("DES-ASC","DES-Same","ASC-Same")) %>%
    mutate(module = md) -> lm_result_list[[x]]
}

lm_result_list %>% 
  do.call(rbind,.)%>% 
  filter(.,module != "grey")-> lm_result_all


saveRDS(lm_result_all,"manuscript/brain/results_use/WGCNA/WGCNA_lm_result_mPFC70.RDS")


