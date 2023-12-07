# libraries 
library(limma)
library(edgeR)
library(WGCNA)
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
library(tidyverse)
grcm38 # mouse genes


#get meta data for 70 min mPFC 
dat70 <- read_csv("manuscript/brain/results_tables/coldata.csv")
colnames(dat70)
dat70 <- dat70 %>% filter(condition != 'control') %>% 
  select(batch, post_idbatch, pre_idbatchcage, Postds, Postrank, Preds, Prerank, wt_d8, CORT = mean_con_ng_ul, SampleID,given1, received1, condition1) %>% 
  mutate(time = 70)

#get meta data for 25hr mPFC 
dat25 <- read_csv("manuscript/brain/results_tables/coldata25hr_use.csv")
colnames(dat25)
dat25 <- dat25 %>% 
  select(batch, post_idbatch, pre_idbatchcage, Postds, Postrank, Preds, Prerank, wt_d8, CORT = mean_con_ng_ul, SampleID = SampleNames,given1, received1, condition1) %>% 
  mutate(time = 25)

dat25$SampleID <- substr(dat25$SampleID, 6,11 )

coldata <- dat70 %>% rbind(dat25)

write.csv(coldata, "manuscript/brain/results_tables/coldata_ALL.csv", row.names = F)


pdata <- coldata  %>% column_to_rownames(., var= "SampleID")
#70 min = 0, 25hr = 1
pdata$time<- ifelse(pdata$time == 70, 0, 1)

# get count data for both time periods
c70<-  read.csv("brain/PFC_counts.csv", row.names = 1)
#trim sample ids
colnames(c70)[c(1:67)] <- substr(colnames(c70)[c(1:67)], 7, 13)


c25 <- read.csv("brain/PFC_25hcounts.csv", row.names = 1)
#trim sample ids
colnames(c25)[c(1:24)] <- substr(colnames(c25)[c(1:24)], 6, 11)
#b1.2.1 outlier DES - remove
c25 <- c25[,c(1:19,21:24)]

#joins counts by gene
c70 <- c70 %>% rownames_to_column(., var = "gene")
c25 <- c25 %>% rownames_to_column(., var = "gene")
counts <-  c70 %>% full_join(c25) %>% column_to_rownames(., var = "gene")


#check before normalizing 
counts<- counts[,rownames(pdata)]
all(rownames(pdata) == colnames(counts))

#normalize and filter with all groups 
dlNorm <- counts[!is.na(rowSums(counts)),]

d = apply(dlNorm, 2, as.numeric)
dim(d)
# [1] 25158    51

d0= DGEList(d, group = pdata$condition1)
dim(d0)
rownames(d0) <- rownames(dlNorm)
d0 <- calcNormFactors(d0)

#won used 10 in 90% of samples for brain paper, which is what tutorial suggest. 
#in liver paper she did 50 in 90% samples 
cutoff <- 10
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
# [1] 13368    51


# Now take out groups that you want
dge.dl$samples$group


pdata %>% 
  dplyr::select(condition1,time) -> var_info  
head(var_info)


dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check


##following Won's code
var_info$condition1 %>%
  factor(.,levels = c("DOM", "DES","ASC", "SUB")) -> group.dl

var_info$time %>% scale() -> time.dl


design.dl <- model.matrix(~ 0 + group.dl*time.dl)
colnames(design.dl) -> mycolnames

v.dl = voom(dge.dl, design.dl, plot = F)


# Many functions expect the matrix to be transposed
datExpr <- t(v.dl$E)
## check rows/cols
nrow(datExpr)
ncol(datExpr)
rownames(datExpr)

# getting trait data
colnames(pdata)
pfc <- pdata[c(1,4:13)]

ifelse(pfc$condition1 == "DOM", 2, pfc$condition1) -> pfc$condition2
ifelse(pfc$condition2 == "DES", 1, pfc$condition2) -> pfc$condition2
ifelse(pfc$condition2 == "ASC", 0, pfc$condition2) -> pfc$condition2
ifelse(pfc$condition2 == "SUB", -1, pfc$condition2) -> pfc$condition2

pfc1 <- pfc %>% 
  dplyr::select(-condition1)

fix <- pfc1 %>% filter(condition2 == -1 & time == 1 )

fix$CORT <- ifelse(is.na(fix$CORT), mean(fix$CORT, na.rm = T), fix$CORT)
fix$wt_d8 <- ifelse(is.na(fix$wt_d8), mean(fix$wt_d8, na.rm = T), fix$wt_d8)
fix$given1 <- ifelse(is.na(fix$given1), mean(fix$given1, na.rm = T), fix$given1)
fix$received1 <- ifelse(is.na(fix$received1), mean(fix$received1, na.rm = T), fix$received1)
pfc1 <- pfc1 %>% full_join(fix) 
pfc1 <- pfc1[c(1:43,45:52),]

#fix weight typo
pfc1$wt_d8 <-  as.character(pfc1$wt_d8)
pfc1$wt_d8 <- gsub("318", '31.8', pfc1$wt_d8)


#make everything numeric 
pfc1[1:11] <- lapply(pfc1[1:11], as.numeric)


sampleTree = hclust(dist(datExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)

par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


# Re-cluster samples with out any samples filtered 
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(pfc1, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(pfc1),
                    main = "AMY dendrogram and trait heatmap(all samples)")

collectGarbage()
#not going to filter because ascenders rreally are not grouping anyways?
saveRDS(datExpr,"manuscript/brain/results/WGCNA_datExpr_ALLSampleswithTime.RDS")


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
     main ="mPFC: Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h

abline(h=0.85,col="red")  #changed to 0.8 (originally 0.9 in the tutorial)
abline(h=0.90,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = "mPFC: Mean connectivity")


text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
print(sft)


# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html
# power 6 since time causing two much difference in groups that power estimate isn't working. 
# I can try consense analysis as well but Hans didn't like that last time. 

cor <- WGCNA::cor
WGCNA_get_net <- function(my_tissue = "ALL",
                          my_power =  6,
                          my_TOMType ="signed", 
                          my_networkType = "signed hybrid"){
  
  x <- readRDS("manuscript/brain/results/WGCNA_datExpr_ALLSampleswithTime.RDS") 
  set.seed(312)
  net = blockwiseModules(x, 
                         power = my_power,
                         TOMType = my_TOMType, 
                         networkType = my_networkType,
                         minModuleSize = 50,
                         reassignThreshold = 0, 
                         mergeCutHeight = 0.25,
                         numericLabels = TRUE, 
                         pamRespectsDendro = FALSE,
                         saveTOMs = FALSE,
                         verbose = 3)
  
  
  saveRDS(net, "manuscript/brain/results/WGCNA_net_mPFC25_Power6_minmod50_ALL.RDS")
  
}

# WGCNA_get_net("mPFC70min", 7, "signed", "signed hybrid")
WGCNA_get_net("mPFC25hr", 6, "signed", "signed hybrid")
# ========================================================================================

net <- readRDS(glue::glue("manuscript/brain/results/WGCNA_net_mPFC25_Power6_minmod50_ALL.RDS"))

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath

dev.off() # make sure you do this before AND after 
png(file = "manuscript/brain/imgs/cluster_dendo_mPFC_Power6_minmod50_ALL.png",
    width=600, height=350)
plotDendroAndColors(net$dendrograms[[2]], mergedColors[net$blockGenes[[2]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "mPFC Power6")

dev.off()

MEs = net$MEs
#22 modules

set.seed(312)
datExpr <- readRDS("manuscript/brain/results/WGCNA_datExpr_ALLSampleswithTime.RDS") 
net <- readRDS("manuscript/brain/results/WGCNA_net_mPFC25_Power6_minmod50_ALL.RDS")
dim(datExpr)

#PART 1
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
moduleNumber = length(unique(moduleColors))
rownames(datExpr)


# getting trait data


# Re-cluster samples with out any samples filtered 
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(pfc1, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(pfc1),
                    main = "mPFC dendrogram and trait heatmap")

datTraits <- pfc1



#PART 2
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

saveRDS(MEs,"manuscript/brain/results/WGCNA_MEs_mPFCALL_Power6.RDS")


#PART3
# sizeGrWindow(10,6)

dev.off() # make sure you do this before AND after

png(file = "manuscript/brain/imgs/module_trait_mPFCALL_Power6.png",
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
               main = "mPFC: Module-trait relationships")
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

saveRDS(wgcna_all, "manuscript/brain/results/WGCNA_WGCNA_MM_GS_mPFCALL_Power6.RDS")



# Module of interest via significance 
my_trait = "Status"
module = "lightgreen"
module = "royalblue"
module = "midnightblue"
module = "blue"
module = 'greenyellow'
module = 'pink'
module = 'purple'
module = 'black'
module = 'salmon'
module = 'lightcyan'
module = 'lighyellow'
module = 'turquoise'
module = 'grey60'
module = 'magenta'

module_list = c("royalblue", 'lightgreen', "midnightblue", "blue", "greenyellow", "pink", "purple", "black", "salmon", "lightcyan", "lightyellow", "turquoise", "grey60", "magenta")
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
my_GS_threshold = 0.0


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

write.csv(hubgenes_df, "manuscript/brain/results_tables/WCGNA_hubgene_list_mPFCALL_Power6.csv", row.names = F)




set.seed(312)

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
          "manuscript/brain/results_tables/WGCNA_all_gos_catogeryBP_mPFCALL_Power6.csv",
          row.names = F)



#######################
# plot
head(coldata)

ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleID") %>%
  pivot_longer(cols = 2:23, names_to = "Module") %>% 
  full_join(coldata) %>% filter(Module != "MEgrey") %>% 
  mutate(status = condition1)


ME_df$Module <- gsub("ME", "", ME_df$Module)


head(ME_df)

#boxplot of eigengenes and status 
source("functions/geom_boxjitter.R")
source("functions/newggtheme.R")

ME_df$time <- factor(ME_df$time, levels = c("70", "25"))

ME_df %>% 
  ggplot(aes(status, value, fill = time, color = time))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+         
  labs(x = "Social condition",
       y = "Module eigengene")+
  facet_wrap(~Module, nrow = 7) +theme_classic()+
  theme(text =element_text(size =10))

cort <- ME_df %>% 
  ggplot(aes(CORT, value, color = status))+
  geom_point()+   stat_smooth(method = "lm", se =F)+
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+         
  labs(x = "CORT",
       y = "Module eigengene")+
  facet_grid(time~Module) + theme_bw()+
  theme(text =element_text(size =10))

ggsave("manuscript/brain/imgs/WGCNA_ALL_cort.png", cort, width = 30, height = 5)
  
ME_df <- ME_df %>% mutate(tot_agg = given1+received1)
ME_df %>% 
  ggplot(aes(tot_agg, value, color = status,time))+
  geom_point()+   stat_smooth(method = "lm", se =F)+
  scale_color_manual(values = viridis::plasma(8))+
  scale_fill_manual(values = viridis::plasma(8))+         
  labs(x = "Rate of Aggression in PreCage",
       y = "Module eigengene")+
  facet_grid(.,cols =Module, rows =time) + theme_classic()+
  theme(text =element_text(size =10))
