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


# Expression values
dlNorm <-  read.csv("brain/PFC_25hcounts.csv", row.names = 1)
#remove zeros
dlNorm <- dlNorm[apply(dlNorm[], 1, function(x) !all(x==0)),]
#trim sample ids
colnames(dlNorm)[c(1:24)] <- substr(colnames(dlNorm)[c(1:24)], 6, 11)

#b1.2.1 outlier DES - remove
dlNorm <- dlNorm[,c(1:19,21:24)]


#Group traits
#Getting metadata ready 
coldata <- read.csv("manuscript/brain/results_tables/coldata25hr_use.csv")
head(coldata)

coldata$SampleNames <-  substr(coldata$SampleNames, 6,11)


coldata <- coldata %>% column_to_rownames(., var = "SampleNames")

coldata
#check before normalizing 
dlNorm<- dlNorm[,rownames(coldata)]
all(rownames(coldata) == colnames(dlNorm))

#normalize and filter with all groups 

dlNorm <- dlNorm[!is.na(rowSums(dlNorm)),]

d = apply(dlNorm, 2, as.numeric)
dim(d)

d0= DGEList(d, group = coldata$condition1)
dim(d0)
rownames(d0) <- rownames(dlNorm)
d0 <- calcNormFactors(d0)

#won used 10 in 90% of samples for brain paper, which is what tutorial suggest. 
#in liver paper she did 50 in 90% samples 
cutoff <- 10
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
# [1] 13066    23


# Now take out groups that you want
dge.dl$samples$group


coldata %>% 
  dplyr::select(condition1) -> var_info  

# row.names <- var_info$SampleID

# row.names(var_info) <- row.names #Assigning row names from as sample names  
head(var_info)

dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check

##following Won's code
var_info$condition1 %>%
  factor(.,levels = c("DOM", "DES","ASC", "SUB")) -> group.dl


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
colnames(coldata)
amy <- coldata[c(2,5,11,13,17,18,21,22)]

ifelse(amy$condition1 == "DOM", 2, amy$condition1) -> amy$condition2
ifelse(amy$condition2 == "DES", 1, amy$condition2) -> amy$condition2
ifelse(amy$condition2 == "ASC", 0, amy$condition2) -> amy$condition2
ifelse(amy$condition2 == "SUB", -1, amy$condition2) -> amy$condition2

amy1 <- amy %>% 
  dplyr::select(-condition1)

fix <- amy1 %>% filter(condition2 == -1)

fix$post_Ncort <- ifelse(is.na(fix$post_Ncort), mean(fix$post_Ncort, na.rm = T), fix$post_Ncort)
fix$wt_d8 <- ifelse(is.na(fix$wt_d8), mean(fix$wt_d8, na.rm = T), fix$wt_d8)
fix$given1 <- ifelse(is.na(fix$given1), mean(fix$given1, na.rm = T), fix$given1)
fix$received1 <- ifelse(is.na(fix$received1), mean(fix$received1, na.rm = T), fix$received1)
amy1 <- amy1 %>% full_join(fix) 
amy1 <- amy1[c(1:15,17:24),]

#make everything numeric 
amy1[1:8] <- lapply(amy1[1:8], as.numeric)

#saving to use for other WGCNA analysis 
write.csv(amy1, "manuscript/brain/results_tables/TraitsforWGCNA_mPFC25.csv", row.names = F)


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
var_info
amy1x <- amy1[c(1,3:9, 11:28),]
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(amy1, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(amy1),
                    main = "AMY dendrogram and trait heatmap(all samples)")

collectGarbage()
#not going to filter because ascenders rreally are not grouping anyways?
# saveRDS(datExpr,"manuscript/brain/results/WGCNA_datExpr25.RDS")
# saveRDS(datExprx,"manuscript/brain/results/WGCNA_datExpr25_outliersRemoved.RDS")


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
     main ="mPFC 25hr: Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h

abline(h=0.85,col="red")  #changed to 0.8 (originally 0.9 in the tutorial)
abline(h=0.90,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = "mPFC 25hr: Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
print(sft)



# TOMType = "signed", networkType = "signed",
# TOMType = "unsigned", networkType = "unsigned",
# TOMType = c("signed", networkType = "signed hybrid")
cor <- WGCNA::cor
WGCNA_get_net <- function(my_tissue = "mPFC25hr",
                          my_power =  4,
                          my_TOMType ="signed", 
                          my_networkType = "signed hybrid"){
  
  x <- readRDS("manuscript/brain/results/WGCNA_datExpr25.RDS") 
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
  
  
  saveRDS(net, "manuscript/brain/results/WGCNA_net_mPFC25_Power4_minmod50.RDS")
  
}

# WGCNA_get_net("mPFC70min", 7, "signed", "signed hybrid")
WGCNA_get_net("mPFC25hr", 4, "signed", "signed hybrid")
# ========================================================================================

net <- readRDS(glue::glue("manuscript/brain/results/WGCNA_net_mPFC25_Power4_minmod50.RDS"))

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath

dev.off() # make sure you do this before AND after 
png(file = "manuscript/brain/imgs/cluster_dendo_mPFC_Power4_minmod50.png",
    width=600, height=350)
plotDendroAndColors(net$dendrograms[[2]], mergedColors[net$blockGenes[[2]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "mPFC 25hr Power4")

dev.off()

MEs = net$MEs


##################################
set.seed(312)
datExpr <- readRDS("manuscript/brain/results/WGCNA_datExpr25.RDS") 
net <- readRDS("manuscript/brain/results/WGCNA_net_mPFC25_Power4_minmod75.RDS")
dim(datExpr)

#PART 1
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
moduleNumber = length(unique(moduleColors))
rownames(datExpr)


# getting trait data
colnames(coldata)
coldata$condition1

p1 <- read_csv("manuscript/brain/results_tables/TraitsforWGCNA_mPFC25.csv")
# Re-cluster samples with out any samples filtered 
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(p1, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(p1),
                    main = "AMY dendrogram and trait heatmap")

datTraits <- p1



#PART 2
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

saveRDS(MEs,"manuscript/brain/results/WGCNA_MEs_mPFC25_Power4.RDS")


#PART3
# sizeGrWindow(10,6)

dev.off() # make sure you do this before AND after

png(file = "manuscript/brain/imgs/module_trait_mPFC25_Power4.png",
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
               main = "mPFC 25hr: Module-trait relationships")
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

saveRDS(wgcna_all, "manuscript/brain/results/WGCNA_WGCNA_MM_GS_mPFC25_Power4.RDS")

##linear models for condition2 

set.seed(312)
datExpr <- readRDS("manuscript/brain/results/WGCNA_datExpr25.RDS") 
net <- readRDS("manuscript/brain/results/WGCNA_net_mPFC25_Power4_minmod75.RDS")
MEs <- readRDS("manuscript/brain/results/WGCNA_MEs_mPFC25_Power4.RDS")

MEs 

MEsx <- MEs[,c(1:15)]

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
moduleNumber = length(unique(moduleColors))
rownames(datExpr)
colnames(MEsx)<- gsub("ME", "", colnames(MEsx))


MEDiss = 1-cor(MEsx)
METree = hclust(as.dist(MEDiss), method = "average")

sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "")

sizeGrWindow(7, 6)
p <- plot(METree, main = "Clustering of module eigengenes", xlab = "", ylab = "",  
          cex.lab=0.65,lwd = 2, rotation = 90) 
p <- plot(METree, main = "", xlab = "", ylab = "",  
          cex.lab=0.65,lwd = 2, rotation = 90) 

par(mfrow = c(1,2), mar = c(5,6,1,6))

dend <- as.dendrogram(hclust(as.dist(MEDiss), method = 'average'))
dend2 <- dendextend::color_labels(dend)
all_dendtree <- dendextend::plot_horiz.dendrogram(dend, side = F, lwd = 2)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes

# get linear model data 
#Getting metadata ready 
coldata <- read.csv("manuscript/brain/results_tables/coldata25hr_use.csv")


#dom 
cd <- coldata %>% filter(condition1 != "SUB")

cd$SampleNames <- substr(cd$SampleNames,6,11)

cd$condition1 <- factor(cd$condition1, levels = c("DOM", "DES", "ASC"))
# cd <- cd %>% rownames_to_column(., var = "SampleID")
cd$pre_id <- str_sub(cd$pre_idbatchcage,1,1)
cd$post_id <- str_sub(cd$post_idbatch, 1,1)

MEs0
orderMEs(MEs0) %>% 
  rownames_to_column("SampleNames") %>%
  left_join(cd) %>% 
  mutate(batch = as.factor(batch)) %>% 
  # mutate_if(is.numeric,scale) %>% 
  relocate(condition1, batch, pre_id, post_id, post_Ncort, Postds, Preds,given1, received1 ) %>% 
  dplyr::select(-SampleNames, -post_idbatch, -mean_con_ng_ul, -pre_idbatch, -pre_idbatchcage,
                -time, -Prerank, -condition, -wt_d4, -wt_d8, -wt_12, -region) %>% na.omit(.)-> ME_df

lm_result_list <- list()

colnames(ME_df)
library(lme4)
library(lmerTest)

for(x in 1:length(MEs0)){
  k = x + 8
  ME_df[,c(1:8,k)] -> df
  md <- gsub("ME","",colnames(df)[9])
  colnames(df)[9] <- "module"
  glmer(module ~ condition1*post_Ncort +(1|batch)+(1|pre_id)+(1|post_id) , data = df) -> mod1
  summary(mod1)
  df
  rbind(summary(lm(df[,9] ~ df[,1]*df[,5]))$coefficients[4,],
        summary(lm(df[,9] ~ df[,1]*df[,5]))$coefficients[5,])%>%
    as.data.frame() %>%
    cbind(key = c("DOM-DES","DOM-ASC")) %>%
    mutate(module = md) -> lm_result_list[[x]]
}

lm_result_list %>% 
  do.call(rbind,.)%>% 
  filter(.,module != "grey") %>% filter(., module != "received1")  %>% format(., scientific = F) -> d

unique(d$module)


d$Estimate <- as.numeric(d$Estimate)
d$`Pr(>|t|)`<- as.numeric(d$`Pr(>|t|)`)
d$`Std. Error` <- as.numeric(d$`Std. Error`)
d$module <- as.factor(d$module)
str(d)

color_above <- d$`Pr(>|t|)` < 0.05
color_below <- d$`Pr(>|t|)` > 0.05

d$color <- as.numeric(d$color)
d$color <- ifelse(d$`Pr(>|t|)`< color_above,d$`Pr(>|t|)`,.055)

dom_dot <- d %>% 
  ggplot(aes(y = module, x = Estimate, color =color)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")+
  geom_point(size = 3)+
  facet_wrap(~key)+
  labs(y="", color = "p-value")+
  xlim(-2,2)+
  scale_color_continuous(low = "lightgrey", high = "lightgray", breaks = c(0.01, 0.02,0.03,0.04, 0.05))+
  geom_errorbar(aes(xmin = Estimate-`Std. Error`, xmax = Estimate+`Std. Error`),width = 0.4)+
  scale_y_discrete(limits = rev(levels(d$module)))+
  theme_bw()+
  theme(text = element_text(size = 15))

dom_dot 

# ggsave("manuscript/brain/imgs/dom70_dot25.png", dom_dot, width = 7, height = 6)


coldata <- read.csv("manuscript/brain/results_tables/coldata25hr_use.csv")
# get rid of controls
coldata <- coldata %>% filter(condition != "control")


cs <- coldata %>% filter(condition1 != "DOM")
cs$SampleNames <- substr(cs$SampleNames, 6,11)

cs$condition1 <- factor(cs$condition1, levels =c("SUB", "DES", "ASC"))

cs$pre_id <- str_sub(cs$pre_idbatchcage,1,1)
cs$post_id <- str_sub(cs$post_idbatch, 1,1)

MEs0
orderMEs(MEs0) %>% 
  rownames_to_column("SampleNames") %>%
  left_join(cs) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  relocate(condition1, batch, pre_id, post_id, post_Ncort, Postds, Preds,given1, received1 ) %>% 
  dplyr::select(-SampleNames, -post_idbatch, -mean_con_ng_ul, -pre_idbatch, -pre_idbatchcage,
                -time, -Prerank, -condition, -wt_d4, -wt_d8, -wt_12, -region) %>% na.omit()-> ME_df

lm_result_list <- list()

ME_df$condition1
library(lme4)
library(lmerTest)

for(x in 1:length(MEs0)){
  k = x + 8
  ME_df[,c(1:8,k)] -> df
  md <- gsub("ME","",colnames(df)[9])
  colnames(df)[9] <- "module"
  lmer(module ~ condition1*post_Ncort +(1|batch)+(1|pre_id)+(1|post_id) , data = df) -> mod1
  summary(mod1)
  df
  rbind(summary(lm(df[,9] ~ df[,1]*df[,5]))$coefficients[4,],
        summary(lm(df[,9] ~ df[,1]*df[,5]))$coefficients[5,])%>%
    as.data.frame() %>%
    cbind(key = c("SUB-DES","SUB-ASC")) %>%
    mutate(module = md) -> lm_result_list[[x]]
}

lm_result_list %>% 
  do.call(rbind,.)%>% 
  filter(.,module != "grey") %>% filter(., module != "received1")  %>% format(., scientific = F) -> d

unique(d$module)


d$Estimate <- as.numeric(d$Estimate)
d$`Pr(>|t|)`<- as.numeric(d$`Pr(>|t|)`)
d$`Std. Error` <- as.numeric(d$`Std. Error`)
d$module <- as.factor(d$module)
str(d)

color_above <- d$`Pr(>|t|)` < 0.05
color_below <- d$`Pr(>|t|)` > 0.05

d$color <- ifelse(d$`Pr(>|t|)`< color_above,d$`Pr(>|t|)`,.055)
d$color <- as.numeric(d$color)

sub_dot <- d %>% 
  ggplot(aes(y = module, x = Estimate, color =color)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")+
  geom_point(size = 3)+
  facet_wrap(~key)+
  labs(y="", color = "p-value")+
  scale_color_continuous(low = "red", high = "lightgray", breaks = c(0.01, 0.02,0.03,0.04, 0.05))+
  geom_errorbar(aes(xmin = Estimate-`Std. Error`, xmax = Estimate+`Std. Error`),width = 0.4)+
  scale_y_discrete(limits = rev(levels(d$module)))+
  theme_bw()+
  theme(text = element_text(size = 15))

sub_dot 

# ggsave("manuscript/brain/imgs/sub70_dot25.png", sub_dot, width = 7, height = 6)



# Module of interest via significance 
my_trait = "Status"
module = "green"
module = "blue"
module = "greenyellow"
module = "pink"
module = 'black'
module = 'turquoise'
module = 'tan'
module = 'pink'
module = 'salmon'
module = 'midnightblue'
module = 'red'
module = 'magenta'
module = 'brown'
module = 'cyan'

module_list = c("black", "blue", "brown", "cyan", "green", "greenyellow", "magenta", "midnightblue", "pink", "purple", "red", "salmon", "tan", "turquoise", "yellow")
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

write.csv(hubgenes_df, "manuscript/brain/results_tables/WCGNA_hubgene_list_mPFC25_Power4.csv", row.names = F)


ADJ=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ, moduleColors) 

Alldegrees1 %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38 %>% dplyr::select(ensgene,symbol)) -> kIN_df

xx <- cbind(ensgene = colnames(datExpr), module = moduleColors) %>% 
  as.tibble()


kIN_df %>% 
  left_join(xx) -> kIN_df

saveRDS(kIN_df,"manuscript/brain/results/kIN/kIN_dataframe_mPFC25_Power4.RDS")


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
          "manuscript/brain/results_tables/WGCNA_all_gos_catogeryBP_mPFC25.csv",
          row.names = F)

