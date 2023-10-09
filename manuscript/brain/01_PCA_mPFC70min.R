library(limma)
library(edgeR)
library(Mus.musculus)
library(DESeq2)
library(limma)
library(tidyverse)
global_size = 10


#Getting metadata ready 
coldata <- read_csv("brain/sample70min_table.csv")
head(coldata)
str(coldata)

#group thing 
coldata$groupEX <- coldata$group
# Normalizing cort data
# df <- transform(df, N = (N - min(N)) / (max(N) - min(N))

coldata <-coldata %>% mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))
coldata <- coldata %>%  dplyr::select(-group, -period)


#getting condition1
table(coldata$condition, coldata$Postrank)
# CDOM, RDOM to Descenders (DOM to SUB)(4->1)
# CSUB, SUB to Ascenders (Sub to DOM)  (1->4)

coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "DOM", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "SUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 3 & coldata$Postrank == 4, "SUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "DES", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "ASC", coldata$condition1)


table(coldata$condition1)
# Second mpfc data data 
pfc_data<- coldata %>% filter(region != "AMY") %>% filter(Prerank !=3) %>% filter(groupEX != "control") %>% filter(Postrank != 3)

row.names <- pfc_data$SampleName

pfc_data1 <- pfc_data%>% dplyr::select(-SampleName)

row.names(pfc_data1) <- row.names #Assigning row names from as sample names  
head(pfc_data1)

# Bring in count data for mpfc
p_countdata <- read_csv("brain/PFC_counts.csv")
p_countdata$X -> nrows
p_count <- p_countdata[,-1]
row.names(p_count) <- nrows
p_count <-as.data.frame(p_count)

#filter
p_count <- p_count[rowSums(p_count >= 20), ]


p_count <- p_count[, rownames(pfc_data1)]
all(rownames(pfc_data1) == colnames(p_count)) #check


d = apply(p_count, 2, as.numeric)
dim(d)

#DGEList from limma
d0<-DGEList(d, group = pfc_data1$condition1)
dim(d0)
rownames(d0) <- rownames(p_count) #redo rownames as get lost when making matrix
d0 <- calcNormFactors(d0) #Normalizing step.


## Removing genes from PCA where average count is below cutoff.
cutoff <- 40
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
# [1] 2658   23 #50
# [1] 2314   23 # 100
# [1] 2404   23 # 75

##making it a factor
pfc_data1$condition1 %>%
  factor(.,levels = c("DES", "DOM","ASC", "SUB")) -> group.dl

## Make the design/contrasts
design.dl <- model.matrix(~ 0 + group.dl)
colnames(design.dl) -> mycolnames

## transform to logcpm counts
v = voom(dge.dl, design.dl, plot = F)


colData <- v$design %>% data.frame()
colData <- pfc_data1 %>% rownames_to_column(var = "SampleName") %>% cbind(colData) 

#get expression levels
rv <- rowVars((v$E)) #measure of variation for each gene
select <- order(rv, decreasing = TRUE)[1:250] #top 250 varying genes
pca1 <- prcomp(t((v$E)[select, ])) #do pca

#get groups
condition3.df <- as.data.frame(colData[, 'condition1', 
                                       drop = FALSE])



#put results into df
d1 <- data.frame(PC1 = pca1$x[, 1], PC2 = pca1$x[, 2], PC3 = pca1$x[, 3], PC4 = pca1$x[, 4], 
                 condition3.df, name = colnames(v$E))


head(d1)

# make scree plot
pv1 <- ((pca1$sdev^2) / (sum(pca1$sdev^2)))*100
barplot(pv1, cex.names=1, xlab=paste("Principal component (PC), 1-", length(pca1$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))

var1 <- ((pca1$sdev[1]^2) / (sum(pca1$sdev^2)))*100
var2 <- ((pca1$sdev[2]^2) / (sum(pca1$sdev^2)))*100
var3 <- ((pca1$sdev[3]^2) / (sum(pca1$sdev^2)))*100
var4 <- ((pca1$sdev[4]^2) / (sum(pca1$sdev^2)))*100

## actually plots 

xlab=paste("PC1, ", round(pv1[1], 2), "%")
ylab=paste("PC2, ", round(pv1[2], 2), "%")
xlab3=paste("PC3, ", round(pv1[3], 2), "%")
ylab4=paste("PC4, ", round(pv1[4], 2), "%")


p1=ggplot(d1, aes(PC1, PC2, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = xlab, y = ylab)+
  scale_color_manual(values = viridis::viridis(6))+
  theme_classic()

p2=ggplot(d1, aes(PC1, PC3, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = xlab, y = xlab3)+
  scale_color_manual(values = viridis::viridis(6))+
  theme_classic()

p3=ggplot(d1, aes(PC2, PC3, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = ylab, y = xlab3)+
  scale_color_manual(values = viridis::viridis(6))+
  theme_classic()

p4=ggplot(d1, aes(PC3, PC4, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = xlab3, y = ylab4)+
  scale_color_manual(values = viridis::viridis(6))+
  theme_classic()

p5=ggplot(d1, aes(PC1, PC4, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = xlab, y = ylab4)+
  scale_color_manual(values = viridis::viridis(6))+
  theme_classic()


p6=ggplot(d1, aes(PC2, PC4, color = condition1, shape = condition1))+
  geom_point(size = 5)+
  scale_shape_manual(values=c(16,17,15,3,7,8,1))+
  labs(x = ylab, y = ylab4)+
  scale_color_manual(values = viridis::viridis(6))+
  theme_classic()

library(gridExtra)
a <- grid.arrange(p1,p2,p3,p4,p5,p6, nrow=2)

ggsave("manuscript/brain/imgs/pca_mPFC_ALL70mn.png", a, width = 12, height = 6)

