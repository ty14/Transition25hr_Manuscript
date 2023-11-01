
# Load the package
library(WGCNA)
library(tidyverse)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)


Exp70 <- readRDS("manuscript/brain/results/WGCNA_datExpr70.RDS") %>% na.omit(.)

Color70 <-readRDS("manuscript/brain/results/kIN/kIN_dataframe_mPFC_Power5.RDS")
Color70 <- Color70 %>% dplyr::select(ensgene, module) %>% 
  dplyr::rename(.,geneID = ensgene)%>% dplyr::rename(.,Module = module) %>% 
  na.omit(.)

Exp25 <- readRDS("manuscript/brain/results/WGCNA_datExpr25.RDS") %>% na.omit()
Color25 <-readRDS("manuscript/brain/results/kIN/kIN_dataframe_mPFC25_Power4.RDS")
Color25 <- Color25 %>% dplyr::select(ensgene, module) %>% dplyr::rename(.,geneID = ensgene) %>% 
  dplyr::rename(.,Module = module) %>% na.omit(.)


Color25 = Color25[Color25$geneID %in% Color70$geneID,]
Color70 = Color70[Color70$geneID %in% Color25$geneID,]
all.equal(Color70$geneID, Color25$geneID)



#match genes in 70 to 25 expression data 
Exp70 <- Exp70[,colnames(Exp70) %in% colnames(Exp25)]
Exp25 <- Exp25[,colnames(Exp25) %in% colnames(Exp70)]
time = match(colnames(Exp70), colnames(Exp25));
table(is.finite(time))

all.equal(colnames(Exp70), colnames(Exp25))








#  Module preservation

setLabels = c("T70", "T25")
multiExpr = list(T70 = list(data = Exp70[,Color70$geneID]), T25 = list(data = Exp25[,Color25$geneID]))
multiColor = list(T70 = Color70$Module, T25= Color25$Module)
nSets = 2

#takes 2.5 hours - make use you save out put. 
# Calculate module preservation statistics
system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = c(1:2),
                          loadPermutedStatistics = FALSE,
                          networkType = "signed hybrid",
                          nPermutations = 1000,
                          verbose = 3)
} );
# Save the results
saveRDS(mp,"manuscript/brain/results/mPFC_modulePreservation.RDS")


