
# Load the package
library(WGCNA)
library(tidyverse)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)


Exp70 <- readRDS("manuscript/brain/results_use/WGCNA/WGCNA_datExpr70NORM_RG.RDS") %>% na.omit(.)

Color70 <-readRDS("manuscript/brain/results_use/kIN/kIN_dataframe_mPFC_power5.RDS")
Color70 <- Color70 %>% dplyr::select(ensgene, module) %>% 
  dplyr::rename(.,geneID = ensgene)%>% dplyr::rename(.,Module = module) %>% 
  na.omit(.)

Exp25 <- readRDS("manuscript/brain/results_use/WGCNA/WGCNA_datExpr25NORM_RG.RDS") %>% na.omit()
Color25 <-readRDS("manuscript/brain/results_use/kIN/kIN_dataframe_mPFC25_power4.RDS")
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
# 
# setLabels = c("T70", "T25")
# multiExpr = list(T70 = list(data = Exp70[,Color70$geneID]), T25 = list(data = Exp25[,Color25$geneID]))
# multiColor = list(T70 = Color70$Module, T25= Color25$Module)
# nSets = 2
# 
# #takes 2.5 hours - make use you save out put. 
# # Calculate module preservation statistics
# system.time( {
#   mp = modulePreservation(multiExpr, multiColor,
#                           referenceNetworks = c(1:2),
#                           loadPermutedStatistics = FALSE,
#                           networkType = "signed hybrid",
#                           nPermutations = 1000,
#                           verbose = 3)
# } );
# # Save the results
# saveRDS(mp,"manuscript/brain/results_use/mPFC_modulePreservation.RDS")

mp <- readRDS("manuscript/brain/results_use/mPFC_modulePreservation.RDS")

ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])
# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )


# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
modColorsx = ifelse(modColors == "blue"|modColors ==  "pink"| modColors == "grey60"|modColors ==  "yellow" |modColors ==  "red"| modColors == "turquoise"| modColors == "grey", modColors, "azure2")

moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColorsx %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5)
png("manuscript/brain/imgs/mPFC_-modulePreservation-Zsummary-medianRank2.png", wi=10, h=5,units = 'in', res = 300)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in c(1:2))
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColorsx[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = .8, offs = 0.135);
    abline(h=0)
    abline(h=2, col = "red", lty = 2)
    abline(h=10, col = "blue", lty = 2)
  }
}
# If plotting into a file, close it
dev.off()




##PLOT2
# Re-initialize module color labels and sizes
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][, 1];
# Exclude improper modules
plotMods = !(modColors %in% c("grey", "gold"));
# Create numeric labels for each module
labs = match(modColors[plotMods], standardColors(50));
# Start the plot: open a suitably sized graphical window and set sectioning and margins. Alternatively,
# plot into a pdf file.
 sizeGrWindow(10, 9);
png("manuscript/brain/imgs/mPFC-PreservationZStatistics2.png", width=10, height=9,units = 'in', res = 150)
par(mfrow = c(4,4))
par(mar = c(3,3,2,2))
par(mgp = c(1.6, 0.4, 0));
for (s in 1:ncol(statsZ))
{
  min = min(statsZ[plotMods, s], na.rm = TRUE);
  max = max(statsZ[plotMods, s], na.rm = TRUE);
  if (min > -max/12) min = -max/12
  plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
       main = colnames(statsZ)[s],
       cex = 2.2,
       ylab = colnames(statsZ)[s], xlab = "Module size", log = "x",
       ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)),
       xlim = c(30, 1000),
       cex.lab = 1.2, cex.axis = 1.2)
  labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = .8, offs = 0.06);
  # text(moduleSizes[-1], statsZ[-c(1:2), s], labels = letter[-c(1:2)], col = "black"); #modColors[-2]);
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
}
# If plotting into a file, close it, otherwise it is unreadable.
dev.off();

data.frame(color = modColors[plotMods], label = labs)


