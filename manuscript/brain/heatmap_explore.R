gh
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

asc.hmx <- as.matrix(asc.hm)
heatmap(asc.hmx,
        Rowv=NA, Colv=NA, col=rev(brewer.pal(9,"RdBu")))

# colnames(asc.hmx) <- substr(colnames(asc.hmx), 1,11)
library(plotly)
library(heatmaply)
heatmaply(th, 
          dendrogram = "none",
          xlab = "", ylab = "", 
          main = "",
          scale = "column",
          margins = c(60,100,40,20),
          grid_color = "white",
          grid_width = 0.00001,
          titleX = FALSE,
          hide_colorbar = TRUE,
          branches_lwd = 0.1,
          # label_names = c("Country", "Feature:", "Value"),
          fontsize_row = 5, fontsize_col = 5,
          # labCol = colnames(mat),
          # labRow = rownames(mat),
          plot_method = c("ggplot"),
          heatmap_layers = theme(axis.line=element_blank(), axis.title = element_text(size =100)))

library(ComplexHeatmap)
th <- t(asc.hmx)
colnames(asc.hmx) <- substr(colnames(asc.hmx),1,12)

Heatmap(asc.hmx, name = "LogFC",cluster_rows = FALSE,
        column_names_rot = 45, show_heatmap_legend = FALSE, rect_gp = gpar(col = "white", lwd = 2))


p <-  Heatmap(asc.hmx, name = "LogFC",cluster_rows = FALSE,
     rect_gp = gpar(col = "white", lwd = 2))
p1
ggsave("manuscript/brain/imgs/ASC_StableHeatmap.png",width =4 , height = 20, dpi = 300)
plot = p1