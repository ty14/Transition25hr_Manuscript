ConvertHumantoMouse <- function(x){
  
  require("biomaRt")
 
  human = useMart(biomart = "genes", dataset = "hsapiens_gene_ensembl",host="http://oct2022.archive.ensembl.org")
  mouse = useMart(biomart = "genes", dataset = "mmusculus_gene_ensembl", hhost="http://oct2022.archive.ensembl.org")
  # mart@host="http://may2012.archive.ensembl.org:80/biomart/martservice"
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  mousex <- unique(genesV2[, 2]) # unique the 2nd column values
  
  human_genes_number <- length(x)
  mouse_genes_number <- length(mousex)
  
  if(human_genes_number != mouse_genes_number){
    genes_not_trans <- setdiff(x, genesV2$HGNC.symbol)
    print("These genes could not be translated:")
    print(genes_not_trans)
    print(paste("A total number of ",length(genes_not_trans),"genes could not be translated!"),sep=" ")
  }else{
    print("All genes were translated successfully!")
  }
  
  return(mousex)
}

ConvertHumantoMouse(agg)

