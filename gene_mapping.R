#GG is the list of gene symbols to be mapped into entrez
mapping_symbol_to_entrez = function(GG){
  geneMapping2= select(org.Hs.eg.db, GG, c("ENTREZID","GENENAME"),"ALIAS")
  geneMapping = c()
  
  mapped_symbols = unique(geneMapping2[,1])
  for(index in mapped_symbols){
    idx = which(geneMapping2[,1] %in% index)
    geneMapping = rbind(geneMapping,geneMapping2[idx[1],])
  }
  
  toRem = which(is.na(geneMapping[,2]))
  if(length(toRem)>0){
    geneMapping = geneMapping[-toRem,]
  }
  
  return(geneMapping)
  
}
