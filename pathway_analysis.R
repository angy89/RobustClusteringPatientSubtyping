
pathway_analysis = function(cluster=RC[[1]]$cluster,red_data){
  LEV = sort(unique(cluster))
  
  topTableRes_clusters =list()
  DF = c()
  for(i in 1:length(LEV)){
    if(i>1)red_data = t(red_data)
    clust_contr=cluster
    clust_contr[clust_contr!=LEV[[i]]] = "control"
    mod <- model.matrix(~-1 + clust_contr)
    rownames(mod) = rownames(red_data)
    colnames(mod) = gsub(pattern = "clust_contr",replacement = "",x = colnames(mod))
    colnames(mod)[1] = paste("Cluster_",colnames(mod)[1],sep="")
    red_data = t(red_data)
    
    fit <- lmFit(red_data, mod)
    contrast = makeContrasts(contrasts = c("Cluster_1-control","Cluster_2-control",
                                           "Cluster_3-control","Cluster_4-control")[i],levels = mod) 
    
    fit2  <- contrasts.fit(fit, contrast)
    fit2 <- eBayes(fit2)
    TT = topTable(fit2,number=nrow(red_data),sort.by = "logFC")
    
    background_entrez =mapping_symbol_to_entrez(rownames(TT))#[,2]
    rownames(background_entrez) = background_entrez$ALIAS
    TT = cbind(TT,background_entrez[rownames(TT),"ENTREZID"])
    toRem = which(is.na(TT$`background_entrez[rownames(TT), "ENTREZID"]`))
    if(length(toRem)>0)TT = TT[-toRem,]
    
    mydf <- data.frame(Entrez=TT$`background_entrez[rownames(TT), "ENTREZID"]`, FC=TT$logFC,adj.P.Val = TT$adj.P.Val, PVal = TT$P.Value)
    mydf <- mydf[mydf$PVal<0.05,]
    mydf$group <- "up"
    mydf$group[mydf$FC < 0] <- "down"
    mydf$othergroup = paste("Cl",i,"-others",sep="")
    
    DF = rbind(DF,mydf)
    topTableRes_clusters[[i]] = TT
    
  }
  

  formula_res <- compareCluster(Entrez~group+othergroup, data=DF, 
                                fun="enrichKEGG",organism = "hsa",use_internal_data = F,
                                pAdjustMethod = "none",universe = TT$`background_entrez[rownames(TT), "ENTREZID"]`,
                                qvalueCutoff = 1)
  
  FR = as.data.frame(formula_res)
  
  
  dp = dotplot(formula_res, x=~group,showCategory=10) + 
    ggplot2::facet_grid(~othergroup) + 
    ggplot2::scale_colour_gradient(low = "black", high = "white") + 
    ggplot2::theme(axis.text.y=ggplot2::element_text(size=ggplot2::rel(0.95)))
  
  return(list(FR= FR,dp = dp,formula_res=formula_res))
}

