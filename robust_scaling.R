robust_scaling = function(dataset,RSC){
  ## Robust scaling
  
  #data centered on median e mean absolute deviation standardization
  marg.center <- apply(dataset, 2, median)
  marg.scale  <- apply(dataset, 2, mad)
  X           <- scale(dataset, center = marg.center, scale = marg.scale)
  ## Spectral decomposition of R
  Eig   <- eigen(RSC, symmetric=TRUE) #autovettori della matrice RSC
  V = Eig$vectors
  Values = Eig$values
  return(list(X=X,V=V,Values=Values))
}