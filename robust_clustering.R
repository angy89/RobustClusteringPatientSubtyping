robust_clustering = function(A,X,V,K,nCores=4){
  ## Execute otrimle for each combination of the parameters
  ans <- list()
  for (i in 1: nrow(A)){
    message('Running case = ', i, ' (out of ', nrow(A), ')')
    ## Project data onto the  reduced coordinate system
    d <- 1:A[i,1]
    Z <- X %*% V[ , d]
    ## Do otrimle clustering
    ## !!! see ncores here !!!
    ans[[i]] <- otrimle(data = Z, G = K, erc = A[i,2], ncores = nCores, monitor = FALSE)
  }
  return(ans)
}