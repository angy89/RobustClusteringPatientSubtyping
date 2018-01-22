RSC_thresholding = function(FLOSSES,R){
  tgrid   <- seq(0.005, 0.995, by=0.01)  
  avgLoss <- colMeans(FLOSSES)
  plot(tgrid, avgLoss, t='b')
  topt  <- tgrid[which.min(avgLoss)]
  RSC   <- R
  RSC[abs(RSC)<=topt] <- 0
  return(RSC)
}