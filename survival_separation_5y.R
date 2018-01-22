require(survival)
require(survRM2)

## Computes L1 based measures of separation between survival curves
## note: higher ==> more separation
##
## Inputs:
##    data = data.frame with samples units on rows and two colums with names
##           'Survival' that is  survivial time
##           'Death'    that is a 0-1 variable indicating death
##
## Outputs: a list with multiple separation measures
##
##
survival.curves.separation <- function(data, cluster, tau){

   ans <- list()
   tmp           <- table(cluster)
   ClusterLabel  <- as.numeric(names(tmp))
   ClusterSize   <- as.numeric(tmp)
   K             <- length(ClusterLabel)
   n             <- sum(ClusterSize)

   ## Set a time grid
   time_grid  <- min(data$Survival) : min(max(data$Survival), tau)
   time_scale <- 1/diff(range(time_grid))

   ## Create estimated survival curve matrix: [time x clusters]
   H <- matrix(NA, nrow=length(time_grid), ncol=K)
   colnames(H) <- paste('Cluster_', ClusterLabel, sep='')

   ## Estimate the KM curve clusterwise on a common support
   for(k in 1:K){
      ## Compute Kaplan-Meier estimator on the kth cluster
      km    <- survfit(Surv(Survival, Death)~1, data = data[cluster==ClusterLabel[k] , ] )
      ## Construct the KM estimator function
      KMfun <- stepfun(x=km$time[-1], y=km$surv)
      H[,k] <- KMfun(time_grid)
   }

   ## construct matrix of pairwise L1 distances
   D <- matrix(0, ncol=K, nrow=K)
   for (i in 1:K){
      for(j in 1:K)
         if(i!=j){
            D[i,j] <- D[j,i] <- sum( abs( H[ , i]  -  H[ , j] ))
         }
   }
   ## Some scaling is given so that these numbers are somewhow interpretable
   ## for the same number of clusters independently of the time interval
   D <- D * time_scale


   ## Metric 1: min pairwise L1 distance
   iut <- which(upper.tri(D, diag=FALSE))
   ans$L1min <- min(D[iut])

   ## Metric 2: compute the summed L1 distance of each
   ## cluster to the nearest one
   diag(D)   <- NA
   ans$L1sum <- sum( D[iut] ) / length(iut)


   return(ans)
}







## Computes the discrepancy between survival curves in temrs of RMST
## (restricted mean survival time).
##
## For all these measures: higher ==> more separation
##
## Inputs:
##    data = data.frame with samples units on rows and two colums with names
##           'Survival' that is  survivial time
##           'Death'    that is a 0-1 variable indicating death
##           tau        truncation time, default 5 years, but if the min is less
##                      this is set to minimum across groups
##
## Outputs: a list with multiple separation measures
##
##
rmst.separation <- function(data, cluster, tau) {
   ans <- list()
   tmp           <- table(cluster)
   ClusterLabel  <- as.numeric(names(tmp))
   ClusterSize   <- as.numeric(tmp)
   K             <- length(ClusterLabel)
   n             <- sum(ClusterSize)


   ## Compute the minimum of the largest observed time in each of the two groups
   max.time <- rep(0, K)
   for (k in 1:K) {
      max.time[k] <- max(data$Survival[cluster == ClusterLabel[k]])
   }
   TAU  <- min(max.time, tau)


   ## Names
   ##    * RMST = restricted mean survival time
   ##    * LER  = life expectancy ratio
   ##
   ## LER: Life Expectancy Ratio  Matrix
   ##      LER[i,j] = LER[j,i] = max{RMST[i] / RMST[j], RMST[j] / RMST[i]}
   ##      note that here we don't have a baseline group so we define the ratio
   ##      always using in the denominator the group that have smaller RMST
   ##
   ## LED: Life Expectancy Difference
   ##    LED[i,j] = LED[j,i] = abs(RMST[i] - RMST[j])
   ##    note that here we don't have a baseline group so we define tha abs difference
   ##
   LER <- LED <-  matrix(0, ncol = K, nrow = K)
   for (i in 1:K) {
      for (j in 1:K)
         if (i != j) {
            ## First select data from  the two groups
            idx <- { cluster == ClusterLabel[i] | cluster == ClusterLabel[j]  }
            x   <- data[idx,]
            ##  Create a 0-1 vector, with gr==1 if cluster==ClusterLabel[i]
            gr0  <- cluster[idx]
            gr   <- ifelse(gr0 == ClusterLabel[i], 1, 0)
            u    <- rmst2(time = x$Survival, status = x$Death, arm = gr, tau = TAU)

            rmst_i <- u$RMST.arm1$rmst[1]
            rmst_j <- u$RMST.arm0$rmst[1]

            LER[i,j]  <- LER[j, i] <- max(rmst_i / rmst_j, rmst_j / rmst_i  )
            LED[i, j] <- LED[j, i] <- abs(rmst_i - rmst_j)
         }
   }

   ## index of the upper triangle
   iut <- which(upper.tri(LER, diag = FALSE))

   ## metric: min of pairwise LER discrepancy
   ans$LERmin <- min(LER[iut])

   ## metric: scaled summed pairwise LER discrepancy
   ans$LERsum <- sum(LER[iut]) / length(iut)

   ## metric: min of pairwise LED discrepancy
   ans$LEDmin <- min(LED[iut])

   ## metric: scaled summed pairwise LED discrepancy
   ans$LEDsum <- sum(LED[iut]) / length(iut)


   return(ans)
}

###

survPval = function(SurvDat,CLg,nYears=5){
  
  fit <- survfit(Surv(Survival, Death) ~ CLg,data = SurvDat,subset = Survival < (365 * nYears))
  suv <- survminer::ggsurvplot(fit, risk.table = TRUE, risk.table.height = 0.5,
                               xlim = c(0,5000), break.time.by = 500, pval = TRUE)
  pVal <- survminer:::surv_pvalue(fit,method = "survdiff",data=SurvDat)
  return(list(fit=fit,suv = suv,pVal=pVal))
  
}
