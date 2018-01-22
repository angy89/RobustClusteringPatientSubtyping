library(survival)
survRes = list()
nYears = 5
tau = 1825
pb = txtProgressBar(min = 0,max = length(RC), style=3)
FinalDataset = rep(0,11 + K)
toRemIdx = c()

for(i in 1:length(RC)){

  if(!is.null(RC[[i]]$cluster) ){# if result is not null we can continue the analyses
    ClustGroup = RC[[i]]$cluster # i-th clustering
    tabClust= table(ClustGroup)  # table
    # if a clustering contain a cluster with few object it is marked bad clustering ("small cluster")
    # RC[[i]]$cpr contains the percentage of patients in each cluster while RC[[i]]$npr contains the percentage of patients into the noise
    # to be marked as good clustering both the regular nor the noise cluster must have enough patients (at least two patients for each cluster)
    if(sum(RC[[i]]$cpr<th)>0 || (RC[[i]]$npr<th && sum(RC[[i]]$cpr) < 0.9999999) ){ 
      if(length(tabClust)>K){
        newRow = c(unlist(A[i,]),rep("-",7),tabClust,i)  
      }else{
        newRow = c(unlist(A[i,]),rep("-",7),0,tabClust,i)
      }
      FinalDataset = rbind(FinalDataset, newRow)
      toRemIdx = c(toRemIdx,"Small Cluster")
      survRes[[i]] = NULL
    }else{
      SurvDat = survivalData
      CLg = ClustGroup
      suvAn = survPval(SurvDat,CLg,nYears)
      pVal = suvAn$pVal
      suv = suvAn$suv
      fit = suvAn$fit
      LNormDist = c(unlist(survival.curves.separation(SurvDat, CLg,tau)),
                    unlist(rmst.separation(SurvDat, CLg,tau)))
      
      survRes[[i]] = list(fit = fit,suv = suv,pVal = pVal,LNormDist=LNormDist)
      
      if(length(tabClust)>K){
        newRow = c(unlist(A[i,]),pVal$pval,LNormDist,tabClust,i)  
      }else{
        newRow = c(unlist(A[i,]),pVal$pval,LNormDist,0,tabClust,i)
      }
      FinalDataset = rbind(FinalDataset, newRow)
      toRemIdx = c(toRemIdx,"Good Clustering")
    }
  }else{
    toRemIdx = c(toRemIdx,"Null Result")
    survRes[[i]] = NULL
    newRow = c(unlist(A[i,]),rep(0,K+8),i)
    FinalDataset = rbind(FinalDataset,newRow)
  }
  
  setTxtProgressBar(pb,i)
}
close(pb)

#FinalDataset = FinalDataset[-1,]
FinalDataset = cbind(FinalDataset,toRemIdx)
colnames(FinalDataset)[3:ncol(FinalDataset)] = c("PVal","L1min","L1sum","LERmin","LERsum","LEDmin","LEDsum","Noise",paste("Cluster",1:K,sep=""),"Comb","ClustEvaluation")
FinalDataset = as.matrix(FinalDataset)
rownames(FinalDataset) = NULL
FinalDataset = FinalDataset[order(unlist(FinalDataset[,3]),decreasing=FALSE),]

goodIdx = which(FinalDataset[,"ClustEvaluation"] == "Good Clustering")
FD =FinalDataset[goodIdx,c(1,2,3,8,10:ncol(FinalDataset))]


goodPvalIDx = which(FD[,"PVal"]<0.05)

if(length(goodPvalIDx)==0){
  FDS = matrix("",nrow=1,ncol=ncol(FD))
  colnames(FDS)=colnames(FD)
  FDS[1,"PVal"]=1
  FDS[1,c("L1min","L1sum","LERmin","LERsum","LEDmin","LEDsum")]=0
  
}else{
  if(length(goodPvalIDx)==1){
    FDS = matrix("",nrow=1,ncol=ncol(FD))
    FDS[1,]=FD[goodPvalIDx,]
    colnames(FDS)=colnames(FD)
  }else{
    FDS = FD[goodPvalIDx,] 
    
  }  
}


finalSurvRes = list(FinalDataset = FinalDataset[-1,],FD = FD[-1,],FDS =FDS[-1,],survRes=survRes)
