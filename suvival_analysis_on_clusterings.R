# RC is the list with all the clustering results of the otrimle function
# A is the matrix that contains the combination of the clustering parameters

# survivalData is the matrix with the survival information. 
#         It has 3 columns: the ID of the patient (PatientID), the days of survival (Survival), and a 0/1 column meaning death or not (Death)

survival_analysis = function(RC,A,survivalData,nYears,th){
  nYears = nYears
  library(survival)
  survRes = list()
  #  th = 0.1
  
  pb = txtProgressBar(min = 0,max = length(RC), style=3)
  FinalDataset = rep(0,7 + K)
  
  toRemIdx = c()
  
  for(i in 1:length(RC)){
    #cat(i,"\n")
    
    if(!is.null(RC[[i]]$cluster) ){
      ClustGroup = RC[[i]]$cluster
      tabClust= table(ClustGroup)
      
      if(sum(RC[[i]]$cpr<th)>0){ #se almeno un cluster ha una percentuale bassa di oggetti lo escludo
        SurvDat = survivalData
        CLg = ClustGroup
        suvAn = survPval(SurvDat,CLg)
        pVal = suvAn$pVal
        suv = suvAn$suv
        
        if(sum(tabClust==1)>0){
          idx = names(tabClust)[which(tabClust==1)]
          remIDx = which(ClustGroup %in% as.numeric(idx))
          LNormDist = unlist(survival.curves.separation(survivalData[-remIDx,], ClustGroup[-remIDx],tau = nYears))
        }else{
          LNormDist = unlist(survival.curves.separation(survivalData, ClustGroup,tau = nYears))
        }
        
        survRes[[i]] = list(suv = suv,pVal = pVal,LNormDist=LNormDist)
        
        if(length(tabClust)>K){
          newRow = c(unlist(A[i,]),pVal$pval,LNormDist,tabClust,i)  
        }else{
          newRow = c(unlist(A[i,]),pVal$pval,LNormDist,0,tabClust,i)
        }
        FinalDataset = rbind(FinalDataset, newRow)
        toRemIdx = c(toRemIdx,"Small Cluster")
        
        #survRes[[i]] = NULL
        #  newRow = c(unlist(A[i,]),rep(0,K+4),i)
        # FinalDataset = rbind(FinalDataset,newRow)
      }
      else{
        SurvDat = survivalData
        CLg = ClustGroup
        suvAn = survPval(SurvDat,CLg)
        pVal = suvAn$pVal
        suv = suvAn$suv
        
        if(sum(tabClust==1)>0){
          idx = which(tabClust==1)
          remIDx = which(ClustGroup==(idx-1))
          LNormDist = unlist(survival.curves.separation(survivalData[-remIDx,], ClustGroup[-remIDx]))
        }else{
          LNormDist = unlist(survival.curves.separation(survivalData, ClustGroup))
        }
        
        survRes[[i]] = list(suv = suv,pVal = pVal,LNormDist=LNormDist)
        
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
      newRow = c(unlist(A[i,]),rep(0,K+4),i)
      FinalDataset = rbind(FinalDataset,newRow)
    }
    
    setTxtProgressBar(pb,i)
  }
  close(pb)
  
  
  FinalDataset = FinalDataset[-1,]
  FinalDataset = cbind(FinalDataset,toRemIdx)
  colnames(FinalDataset)[3:ncol(FinalDataset)] = c("PVal","L1min","L1sum","Noise",paste("Cluster",1:K,sep=""),"Comb","ClustEvaluation")
  FinalDataset = as.matrix(FinalDataset)
  rownames(FinalDataset) = NULL
  FinalDataset = FinalDataset[order(unlist(FinalDataset[,3]),decreasing=FALSE),]
  
  goodIdx = which(FinalDataset[,"ClustEvaluation"] == "Good Clustering")
  FD =FinalDataset[goodIdx,]
  View(FD)
  
  goodPvalIDx = which(FD[,"PVal"]<0.05)
  FDS = FD[goodPvalIDx,] 
  View(FDS)
  
  return(list(FinalDataset = FinalDataset,FD = FD,FDS =FDS,survRes=survRes))
}