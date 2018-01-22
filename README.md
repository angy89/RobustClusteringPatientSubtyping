# RobustClusteringPatientSubtyping

This is the R implementation of a robust clustering of noisy high-dimensional gene
expression data technique for patients subtyping

### Demo
```R
library(otrimle) 
library(survRM2)
library(survminer)
library(survival)
require(org.Hs.eg.db)
library(limma)
require(clusterProfiler)

source("robust_clustering.R")
source("robust_scaling.R")
source("RSC_thresholding.R")
source("survival_separation_5y.R")
source("suvival_analysis_on_clusterings.R")
source("gene_mapping.R")
source("pathway_analysis.R")

# Example with the LUNG DATASET
nYears = 5
tau = 1825

# Read survival data
survivalData = read.table("data/LUNG_Survival.txt",sep="\t",header=TRUE)

#load gene expression dataset and the already computed robust correlation matrix
load("data/LUNG_rob_cor.RData")

# remove patients with survival longer than 5 years
toRem = which(survivalData$Survival>1825)
if(length(toRem)>0) survivalData =  survivalData[-toRem,]
if(length(toRem)>0) red_data = red_data[-toRem,]

## Thresholding
RSC=  RSC_thresholding(FLOSSES,R)

#Compute spectral decomposition 
scaled = robust_scaling(dataset = red_data,RSC = RSC)

#number of clusters as in the SNF paper
K = 4 

# clusters with less than th patients are removed
th = (nrow(survivalData) * 3/100)/100

#Compute all combination of D and ERC
A   <- expand.grid(nProjections = 11, gamma=10) #Optimal values obtained during the analysis
 
RC = robust_clustering(A,X=scaled$X,V=scaled$V,K,nCores=4)

source("good_surv_script.R")  

#plot survival results
finalSurvRes$survRes

#clustering information
View(finalSurvRes$FDS)

FR = pathway_analysis(cluster=RC[[1]]$cluster,red_data)
View(FR$FR)

#Dotplot of relevant pathways
FR$dp  

```R
