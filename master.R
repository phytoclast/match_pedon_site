library(soilDB)
library(aqp)#load before dplyr
library(stringr)
library(BiodiversityR)
library(cluster)
library(ape)
library(dendextend)
library(dplyr)
library(dynamicTreeCut)
library(rpart)
library(rpart.plot)
library(goeveg)
library(proxy)
library(foreign)
library(optpart)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load Veg tables ---- 
source('processplot.R') 
#load clustering functions ----
source('clusterfunctions.R') 


if (T){
  a1 <- 'bray-flex25'
  k=8
  d <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  t1 <- flexbeta(d, beta= -0.25)
  makeplot(a1,d,t1,k)
}
timeA = Sys.time()
indob <- indanalysis(plotdata)
Sys.time() - timeA  

ind.table <- indob[[1]]
dni.table <- indob[[2]]
clu.table <- indob[[3]]
ulc.table <- indob[[4]]
clind.table <- indob[[5]]
dnilc.table <- indob[[6]]
sil.table <- indob[[7]]
lis.table <- indob[[8]]
Sys.time() - timeA  

