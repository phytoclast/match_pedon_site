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

#sil analysis by branch
d <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
t <- flexbeta(d, beta= -0.15)
k = 2
for(k in 2:10){#k=2
  labels = names(d)
  sil <- (t %>% cutree(k=k) %>% silhouette(d)) 
  sil <- as.data.frame(cbind(labels = labels, cluster=sil[,1], sil0=sil[,3]))
  sil$sil0 <- as.numeric(sil$sil0)
  sil <- sil %>% group_by(cluster) %>% mutate(sil.mean = mean(sil0))
  sil0 <- sil[,c('cluster','sil.mean')]
  colnames(sil0) <- c(paste0('cluster.',k),paste0('sil.',k))
  if(k==2){sil.table = cbind(labels = labels, sil0)}else{sil.table = cbind(sil.table,sil0)}
}
sil.table <- sil.table[c(1,
                         unique(floor(2:ncol(sil.table)/2)*2)
                         ,
                         unique(floor(2:ncol(sil.table)/2)*2)+1
)]

if (T){
  a <- 'flex20'
  k=4
  d <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  t <- flexbeta(d, beta= -0.20)
  makeplot(a,d,t,k)
}
groups = cutreeHybrid(t, minClusterSize = 3, distM = as.matrix(d))$labels


#dynamicTreeCut

deepSplit <- 1
maxCoreScatter <- 0.61  
minGap <- (1 - maxCoreScatter) * 3/4
for(i in 1:400){#i=2
  maxCoreScatter <- i/400  
  minGap <- (1 - maxCoreScatter) * 3/4
  
groups <- cutreeDynamic(t, minClusterSize=3, method = 'hybrid', distM=as.matrix(d), 
              deepSplit=deepSplit, maxCoreScatter=maxCoreScatter, minGap=minGap, maxAbsCoreScatter=NULL, minAbsGap=NULL)
nclust = length(unique(groups))
maxcore.table0 = as.data.frame(cbind(maxcor=maxCoreScatter, nclust)) 
if(i==1){maxcore.table <- maxcore.table0}else{maxcore.table <- rbind(maxcore.table, maxcore.table0)}
}
maxcore.table <-  maxcore.table %>% group_by(nclust) %>% summarise(maxcor = mean(maxcor))
