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


#group parameters ----
beta= -0.20
k = 6
d <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
t <- d %>% flexbeta(beta= beta)
tward <- agnes(d, method = 'ward') %>% as.hclust()
dk <- dynamicK(t, d)
dkward <- dynamicK(tward, d)
groups <- cutree(t, k = k)
maxCoreScatter <- dk[dk$nclust %in% k,]$maxcor 
minGap <- (1 - maxCoreScatter) * 3/4
dyngroups <- 
  cutreeDynamic(tbraydynam, minClusterSize=1, method = 'hybrid', distM=as.matrix(distbray),deepSplit=1, maxCoreScatter=maxCoreScatter, minGap=minGap, maxAbsCoreScatter=NULL, minAbsGap=NULL)
groups=dyngroups

source('groupplotsummary.R') 


Com.Structure[order(as.numeric(as.character(Com.Structure$cluster))),c("cluster", "association", "WetStructure")]




















if (T){
  a1 <- 'bray-flex20'
  k=6
  d <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  t1 <- flexbeta(d, beta= -0.20)
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


beta= -0.20
d <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
t <- d %>% flexbeta(beta= beta)
tward <- agnes(d, method = 'ward') %>% as.hclust()
dk <- dynamicK(t, d)
dkward <- dynamicK(tward, d)
timeA = Sys.time()
indob <- indanalysisdynamicward(plotdata,dk,dkward, beta)
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
  k=5
  d <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  t <- flexbeta(d, beta= -0.20)
  makeplot(a,d,t,k)
}
if (T){
  a <- 'flex20dynamic'
  k=6
  d <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  t <- flexbeta(d, beta= -0.20)
  makeplotdynamic(a,d,t,k)
}
if (T){
  a <- 'warddynamic'
  k=6
  d <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  t <- agnes(d, method='ward')
  makeplotdynamic(a,d,t,k)
}
groups = cutreeHybrid(t, minClusterSize = 3, distM = as.matrix(d))$labels

dk=dynamicK(t,d)



agnes(d, method = 'average')





cor(cophenetic(agnes(d, method = 'average')), cophenetic(agnes(d, method = 'average')))
cor(cophenetic(flexbeta(d, beta= -0)), cophenetic(agnes(d, method = 'average')))
cor(cophenetic(flexbeta(d, beta= -0.05)), cophenetic(agnes(d, method = 'average')))
cor(cophenetic(flexbeta(d, beta= -0.10)), cophenetic(agnes(d, method = 'average')))
cor(cophenetic(flexbeta(d, beta= -0.15)), cophenetic(agnes(d, method = 'average')))
cor(cophenetic(flexbeta(d, beta= -0.20)), cophenetic(agnes(d, method = 'average')))
cor(cophenetic(flexbeta(d, beta= -0.25)), cophenetic(agnes(d, method = 'average')))
cor(cophenetic(flexbeta(d, beta= -0.35)), cophenetic(agnes(d, method = 'average')))
cor(cophenetic(flexbeta(d, beta= -0.40)), cophenetic(agnes(d, method = 'average')))

cor(cophenetic(agnes(d, method = 'average')), cophenetic(agnes(d, method = 'ward')))
cor(cophenetic(flexbeta(d, beta= -0)), cophenetic(agnes(d, method = 'ward')))
cor(cophenetic(flexbeta(d, beta= -0.05)), cophenetic(agnes(d, method = 'ward')))
cor(cophenetic(flexbeta(d, beta= -0.10)), cophenetic(agnes(d, method = 'ward')))
cor(cophenetic(flexbeta(d, beta= -0.15)), cophenetic(agnes(d, method = 'ward')))
cor(cophenetic(flexbeta(d, beta= -0.20)), cophenetic(agnes(d, method = 'ward')))
cor(cophenetic(flexbeta(d, beta= -0.225)), cophenetic(agnes(d, method = 'ward')))
cor(cophenetic(flexbeta(d, beta= -0.25)), cophenetic(agnes(d, method = 'ward')))
cor(cophenetic(flexbeta(d, beta= -0.275)), cophenetic(agnes(d, method = 'ward')))
cor(cophenetic(flexbeta(d, beta= -0.30)), cophenetic(agnes(d, method = 'ward')))
cor(cophenetic(flexbeta(d, beta= -0.35)), cophenetic(agnes(d, method = 'ward')))
cor(cophenetic(flexbeta(d, beta= -0.40)), cophenetic(agnes(d, method = 'ward')))
cor(cophenetic(agnes(d, method = 'ward')), cophenetic(agnes(d, method = 'ward')))

