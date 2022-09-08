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
library(dendsort)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load Veg tables ---- 
source('processplot.R') 
#load clustering functions ----
source('clusterfunctions.R') 


#run seperate open vs forest analyses

newnames <- paste0('t',str_pad(round(overstorycover$overstorycover, 0), 2,'left',0),'s',str_pad(round(overstorycover$shrubcover, 0), 2,'left',0), overstorycover$soilplot)

plotdatax <- plotdata
rownames(plotdatax) <- newnames
shrubplots <- subset(overstorycover, shrubcover >= 10 & overstorycover < 10)$soilplot
plotdata.shrub <- subset(plotdata, rownames(plotdata) %in%  shrubplots )
openplots <- subset(overstorycover, shrubcover < 10 & overstorycover < 10)$soilplot
plotdata.open <- subset(plotdata, rownames(plotdata) %in%  openplots )
openplots <- subset(overstorycover, shrubcover < 10 & overstorycover < 10)$soilplot
plotdata.open <- subset(plotdata, rownames(plotdata) %in%  c(openplots,shrubplots) )
forestplots <- subset(overstorycover, overstorycover >= 10)$soilplot
plotdata.forest <- subset(plotdata, rownames(plotdata) %in%  forestplots )

plotdata.open.normal <- cleanplotdata(plotdata.open)
plotdata.forest.normal <- cleanplotdata(plotdata.forest)

timeA = Sys.time()
indob <- indanalysis2(plotdata.open.normal)
Sys.time() - timeA  

open.ind.table <- indob[[1]]
open.dni.table <- indob[[2]]
open.weak.table <- indob[[3]]
open.kaew.table <- indob[[4]]
Sys.time() - timeA  

timeA = Sys.time()
indob <- indanalysis2(plotdata.forest.normal)
Sys.time() - timeA  

forest.ind.table <- indob[[1]]
forest.dni.table <- indob[[2]]
forest.weak.table <- indob[[3]]
forest.kaew.table <- indob[[4]]
Sys.time() - timeA  


#forest ----
p1=plotdata.forest.normal
d <- vegdist(p1, method='bray', binary=FALSE, na.rm=T)
k = 5
t  <- d %>% flexbeta(beta= -0.15) %>% as.hclust() %>% dendsort()
groups <- cutree(t, k = k)
groups <- grouporder(t, groups)
if (T){
  a <- 'forest'
  makeplotgroup(a,d,t,groups)
}

indicators <- indgroup(p1, groups, F)
indicators2 <- indgroup2(p1, groups, F)
iplot <- indicators2[[1]]
gplot <- indicators2[[2]]
aplot <- indicators2[[3]]

source('groupplotsummary.R') 
source('USNVC_compare_specieslists_loop_by_cluster.R') 

Com.Structure[order(as.numeric(as.character(Com.Structure$cluster))),c("cluster", "association", "WetStructure")]
plotassociations[order(as.numeric(as.character(plotassociations$clust))),c("clust", "scientificname")]


#open ----
p1=plotdata.open.normal
d <- vegdist(p1, method='kulczynski', binary=FALSE, na.rm=T)
k = 5
t  <- d %>% agnes(method = 'ward') %>% as.hclust() %>% dendsort()
groups <- cutree(t, k = k)
groups <- grouporder(t, groups)
if (T){
  a <- 'open'
  makeplotgroup(a,d,t,groups)
}

indicators <- indgroup(p1, groups, F)
indicators2 <- indgroup2(p1, groups, F)
iplot <- indicators2[[1]]
gplot <- indicators2[[2]]
aplot <- indicators2[[3]]

source('groupplotsummary.R') 
source('USNVC_compare_specieslists_loop_by_cluster.R') 

Com.Structure[order(as.numeric(as.character(Com.Structure$cluster))),c("cluster", "association", "WetStructure")]
plotassociations[order(as.numeric(as.character(plotassociations$clust))),c("clust", "scientificname")]



#group parameters ----
# beta= -0.25
# k = 4
# d <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
# tbeta <- d %>% flexbeta(beta= beta) %>% dendsort()
# tward <- agnes(d, method = 'ward') %>% as.hclust() %>% dendsort()

# dk <- dynamicK(t, d)
# dkward <- dynamicK(tward, d)
# groups <- cutree(t, k = k)
# maxCoreScatter <- dk[dk$nclust %in% k,]$maxcor 
# minGap <- (1 - maxCoreScatter) * 3/4
# dyngroups <- 
#   cutreeDynamic(t, minClusterSize=1, method = 'hybrid', distM=as.matrix(d),deepSplit=1, maxCoreScatter=maxCoreScatter, minGap=minGap, maxAbsCoreScatter=NULL, minAbsGap=NULL)
# groups=dyngroups
# groups <- cutree(t, k = k)
# groups <- grouporder(t, groups)
# grouporder <- function(t,groups){
# soilplot <- names(d)
# clust <- unname(groups)
# groupdf <- as.data.frame(cbind(soilplot, clust))
# groupdf$clust <- (as.numeric(as.character(groupdf$clust)))
# torder <- as.data.frame(cbind(trow=t$order))
# torder$torder <- row(torder)[,1]
# newlabels <- as.data.frame(cbind(labels=t$labels))
# newlabels$trow <- row(newlabels)[,1]
# newlabels <- merge(newlabels, groupdf, by.x='labels', by.y ='soilplot')
# newlabels <- merge(newlabels, torder, by='trow')
# 
# grouporder <- newlabels %>% group_by(clust)  %>%  summarise(thisorder = min(torder))
# grouporder$newgroup <- order(order(grouporder$thisorder))
# newlabels <- merge(newlabels, grouporder, by='clust')
# newgroups <- newlabels$newgroup
# names(newgroups) <- newlabels$labels
# return(newgroups)
# }
# newlabels$newlabels <- paste(newlabels$clust, newlabels$newlabels)
# newlabels <- newlabels[order(newlabels$row),1]
# newtree <- t
# newtree$labels <- newlabels


beta= -0.25
k = 3
d <- vegdist(p.normal, method='bray', binary=FALSE, na.rm=T)
tbeta <- d %>% flexbeta(beta= beta) %>% dendsort()
tward <- agnes(d, method = 'ward') %>% as.hclust() %>% dendsort()

t <- tbeta
groups <- cutree(t, k = k)
groups <- grouporder(t, groups)
if (T){
  a <- 'flex'
  makeplotgroup(a,d,t,groups)
}

t <- tward
groups <- cutree(t, k = k)
groups <- grouporder(t, groups)
if (T){
  a <- 'ward'
  makeplotgroup(a,d,t,groups)
}

distsim <- as.dist(simil(p.normal,method='Simpson'))
distkulc <- vegdist(p.normal, method='kulczynski', binary=FALSE, na.rm=T)
distbray <- vegdist(p.normal, method='bray', binary=FALSE, na.rm=T)
tbrayagnes <- distbray %>% agnes(method = 'average') %>% as.hclust() %>% dendsort()
tsimpagnes <- distsim %>% agnes(method = 'average') %>% as.hclust() %>% dendsort()
tkulcagnes <- distkulc %>% agnes(method = 'average')%>% as.hclust() %>% dendsort()
tbrayflex25 <- distbray %>% flexbeta(beta= -0.25) %>% as.hclust() %>% dendsort()
tbrayflex15 <- distbray %>% flexbeta(beta= -0.15) %>% as.hclust() %>% dendsort()
tbrayward <- distbray %>% agnes(method = 'ward') %>% as.hclust() %>% dendsort()
tbraydiana <- distbray %>% diana() %>% as.hclust() %>% dendsort()



if (T){
  a <- 'kulcagnes'
  t <- tkulcagnes 
  d <- distkulc
  groups <- cutree(t, k = k)
  groups <- grouporder(t, groups)
  makeplotgroup(a,d,t,groups)}

if (T){
  a <- 'simpagnes'
  t <- tsimpagnes 
  d <- distsim
  groups <- cutree(t, k = k)
  groups <- grouporder(t, groups)
  makeplotgroup(a,d,t,groups)
}

if (T){
  a <- 'brayagnes'
  t <- tbrayagnes
  d <- distbray
  groups <- cutree(t, k = k)
  groups <- grouporder(t, groups)
  
  makeplotgroup(a,d,t,groups)}

if (T){
  a <- 'brayflex25'
  t <- tbrayflex25
  d <- distbray
  groups <- cutree(t, k = k)
  groups <- grouporder(t, groups)
  makeplotgroup(a,d,t,groups)}

if (T){
  k=3
  a <- 'brayward'
  t <- tbrayward
  d <- distbray
  groups <- cutree(t, k = k)
  groups <- grouporder(t, groups)
  makeplotgroup(a,d,t,groups)
  }

if (F){
  k=6
  a <- 'brayflex15'
  t <- tbrayflex15
  d <- distbray
  groups <- cutree(t, k = k)
  groups <- grouporder(t, groups)
  makeplotgroup(a,d,t,groups)
}

if (F){
  k=6
  a <- 'braydiana'
  t <- tbraydiana
  d <- distbray
  groups <- cutree(t, k = k)
  groups <- grouporder(t, groups)
  makeplotgroup(a,d,t,groups)
}



  
  
indicators <- indgroup(p1, groups, F)
indicators2 <- indgroup2(p1, groups, F)
iplot <- indicators2[[1]]
gplot <- indicators2[[2]]
aplot <- indicators2[[3]]

#dclust <- clustvar(d, groups)


source('groupplotsummary.R') 
source('USNVC_compare_specieslists_loop_by_cluster.R') 

Com.Structure[order(as.numeric(as.character(Com.Structure$cluster))),c("cluster", "association", "WetStructure")]
plotassociations[order(as.numeric(as.character(plotassociations$clust))),c("clust", "scientificname")]


timeA = Sys.time()
indob <- indanalysis2(p.normal)
Sys.time() - timeA  

ind.table <- indob[[1]]
dni.table <- indob[[2]]
weak.table <- indob[[3]]
kaew.table <- indob[[4]]
Sys.time() - timeA  


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
  k=5
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

