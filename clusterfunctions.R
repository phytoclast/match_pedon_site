
grouporder <- function(t,groups){
  t=as.hclust(t)
  soilplot <- names(groups)
  clust <- unname(groups)
  groupdf <- as.data.frame(cbind(soilplot, clust))
  groupdf$clust <- (as.numeric(as.character(groupdf$clust)))
  torder <- as.data.frame(cbind(trow=t$order))
  torder$torder <- row(torder)[,1]
  newlabels <- as.data.frame(cbind(labels=t$labels))
  newlabels$trow <- row(newlabels)[,1]
  newlabels <- merge(newlabels, groupdf, by.x='labels', by.y ='soilplot')
  newlabels <- merge(newlabels, torder, by='trow')
  
  grouporder <- newlabels %>% group_by(clust)  %>%  summarise(thisorder = min(torder))
  grouporder$newgroup <- order(order(grouporder$thisorder, decreasing=T))
  newlabels <- merge(newlabels, grouporder, by='clust')
  newlabels <- newlabels[order(newlabels$trow),]
  newgroups <- newlabels$newgroup
  names(newgroups) <- newlabels$labels
  return(newgroups)
}

dynamicK <- function(t,d){
  deepSplit <- 1
  maxCoreScatter <- .5  
  minGap <- (1 - maxCoreScatter) * 3/4
  for(i in 1:400){#i=2
    maxCoreScatter <- i/400  
    minGap <- (1 - maxCoreScatter) * 3/4
    
    groups <- cutreeDynamic(t, minClusterSize=1, method = 'hybrid', distM=as.matrix(d), 
                            deepSplit=deepSplit, maxCoreScatter=maxCoreScatter, minGap=minGap, maxAbsCoreScatter=NULL, minAbsGap=NULL)
    nclust = length(unique(groups))
    maxcore.table0 = as.data.frame(cbind(maxcor=maxCoreScatter, nclust)) 
    if(i==1){maxcore.table <- maxcore.table0}else{maxcore.table <- rbind(maxcore.table, maxcore.table0)}
  }
  maxcore.table <-  maxcore.table %>% group_by(nclust) %>% summarise(maxcor = mean(maxcor))
  return(maxcore.table)}

makeplot <- function(a,d,t,k){
  filename <- paste0('output/Vegplot_',a,'.png')
  t <- as.hclust(t)
  #make cuts and reformat dendrogram
  ngroups=k
  groups <- cutree(t, k = ngroups)
  soilplot <- names(d)
  clust <- unname(groups)
  groupdf <- as.data.frame(cbind(soilplot, clust))
  groupdf$clust <- (as.numeric(as.character(groupdf$clust)))
  maxcluster <- max(groupdf$clust)
  numberzeros <- nrow(groupdf[(groupdf$clust == 0),])
  whichrecords <- which(groupdf$clust == 0)
  if (nrow(groupdf[groupdf$clust == 0,]) != 0){
    for (i in 1:numberzeros){ #assign all zero clusters to unique cluster number.
      groupdf[whichrecords[i],]$clust <- maxcluster+i}}
  
  newlabels <- t$labels
  newlabels <- as.data.frame(newlabels)
  newlabels$row <- row(newlabels)[,1]
  newlabels <- merge(newlabels, groupdf, by.x='newlabels', by.y ='soilplot')
  newlabels$newlabels <- paste(newlabels$clust, newlabels$newlabels)
  newlabels <- newlabels[order(newlabels$row),1]
  newtree <- t
  newtree$labels <- newlabels
  
  dend1 <- color_branches(as.dendrogram(as.hclust(newtree)), clusters = groups[order.dendrogram(as.dendrogram(t))])
  dend1 <- color_labels(dend1, col = get_leaves_branches_col(dend1))
  
  #output file
  
  w <- 800
  h <- nrow(plotdata)*12+80
  u <- 12
  png(filename=filename,width = w, height = h, units = "px", pointsize = u)
  
  par(mar = c(2,0,1,13))
  plot(dend1, horiz = TRUE, main=paste('floristic simularity', a,'method of', 'vegetation'), font=1, cex=0.84)
  #rect.dendrogram(dend1, k = ngroups, horiz = TRUE)
  dev.off()
  
}

makeplotgroup <- function(a,d,t,groups){
  filename <- paste0('output/Vegplot_',a,'.png')
  t <- as.hclust(t)
  #make cuts and reformat dendrogram
  groups <- groups
  
  soilplot <- names(d)
  clust <- unname(groups)
  groupdf <- as.data.frame(cbind(soilplot, clust))
  groupdf$clust <- (as.numeric(as.character(groupdf$clust)))
  maxcluster <- max(groupdf$clust)
  numberzeros <- nrow(groupdf[(groupdf$clust == 0),])
  whichrecords <- which(groupdf$clust == 0)
  if (nrow(groupdf[groupdf$clust == 0,]) != 0){
    for (i in 1:numberzeros){ #assign all zero clusters to unique cluster number.
      groupdf[whichrecords[i],]$clust <- maxcluster+i}}
  
  newlabels <- t$labels
  newlabels <- as.data.frame(newlabels)
  newlabels$row <- row(newlabels)[,1]
  newlabels <- merge(newlabels, groupdf, by.x='newlabels', by.y ='soilplot')
  newlabels$newlabels <- paste(newlabels$clust, newlabels$newlabels)
  newlabels <- newlabels[order(newlabels$row),1]
  newtree <- t
  newtree$labels <- newlabels
  
  dend1 <- color_branches(as.dendrogram(as.hclust(newtree)), clusters = groups[order.dendrogram(as.dendrogram(t))])
  dend1 <- color_labels(dend1, col = get_leaves_branches_col(dend1))
  
  #output file
  
  w <- 800
  h <- nrow(plotdata)*12+80
  u <- 12
  png(filename=filename,width = w, height = h, units = "px", pointsize = u)
  
  par(mar = c(2,0,1,13))
  plot(dend1, horiz = TRUE, main=paste('floristic simularity', a,'method of', 'vegetation'), font=1, cex=0.84)
  #rect.dendrogram(dend1, k = ngroups, horiz = TRUE)
  dev.off()
  
}

makeplotdynamic <- function(a,d,t,k){#k=2
  filename <- paste0('output/Vegplot_',a,'.png')
  t <- as.hclust(t)
  #make cuts and reformat dendrogram
  ngroups=k
  # groups <- cutree(t, k = ngroups)
  # groups <- kmeans(d, centers = ngroups)$cluster
  # groups = optind$clustering
  dk= dynamicK(t,d)
  
  deepSplit <- 1
  maxCoreScatter <- dk[dk$nclust %in% k,]$maxcor  
  minGap <- (1 - maxCoreScatter) * 3/4
  
  groups <- cutreeDynamic(t, minClusterSize=1, method="hybrid", distM=as.matrix(d), 
                          deepSplit=deepSplit, maxCoreScatter=maxCoreScatter, minGap=minGap, maxAbsCoreScatter=NULL, minAbsGap=NULL, pamRespectsDendro=T)
  
  soilplot <- names(d)
  clust <- unname(groups)
  groupdf <- as.data.frame(cbind(soilplot, clust))
  groupdf$clust <- (as.numeric(as.character(groupdf$clust)))
  maxcluster <- max(groupdf$clust)
  numberzeros <- nrow(groupdf[(groupdf$clust == 0),])
  whichrecords <- which(groupdf$clust == 0)
  if (nrow(groupdf[groupdf$clust == 0,]) != 0){
    for (i in 1:numberzeros){ #assign all zero clusters to unique cluster number.
      groupdf[whichrecords[i],]$clust <- maxcluster+i}}
  
  newlabels <- t$labels
  newlabels <- as.data.frame(newlabels)
  newlabels$row <- row(newlabels)
  newlabels <- merge(newlabels, groupdf, by.x='newlabels', by.y ='soilplot')
  newlabels$newlabels <- paste(newlabels$clust, newlabels$newlabels)
  newlabels <- newlabels[order(newlabels$row),1]
  newtree <- t
  newtree$labels <- newlabels
  
  dend1 <- color_branches(as.dendrogram(as.hclust(newtree)), clusters = groups[order.dendrogram(as.dendrogram(t))])
  dend1 <- color_labels(dend1, col = get_leaves_branches_col(dend1))
  
  #output file
  
  w <- 800
  h <- nrow(plotdata)*12+80
  u <- 12
  png(filename=filename,width = w, height = h, units = "px", pointsize = u)
  
  par(mar = c(2,0,1,13))
  plot(dend1, horiz = TRUE, main=paste('floristic simularity', a,'method of', 'vegetation'), font=1, cex=0.84)
  #rect.dendrogram(dend1, k = ngroups, horiz = TRUE)
  dev.off()
  
}

indanalysis <- function(plotdata){
plotdata.total <- apply(plotdata, MARGIN = 2, FUN = 'sum')
plotdata.total <- as.data.frame(cbind(name=names(plotdata.total),total=plotdata.total))
removetaxon <- plotdata.total[plotdata.total$total %in% 0,]$name
plotdata1 <- plotdata[,!colnames(plotdata) %in% removetaxon]

distbray <- vegdist(plotdata1, method='bray', binary=FALSE, na.rm=T)# dbray <- as.data.frame(as.matrix(distbray))
distjac <- vegdist(plotdata1, method='jaccard', binary=FALSE, na.rm=T)
distsim <- as.dist(simil(plotdata1,method='Simpson'))
distkulc <- vegdist(plotdata1, method='kulczynski', binary=FALSE, na.rm=T)
tbrayagnes <- distbray %>% agnes(method = 'average')
tbrayflex05 <- distbray %>% flexbeta(beta= -0.05)
tbrayflex10 <- distbray %>% flexbeta(beta= -0.10)
tbrayflex15 <- distbray %>% flexbeta(beta= -0.15)
tbrayflex20 <- distbray %>% flexbeta(beta= -0.20)
tbrayflex25 <- distbray %>% flexbeta(beta= -0.25)
tbrayflex30 <- distbray %>% flexbeta(beta= -0.30)
tbrayflex35 <- distbray %>% flexbeta(beta= -0.35)
tbrayward <- distbray %>% agnes(method = 'ward')
tbraydiana <- distbray %>% diana 
tsimpagnes <- distsim %>% agnes(method = 'average') 
tjacagnes <- distjac %>% agnes(method = 'average') 
tkulcagnes <- distkulc %>% agnes(method = 'average')
tkulcward <- distkulc %>% agnes(method = 'ward')



k <- 2
klevel <- 0

sil.upgma <- 0
sil.flex05 <- 0
sil.flex10 <- 0
sil.flex15 <- 0
sil.flex20 <- 0
sil.flex25 <- 0
sil.flex30 <- 0
sil.flex35 <- 0
sil.ward <- 0
sil.jac <- 0
sil.sim <- 0
sil.diana <- 0
sil.kmeans <- 0
sil.kulc <- 0
sil.kward <- 0


ind.upgma <- 0
ind.flex05 <- 0
ind.flex10 <- 0
ind.flex15 <- 0
ind.flex20 <- 0
ind.flex25 <- 0
ind.flex30 <- 0
ind.flex35 <- 0
ind.ward <- 0

clu.upgma <- 0
clu.flex05 <- 0
clu.flex10 <- 0
clu.flex15 <- 0
clu.flex20 <- 0
clu.flex25 <- 0
clu.flex30 <- 0
clu.flex35 <- 0
clu.ward <- 0

clind.upgma <- 0
clind.flex05 <- 0
clind.flex10 <- 0
clind.flex15 <- 0
clind.flex20 <- 0
clind.flex25 <- 0
clind.flex30 <- 0
clind.flex35 <- 0
clind.ward <- 0

timeA = Sys.time()
for (k in 2:10){
  
   ind.upgma.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayagnes)))
  ind.flex05.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex05)))
  ind.flex10.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex10)))
  ind.flex15.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex15)))
  ind.flex20.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex20)))
  ind.flex25.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex25)))
  ind.flex30.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex30)))
  ind.flex35.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex35)))
  ind.ward.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayward)))
  
  ind.upgma.1 <- ind.upgma.0$sig.spc
  ind.flex05.1 <- ind.flex05.0$sig.spc
  ind.flex10.1 <- ind.flex10.0$sig.spc
  ind.flex15.1 <- ind.flex15.0$sig.spc
  ind.flex20.1 <- ind.flex20.0$sig.spc
  ind.flex25.1 <- ind.flex25.0$sig.spc
  ind.flex30.1 <- ind.flex30.0$sig.spc
  ind.flex35.1 <- ind.flex35.0$sig.spc
  ind.ward.1 <- ind.ward.0$sig.spc
  
  clu.upgma.1 <- ind.upgma.0$sig.clust
  clu.flex05.1 <- ind.flex05.0$sig.clust
  clu.flex10.1 <- ind.flex10.0$sig.clust
  clu.flex15.1 <- ind.flex15.0$sig.clust
  clu.flex20.1 <- ind.flex20.0$sig.clust
  clu.flex25.1 <- ind.flex25.0$sig.clust
  clu.flex30.1 <- ind.flex30.0$sig.clust
  clu.flex35.1 <- ind.flex35.0$sig.clust
  clu.ward.1 <- ind.ward.0$sig.clust

  clind.upgma.1 <- ind.upgma.1*clu.upgma.1/k
  clind.flex05.1 <- ind.flex05.1*clu.flex05.1/k
  clind.flex10.1 <- ind.flex10.1*clu.flex10.1/k
  clind.flex15.1 <- ind.flex15.1*clu.flex15.1/k
  clind.flex20.1 <- ind.flex20.1*clu.flex20.1/k
  clind.flex25.1 <- ind.flex25.1*clu.flex25.1/k
  clind.flex30.1 <- ind.flex30.1*clu.flex30.1/k
  clind.flex35.1 <- ind.flex35.1*clu.flex35.1/k
  clind.ward.1 <- ind.ward.1*clu.ward.1/k
  
   sil.upgma.1 <- (tbrayagnes %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.flex05.1 <- (tbrayflex05 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.flex10.1 <- (tbrayflex10 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.flex15.1 <- (tbrayflex15 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.flex20.1 <- (tbrayflex20 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.flex25.1 <- (tbrayflex25 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.flex30.1 <- (tbrayflex30 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.flex35.1 <- (tbrayflex35 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.ward.1 <- (tbrayward %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  
  sil.jac.1 <- (tjacagnes %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.sim.1 <- (tsimpagnes %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.diana.1 <- (tbraydiana %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.kmeans.1 <- (kmeans(distbray, centers = k)$cluster %>% silhouette(distbray))[,3] %>% mean
  sil.kulc.1 <- (tkulcagnes %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.kward.1 <- (tkulcward %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  
klevel <- c(klevel, k)
  ind.upgma <- c(ind.upgma, ind.upgma.1)
  ind.flex05 <- c(ind.flex05, ind.flex05.1)
  ind.flex10 <- c(ind.flex10, ind.flex10.1)
  ind.flex15 <- c(ind.flex15, ind.flex15.1)
  ind.flex20 <- c(ind.flex20, ind.flex20.1)
  ind.flex25 <- c(ind.flex25, ind.flex25.1)
  ind.flex30 <- c(ind.flex30, ind.flex30.1)
  ind.flex35 <- c(ind.flex35, ind.flex35.1)
  ind.ward <- c(ind.ward, ind.ward.1)
  
  clu.upgma <- c(clu.upgma, clu.upgma.1)
  clu.flex05 <- c(clu.flex05, clu.flex05.1)
  clu.flex10 <- c(clu.flex10, clu.flex10.1)
  clu.flex15 <- c(clu.flex15, clu.flex15.1)
  clu.flex20 <- c(clu.flex20, clu.flex20.1)
  clu.flex25 <- c(clu.flex25, clu.flex25.1)
  clu.flex30 <- c(clu.flex30, clu.flex30.1)
  clu.flex35 <- c(clu.flex35, clu.flex35.1)
  clu.ward <- c(clu.ward, clu.ward.1)
  
  clind.upgma <- c(clind.upgma, clind.upgma.1)
  clind.flex05 <- c(clind.flex05, clind.flex05.1)
  clind.flex10 <- c(clind.flex10, clind.flex10.1)
  clind.flex15 <- c(clind.flex15, clind.flex15.1)
  clind.flex20 <- c(clind.flex20, clind.flex20.1)
  clind.flex25 <- c(clind.flex25, clind.flex25.1)
  clind.flex30 <- c(clind.flex30, clind.flex30.1)
  clind.flex35 <- c(clind.flex35, clind.flex35.1)
  clind.ward <- c(clind.ward, clind.ward.1)
  
  sil.upgma <- c(sil.upgma, sil.upgma.1)
  sil.flex05 <- c(sil.flex05, sil.flex05.1)
  sil.flex10 <- c(sil.flex10, sil.flex10.1)
  sil.flex15 <- c(sil.flex15, sil.flex15.1)
  sil.flex20 <- c(sil.flex20, sil.flex20.1)
  sil.flex25 <- c(sil.flex25, sil.flex25.1)
  sil.flex30 <- c(sil.flex30, sil.flex30.1)
  sil.flex35 <- c(sil.flex35, sil.flex35.1)
  sil.ward <- c(sil.ward, sil.ward.1)
  
  sil.jac <- c(sil.jac, sil.jac.1)
  sil.sim <- c(sil.sim, sil.sim.1)
  sil.diana <- c(sil.diana, sil.diana.1)
  sil.kmeans <- c(sil.kmeans, sil.kmeans.1)
  sil.kulc <- c(sil.kulc, sil.kulc.1)
  sil.kward <- c(sil.kward, sil.kward.1)
  
}
Sys.time() - timeA  

sil.table <- as.data.frame(cbind(klevel,sil.upgma,sil.flex05,sil.flex10,sil.flex15,sil.flex20,sil.flex25,sil.flex30,sil.flex35,sil.ward, 
                                 sil.jac,  sil.sim, sil.diana, sil.kmeans, sil.kulc, sil.kward))
sil.table <- sil.table[-1,]
lis.table <- sil.table[,-1] %>% t() %>% as.data.frame()
lis.table$s2to8 <- apply(lis.table[,2:7], MARGIN=1, FUN = 'mean')
lis.table <- lis.table %>% mutate(s2to8 = apply(lis.table[,1:7], MARGIN=1, FUN = 'mean'))
sil.table <- sil.table %>% mutate(mean = apply(sil.table[,2:16], MARGIN=1, FUN = 'mean'))

ind.table <- as.data.frame(cbind(klevel,ind.upgma,ind.flex05,ind.flex10,ind.flex15,ind.flex20,ind.flex25,ind.flex30,ind.flex35,ind.ward))
ind.table <- ind.table[-1,]
dni.table <- ind.table[,-1] %>% t() %>% as.data.frame()
dni.table <- dni.table %>% mutate(s2to8 = apply(dni.table[,1:7], MARGIN=1, FUN = 'mean'))

clu.table <- as.data.frame(cbind(klevel,clu.upgma,clu.flex05,clu.flex10,clu.flex15,clu.flex20,clu.flex25,clu.flex30,clu.flex35,clu.ward))
clu.table <- clu.table[-1,]
ulc.table <- clu.table[,-1] %>% t() %>% as.data.frame()
ulc.table <- ulc.table %>% mutate(s2to8 = apply(ulc.table[,1:7], MARGIN=1, FUN = 'mean'))

clind.table <- as.data.frame(cbind(klevel,clind.upgma,clind.flex05,clind.flex10,clind.flex15,clind.flex20,clind.flex25,clind.flex30,clind.flex35,clind.ward))
clind.table <- clind.table[-1,]
dnilc.table <- clind.table[,-1] %>% t() %>% as.data.frame()
dnilc.table <- dnilc.table %>% mutate(s2to8 = apply(dnilc.table[,1:7], MARGIN=1, FUN = 'mean'))
clind.table <- clind.table %>% mutate(mean = apply(clind.table[,2:10], MARGIN=1, FUN = 'mean'))


listofoutput <- list(ind.table, dni.table, clu.table, ulc.table, clind.table, dnilc.table, sil.table, lis.table)
return(listofoutput)}

indanalysisdynamic <- function(plotdata,dk,beta){
  plotdata.total <- apply(plotdata, MARGIN = 2, FUN = 'sum')
  plotdata.total <- as.data.frame(cbind(name=names(plotdata.total),total=plotdata.total))
  removetaxon <- plotdata.total[plotdata.total$total %in% 0,]$name
  plotdata1 <- plotdata[,!colnames(plotdata) %in% removetaxon]
  
  distbray <- vegdist(plotdata1, method='bray', binary=FALSE, na.rm=T)# dbray <- as.data.frame(as.matrix(distbray))
  distjac <- vegdist(plotdata1, method='jaccard', binary=FALSE, na.rm=T)
  distsim <- as.dist(simil(plotdata1,method='Simpson'))
  distkulc <- vegdist(plotdata1, method='kulczynski', binary=FALSE, na.rm=T)
  tbrayagnes <- distbray %>% agnes(method = 'average')
  tbrayflex05 <- distbray %>% flexbeta(beta= -0.05)
  tbrayflex10 <- distbray %>% flexbeta(beta= -0.10)
  tbrayflex15 <- distbray %>% flexbeta(beta= -0.15)
  tbrayflex20 <- distbray %>% flexbeta(beta= -0.20)
  tbrayflex25 <- distbray %>% flexbeta(beta= -0.25)
  tbrayflex30 <- distbray %>% flexbeta(beta= -0.30)
  tbrayflex35 <- distbray %>% flexbeta(beta= -0.35)
  tbrayward <- distbray %>% agnes(method = 'ward')
  tbraydiana <- distbray %>% diana 
  tsimpagnes <- distsim %>% agnes(method = 'average') 
  tjacagnes <- distjac %>% agnes(method = 'average') 
  tkulcagnes <- distkulc %>% agnes(method = 'average')
  tkulcward <- distkulc %>% agnes(method = 'ward')
  #dynamic
  tbraydynam <- distbray %>% flexbeta(beta= beta)
  # dk <- dynamicK(tbraydynam, distbray)
  
  
  k <- 2
  klevel <- 0
  
  sil.upgma <- 0
  sil.flex05 <- 0
  sil.flex10 <- 0
  sil.flex15 <- 0
  sil.flex20 <- 0
  sil.flex25 <- 0
  sil.flex30 <- 0
  sil.flex35 <- 0
  sil.ward <- 0
  sil.jac <- 0
  sil.sim <- 0
  sil.diana <- 0
  sil.kmeans <- 0
  sil.kulc <- 0
  sil.kward <- 0
  
  
  ind.upgma <- 0
  ind.flex05 <- 0
  ind.flex10 <- 0
  ind.flex15 <- 0
  ind.flex20 <- 0
  ind.flex25 <- 0
  ind.flex30 <- 0
  ind.flex35 <- 0
  ind.ward <- 0
  
  clu.upgma <- 0
  clu.flex05 <- 0
  clu.flex10 <- 0
  clu.flex15 <- 0
  clu.flex20 <- 0
  clu.flex25 <- 0
  clu.flex30 <- 0
  clu.flex35 <- 0
  clu.ward <- 0
  
  clind.upgma <- 0
  clind.flex05 <- 0
  clind.flex10 <- 0
  clind.flex15 <- 0
  clind.flex20 <- 0
  clind.flex25 <- 0
  clind.flex30 <- 0
  clind.flex35 <- 0
  clind.ward <- 0
  
  ind.dynamic <- 0
  clu.dynamic <- 0
  clind.dynamic <- 0
  sil.dynamic <- 0
  ind.dynamicward <- 0
  clu.dynamicward <- 0
  clind.dynamicward <- 0
  sil.dynamicward <- 0
  
  
  timeA = Sys.time()
  for (k in 2:10){
    
    ind.upgma.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayagnes)))
    ind.flex05.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex05)))
    ind.flex10.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex10)))
    ind.flex15.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex15)))
    ind.flex20.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex20)))
    ind.flex25.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex25)))
    ind.flex30.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex30)))
    ind.flex35.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex35)))
    ind.ward.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayward)))
    
    ind.upgma.1 <- ind.upgma.0$sig.spc
    ind.flex05.1 <- ind.flex05.0$sig.spc
    ind.flex10.1 <- ind.flex10.0$sig.spc
    ind.flex15.1 <- ind.flex15.0$sig.spc
    ind.flex20.1 <- ind.flex20.0$sig.spc
    ind.flex25.1 <- ind.flex25.0$sig.spc
    ind.flex30.1 <- ind.flex30.0$sig.spc
    ind.flex35.1 <- ind.flex35.0$sig.spc
    ind.ward.1 <- ind.ward.0$sig.spc
    
    clu.upgma.1 <- ind.upgma.0$sig.clust
    clu.flex05.1 <- ind.flex05.0$sig.clust
    clu.flex10.1 <- ind.flex10.0$sig.clust
    clu.flex15.1 <- ind.flex15.0$sig.clust
    clu.flex20.1 <- ind.flex20.0$sig.clust
    clu.flex25.1 <- ind.flex25.0$sig.clust
    clu.flex30.1 <- ind.flex30.0$sig.clust
    clu.flex35.1 <- ind.flex35.0$sig.clust
    clu.ward.1 <- ind.ward.0$sig.clust
    
    clind.upgma.1 <- ind.upgma.1*clu.upgma.1/k
    clind.flex05.1 <- ind.flex05.1*clu.flex05.1/k
    clind.flex10.1 <- ind.flex10.1*clu.flex10.1/k
    clind.flex15.1 <- ind.flex15.1*clu.flex15.1/k
    clind.flex20.1 <- ind.flex20.1*clu.flex20.1/k
    clind.flex25.1 <- ind.flex25.1*clu.flex25.1/k
    clind.flex30.1 <- ind.flex30.1*clu.flex30.1/k
    clind.flex35.1 <- ind.flex35.1*clu.flex35.1/k
    clind.ward.1 <- ind.ward.1*clu.ward.1/k
    
    sil.upgma.1 <- (tbrayagnes %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.flex05.1 <- (tbrayflex05 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.flex10.1 <- (tbrayflex10 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.flex15.1 <- (tbrayflex15 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.flex20.1 <- (tbrayflex20 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.flex25.1 <- (tbrayflex25 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.flex30.1 <- (tbrayflex30 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.flex35.1 <- (tbrayflex35 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.ward.1 <- (tbrayward %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    
    sil.jac.1 <- (tjacagnes %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.sim.1 <- (tsimpagnes %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.diana.1 <- (tbraydiana %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.kmeans.1 <- (kmeans(distbray, centers = k)$cluster %>% silhouette(distbray))[,3] %>% mean
    sil.kulc.1 <- (tkulcagnes %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.kward.1 <- (tkulcward %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    
    klevel <- c(klevel, k)
    ind.upgma <- c(ind.upgma, ind.upgma.1)
    ind.flex05 <- c(ind.flex05, ind.flex05.1)
    ind.flex10 <- c(ind.flex10, ind.flex10.1)
    ind.flex15 <- c(ind.flex15, ind.flex15.1)
    ind.flex20 <- c(ind.flex20, ind.flex20.1)
    ind.flex25 <- c(ind.flex25, ind.flex25.1)
    ind.flex30 <- c(ind.flex30, ind.flex30.1)
    ind.flex35 <- c(ind.flex35, ind.flex35.1)
    ind.ward <- c(ind.ward, ind.ward.1)
    
    clu.upgma <- c(clu.upgma, clu.upgma.1)
    clu.flex05 <- c(clu.flex05, clu.flex05.1)
    clu.flex10 <- c(clu.flex10, clu.flex10.1)
    clu.flex15 <- c(clu.flex15, clu.flex15.1)
    clu.flex20 <- c(clu.flex20, clu.flex20.1)
    clu.flex25 <- c(clu.flex25, clu.flex25.1)
    clu.flex30 <- c(clu.flex30, clu.flex30.1)
    clu.flex35 <- c(clu.flex35, clu.flex35.1)
    clu.ward <- c(clu.ward, clu.ward.1)
    
    clind.upgma <- c(clind.upgma, clind.upgma.1)
    clind.flex05 <- c(clind.flex05, clind.flex05.1)
    clind.flex10 <- c(clind.flex10, clind.flex10.1)
    clind.flex15 <- c(clind.flex15, clind.flex15.1)
    clind.flex20 <- c(clind.flex20, clind.flex20.1)
    clind.flex25 <- c(clind.flex25, clind.flex25.1)
    clind.flex30 <- c(clind.flex30, clind.flex30.1)
    clind.flex35 <- c(clind.flex35, clind.flex35.1)
    clind.ward <- c(clind.ward, clind.ward.1)
    
    sil.upgma <- c(sil.upgma, sil.upgma.1)
    sil.flex05 <- c(sil.flex05, sil.flex05.1)
    sil.flex10 <- c(sil.flex10, sil.flex10.1)
    sil.flex15 <- c(sil.flex15, sil.flex15.1)
    sil.flex20 <- c(sil.flex20, sil.flex20.1)
    sil.flex25 <- c(sil.flex25, sil.flex25.1)
    sil.flex30 <- c(sil.flex30, sil.flex30.1)
    sil.flex35 <- c(sil.flex35, sil.flex35.1)
    sil.ward <- c(sil.ward, sil.ward.1)
    
    sil.jac <- c(sil.jac, sil.jac.1)
    sil.sim <- c(sil.sim, sil.sim.1)
    sil.diana <- c(sil.diana, sil.diana.1)
    sil.kmeans <- c(sil.kmeans, sil.kmeans.1)
    sil.kulc <- c(sil.kulc, sil.kulc.1)
    sil.kward <- c(sil.kward, sil.kward.1)
    
    #dynamic
    strd0 <- stride(k, as.hclust(tbraydynam))
    maxCoreScatter <- dk[dk$nclust %in% k,]$maxcor 
    minGap <- (1 - maxCoreScatter) * 3/4
    dyngroups <- 
      cutreeDynamic(tbraydynam, minClusterSize=1, method = 'hybrid', distM=as.matrix(distbray),deepSplit=1, maxCoreScatter=maxCoreScatter, minGap=minGap, maxAbsCoreScatter=NULL, minAbsGap=NULL)
    strd0.clustering <- strd0$clustering
    strd0.clustering[,1] <- dyngroups
    strd0$clustering <- strd0.clustering
    ind.dynamic.0 <- optimclass(plotdata1, strd0)
    ind.dynamic.1 <- ind.dynamic.0$sig.spc
    clu.dynamic.1 <- ind.dynamic.0$sig.clust
    clind.dynamic.1 <- ind.dynamic.1*clu.dynamic.1/k
    sil.dynamic.1 <- (dyngroups %>% silhouette(distbray))[,3] %>% mean
    
    ind.dynamic <- c(ind.dynamic, ind.dynamic.1)
    clu.dynamic <- c(clu.dynamic, clu.dynamic.1)
    clind.dynamic <- c(clind.dynamic, clind.dynamic.1)
    sil.dynamic <- c(sil.dynamic, sil.dynamic.1)
    

    
  }
  Sys.time() - timeA  
  
  sil.table <- as.data.frame(cbind(klevel,sil.upgma,sil.flex05,sil.flex10,sil.flex15,sil.flex20,sil.flex25,sil.flex30,sil.flex35,sil.ward, sil.dynamic,
                                   sil.jac,  sil.sim, sil.diana, sil.kmeans, sil.kulc, sil.kward))
  sil.table <- sil.table[-1,]
  lis.table <- sil.table[,-1] %>% t() %>% as.data.frame()
  lis.table$s2to8 <- apply(lis.table[,2:7], MARGIN=1, FUN = 'mean')
  lis.table <- lis.table %>% mutate(s2to8 = apply(lis.table[,1:7], MARGIN=1, FUN = 'mean'))
  sil.table <- sil.table %>% mutate(mean = apply(sil.table[,2:16], MARGIN=1, FUN = 'mean'))
  
  ind.table <- as.data.frame(cbind(klevel,ind.upgma,ind.flex05,ind.flex10,ind.flex15,ind.flex20,ind.flex25,ind.flex30,ind.flex35,ind.ward, ind.dynamic))
  ind.table <- ind.table[-1,]
  dni.table <- ind.table[,-1] %>% t() %>% as.data.frame()
  dni.table <- dni.table %>% mutate(s2to8 = apply(dni.table[,1:7], MARGIN=1, FUN = 'mean'))
  
  clu.table <- as.data.frame(cbind(klevel,clu.upgma,clu.flex05,clu.flex10,clu.flex15,clu.flex20,clu.flex25,clu.flex30,clu.flex35,clu.ward, clu.dynamic))
  clu.table <- clu.table[-1,]
  ulc.table <- clu.table[,-1] %>% t() %>% as.data.frame()
  ulc.table <- ulc.table %>% mutate(s2to8 = apply(ulc.table[,1:7], MARGIN=1, FUN = 'mean'))
  
  clind.table <- as.data.frame(cbind(klevel,clind.upgma,clind.flex05,clind.flex10,clind.flex15,clind.flex20,clind.flex25,clind.flex30,clind.flex35,clind.ward, clind.dynamic))
  clind.table <- clind.table[-1,]
  dnilc.table <- clind.table[,-1] %>% t() %>% as.data.frame()
  dnilc.table <- dnilc.table %>% mutate(s2to8 = apply(dnilc.table[,1:7], MARGIN=1, FUN = 'mean'))
  clind.table <- clind.table %>% mutate(mean = apply(clind.table[,2:10], MARGIN=1, FUN = 'mean'))
  
  
  listofoutput <- list(ind.table, dni.table, clu.table, ulc.table, clind.table, dnilc.table, sil.table, lis.table)
  return(listofoutput)}

indanalysisdynamicward <- function(plotdata,dk,dkward,beta){
  plotdata.total <- apply(plotdata, MARGIN = 2, FUN = 'sum')
  plotdata.total <- as.data.frame(cbind(name=names(plotdata.total),total=plotdata.total))
  removetaxon <- plotdata.total[plotdata.total$total %in% 0,]$name
  plotdata1 <- plotdata[,!colnames(plotdata) %in% removetaxon]
  
  distbray <- vegdist(plotdata1, method='bray', binary=FALSE, na.rm=T)# dbray <- as.data.frame(as.matrix(distbray))
  distjac <- vegdist(plotdata1, method='jaccard', binary=FALSE, na.rm=T)
  distsim <- as.dist(simil(plotdata1,method='Simpson'))
  distkulc <- vegdist(plotdata1, method='kulczynski', binary=FALSE, na.rm=T)
  tbrayagnes <- distbray %>% agnes(method = 'average')
  tbrayflex05 <- distbray %>% flexbeta(beta= -0.05)
  tbrayflex10 <- distbray %>% flexbeta(beta= -0.10)
  tbrayflex15 <- distbray %>% flexbeta(beta= -0.15)
  tbrayflex20 <- distbray %>% flexbeta(beta= -0.20)
  tbrayflex25 <- distbray %>% flexbeta(beta= -0.25)
  tbrayflex30 <- distbray %>% flexbeta(beta= -0.30)
  tbrayflex35 <- distbray %>% flexbeta(beta= -0.35)
  tbrayward <- distbray %>% agnes(method = 'ward')
  tbraydiana <- distbray %>% diana 
  tsimpagnes <- distsim %>% agnes(method = 'average') 
  tjacagnes <- distjac %>% agnes(method = 'average') 
  tkulcagnes <- distkulc %>% agnes(method = 'average')
  tkulcward <- distkulc %>% agnes(method = 'ward')
  #dynamic
  tbraydynam <- distbray %>% flexbeta(beta= beta)
  # dk <- dynamicK(tbraydynam, distbray)
  # dkward <- dynamicK(as.hclust(tbrayward), distbray)
  
  
  k <- 2
  klevel <- 0
  
  sil.upgma <- 0
  sil.flex05 <- 0
  sil.flex10 <- 0
  sil.flex15 <- 0
  sil.flex20 <- 0
  sil.flex25 <- 0
  sil.flex30 <- 0
  sil.flex35 <- 0
  sil.ward <- 0
  sil.jac <- 0
  sil.sim <- 0
  sil.diana <- 0
  sil.kmeans <- 0
  sil.kulc <- 0
  sil.kward <- 0
  
  
  ind.upgma <- 0
  ind.flex05 <- 0
  ind.flex10 <- 0
  ind.flex15 <- 0
  ind.flex20 <- 0
  ind.flex25 <- 0
  ind.flex30 <- 0
  ind.flex35 <- 0
  ind.ward <- 0
  
  clu.upgma <- 0
  clu.flex05 <- 0
  clu.flex10 <- 0
  clu.flex15 <- 0
  clu.flex20 <- 0
  clu.flex25 <- 0
  clu.flex30 <- 0
  clu.flex35 <- 0
  clu.ward <- 0
  
  clind.upgma <- 0
  clind.flex05 <- 0
  clind.flex10 <- 0
  clind.flex15 <- 0
  clind.flex20 <- 0
  clind.flex25 <- 0
  clind.flex30 <- 0
  clind.flex35 <- 0
  clind.ward <- 0
  
  ind.dynamic <- 0
  clu.dynamic <- 0
  clind.dynamic <- 0
  sil.dynamic <- 0
  ind.dynamicward <- 0
  clu.dynamicward <- 0
  clind.dynamicward <- 0
  sil.dynamicward <- 0
  
  
  timeA = Sys.time()
  for (k in 2:10){
    
    ind.upgma.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayagnes)))
    ind.flex05.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex05)))
    ind.flex10.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex10)))
    ind.flex15.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex15)))
    ind.flex20.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex20)))
    ind.flex25.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex25)))
    ind.flex30.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex30)))
    ind.flex35.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex35)))
    ind.ward.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayward)))
    
    ind.upgma.1 <- ind.upgma.0$sig.spc
    ind.flex05.1 <- ind.flex05.0$sig.spc
    ind.flex10.1 <- ind.flex10.0$sig.spc
    ind.flex15.1 <- ind.flex15.0$sig.spc
    ind.flex20.1 <- ind.flex20.0$sig.spc
    ind.flex25.1 <- ind.flex25.0$sig.spc
    ind.flex30.1 <- ind.flex30.0$sig.spc
    ind.flex35.1 <- ind.flex35.0$sig.spc
    ind.ward.1 <- ind.ward.0$sig.spc
    
    clu.upgma.1 <- ind.upgma.0$sig.clust
    clu.flex05.1 <- ind.flex05.0$sig.clust
    clu.flex10.1 <- ind.flex10.0$sig.clust
    clu.flex15.1 <- ind.flex15.0$sig.clust
    clu.flex20.1 <- ind.flex20.0$sig.clust
    clu.flex25.1 <- ind.flex25.0$sig.clust
    clu.flex30.1 <- ind.flex30.0$sig.clust
    clu.flex35.1 <- ind.flex35.0$sig.clust
    clu.ward.1 <- ind.ward.0$sig.clust
    
    clind.upgma.1 <- ind.upgma.1*clu.upgma.1/k
    clind.flex05.1 <- ind.flex05.1*clu.flex05.1/k
    clind.flex10.1 <- ind.flex10.1*clu.flex10.1/k
    clind.flex15.1 <- ind.flex15.1*clu.flex15.1/k
    clind.flex20.1 <- ind.flex20.1*clu.flex20.1/k
    clind.flex25.1 <- ind.flex25.1*clu.flex25.1/k
    clind.flex30.1 <- ind.flex30.1*clu.flex30.1/k
    clind.flex35.1 <- ind.flex35.1*clu.flex35.1/k
    clind.ward.1 <- ind.ward.1*clu.ward.1/k
    
    sil.upgma.1 <- (tbrayagnes %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.flex05.1 <- (tbrayflex05 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.flex10.1 <- (tbrayflex10 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.flex15.1 <- (tbrayflex15 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.flex20.1 <- (tbrayflex20 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.flex25.1 <- (tbrayflex25 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.flex30.1 <- (tbrayflex30 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.flex35.1 <- (tbrayflex35 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.ward.1 <- (tbrayward %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    
    sil.jac.1 <- (tjacagnes %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.sim.1 <- (tsimpagnes %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.diana.1 <- (tbraydiana %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.kmeans.1 <- (kmeans(distbray, centers = k)$cluster %>% silhouette(distbray))[,3] %>% mean
    sil.kulc.1 <- (tkulcagnes %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.kward.1 <- (tkulcward %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    
    klevel <- c(klevel, k)
    ind.upgma <- c(ind.upgma, ind.upgma.1)
    ind.flex05 <- c(ind.flex05, ind.flex05.1)
    ind.flex10 <- c(ind.flex10, ind.flex10.1)
    ind.flex15 <- c(ind.flex15, ind.flex15.1)
    ind.flex20 <- c(ind.flex20, ind.flex20.1)
    ind.flex25 <- c(ind.flex25, ind.flex25.1)
    ind.flex30 <- c(ind.flex30, ind.flex30.1)
    ind.flex35 <- c(ind.flex35, ind.flex35.1)
    ind.ward <- c(ind.ward, ind.ward.1)
    
    clu.upgma <- c(clu.upgma, clu.upgma.1)
    clu.flex05 <- c(clu.flex05, clu.flex05.1)
    clu.flex10 <- c(clu.flex10, clu.flex10.1)
    clu.flex15 <- c(clu.flex15, clu.flex15.1)
    clu.flex20 <- c(clu.flex20, clu.flex20.1)
    clu.flex25 <- c(clu.flex25, clu.flex25.1)
    clu.flex30 <- c(clu.flex30, clu.flex30.1)
    clu.flex35 <- c(clu.flex35, clu.flex35.1)
    clu.ward <- c(clu.ward, clu.ward.1)
    
    clind.upgma <- c(clind.upgma, clind.upgma.1)
    clind.flex05 <- c(clind.flex05, clind.flex05.1)
    clind.flex10 <- c(clind.flex10, clind.flex10.1)
    clind.flex15 <- c(clind.flex15, clind.flex15.1)
    clind.flex20 <- c(clind.flex20, clind.flex20.1)
    clind.flex25 <- c(clind.flex25, clind.flex25.1)
    clind.flex30 <- c(clind.flex30, clind.flex30.1)
    clind.flex35 <- c(clind.flex35, clind.flex35.1)
    clind.ward <- c(clind.ward, clind.ward.1)
    
    sil.upgma <- c(sil.upgma, sil.upgma.1)
    sil.flex05 <- c(sil.flex05, sil.flex05.1)
    sil.flex10 <- c(sil.flex10, sil.flex10.1)
    sil.flex15 <- c(sil.flex15, sil.flex15.1)
    sil.flex20 <- c(sil.flex20, sil.flex20.1)
    sil.flex25 <- c(sil.flex25, sil.flex25.1)
    sil.flex30 <- c(sil.flex30, sil.flex30.1)
    sil.flex35 <- c(sil.flex35, sil.flex35.1)
    sil.ward <- c(sil.ward, sil.ward.1)
    
    sil.jac <- c(sil.jac, sil.jac.1)
    sil.sim <- c(sil.sim, sil.sim.1)
    sil.diana <- c(sil.diana, sil.diana.1)
    sil.kmeans <- c(sil.kmeans, sil.kmeans.1)
    sil.kulc <- c(sil.kulc, sil.kulc.1)
    sil.kward <- c(sil.kward, sil.kward.1)
    
    #dynamic
    strd0 <- stride(k, as.hclust(tbraydynam))
    maxCoreScatter <- dk[dk$nclust %in% k,]$maxcor 
    minGap <- (1 - maxCoreScatter) * 3/4
    dyngroups <- 
      cutreeDynamic(tbraydynam, minClusterSize=1, method = 'hybrid', distM=as.matrix(distbray),deepSplit=1, maxCoreScatter=maxCoreScatter, minGap=minGap, maxAbsCoreScatter=NULL, minAbsGap=NULL)
    strd0.clustering <- strd0$clustering
    strd0.clustering[,1] <- dyngroups
    strd0$clustering <- strd0.clustering
    ind.dynamic.0 <- optimclass(plotdata1, strd0)
    ind.dynamic.1 <- ind.dynamic.0$sig.spc
    clu.dynamic.1 <- ind.dynamic.0$sig.clust
    clind.dynamic.1 <- ind.dynamic.1*clu.dynamic.1/k
    sil.dynamic.1 <- (dyngroups %>% silhouette(distbray))[,3] %>% mean
    
    ind.dynamic <- c(ind.dynamic, ind.dynamic.1)
    clu.dynamic <- c(clu.dynamic, clu.dynamic.1)
    clind.dynamic <- c(clind.dynamic, clind.dynamic.1)
    sil.dynamic <- c(sil.dynamic, sil.dynamic.1)
    
    #dynamic ward
    strdw0 <- stride(k, as.hclust(tbrayward))
    maxCoreScatter <- dkward[dkward$nclust %in% k,]$maxcor 
    minGap <- (1 - maxCoreScatter) * 3/4
    dyngroupw <- 
      cutreeDynamic(as.hclust(tbrayward), minClusterSize=1, method = 'hybrid', distM=as.matrix(distbray),deepSplit=1, maxCoreScatter=maxCoreScatter, minGap=minGap, maxAbsCoreScatter=NULL, minAbsGap=NULL)
    strdw0.clustering <- strdw0$clustering
    strdw0.clustering[,1] <- dyngroupw
    strdw0$clustering <- strdw0.clustering
    ind.dynamicward.0 <- optimclass(plotdata1, strdw0)
    ind.dynamicward.1 <- ind.dynamicward.0$sig.spc
    clu.dynamicward.1 <- ind.dynamicward.0$sig.clust
    clind.dynamicward.1 <- ind.dynamicward.1*clu.dynamicward.1/k
    sil.dynamicward.1 <- (dyngroupw %>% silhouette(distbray))[,3] %>% mean
    
    ind.dynamicward <- c(ind.dynamicward, ind.dynamicward.1)
    clu.dynamicward <- c(clu.dynamicward, clu.dynamicward.1)
    clind.dynamicward <- c(clind.dynamicward, clind.dynamicward.1)
    sil.dynamicward <- c(sil.dynamicward, sil.dynamicward.1)
    
    
  }
  Sys.time() - timeA  
  
  sil.table <- as.data.frame(cbind(klevel,sil.upgma,sil.flex05,sil.flex10,sil.flex15,sil.flex20,sil.flex25,sil.flex30,sil.flex35,sil.ward, sil.dynamic,sil.dynamicward,
                                   sil.jac,  sil.sim, sil.diana, sil.kmeans, sil.kulc, sil.kward))
  sil.table <- sil.table[-1,]
  lis.table <- sil.table[,-1] %>% t() %>% as.data.frame()
  lis.table$s2to8 <- apply(lis.table[,2:7], MARGIN=1, FUN = 'mean')
  lis.table <- lis.table %>% mutate(s2to8 = apply(lis.table[,1:7], MARGIN=1, FUN = 'mean'))
  sil.table <- sil.table %>% mutate(mean = apply(sil.table[,2:16], MARGIN=1, FUN = 'mean'))
  
  ind.table <- as.data.frame(cbind(klevel,ind.upgma,ind.flex05,ind.flex10,ind.flex15,ind.flex20,ind.flex25,ind.flex30,ind.flex35,ind.ward, ind.dynamic, ind.dynamicward))
  ind.table <- ind.table[-1,]
  dni.table <- ind.table[,-1] %>% t() %>% as.data.frame()
  dni.table <- dni.table %>% mutate(s2to8 = apply(dni.table[,1:7], MARGIN=1, FUN = 'mean'))
  
  clu.table <- as.data.frame(cbind(klevel,clu.upgma,clu.flex05,clu.flex10,clu.flex15,clu.flex20,clu.flex25,clu.flex30,clu.flex35,clu.ward, clu.dynamic,clu.dynamicward))
  clu.table <- clu.table[-1,]
  ulc.table <- clu.table[,-1] %>% t() %>% as.data.frame()
  ulc.table <- ulc.table %>% mutate(s2to8 = apply(ulc.table[,1:7], MARGIN=1, FUN = 'mean'))
  
  clind.table <- as.data.frame(cbind(klevel,clind.upgma,clind.flex05,clind.flex10,clind.flex15,clind.flex20,clind.flex25,clind.flex30,clind.flex35,clind.ward, clind.dynamic,clind.dynamicward))
  clind.table <- clind.table[-1,]
  dnilc.table <- clind.table[,-1] %>% t() %>% as.data.frame()
  dnilc.table <- dnilc.table %>% mutate(s2to8 = apply(dnilc.table[,1:7], MARGIN=1, FUN = 'mean'))
  clind.table <- clind.table %>% mutate(mean = apply(clind.table[,2:10], MARGIN=1, FUN = 'mean'))
  
  
  listofoutput <- list(ind.table, dni.table, clu.table, ulc.table, clind.table, dnilc.table, sil.table, lis.table)
  return(listofoutput)}

indgroup <- function(plotdata, groups, qualitative){#qualitative=F
  p.rowmax <- apply(plotdata, MARGIN = 1, FUN=max)
  p.normal <- plotdata/p.rowmax
  if(qualitative){p.normal <- as.data.frame((p.normal>0)*1)}
  
  plotdata1 <- p.normal %>% cbind(groups=groups)
  # plotdata1 <- t(plotdata1)
  plotdata.group <- plotdata1 %>% group_by(groups) %>% summarise(across(colnames(plotdata)[1:ncol(plotdata1)-1], mean))
  # plotdata.group.total <- (apply(plotdata.group, MARGIN = 2, FUN = sum))
  plotdata.group.total <- plotdata.group %>% mutate(across(colnames(plotdata.group)[2:ncol(plotdata.group)], sum))
  
  plotdata.affinity <- plotdata.group/plotdata.group.total
  plotdata.indicator <- (plotdata.affinity*plotdata.group)^0.5
  plotdata.indicator$groups <- plotdata.group$groups
  plotdata.indicator.total <- as.data.frame(cbind(groups = plotdata.group$groups,
                                                  total = apply(plotdata.indicator[,2:ncol(plotdata.indicator)], MARGIN = 1, FUN = sum),
                                                  maxval = apply(plotdata.indicator[,2:ncol(plotdata.indicator)], MARGIN = 1, FUN = max)
  ))
  return(plotdata.indicator.total)
}

indgroup2 <- function(plotdata, groups, qualitative){#qualitative=F
  p.rowmax <- apply(plotdata, MARGIN = 1, FUN=max)
  p.normal <- plotdata/p.rowmax
  if(qualitative){p.normal <- as.data.frame((p.normal>0)*1)}
  
  plotdata1 <- p.normal %>% cbind(groups=groups)
  # plotdata1 <- t(plotdata1)
  plotdata.group <- plotdata1 %>% group_by(groups) %>% summarise(across(colnames(plotdata)[1:ncol(plotdata1)-1], mean))
  # plotdata.group.total <- (apply(plotdata.group, MARGIN = 2, FUN = sum))
  plotdata.group.total <- plotdata.group %>% mutate(across(colnames(plotdata.group)[2:ncol(plotdata.group)], sum))
  
  plotdata.affinity <- plotdata.group/plotdata.group.total
  plotdata.indicator <- (plotdata.affinity*plotdata.group)^0.5
  plotdata.indicator$groups <- plotdata.group$groups
  iplot <- t(plotdata.indicator[,-1]) %>% as.data.frame()
  gplot <- t(plotdata.group[,-1]) %>% as.data.frame()
  aplot <- t(plotdata.affinity[,-1]) %>% as.data.frame()
  indlist <- list(iplot, gplot, aplot)
  return(indlist)
}

clustvar <- function(d, groups){
  df <- as.matrix(d)*-1+1
  grps <- sort(unique(groups))
  for(i in 1:length(grps)){#i=3
    cluster = grps[i]
    grpclust <- df[which(groups %in% i),which(groups %in% i)]
    dmean <- mean(grpclust)
    if(length(grpclust)>1){
      dmin <- min(apply(grpclust, MARGIN=2, FUN=mean))
    }else{dmax = max(grpclust)}
       dclust0 <- as.data.frame(cbind(cluster,dmean, dmin))
    if(i==1){dclust=dclust0}else{dclust=rbind(dclust,dclust0)}
  }
  return(dclust)
}

indanalysis2 <- function(plotdata){
  plotdata.total <- apply(plotdata, MARGIN = 2, FUN = 'sum')
  plotdata.total <- as.data.frame(cbind(name=names(plotdata.total),total=plotdata.total))
  removetaxon <- plotdata.total[plotdata.total$total %in% 0,]$name
  plotdata1 <- plotdata[,!colnames(plotdata) %in% removetaxon]
   
  distbray <- vegdist(plotdata1, method='bray', binary=FALSE, na.rm=T)# dbray <- as.data.frame(as.matrix(distbray))
  distjac <- vegdist(plotdata1, method='jaccard', binary=FALSE, na.rm=T)
  distsim <- as.dist(simil(plotdata1,method='Simpson'))
  distkulc <- vegdist(plotdata1, method='kulczynski', binary=FALSE, na.rm=T)
  tbrayagnes <- distbray %>% agnes(method = 'average')
  tbrayflex05 <- distbray %>% flexbeta(beta= -0.05)
  tbrayflex10 <- distbray %>% flexbeta(beta= -0.10)
  tbrayflex15 <- distbray %>% flexbeta(beta= -0.15)
  tbrayflex20 <- distbray %>% flexbeta(beta= -0.20)
  tbrayflex25 <- distbray %>% flexbeta(beta= -0.25)
  tbrayflex30 <- distbray %>% flexbeta(beta= -0.30)
  tbrayflex35 <- distbray %>% flexbeta(beta= -0.35)
  tbrayward <- distbray %>% agnes(method = 'ward')
  tbraydiana <- distbray %>% diana 
  tsimpagnes <- distsim %>% agnes(method = 'average') 
  tjacagnes <- distjac %>% agnes(method = 'average') 
  tkulcagnes <- distkulc %>% agnes(method = 'average')
  tkulcward <- distkulc %>% agnes(method = 'ward')
  
  klevel <- 0
  
  ind.upgma <- 0
  ind.flex05 <- 0
  ind.flex10 <- 0
  ind.flex15 <- 0
  ind.flex20 <- 0
  ind.flex25 <- 0
  ind.flex30 <- 0
  ind.flex35 <- 0
  ind.ward <- 0
  
  ind.jac <- 0
  ind.sim <- 0
  ind.diana <- 0
  ind.kmeans <- 0
  ind.kulc <- 0
  ind.kward <- 0
  
  
  weak.upgma <- 0
  weak.flex05 <- 0
  weak.flex10 <- 0
  weak.flex15 <- 0
  weak.flex20 <- 0
  weak.flex25 <- 0
  weak.flex30 <- 0
  weak.flex35 <- 0
  weak.ward <- 0
  
  weak.jac <- 0
  weak.sim <- 0
  weak.diana <- 0
  weak.kmeans <- 0
  weak.kulc <- 0
  weak.kward <- 0
  
  for (k in 2:min(nrow(plotdata1)-1,15)){
    
    ind.upgma.0 <- indgroup(plotdata1, cutree(tbrayagnes, k = k), F)
    ind.flex05.0 <- indgroup(plotdata1, cutree(tbrayflex05, k = k), F)
    ind.flex10.0 <- indgroup(plotdata1, cutree(tbrayflex10, k = k), F)
    ind.flex15.0 <- indgroup(plotdata1, cutree(tbrayflex15, k = k), F)
    ind.flex20.0 <- indgroup(plotdata1, cutree(tbrayflex20, k = k), F)
    ind.flex25.0 <- indgroup(plotdata1, cutree(tbrayflex25, k = k), F)
    ind.flex30.0 <- indgroup(plotdata1, cutree(tbrayflex30, k = k), F)
    ind.flex35.0 <- indgroup(plotdata1, cutree(tbrayflex35, k = k), F)
    ind.ward.0 <- indgroup(plotdata1, cutree(tbrayward, k = k), F)
    
    ind.jac.0 <- indgroup(plotdata1, cutree(tjacagnes, k = k), F)
    ind.sim.0 <- indgroup(plotdata1, cutree(tsimpagnes, k = k), F)
    ind.diana.0 <- indgroup(plotdata1, cutree(tbraydiana, k = k), F)
    ind.kmeans.0 <- indgroup(plotdata1, kmeans(distbray, centers = k)$cluster, F)
    ind.kulc.0 <- indgroup(plotdata1, cutree(tkulcagnes, k = k), F)
    ind.kward.0 <- indgroup(plotdata1, cutree(tkulcward, k = k), F)
    
    ind.upgma.1 <- ind.upgma.0[,'total'] %>% mean()
    ind.flex05.1 <- ind.flex05.0[,'total'] %>% mean()
    ind.flex10.1 <- ind.flex10.0[,'total'] %>% mean()
    ind.flex15.1 <- ind.flex15.0[,'total'] %>% mean()
    ind.flex20.1 <- ind.flex20.0[,'total'] %>% mean()
    ind.flex25.1 <- ind.flex25.0[,'total'] %>% mean()
    ind.flex30.1 <- ind.flex30.0[,'total'] %>% mean()
    ind.flex35.1 <- ind.flex35.0[,'total'] %>% mean()
    ind.ward.1 <- ind.ward.0[,'total'] %>% mean()
    
    ind.jac.1 <-  ind.jac.0[,'total'] %>% mean()
    ind.sim.1 <- ind.sim.0[,'total'] %>% mean()
    ind.diana.1 <- ind.diana.0[,'total'] %>% mean()
    ind.kmeans.1 <- ind.kmeans.0[,'total'] %>% mean()
    ind.kulc.1 <- ind.kulc.0[,'total'] %>% mean()
    ind.kward.1 <- ind.kward.0[,'total'] %>% mean()
    
    weak.upgma.1 <- ind.upgma.0[,'maxval'] %>% min()
    weak.flex05.1 <- ind.flex05.0[,'maxval'] %>% min()
    weak.flex10.1 <- ind.flex10.0[,'maxval'] %>% min()
    weak.flex15.1 <- ind.flex15.0[,'maxval'] %>% min()
    weak.flex20.1 <- ind.flex20.0[,'maxval'] %>% min()
    weak.flex25.1 <- ind.flex25.0[,'maxval'] %>% min()
    weak.flex30.1 <- ind.flex30.0[,'maxval'] %>% min()
    weak.flex35.1 <- ind.flex35.0[,'maxval'] %>% min()
    weak.ward.1 <- ind.ward.0[,'maxval'] %>% min()
    
    weak.jac.1 <-  ind.jac.0[,'maxval'] %>% min()
    weak.sim.1 <- ind.sim.0[,'maxval'] %>% min()
    weak.diana.1 <- ind.diana.0[,'maxval'] %>% min()
    weak.kmeans.1 <- ind.kmeans.0[,'maxval'] %>% min()
    weak.kulc.1 <- ind.kulc.0[,'maxval'] %>% min()
    weak.kward.1 <- ind.kward.0[,'maxval'] %>% min()
    
    
    klevel <- c(klevel, k)
    
    ind.upgma <- c(ind.upgma, ind.upgma.1)
    ind.flex05 <- c(ind.flex05, ind.flex05.1)
    ind.flex10 <- c(ind.flex10, ind.flex10.1)
    ind.flex15 <- c(ind.flex15, ind.flex15.1)
    ind.flex20 <- c(ind.flex20, ind.flex20.1)
    ind.flex25 <- c(ind.flex25, ind.flex25.1)
    ind.flex30 <- c(ind.flex30, ind.flex30.1)
    ind.flex35 <- c(ind.flex35, ind.flex35.1)
    ind.ward <- c(ind.ward, ind.ward.1)
    
    ind.jac <- c(ind.jac, ind.jac.1)
    ind.sim <- c(ind.sim, ind.sim.1)
    ind.diana <- c(ind.diana, ind.diana.1)
    ind.kmeans <- c(ind.kmeans, ind.kmeans.1)
    ind.kulc <- c(ind.kulc, ind.kulc.1)
    ind.kward <- c(ind.kward, ind.kward.1)
    
    weak.upgma <- c(weak.upgma, weak.upgma.1)
    weak.flex05 <- c(weak.flex05, weak.flex05.1)
    weak.flex10 <- c(weak.flex10, weak.flex10.1)
    weak.flex15 <- c(weak.flex15, weak.flex15.1)
    weak.flex20 <- c(weak.flex20, weak.flex20.1)
    weak.flex25 <- c(weak.flex25, weak.flex25.1)
    weak.flex30 <- c(weak.flex30, weak.flex30.1)
    weak.flex35 <- c(weak.flex35, weak.flex35.1)
    weak.ward <- c(weak.ward, weak.ward.1)
    
    weak.jac <- c(weak.jac, weak.jac.1)
    weak.sim <- c(weak.sim, weak.sim.1)
    weak.diana <- c(weak.diana, weak.diana.1)
    weak.kmeans <- c(weak.kmeans, weak.kmeans.1)
    weak.kulc <- c(weak.kulc, weak.kulc.1)
    weak.kward <- c(weak.kward, weak.kward.1)
    
  }    
  
  ind.table <- as.data.frame(cbind(klevel,ind.upgma,ind.flex05,ind.flex10,ind.flex15,ind.flex20,ind.flex25,ind.flex30,ind.flex35,ind.ward,ind.jac,ind.sim,ind.diana,ind.kmeans,ind.kulc,ind.kward))
  ind.table <- ind.table[-1,]
  
  weak.table <- as.data.frame(cbind(klevel,weak.upgma,weak.flex05,weak.flex10,weak.flex15,weak.flex20,weak.flex25,weak.flex30,weak.flex35,weak.ward,weak.jac,weak.sim,weak.diana,weak.kmeans,weak.kulc,weak.kward))
  weak.table <- weak.table[-1,]
  
  dni.table <- ind.table[,-1] %>% t() %>% as.data.frame()
  dni.table <- dni.table %>% mutate(s2to6 = apply(dni.table[,1:5], MARGIN=1, FUN = 'mean'))
  
  kaew.table <- weak.table[,-1] %>% t() %>% as.data.frame()
  kaew.table <- kaew.table %>% mutate(s2to6 = apply(kaew.table[,1:5], MARGIN=1, FUN = 'mean'))
  
  ind.table <- ((ind.table[,c(2:ncol(ind.table))]) -
                          (apply(ind.table[,c(2:ncol(ind.table))], MARGIN=2, FUN = 'mean')+0.001))/
    (apply(ind.table[,c(2:ncol(ind.table))], MARGIN=2, FUN = 'sd')+0.001)
  
  
  listofoutput <- list(ind.table, dni.table, weak.table, kaew.table)
  return(listofoutput)}




adjustcover <- function(cover.est, cover.total){
  #this function adjusts ocular estimates of individual species to be coherent with the ocular estimate of total cover for a stratum.
  #cover.est = c(10,50,80) #estimated cover for individual taxa in stratum
  #cover.total = 80 #estimated cover for whole stratum
  cover.est = cover.est/100
  cover.total = cover.total/100
  cover.agg1 = 1-10^(sum(log10(1-cover.est)))
  cover.fac1 = (cover.total/cover.agg1)^1.5 #first pass makes a linear adjustment so that relative cover is consistent with field estimate.
  cover.adj1 = cover.est*cover.fac1
  cover.agg2 = 1-10^(sum(log10(1-cover.adj1)))
  cover.fac2 = (log10(1-cover.total)/log10(1-cover.agg2))#second pass fine tunes adjusted cover so that aggregate cover matches ocular total cover.
  cover.agg3 = 1-10^(cover.fac2*log10(1-cover.adj1))
  cover.adj  = cover.agg3*100
  #cover.agg.test <- 100*(1-10^(sum(log10(1-cover.adj/100)))) #produces aggregate of stratum after adjustment; should be the same as cover.total
  return(cover.adj)
}

cleanplotdata <- function(plotdata){
  plotdata.total <- apply(plotdata, MARGIN = 2, FUN = 'sum')
  plotdata.total <- as.data.frame(cbind(name=names(plotdata.total),total=plotdata.total))
  removetaxon <- plotdata.total[plotdata.total$total %in% 0,]$name
  plotdata1 <- plotdata[,!colnames(plotdata) %in% removetaxon]
  p.rowmax <- apply(plotdata1, MARGIN = 1, FUN=max)
  p.normal <- plotdata1/p.rowmax
  return(p.normal)
}