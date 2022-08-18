


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

makeplotdynamic <- function(a,d,t,k){
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
                          deepSplit=deepSplit, maxCoreScatter=maxCoreScatter, minGap=minGap, maxAbsCoreScatter=NULL, minAbsGap=NULL)
  
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
