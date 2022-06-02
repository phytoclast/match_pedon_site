setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(soilDB)
library(aqp)#load before dplr
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


#----
Com.Sp.mean <- readRDS('data/Com.Sp.mean.RDS')
plotdata <- makecommunitydataset(Com.Sp.mean, row = 'soilplot', column = 'Species', value = 'sqrttotal', drop = TRUE)


distbray <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
distjac <- vegdist(plotdata, method='jaccard', binary=FALSE, na.rm=T)

#validity using silhouette index
sil.table <-'x'
#----
silanalysis <- function(input){
  distbray <- vegdist(input, method='bray', binary=FALSE, na.rm=T)
  distjac <- vegdist(input, method='jaccard', binary=FALSE, na.rm=T)
  distsim <- as.dist(simil(input,method='Simpson'))
  
  maxcluster <- min(20, nrow(input)-1)
  k <- 2
  klevel <- 0
  sil.bray <- 0
  sil.jac <- 0
  sil.sim <- 0
  sil.ward <- 0
  sil.diana <- 0
  sil.kmeans <- 0
  sil.single <- 0
  sil.complete <- 0
  sil.wardeuc <- 0
  sil.kmeanseuc <- 0

  for (k in 2:20){
    sil.bray1 <- (distbray %>% agnes(method = 'average') %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.jac1 <- (distjac %>% agnes(method = 'average') %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.sim1 <- (distsim %>% agnes(method = 'average') %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.ward1 <- (distbray %>% agnes(method = 'ward') %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.diana1 <- (distbray %>% diana %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.kmeans1 <- (kmeans(distbray, centers = k)$cluster %>% silhouette(distbray))[,3] %>% mean
    sil.single1 <- (distbray %>% agnes(method = 'single') %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.complete1 <- (distbray %>% agnes(method = 'complete') %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.wardeuc1 <- (plotdata %>% agnes(method = 'ward') %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
    sil.kmeanseuc1 <- (kmeans(plotdata, centers = k)$cluster %>% silhouette(distbray))[,3] %>% mean

    
    klevel <- c(klevel, k)
    sil.bray <- c(sil.bray, sil.bray1)
    sil.jac <- c(sil.jac, sil.jac1)
    sil.sim <- c(sil.sim, sil.sim1)
    sil.ward <- c(sil.ward, sil.ward1)
    sil.diana <- c(sil.diana, sil.diana1)
    sil.kmeans <- c(sil.kmeans, sil.kmeans1)
    sil.single <- c(sil.single, sil.single1)
    sil.complete <- c(sil.complete, sil.complete1)
    sil.wardeuc <- c(sil.wardeuc, sil.wardeuc1)
    sil.kmeanseuc <- c(sil.kmeanseuc, sil.kmeanseuc1)
  }
  sil.table <- as.data.frame(cbind(klevel,sil.bray,sil.jac,sil.sim,sil.ward,sil.diana,sil.kmeans,sil.single,sil.complete))
  sil.table <- sil.table[-1,]
  sil.table <<- sil.table
  return(sil.table)}
#----

silanalysis(plotdata)
#analysis method
makeplot <- function(amethod,jacdist,jactree, soilgroup,k){
  filename <- paste0('output/Soils_',soilgroup,"_",amethod,'.png')
  
  #make cuts and reformat dendrogram
  ngroups=k
  groups <- cutree(jactree, k = ngroups)
  
  soilplot <- names(groups)
  clust <- unname(groups)
  groupdf <- as.data.frame(cbind(soilplot, clust))
  groupdf$clust <- (as.numeric(as.character(groupdf$clust)))
  maxcluster <- max(groupdf$clust)
  numberzeros <- nrow(groupdf[(groupdf$clust == 0),])
  whichrecords <- which(groupdf$clust == 0)
  if (nrow(groupdf[groupdf$clust == 0,]) != 0){
    for (i in 1:numberzeros){ #assign all zero clusters to unique cluster number.
      groupdf[whichrecords[i],]$clust <- maxcluster+i}}
  
  newlabels <- jactree$order.lab
  newlabels <- as.data.frame(newlabels)
  newlabels$row <- row(newlabels)
  newlabels <- merge(newlabels, groupdf, by.x='newlabels', by.y ='soilplot')
  newlabels$newlabels <- paste(newlabels$clust, newlabels$newlabels)
  newlabels <- newlabels[order(newlabels$row),1]
  newtree <- jactree
  newtree$order.lab <- newlabels
  
  dend1 <- color_branches(as.hclust(newtree), k = ngroups)
  dend1 <- color_labels(dend1, k = ngroups)
  
  #output file
  
  w <- 800
  h <- nrow(jacdist)*12+80
  u <- 12
  png(filename=filename,width = w, height = h, units = "px", pointsize = u)
  
  par(mar = c(2,0,1,13))
  plot(dend1, horiz = TRUE, main=paste('floristic simularity', amethod,' method of', soilgroup, 'soils'), font=1, cex=0.84)
  dev.off()
  
}

amethod <- 'bray-agnes' 
k=16
if (T){
  amethod <- 'bray-agnes' 
  k=8
  jacdist <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
  jactree <- agnes(jacdist, method='average')
  makeplot(amethod,jacdist,jactree,soilgroup,k)
}
if (T){
  amethod <- 'bray-single' 
  jacdist <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
  jactree <- agnes(jacdist, method='single')
  makeplot(amethod,jacdist,jactree,soilgroup,k)
}
if (T){
  amethod <- 'bray-complete' 
  jacdist <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
  jactree <- agnes(jacdist, method='complete')
  makeplot(amethod,jacdist,jactree,soilgroup,k)
}
if (T){
  amethod <- 'bray-diana' 
  jacdist <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
  jactree <- diana(jacdist)
  makeplot(amethod,jacdist,jactree,soilgroup,k)
}
if (T){
  amethod <- 'bray-ward'
  k=8
  jacdist <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
  jactree <- agnes(jacdist, method='ward')
  makeplot(amethod,jacdist,jactree,soilgroup,k)
}
if (T){
  amethod <- 'jaccard-agnes' 
  jacdist <- vegdist(plotdata, method='jaccard', binary=FALSE, na.rm=T)
  jactree <- agnes(jacdist, method='average')
  makeplot(amethod,jacdist,jactree,soilgroup,k)
}
if (T){
  amethod <- 'kulczynski-agnes' 
  jacdist <- vegdist(plotdata, method='kulczynski', binary=FALSE, na.rm=T)
  jactree <- agnes(jacdist, method='average')
  makeplot(amethod,jacdist,jactree,soilgroup,k)
}
if (T){
  amethod <- 'kulczynski-ward' 
  jacdist <- vegdist(plotdata, method='kulczynski', binary=FALSE, na.rm=T)
  jactree <- agnes(jacdist, method='ward')
  makeplot(amethod,jacdist,jactree,soilgroup,k)
}
#----
#group dominant and indicator species
k=8
d <- ((vegdist(plotdata, method='kulczynski', binary=FALSE, na.rm=T)))
t <- agnes(d, method='ward')
groups <- cutree(t, k = k)

soilplot <- names(groups)
clust <- unname(groups)
groupdf <- as.data.frame(cbind(soilplot, clust))
groupdf$clust <- (as.numeric(as.character(groupdf$clust)))
maxcluster <- max(groupdf$clust)
numberzeros <- nrow(groupdf[(groupdf$clust == 0),])
whichrecords <- which(groupdf$clust == 0)
if (nrow(groupdf[groupdf$clust == 0,]) != 0){
  for (i in 1:numberzeros){ #assign all zero clusters to unique cluster number.
    groupdf[whichrecords[i],]$clust <- maxcluster+i}}



Com.Sp.groups <- merge(groupdf,  Com.Sp.mean, by='soilplot', all.x=TRUE, all.y = TRUE)

#----
#average spp by cluster

Com.Sp.groups.sum <- aggregate(Com.Sp.groups[,c('Field', 'Shrub', 'Subcanopy','Tree', 'Total')],
                               by=list(Com.Sp.groups$clust, Com.Sp.groups$Species), FUN='sum')
colnames(Com.Sp.groups.sum) <- c('cluster', 'taxon', 'Field', 'Shrub', 'Subcanopy','Tree', 'Total')
Com.Sp.groups.count <- aggregate(unique(Com.Sp.groups[c('clust', 'soilplot')])$soilplot, 
                                 by=list(unique(Com.Sp.groups[c('clust', 'soilplot')])$clust), FUN='length')
colnames(Com.Sp.groups.count) <- c('cluster', 'count')
Com.Sp.groups.mean <- merge(Com.Sp.groups.sum, Com.Sp.groups.count, by = 'cluster')
Com.Sp.groups.mean[,c('Field', 'Shrub', 'Subcanopy','Tree', 'Total')] <- Com.Sp.groups.mean[,c('Field', 'Shrub', 'Subcanopy','Tree', 'Total')]/Com.Sp.groups.mean$count
rm(Com.Sp.groups.sum, Com.Sp.groups.count)

#frequency spp by cluster
Com.Sp.prefreq <- Com.Sp.groups
Com.Sp.prefreq$Total <- ifelse(Com.Sp.prefreq$Total >0, 1,0)
Com.Sp.freq.sum <- aggregate(Com.Sp.prefreq$Total,
                             by=list(Com.Sp.prefreq$clust, Com.Sp.prefreq$Species), FUN='sum')
colnames(Com.Sp.freq.sum) <- c('cluster', 'taxon', 'freq')
Com.Sp.groups.count <- aggregate(unique(Com.Sp.prefreq[c('clust', 'soilplot')])$soilplot, 
                                 by=list(unique(Com.Sp.prefreq[c('clust', 'soilplot')])$clust), FUN='length')
colnames(Com.Sp.groups.count) <- c('cluster', 'count')
Com.Sp.groups.freq <- merge(Com.Sp.freq.sum, Com.Sp.groups.count, by = 'cluster')
Com.Sp.groups.freq$freq <- Com.Sp.groups.freq$freq/Com.Sp.groups.freq$count*100
Com.Sp.groups.mean <- merge(Com.Sp.groups.mean, Com.Sp.groups.freq[,c('cluster', 'taxon', 'freq')], by = c('cluster', 'taxon'))
Com.Sp.groups.mean$freqcover <- (Com.Sp.groups.mean$Total+Com.Sp.groups.mean$freq*3)/4
rm(Com.Sp.freq.sum, Com.Sp.groups.count)

#Affinity analysis
Com.Sp.groups.spptotals <- aggregate(Com.Sp.groups.mean[,c('freqcover')],
                                     by=list(Com.Sp.groups.mean$taxon), FUN='sum')
colnames(Com.Sp.groups.spptotals) <- c('taxon', 'Totalfreqcover')
Com.Sp.groups.mean <- merge(Com.Sp.groups.mean, Com.Sp.groups.spptotals, by = 'taxon')
Com.Sp.groups.mean$affinity <- (Com.Sp.groups.mean$freq*3 + Com.Sp.groups.mean$freqcover/(Com.Sp.groups.mean$Totalfreqcover+0.0001)*100)/4

#calculate dominant stratum by species in group.
Com.Sp.groups.mean$stratum <- round(
  ( (Com.Sp.groups.mean$Tree)^2*3+
      (Com.Sp.groups.mean$Subcanopy)^2*2.5+
      (Com.Sp.groups.mean$Shrub)^2*2+
      (Com.Sp.groups.mean$Field)^2*1
  )/
    ( (Com.Sp.groups.mean$Tree)^2+
        (Com.Sp.groups.mean$Subcanopy)^2+
        (Com.Sp.groups.mean$Shrub)^2+
        (Com.Sp.groups.mean$Field)^2
      +0.000000001),0)

#Classify Structure
Com.Sp.Agg.groups <- merge(groupdf,  Com.Sp.Agg, by='soilplot', all.x=TRUE, all.y = TRUE)
Com.Sp.Agg.groups.sum <- aggregate(Com.Sp.Agg.groups$Total, by=list(Com.Sp.Agg.groups$clust, Com.Sp.Agg.groups$Simple), FUN='sum')
colnames(Com.Sp.Agg.groups.sum) <- c('cluster', 'simple', 'Total')
Com.Sp.Agg.count <- aggregate(unique(Com.Sp.Agg.groups[c('clust', 'soilplot')])$soilplot, 
                              by=list(unique(Com.Sp.Agg.groups[c('clust', 'soilplot')])$clust), FUN='length')
colnames(Com.Sp.Agg.count) <- c('cluster', 'count')
Com.Sp.Agg.groups.mean <- merge(Com.Sp.Agg.groups.sum, Com.Sp.Agg.count, by = 'cluster')
Com.Sp.Agg.groups.mean$Total <- Com.Sp.Agg.groups.mean$Total/(Com.Sp.Agg.groups.mean$count+0.0001)
rm(Com.Sp.Agg.groups.sum, Com.Sp.Agg.count)
listofsimple <- c("DeciduousShrub", "DeciduousTree", "EvergreenShrub", "EvergreenTree", "Forb", "Graminoid", "Nonvascular")
cluster = unique(Com.Sp.Agg.groups.mean$cluster)
Com.Structure <- as.data.frame(cluster)
for (i in 1:length(listofsimple)){
  x <- Com.Sp.Agg.groups.mean[Com.Sp.Agg.groups.mean$simple %in% listofsimple[i],c('cluster','Total')]
  
  Com.Structure <- merge(Com.Structure,x, by='cluster', all.x = TRUE)
  Com.Structure$Total <- ifelse(is.na(Com.Structure$Total),0,Com.Structure$Total)
  colnames(Com.Structure)[colnames(Com.Structure)=="Total"] <- listofsimple[i]
}

Com.Structure$Structure <-ifelse(Com.Structure$EvergreenTree + Com.Structure$DeciduousTree >= 60,'Forestland',
                                 ifelse(Com.Structure$EvergreenTree + Com.Structure$DeciduousTree >= 25,'Woodland',
                                        ifelse(Com.Structure$EvergreenShrub + Com.Structure$DeciduousShrub >= 25,'Shrubland',
                                               ifelse(Com.Structure$Forb + Com.Structure$Graminoid >= 25,'Herbland',
                                                      ifelse(Com.Structure$Nonvascular >= 25,'Mossland',
                                                             ifelse(Com.Structure$EvergreenTree + Com.Structure$DeciduousTree +
                                                                      Com.Structure$EvergreenShrub + Com.Structure$DeciduousShrub +
                                                                      Com.Structure$Forb + Com.Structure$Graminoid + 
                                                                      Com.Structure$Nonvascular >= 1,'Sparse','Barren'))))))
Com.Structure$Structure <- 
  ifelse(Com.Structure$Structure %in% 'Forestland' &
           Com.Structure$EvergreenTree/(Com.Structure$EvergreenTree+Com.Structure$DeciduousTree)*100 >=75, 'Evergreen Forest',
         ifelse(Com.Structure$Structure %in% 'Forestland' &
                  Com.Structure$EvergreenTree/(Com.Structure$EvergreenTree+Com.Structure$DeciduousTree)*100 >=25, 'Mixed Forest',
                ifelse(Com.Structure$Structure %in% 'Forestland', 'Deciduous Forest',Com.Structure$Structure)))

Com.Structure$Structure <- 
  ifelse(Com.Structure$Structure %in% 'Shrubland' &
           Com.Structure$EvergreenShrub/(Com.Structure$EvergreenShrub+Com.Structure$DeciduousShrub)*100 >=50, 'Evergreen Shrubland',
         ifelse(Com.Structure$Structure %in% 'Shrubland', 'Deciduous Shrubland',Com.Structure$Structure))

Com.Structure$Structure <- 
  ifelse(Com.Structure$Structure %in% 'Herbland' &
           Com.Structure$Graminoid/(Com.Structure$Graminoid+Com.Structure$Forb)*100 >=50, 'Grassland',
         ifelse(Com.Structure$Structure %in% 'Herbland', 'Meadow',Com.Structure$Structure))

Com.Structure$Structure <- 
  ifelse(Com.Structure$Structure %in% 'Grassland' &
           Com.Structure$EvergreenTree+Com.Structure$DeciduousTree >= 10, 'Treed Grassland',
         ifelse(Com.Structure$Structure %in% 'Grassland' &
                  Com.Structure$EvergreenShrub+Com.Structure$DeciduousShrub >= 10, 'Shrubby Grassland',Com.Structure$Structure))
Com.Structure$Structure <- 
  ifelse(Com.Structure$Structure %in% 'Meadow' &
           Com.Structure$EvergreenTree+Com.Structure$DeciduousTree >= 10, 'Treed Meadow',
         ifelse(Com.Structure$Structure %in% 'Meadow' &
                  Com.Structure$EvergreenShrub+Com.Structure$DeciduousShrub >= 10, 'Shrubby Meadow',Com.Structure$Structure))
#group wetness
Com.Sp.Wet.groups <- merge(groupdf,  Com.Sp.Agg.wet, by='soilplot', all.x=TRUE, all.y = TRUE)
Com.Sp.Wet.groups <- Com.Sp.Wet.groups[Com.Sp.Wet.groups$Simple %in% 'Wet',]
Com.Sp.Wet.groups <- aggregate(Com.Sp.Wet.groups$Total, by=list(Com.Sp.Wet.groups$clust), FUN='mean')
colnames(Com.Sp.Wet.groups) <- c('cluster','Wetness')
Com.Structure <- merge(Com.Structure, Com.Sp.Wet.groups, by='cluster', all.x = TRUE)
Com.Structure$WetStructure <- 
  ifelse(Com.Structure$Wetness >= 50, paste('Wet', Com.Structure$Structure), Com.Structure$Structure)

#----
#rank
Com.Sp.groups.mean<- merge(Com.Sp.groups.mean, Com.Structure[,c('cluster','Structure')], by='cluster', all.x = TRUE)
Com.Sp.groups.mean$overunder <- ifelse(Com.Sp.groups.mean$stratum == 3, 1,0)
#rank
Com.Sp.groups.mean <- 
  Com.Sp.groups.mean %>%
  group_by(cluster) %>%
  mutate(ranks = order(order(Total, decreasing=TRUE)))

Com.Sp.groups.mean <- 
  Com.Sp.groups.mean %>%
  group_by(cluster) %>%
  mutate(freqranks = order(order(freqcover, decreasing=TRUE)))

Com.Sp.groups.mean <- 
  Com.Sp.groups.mean %>%
  group_by(cluster) %>%
  mutate(affranks = order(order(affinity, decreasing=TRUE)))

Com.Sp.groups.mean <- 
  Com.Sp.groups.mean %>%
  group_by(cluster, stratum) %>%
  mutate(subranks = order(order(freqcover, decreasing=TRUE)))

Com.Sp.groups.mean <- 
  Com.Sp.groups.mean %>%
  group_by(cluster, overunder) %>%
  mutate(overunderranks = order(order(freqcover, decreasing=TRUE)))

Com.rank <- subset(Com.Sp.groups.mean, !stratum %in% 0)




Com.rankA <- subset(Com.rank,
                    ((grepl('Forest',Structure))&
                       (((freqranks <= 2 | affranks <= 1) & stratum %in% c(1,2,3))|(overunderranks <= 1))#forced understory
                     #(((freqranks <= 3 | affranks <= 1) & stratum %in% c(1,2,3))|(subranks <= 1 & stratum %in% c(2)))#can be all overstory
                    )|
                      ((grepl('Woodland',Structure))&
                         (((freqranks <= 2 | affranks <= 1) & stratum %in% c(1,2,3))|(overunderranks <= 1))
                      )|
                      ((grepl('Shrubland',Structure))&
                         (((freqranks <= 3 | affranks <= 1) & stratum %in% c(1,2))|(subranks <= 1 & stratum %in% c(2)))
                      )|
                      ((Structure %in% c('Treed Grassland', 'Treed Meadow'))&
                         (((freqranks <= 2 | affranks <= 1) & stratum %in% c(1,2,3))|(subranks <= 1 & stratum %in% c(1,3)))
                      )|
                      ((Structure %in% c('Shrubby Grassland', 'Shrubby Meadow'))&
                         (((freqranks <= 2 | affranks <= 1) & stratum %in% c(1,2))|(subranks <= 1 & stratum %in% c(1,2)))
                      )|
                      ((Structure %in% c('Grassland', 'Meadow'))&
                         (((freqranks <= 3 | affranks <= 1) & stratum %in% c(1))|(subranks <= 1 & stratum %in% c(1)))
                      ))

Com.rankB <- subset(Com.rank,
                    affranks <= 2|(subranks <= 1 & freqranks <= 5)|(subranks <= 1 & affranks <= 5))

Com.Ass <- Com.rankA
Com.Ass <- Com.Ass[,c('cluster', 'taxon', 'stratum', 'freqcover')]
Com.Ass <- 
  Com.Ass %>%
  group_by(cluster, stratum) %>%
  mutate(ranks = order(order(freqcover, decreasing=FALSE)))
Com.Ass$ranks <- Com.Ass$stratum*10+Com.Ass$ranks 
Com.Ass <- 
  Com.Ass %>%
  group_by(cluster) %>%
  mutate(ranks = order(order(ranks, decreasing=TRUE)))

Com.Structure$association <- ""

Com.Ass$taxon <- as.character(Com.Ass$taxon)
nclust <- unique(Com.Ass$cluster)
for (i in 1:length(nclust)){
  Com.B <- subset(Com.Ass, cluster %in% nclust[i])
  nrank <- length(unique(Com.B$ranks))
  assname <- ""
  for (j in 1:nrank){
    assname <- ifelse(j == 1,Com.B[Com.B$ranks %in% j,]$taxon,
                      ifelse(Com.B[Com.B$ranks %in% j,]$stratum == Com.B[Com.B$ranks %in% (j-1),]$stratum,
                             paste0(assname, '-',Com.B[Com.B$ranks %in% j,]$taxon),paste0(assname, '/',Com.B[Com.B$ranks %in% j,]$taxon)))
    
  }
  Com.Structure[Com.Structure$cluster %in% nclust[i],]$association <- assname
}
Com.Structure[order(as.numeric(as.character(Com.Structure$cluster))),c("cluster", "association", "WetStructure")]
###end associated spp ###
sil <- (distbray %>% agnes(method = 'ward') %>% cutree(k=10) %>% silhouette(distbray))[,] %>% as.data.frame() 
sil.summary <- aggregate(sil[,c('sil_width')], by=list(cluster = sil$cluster ), FUN='mean') %>%  `colnames<-`(c('cluster', 'sil_width'))
