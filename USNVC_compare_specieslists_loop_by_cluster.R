betasim <- function(p){
  d <- matrix(1, nrow = nrow(p), ncol = nrow(p))
  rownames(d) <- rownames(p)
  colnames(d) <- rownames(p)
  for(j in 1:nrow(p)){
    for(k in 1:nrow(p)){
      d[j,k] <- 1-sum((p[j,]*p[k,])^0.5)/sqrt(sum(p[j,])*sum(p[k,]))
    }}
  d<-as.dist(d)
}

ecoregion <-  read.delim('data/USNVC/d_usfs_ecoregion2007.txt', encoding = 'UTF-8', na.strings = '', stringsAsFactors=FALSE)
vegecoregion <-  read.delim('data/USNVC/UnitXEcoregionUsfs2007.txt', encoding = 'UTF-8', na.strings = '', stringsAsFactors=FALSE)
vegecoregion <- merge(vegecoregion, ecoregion, by='usfs_ecoregion_2007_id')

unit <-  read.delim('data/USNVC/unit.txt', encoding = 'UTF-8', na.strings = '', stringsAsFactors=FALSE)
states <-read.delim('data/USNVC/d_subnation.txt')
vegstates <-read.delim('data/USNVC/UnitXSubnation.txt')
vegstates <- merge(vegstates, states, by='subnation_id')

USNVClist <- readRDS('data/USNVC/USNVClist.RDS')
USNVClist <- subset(USNVClist, !grepl('\\.', acctaxon))



plotgroup <- subset(Com.Sp.groups.mean, select=c(cluster, taxon, freqcover))
plotgroupsum <- plotgroup %>% group_by(cluster, taxon) %>% summarise(sum = sum(freqcover))
plotgroupsum <- plotgroupsum %>% group_by(cluster) %>% mutate(max = max(sum), Imp = sum/max)


states <- c('MI','IN','OH')
core <- c('MI')
ecoregion <- c('222')
level <- 'Association'
level <- unique(unit[unit$hierarchylevel %in% level,'element_global_id'])
states <- unique(vegstates[vegstates$subnation_code %in% states,'element_global_id'])
core <- unique(vegstates[vegstates$subnation_code %in% core,'element_global_id'])
ecoregion <- unique(vegecoregion[vegecoregion$usfs_ecoregion_2007_concat_cd %in% ecoregion,'element_global_id'])

plotassociations <- as.data.frame(lapply(as.data.frame(cbind(clust = 'x', 'element_global_id'=0, 'scientificname'='x')), as.character), stringsAsFactors=FALSE)

for(i in 1:k){ #i=1
  g <- subset(plotgroupsum, cluster %in% i)
  g$Imp <- g$Imp^1 #aggregate data should not be square rooted
  
  gtotal <- sum(g$Imp)
  vegtotal <- aggregate(USNVClist$x, by=list(element_global_id = USNVClist$element_global_id, scientificname = USNVClist$scientificname), FUN='sum')
  names(vegtotal)[names(vegtotal)=='x'] <-'vegtotal'
  gmerge <- merge(g, USNVClist, by.x = 'taxon', by.y = 'acctaxon')
  gmerge$intersect <- (gmerge$Imp*gmerge$x)^0.5
  gintersect <- aggregate(gmerge$intersect, by=list(element_global_id = gmerge$element_global_id, scientificname = gmerge$scientificname), FUN='sum')
  names(gintersect)[names(gintersect)=='x'] <-'intersect'
  g <- merge(gintersect, vegtotal, by=c('element_global_id', 'scientificname'))
  g$affinity <- g$intersect/(g$vegtotal*gtotal)^0.5*100
  g$state <- 'no'
  if(nrow(g[g$element_global_id %in% states,])>0){
    g[g$element_global_id %in% states,]$state <- 'near'}
  if(nrow(g[g$element_global_id %in% core,])>0){
    g[g$element_global_id %in% core,]$state <- 'yes'}
  g$ecoregion <- 'no'
  if(nrow(g[g$element_global_id %in% ecoregion,])>0){
    g[g$element_global_id %in% ecoregion,]$ecoregion <- 'yes'}
  g$level <- 'no'
  g[g$element_global_id %in% level,]$level <- 'yes'
  g <- subset(g, level == 'yes')
  g$best <- g$affinity
  g[g$ecoregion == 'no',]$best <- g[g$ecoregion == 'no',]$best*0.75
  g[g$element_global_id %in% states,]$best <- g[g$element_global_id %in% states,]$best*1/0.75
  g[g$element_global_id %in% core,]$best <- g[g$element_global_id %in% core,]$best*1/0.75
  g$best <- g$best/max(g$best)*100
  g <- g[,!colnames(g)%in% c('intersect','vegtotal')]
  rm(gmerge, vegtotal, gintersect)
  g <- subset(g, best >= 0 & level == 'yes')# & state == 'yes')
  g <- g[order(g$best, decreasing = TRUE),]
  plotassociations1 <- as.data.frame(lapply(as.data.frame(cbind(clust = i,g[1,1:2])), as.character), stringsAsFactors=FALSE)
  
  plotassociations <- rbind(plotassociations,plotassociations1)
}
rm(plotassociations1)
plotassociations <- plotassociations[-1,]

write.csv(plotassociations, 'output/allclustassociations.csv', row.names = F)