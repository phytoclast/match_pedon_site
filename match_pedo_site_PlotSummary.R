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
library(isopam)


#----
NASISPEDONS <- readRDS("data/NASISPEDONS.RDS")
fp <- readRDS('data/fp.RDS')
names(NASISPEDONS)[names(NASISPEDONS) == "x_std"] <- 'Std.Longitude'
names(NASISPEDONS)[names(NASISPEDONS) == "y_std"] <- 'Std.Latitude'
names(NASISPEDONS)[names(NASISPEDONS) == "taxonname"] <- 'Current.Taxon.Name'
NASISPEDONS$Current.Taxonomic.Class <- paste(NASISPEDONS$taxpartsize, NASISPEDONS$taxtempregime, NASISPEDONS$taxsubgrp, sep = ', ')
NASISPEDONS$Current.Taxonomic.Class <- str_replace(NASISPEDONS$Current.Taxonomic.Class, 'NA, ','')
NASISPEDONS$Current.Taxonomic.Class <- str_replace(NASISPEDONS$Current.Taxonomic.Class, 'not used, ','')
#write.csv(NASISPEDONS, 'data/NASISPEDON.csv')
names(NASISPEDONS)[names(NASISPEDONS) == "pedon_id"] <- 'User.Site.ID'

pedonoverride <- read.delim("data/pedonoverride.txt")
s <- read.delim("data/s.txt")

List_Habits <- read.delim("data/List_Habits.txt", na.strings = '')
obsspp <- read.delim("data/Observed_Species.txt", encoding = 'UTF-8', na.strings = '')
obsspp[obsspp$AcTaxon == 'Phalaris' & !is.na(obsspp$AcTaxon) ,]$AcTaxon <- 'Phalaris arundinacea'
obsspp <- obsspp[!grepl("\\?", obsspp$AcTaxon) | is.na(obsspp$AcTaxon) ,]
obsspp <- subset(obsspp, !is.na(specific_epithet) | AcTaxon %in% c('Sphagnum', 'Chara'))

#listspp <- read.delim("data/List_Species2011.txt", encoding = 'UTF-8', na.strings = '')
listspp <- readRDS("data/listspp.RDS")
#saveRDS(listspp, 'data/listspp.RDS')
obsspp <- subset(obsspp, !substr(AcTaxon,1,1) %in% '-'& !AcTaxon %in% '' & !is.na(Habit))
obsspp <- merge(obsspp, List_Habits[,c('Form','Simple')], by.x = 'Habit', by.y = 'Form', all.x = TRUE)
obs <- read.delim("data/Sites.txt")
obs <- subset(obs,Observer_Code %in% 'GRR.GJS' & Year >=2011 & !Observation_Type %in% c('Bogus', 'Floristics')) #filter out bogus points
unique(obs$Observation_Type)

obsspp <- merge(obs[,c('Observation_ID','Observation_Label')],obsspp, by='Observation_ID')
obsspp <- subset(obsspp, Field+Shrub+Subcanopy+Tree > 0)
if (F){ #if true, remove ambiguous taxa
  obsspp <- subset(obsspp, grepl(' ', AcTaxon) & !grepl('\\?', AcTaxon))}
if (F){ #if true, remove bad invasives
  obsspp <- subset(obsspp, !grepl('Phalaris', AcTaxon) | grepl('Rosa multiflora', AcTaxon))}
#VEGOBS <- read.delim("data/VEGOBS.txt")
VEGOBS <- subset(obs, !(Latitude == 0 & Longitude == 0) & Year > 1990 & Mon > 0, select = c("Observation_ID", "Observation_Label", "Observation_Type","Latitude","Longitude","Year","Mon","Day","State","County", "Soil.Series" ))
VEGOBS$Soil.Series <- str_replace_all(VEGOBS$Soil.Series, ',', ' ')
VEGOBS$Soil.Series <- str_replace_all(VEGOBS$Soil.Series, '\\?', ' ')
VEGOBS$Soil.Series <- str_replace_all(VEGOBS$Soil.Series, '/', ' ')
VEGOBS$pedon <- ""
VEGOBS$taxonname <- ""
VEGOBS$taxonclass <- ""
VEGOBS$pedondist <- NA
VEGOBS$pedondate <- 0
n <- nrow(VEGOBS)
for (i in 1:n){
  
  
  NASISPEDONS$distance <- (((VEGOBS[i,]$Latitude - NASISPEDONS$Std.Latitude)/360*40041.47*1000)^2 +
                             ((VEGOBS[i,]$Longitude - NASISPEDONS$Std.Longitude)/360*40041.47*1000*cos(VEGOBS[i,]$Latitude/2/360*2*3.141592))^2)^0.5
  NASISPEDONS$distance2 <- ifelse(!is.na(VEGOBS[i,]$Soil.Series)&!is.na(NASISPEDONS$Current.Taxon.Name)&
                                    (str_split_fixed(VEGOBS[i,]$Soil.Series, " ",2)[,1])==
                                    (str_split_fixed(NASISPEDONS$Current.Taxon.Name, " ",2)[,1]), 0.1,1)
  #NASISPEDONS$distance <- NASISPEDONS$distance*NASISPEDONS$distance2
  NASISPEDONS$distance <- NASISPEDONS$distance+(NASISPEDONS$distance2-1)*222.2
  mindist <- min(NASISPEDONS$distance, na.rm = TRUE)
  if(!is.na(subset(NASISPEDONS, NASISPEDONS$User.Site.ID == VEGOBS[i,]$Observation_Label )[1,]$User.Site.ID)){
    VEGOBS[i,]$pedondist <- -1000
    VEGOBS[i,]$pedon <- as.character(subset(NASISPEDONS, NASISPEDONS$User.Site.ID == VEGOBS[i,]$Observation_Label )[1,]$User.Site.ID)
    VEGOBS[i,]$taxonname <- as.character(subset(NASISPEDONS, NASISPEDONS$User.Site.ID == VEGOBS[i,]$Observation_Label )[1,]$Current.Taxon.Name)
    VEGOBS[i,]$taxonclass <- as.character(subset(NASISPEDONS, NASISPEDONS$User.Site.ID == VEGOBS[i,]$Observation_Label )[1,]$Current.Taxonomic.Class)
    VEGOBS[i,]$pedondate <- as.numeric(substr(as.character(subset(NASISPEDONS, NASISPEDONS$User.Site.ID == VEGOBS[i,]$Observation_Label )[1,]$obs_date), 1,4))
  }else{
    VEGOBS[i,]$pedondist <- round(mindist,1)
    VEGOBS[i,]$pedon <- as.character(subset(NASISPEDONS, NASISPEDONS$distance == mindist)[1,]$User.Site.ID)
    VEGOBS[i,]$taxonname <- as.character(subset(NASISPEDONS, NASISPEDONS$distance == mindist)[1,]$Current.Taxon.Name)
    VEGOBS[i,]$taxonclass <- as.character(subset(NASISPEDONS, NASISPEDONS$distance == mindist)[1,]$Current.Taxonomic.Class)
    VEGOBS[i,]$pedondate <- as.numeric(substr(as.character(subset(NASISPEDONS, NASISPEDONS$distance == mindist)[1,]$obs_date), 1,4))
    #VEGOBS[i,]$pedondate <- as.numeric(substr(as.character(NASISPEDONS[NASISPEDONS$distance == mindist,]$obs_date[1]), 1,4))
  }
}
#VEGOBS <- merge(VEGOBS, NASISPEDONS[,c('User.Site.ID', 'Current.Taxon.Name', 'Current.Taxonomic.Class')], by.x = 'Observation_Label', by.y = 'User.Site.ID', all.x = TRUE)


mu <- readRDS(file='data/mu.RDS')

VEGOBS_mukeys <- foreign::read.dbf("data/sitemukey.dbf") 
VEGOBS_soilnames <- merge(VEGOBS_mukeys[,c('Observat_3','RASTERVALU')], mu[,c('lmapunitiid', 'muname')], by.x='RASTERVALU', by.y= 'lmapunitiid')
VEGOBS <- merge(VEGOBS, VEGOBS_soilnames[,c('Observat_3', 'muname')], by.x='Observation_Label', by.y= 'Observat_3')

VEGOBS$eval <- "dump"
VEGOBS[VEGOBS$pedondist < 50,]$eval <- "keep1" 
VEGOBS[(str_split_fixed(VEGOBS$muname, " ",2)[,1])==VEGOBS$taxonname & VEGOBS$pedondist < 1000 & VEGOBS$eval == 'dump', ]$eval <- "keep2"
VEGOBS[VEGOBS$pedondate==VEGOBS$Year & VEGOBS$pedondist < 200 & VEGOBS$eval == 'dump', ]$eval <- "keep3"
#VEGOBS[VEGOBS$Observation_Label==VEGOBS$pedon & VEGOBS$eval == 'dump', ]$eval <- "keep4"
VEGOBS[VEGOBS$eval == 'dump',]$taxonname <- ""
VEGOBS[VEGOBS$eval == 'dump',]$taxonclass <- ""
VEGOBS[VEGOBS$eval == 'dump',]$pedon <- ""
VEGOBS$Soil <- VEGOBS$taxonname
for (i in 1:nrow(pedonoverride)){ #replace known soils that disagree with map unit and don't have pedon records
  VEGOBS[VEGOBS$Observation_Label %in% pedonoverride$sitelabel[i],]$Soil <- as.character(pedonoverride$soil[i])}
VEGOBS[VEGOBS$Soil %in% '',]$Soil <- str_split_fixed(VEGOBS[VEGOBS$Soil %in% '',]$muname, " ",2)[,1]
write.csv(VEGOBS, 'output/VEGOBS.csv', row.names = FALSE, na = "")


#----
#narrow to soil series
soilgroup <- 'all'
ngroups <- 18
if (T){
  remove <- c('2016MI037002',
              'S12060501', 
              'S12060502',
              'S12060503',
              'S12071304',
              's20190701.02'
  )
  add <- 'S12062503'
  sortsoils <- unique(subset(s, (T150_OM >= 10 | grepl('histic',taxsubgrp) |grepl('histosols',taxorder)) & !grepl('Dysic',taxclname), select = 'compname'))[,1]
  VEGOBS <- subset(VEGOBS,Soil %in% sortsoils & !Observation_Label %in% remove |Observation_Label %in% add )
  ngroups <- 8
  soilgroup <- 'euic_mucks'}


#----
#observed species means by plot

Com.Sp.sum<-aggregate(obsspp[,c('Field', 'Shrub', 'Subcanopy', 'Tree')], by=list(obsspp$Observation_Label, obsspp$AcTaxon, obsspp$Simple), FUN=sum) #sum within plot
colnames(Com.Sp.sum)<-c('Observation_Label', 'Species', 'Simple', 'Field', 'Shrub', 'Subcanopy', 'Tree') #restore column names
#freq
Com.Sp.freq<-aggregate(obsspp[,c('AcTaxon')], by=list(obsspp$Observation_Label, obsspp$AcTaxon), FUN=length) #frequency within plot
colnames(Com.Sp.freq)<- c('Observation_Label', 'Species', 'freq')
Com.max.freq<-aggregate(Com.Sp.freq[,c('freq')], by=list(Com.Sp.freq$Observation_Label), FUN=max) #freq within plot
colnames(Com.max.freq)<- c('Observation_Label', 'mfreq')
Com.max.freq$mfreq<-ifelse(Com.max.freq$mfreq>4,4,Com.max.freq$mfreq)#effectively ensureing values do not exceed 4. Species listed 5 times might occur if surveyer was unaware of species already counted in subplots, but this only adds a trace amount.
Com.Sp.mean<-merge(Com.Sp.sum, Com.max.freq[,c("Observation_Label","mfreq")], by="Observation_Label")
Com.Sp.mean$Field<-Com.Sp.mean$Field/Com.Sp.mean$mfreq
Com.Sp.mean$Shrub<-Com.Sp.mean$Shrub/Com.Sp.mean$mfreq
Com.Sp.mean$Subcanopy<-Com.Sp.mean$Subcanopy/Com.Sp.mean$mfreq
Com.Sp.mean$Tree<-Com.Sp.mean$Tree/Com.Sp.mean$mfreq
rm(Com.max.freq)
#heights
Com.Sp.hts<-aggregate(obsspp[,c('Fmin', 'Fmax', 'Smin', 'Smax', 'SCmin', 'SCmax','Tmin', 'Tmax')], by=list(obsspp$Observation_Label, obsspp$AcTaxon, obsspp$Simple), FUN=mean, na.action = na.omit) #sum within plot
colnames(Com.Sp.hts)<-c('Observation_Label', 'Species', 'Simple', 'Fmin', 'Fmax', 'Smin', 'Smax', 'SCmin', 'SCmax','Tmin', 'Tmax') #restore column names
Com.Sp.mean <- merge(Com.Sp.mean,Com.Sp.hts, by=c('Observation_Label', 'Species', 'Simple'), all.x = T)

#ensure not to exceed 100%
Com.Sp.mean$Field <- ifelse(Com.Sp.mean$Field > 100,100,Com.Sp.mean$Field)
Com.Sp.mean$Shrub <- ifelse(Com.Sp.mean$Shrub > 100,100,Com.Sp.mean$Shrub)
Com.Sp.mean$Subcanopy <- ifelse(Com.Sp.mean$Subcanopy > 100,100,Com.Sp.mean$Subcanopy)
Com.Sp.mean$Tree <- ifelse(Com.Sp.mean$Tree > 100,100,Com.Sp.mean$Tree)
#average overstory and understory
Com.Sp.mean$Total <- 100*(1-10^(apply(log10(1-(Com.Sp.mean[,c('Field', 'Shrub', 'Subcanopy', 'Tree')]/100.001)), MARGIN = 1, FUN='sum')))
Com.Sp.mean <- merge(Com.Sp.mean, VEGOBS[,c('Observation_Label', 'Soil')], by='Observation_Label')
Com.Sp.mean <-subset(Com.Sp.mean, !substr(Species,1,1) %in% '-'& !Species %in% '')
Com.Sp.mean$soilplot <- paste(Com.Sp.mean$Soil , Com.Sp.mean$Observation_Label)
Com.Sp.mean$soilplot <- str_replace_all(Com.Sp.mean$soilplot, ' ', '.')
Com.Sp.mean$soilplot <- str_replace_all(Com.Sp.mean$soilplot, '-', '.')
Com.Sp.mean$soilplot <- str_replace_all(Com.Sp.mean$soilplot, ',', '.')
Com.Sp.mean$soilplot <- str_replace_all(Com.Sp.mean$soilplot, ':', '.')
Com.Sp.mean$soilplot <- str_replace_all(Com.Sp.mean$soilplot, ';', '.')
Com.Sp.mean$soilplot <- str_replace_all(Com.Sp.mean$soilplot, '&', '.')
#----
#need formula for aggregating within strata... 
Com.Sp.preagg <- Com.Sp.mean

Com.Sp.preagg[,c('Field', 'Shrub', 'Subcanopy', 'Tree')] <- log10(1-(Com.Sp.preagg[,c('Field', 'Shrub', 'Subcanopy', 'Tree')]/100.001))
Com.Sp.Agg <- aggregate(Com.Sp.preagg[,c('Field', 'Shrub', 'Subcanopy', 'Tree')], by=list(Com.Sp.preagg$soilplot, Com.Sp.preagg$Simple), FUN = 'sum')
colnames(Com.Sp.Agg) <- c('soilplot', 'Simple', 'Field', 'Shrub', 'Subcanopy', 'Tree')
Com.Sp.Agg$Total <- Com.Sp.Agg$Field + Com.Sp.Agg$Shrub
Com.Sp.AggCanopy <- subset(Com.Sp.Agg,  Simple %in% c('Deciduous','Evergreen'))
Com.Sp.AggCanopy$Simple <- paste0(Com.Sp.AggCanopy$Simple, 'Tree')
Com.Sp.AggCanopy$Total <- Com.Sp.AggCanopy$Subcanopy + Com.Sp.AggCanopy$Tree
Com.Sp.Agg$Simple <-as.character(Com.Sp.Agg$Simple)
Com.Sp.Agg[Com.Sp.Agg$Simple %in% c('Evergreen', 'Deciduous'),]$Simple <- 
  paste0(Com.Sp.Agg[Com.Sp.Agg$Simple %in% c('Evergreen', 'Deciduous'),]$Simple, 'Shrub')
Com.Sp.Agg <- rbind(Com.Sp.Agg, Com.Sp.AggCanopy)
rm(Com.Sp.preagg, Com.Sp.AggCanopy)
Com.Sp.Agg$Total <- (10^(Com.Sp.Agg$Total)*-1+1)*100
#----
#wetland indicator status 
Com.Sp.wet <- merge(Com.Sp.mean,listspp[!is.na(listspp$Wetness),c('Taxon', 'Wetness')], by.x = 'Species', by.y = 'Taxon')
Com.Sp.wet$wettotal <- Com.Sp.wet$Total * Com.Sp.wet$Wetness
Com.Sp.wet$drytotal <- Com.Sp.wet$Total * (100-Com.Sp.wet$Wetness)
Com.Sp.wet.agg <- aggregate(Com.Sp.wet[,c('wettotal', 'drytotal', 'Total')], by = list(Com.Sp.wet$soilplot), FUN = 'sum')
colnames(Com.Sp.wet.agg) <- c('soilplot', 'wettotal', 'drytotal', 'Total')
Com.Sp.wet.agg$Wet <- Com.Sp.wet.agg$wettotal/Com.Sp.wet.agg$Total
Com.Sp.wet.agg$Dry <- Com.Sp.wet.agg$drytotal/Com.Sp.wet.agg$Total
Com.Sp.wet.agg$Total <- Com.Sp.wet.agg$Wet
Com.Sp.wet.agg$Simple <- 'Wet'
Com.Sp.wet.agg2 <- Com.Sp.wet.agg
Com.Sp.wet.agg2$Total <- Com.Sp.wet.agg$Dry
Com.Sp.wet.agg2$Simple <- 'Dry'
Com.Sp.wet.agg <- rbind(Com.Sp.wet.agg, Com.Sp.wet.agg2)
Com.Sp.wet.agg <- Com.Sp.wet.agg[,c('soilplot', 'Simple', 'Total')]

Com.Sp.Agg.wet <- Com.Sp.Agg[,c('soilplot', 'Simple', 'Total')]
Com.Sp.Agg.wet <- rbind(Com.Sp.Agg.wet, Com.Sp.wet.agg)
rm(Com.Sp.wet.agg2,Com.Sp.wet.agg)
#----
#cluster analysis
Com.Sp.mean$sqrttotal <- sqrt(Com.Sp.mean$Total)
#saveRDS(Com.Sp.mean, 'C:/workspace2/USNVC/data/plotdata.RDS')

plotdata <- makecommunitydataset(Com.Sp.mean, row = 'soilplot', column = 'Species', value = 'sqrttotal', drop = TRUE)
Com.Sp.Agg$sqrttotal <- Com.Sp.Agg$Total^0.5
plotinputs2 <- makecommunitydataset(Com.Sp.Agg.wet, row = 'soilplot', column = 'Simple', value = 'Total', drop = TRUE)

#evaluating methods 

distbray <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
distjac <- vegdist(plotdata, method='jaccard', binary=FALSE, na.rm=T)
distkulc <- vegdist(plotdata, method='kulczynski', binary=FALSE, na.rm=T)

#validity using silhouette index
sil.table <-'x'
#----
silanalysis <- function(input){
  distbray <- vegdist(input, method='bray', binary=FALSE, na.rm=T)
  distjac <- vegdist(input, method='jaccard', binary=FALSE, na.rm=T)
  distsim <- as.dist(simil(input,method='Simpson'))
  distkulc <- vegdist(plotdata, method='kulczynski', binary=FALSE, na.rm=T)
  distbraykulc <- (distbray+distkulc)/2
  dward <- agnes(distkulc, method='ward') %>% cophenetic()
  dward <- dward/mean(dward)
  dbyrid <- (dward^(1/2)+distkulc^2*0.5)/1.5
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
  sil.kulc <- 0
  sil.kulcward <- 0
  sil.kulchybrid <- 0

  for (k in 2:24){
    sil.bray1 <- (distbray %>% agnes(method = 'average') %>% cutree(k=k) %>% silhouette(distbraykulc))[,3]%>% mean
    sil.jac1 <- (distjac %>% agnes(method = 'average') %>% cutree(k=k) %>% silhouette(distbraykulc))[,3]%>% mean
    sil.sim1 <- (distsim %>% agnes(method = 'average') %>% cutree(k=k) %>% silhouette(distbraykulc))[,3]%>% mean
    sil.ward1 <- (distbray %>% agnes(method = 'ward') %>% cutree(k=k) %>% silhouette(distbraykulc))[,3]%>% mean
    sil.diana1 <- (distkulc %>% diana %>% cutree(k=k) %>% silhouette(distbraykulc))[,3]%>% mean
    sil.kmeans1 <- (kmeans(distbray, centers = k)$cluster %>% silhouette(distbraykulc))[,3] %>% mean
    sil.kulc1 <- (distkulc %>% agnes(method = 'average') %>% cutree(k=k) %>% silhouette(distbraykulc))[,3]%>% mean
    sil.kulcward1 <- (distkulc %>% agnes(method = 'ward') %>% cutree(k=k) %>% silhouette(distbraykulc))[,3]%>% mean
    sil.kulchybrid1 <- (dbyrid %>% agnes(method = 'average') %>% cutree(k=k) %>% silhouette(distbraykulc))[,3]%>% mean
    
    
    klevel <- c(klevel, k)
    sil.bray <- c(sil.bray, sil.bray1)
    sil.jac <- c(sil.jac, sil.jac1)
    sil.sim <- c(sil.sim, sil.sim1)
    sil.ward <- c(sil.ward, sil.ward1)
    sil.diana <- c(sil.diana, sil.diana1)
    sil.kmeans <- c(sil.kmeans, sil.kmeans1)
    sil.kulc <- c(sil.kulc, sil.kulc1)
    sil.kulcward <- c(sil.kulcward, sil.kulcward1)
    sil.kulchybrid <- c(sil.kulchybrid, sil.kulchybrid1)
    
  }
  sil.table <- as.data.frame(cbind(klevel,sil.bray,sil.jac,sil.sim,sil.ward,sil.diana,sil.kmeans,sil.kulc,sil.kulcward,sil.kulchybrid))
  sil.table <- sil.table[-1,]
  sil.table <<- sil.table
  return(sil.table)}
#----
silanalysis(plotdata)

#plotting function
makeplot <- function(a,d,t, soilgroup,k){
  filename <- paste0('output/Soils_',soilgroup,"_",a,'.png')
  t <- as.hclust(t)
  #make cuts and reformat dendrogram
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
  
  newlabels <- t$labels
  newlabels <- as.data.frame(newlabels)
  newlabels$row <- row(newlabels)
  newlabels <- merge(newlabels, groupdf, by.x='newlabels', by.y ='soilplot')
  newlabels$newlabels <- paste(newlabels$clust, newlabels$newlabels)
  newlabels <- newlabels[order(newlabels$row),1]
  newtree <- t
  newtree$labels <- newlabels
  
  dend1 <- color_branches(as.hclust(newtree), k = k)
  dend1 <- color_labels(dend1, k = k)
  
  #output file
  
  w <- 800
  h <- nrow(d)*12+80
  u <- 12
  png(filename=filename,width = w, height = h, units = "px", pointsize = u)
  
  par(mar = c(2,0,1,13))
  plot(dend1, horiz = TRUE, main=paste('floristic simularity', a,' method of', soilgroup, 'soils'), font=1, cex=0.84)
  dev.off()
  
}

amethod <- 'bray-agnes' 
k=16
if (T){
  a <- 'kulczynski-agnes' 
  k=8
  d <- vegdist(plotdata, method='kulczynski', binary=FALSE, na.rm=T)
  t <- agnes(d, method='average')
  makeplot(a,d,t,soilgroup,k)
}
if (T){
  amethod <- 'kulczynski-ward' 
  k=8
  d <- vegdist(plotdata, method='kulczynski', binary=FALSE, na.rm=T)
  t <- agnes(d, method='ward')
  makeplot(amethod,d,t,soilgroup,k)
}

if (T){
  amethod <- 'kulczynski-hybrid' 
  k=8
  d <- vegdist(plotdata, method='kulczynski', binary=FALSE, na.rm=T)
  dward <- agnes(d, method='ward') %>% cophenetic()
  dward <- dward/mean(dward)
  dbyrid <- (dward^(1/2)+d^2*0.5)/1.5
  t <- agnes(dbyrid, method='average')
  makeplot(amethod,d,t,soilgroup,k)
}

if (T){
  amethod <- 'bray-agnes' 
  k=8
  d <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
  t <- agnes(d, method='average')
  makeplot(amethod,d,t,soilgroup,k)
}
if (T){
  amethod <- 'bray-single' 
  d <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
  t <- agnes(d, method='single')
  makeplot(amethod,d,t,soilgroup,k)
}
if (T){
  amethod <- 'bray-complete' 
  d <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
  t <- agnes(d, method='complete')
  makeplot(amethod,d,t,soilgroup,k)
}
if (T){
  amethod <- 'bray-diana' 
  d <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
  t <- diana(d)
  makeplot(amethod,d,t,soilgroup,k)
}

if (T){
  amethod <- 'bray-ward'
  k=4
  d <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
  t <- agnes(d, method='ward')
  makeplot(amethod,d,t,soilgroup,k)
  amethod <- 'bray-ward-binary'
  k=4
  d <- vegdist(plotdata, method='bray', binary=TRUE, na.rm=T)
  t <- agnes(d, method='ward')
  makeplot(amethod,d,t,soilgroup,k)
}
if (T){
  amethod <- 'jaccard-agnes' 
  d <- vegdist(plotdata, method='jaccard', binary=FALSE, na.rm=T)
  t <- agnes(d, method='average')
  makeplot(amethod,d,t,soilgroup,k)
}

if (T){
  amethod <- 'isopam'
  d <- vegdist(plotdata, method='kulczynski', binary=FALSE, na.rm=T)
  k=7
  tpam <- isopam(plotdata, distance = 'kulczynski', stopat = c(1,7))
  t <- tpam$dendro
  makeplot(amethod,d,t,soilgroup,k)
  pamtab <- isotab(tpam, level = 3)
  pamtab <- pamtab$tab
  pamtab$x1 <- as.numeric(gsub("[^0-9.-]", "", (pamtab$`1`)))
  pamtab$x2 <- as.numeric(gsub("[^0-9.-]", "", (pamtab$`2`)))
  pamtab$dif1_2 <- pamtab$x1-pamtab$x2
  pamtab$x3.1 <- as.numeric(gsub("[^0-9.-]", "", (pamtab$`3.1`)))
  pamtab$x3.2 <- as.numeric(gsub("[^0-9.-]", "", (pamtab$`3.2`)))
  pamtab$dif3.1_3.2 <- pamtab$x3.1-pamtab$x3.2
  pamtab$dif2_3.1 <- pamtab$x2-pamtab$x3.1
  coph <- cophenetic(t)
  tward <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T) %>% agnes(method = 'ward') 
  dward <- tward %>% cophenetic()
  d2 <- (coph/mean(coph)*2 + dward/mean(dward))/3
  t2 <- agnes(d2, method = 'average')
  makeplot('kulczynski-isopam-ward-hybrid',d2,t2,soilgroup,k)
  
  w <- 900
  h <- nrow(d)*15+80
  u <- 18
  png(filename='output/duelingdendro.png',width = w, height = h, units = "px", pointsize = u)
  
  par(mar = c(2,0,1,13))
  sharpshootR::dueling.dendrograms(as.phylo(as.hclust(tward)), as.phylo(as.hclust(t2)))
  dev.off()
  
}



distbray <- vegdist(plotdata, method='bray', binary=F, na.rm=T)
distjac <- vegdist(plotdata, method='jaccard', binary=F, na.rm=T)
distsim <- as.dist(simil(plotdata,method='Simpson'))
distkulc <- vegdist(plotdata, method='kulczynski', binary=F, na.rm=T)
distbraykulc <- (distbray+distkulc)/2
tbrayagnes <- distbray %>% agnes(method = 'average')
tjacagnes <- distjac %>% agnes(method = 'average')
tjacward <- distjac %>% agnes(method = 'ward')
tbrayward <- distbray %>% agnes(method = 'ward')
tbraydiana <- distbray %>% diana 
tsimpagnes <- distsim %>% agnes(method = 'average') 
tkulcagnes <- distkulc %>% agnes(method = 'average')
tkulcward <- distkulc %>% agnes(method = 'ward')
tkulcbrayward <- distbraykulc %>% agnes(method = 'ward')
cor(cophenetic(tbrayagnes), cophenetic(tpam$dendro))
cor(cophenetic(tjacagnes), cophenetic(tpam$dendro))
cor(cophenetic(tjacward), cophenetic(tpam$dendro))
cor(cophenetic(tbraydiana), cophenetic(tpam$dendro))
cor(cophenetic(tbrayward), cophenetic(tpam$dendro))
cor(cophenetic(tsimpagnes), cophenetic(tpam$dendro))
cor(cophenetic(tkulcagnes), cophenetic(tpam$dendro))
cor(cophenetic(tkulcward), cophenetic(tpam$dendro))
cor(cophenetic(tkulcbrayward), cophenetic(tpam$dendro))
#----
#group dominant and indicator species
k=8
d <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
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


#average heights by cluster

Com.Sp.groups.hts <- aggregate(Com.Sp.groups[,c('Fmin', 'Fmax', 'Smin', 'Smax', 'SCmin', 'SCmax','Tmin', 'Tmax')],
                               by=list(Com.Sp.groups$clust, Com.Sp.groups$Species), FUN='mean')
colnames(Com.Sp.groups.hts) <- c('cluster', 'taxon', 'Fmin', 'Fmax', 'Smin', 'Smax', 'SCmin', 'SCmax','Tmin', 'Tmax')

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

#percentile spp by cluster
Com.Sp.groups <- subset(Com.Sp.groups, !is.na(Field) & !is.na(Shrub) & !is.na(Subcanopy) & !is.na(Tree))#remove nulls
Com.Sp.groups.pctl <- ddply(Com.Sp.groups, c('clust', 'Species'), summarise,
                            f05 = quantile(Field, 0.05),
                            f75 = quantile(Field, 0.75),
                            s05 = quantile(Shrub, 0.05),
                            s75 = quantile(Shrub, 0.75),
                            sc05 = quantile(Subcanopy, 0.05),
                            sc75 = quantile(Subcanopy, 0.75),
                            t05 = quantile(Tree, 0.05),
                            t75 = quantile(Tree, 0.75)
                            )
Com.Sp.groups.pctl <- merge(Com.Sp.groups.pctl, Com.Sp.groups.freq, by.x=c('clust', 'Species'), by.y=c('cluster', 'taxon'), all.x = T)
Com.Sp.groups.pctl$f05 <- ifelse(Com.Sp.groups.pctl$freq < 75, 0,
                                 Com.Sp.groups.pctl$f05)
Com.Sp.groups.pctl$s05 <- ifelse(Com.Sp.groups.pctl$freq < 75, 0,
                                 Com.Sp.groups.pctl$s05)
Com.Sp.groups.pctl$sc05 <- ifelse(Com.Sp.groups.pctl$freq < 75, 0,
                                 Com.Sp.groups.pctl$sc05)
Com.Sp.groups.pctl$t05 <- ifelse(Com.Sp.groups.pctl$freq < 75, 0,
                                 Com.Sp.groups.pctl$t05)

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

#saveRDS(Com.Sp.mean, 'C:/workspace2/USNVC/data/Com.Sp.mean.RDS')
#saveRDS(Com.Sp.mean, 'C:/workspace2/USNVC/data/allplotdata.RDS')
