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

obsover <- obsspp
obsover$AcTaxon <- paste('over.',obsover$AcTaxon)
obsover$Field <- obsover$Field/1000
obsover$Shrub <- obsover$Shrub/1000
obsunder <- obsspp
obsunder$AcTaxon <- paste('under.',obsunder$AcTaxon)
obsunder$Subcanopy <- obsunder$Subcanopy/1000
obsunder$Tree <- obsunder$Tree/1000
obsspp <- rbind(obsspp, obsover,obsunder)
obsspp <- subset(obsspp, Field+Shrub+Subcanopy+Tree > 0)



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

#generic structure
Com.Sp.agg <- aggregate(log10(1-(Com.Sp.mean[,c('Field', 'Shrub', 'Subcanopy', 'Tree')]/100.001)), by=list(soilplot = Com.Sp.mean$soilplot),  FUN='sum')
Com.Sp.agg[,c('Field', 'Shrub', 'Subcanopy', 'Tree')] <- 100*(1-10^(Com.Sp.agg[,c('Field', 'Shrub', 'Subcanopy', 'Tree')]))                 
Com.Sp.agg$over <- 100*(1-10^(apply(log10(1-(Com.Sp.agg[,c('Subcanopy', 'Tree')]/100.001)), MARGIN = 1, FUN='sum')))
Com.Sp.agg$under <- 100*(1-10^(apply(log10(1-(Com.Sp.agg[,c('Field', 'Shrub')]/100.001)), MARGIN = 1, FUN='sum')))
rownames(Com.Sp.agg) <- Com.Sp.agg$soilplot
Com.Sp.agg <- Com.Sp.agg[,-1]

#----
#cluster analysis
Com.Sp.mean$sqrttotal <- sqrt(Com.Sp.mean$Total)
#saveRDS(Com.Sp.mean, 'C:/workspace2/USNVC/data/plotdata.RDS')
overlist <- grepl('over.', Com.Sp.mean$Species)
Com.Sp.mean.over <- subset(Com.Sp.mean, overlist)
underlist <- grepl('under.', Com.Sp.mean$Species)
Com.Sp.mean.under <- subset(Com.Sp.mean, underlist)
Com.Sp.mean.structure <-  subset(Com.Sp.mean, underlist|overlist)
Com.Sp.mean <-  subset(Com.Sp.mean, !underlist & !overlist)
plotdata.structure <-  makecommunitydataset(Com.Sp.mean.structure, row = 'soilplot', column = 'Species', value = 'sqrttotal', drop = TRUE)
plotdata.over <- makecommunitydataset(Com.Sp.mean.over, row = 'soilplot', column = 'Species', value = 'sqrttotal', drop = TRUE)
plotdata.under <- makecommunitydataset(Com.Sp.mean.under, row = 'soilplot', column = 'Species', value = 'sqrttotal', drop = TRUE)
plotdata <-  makecommunitydataset(Com.Sp.mean, row = 'soilplot', column = 'Species', value = 'sqrttotal', drop = TRUE)

#generic structure
Com.Sp.agg <- aggregate(log10(1-(Com.Sp.mean[,c('Field', 'Shrub', 'Subcanopy', 'Tree')]/100.001)), by=list(soilplot = Com.Sp.mean$soilplot),  FUN='sum')
Com.Sp.agg[,c('Field', 'Shrub', 'Subcanopy', 'Tree')] <- 100*(1-10^(Com.Sp.agg[,c('Field', 'Shrub', 'Subcanopy', 'Tree')]))                 
Com.Sp.agg$over <- 100*(1-10^(apply(log10(1-(Com.Sp.agg[,c('Subcanopy', 'Tree')]/100.001)), MARGIN = 1, FUN='sum')))
Com.Sp.agg$under <- 100*(1-10^(apply(log10(1-(Com.Sp.agg[,c('Field', 'Shrub')]/100.001)), MARGIN = 1, FUN='sum')))
rownames(Com.Sp.agg) <- Com.Sp.agg$soilplot
Com.Sp.agg <- Com.Sp.agg[,-1]

#evaluating methods 

#----

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

if (T){
  k=9
  d.under <- vegdist(plotdata.under, method='bray', binary=FALSE, na.rm=T)
  t.under <- agnes(d.under, method = 'ward')
  makeplot('bray-ward-understory',d.under,t.under,soilgroup,k)
  
  d.over <- vegdist(plotdata.over, method='bray', binary=FALSE, na.rm=T)
  t.over <- agnes(d.over, method = 'ward')
  makeplot('bray-ward-overstory',d.over,t.over,soilgroup,k)
  
  d.structure <- vegdist(plotdata.structure, method='bray', binary=FALSE, na.rm=T)
  t.structure <- agnes(d.structure, method = 'ward')
  makeplot('bray-ward-structure',d.structure,t.structure,soilgroup,k)
  
  d <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
  t <- agnes(d, method = 'ward')
  makeplot('bray-ward-nostructure',d,t,soilgroup,k)
  
  d.genstr <- vegdist(Com.Sp.agg, method = 'euclidean', binary=FALSE, na.rm=T)
  t.genstr <- agnes(d.genstr, method = 'ward')
  makeplot('bray-ward-genericstructure',d.genstr,t.genstr,soilgroup,k)
  
  d.compound <- (d.structure/mean(d.structure)*2 + d.genstr/mean(d.genstr))/3
  t.compound <- agnes(d.compound, method = 'ward')
  makeplot('bray-ward-compound',d.compound,t.compound,soilgroup,k)

  
  w <- 900
  h <- nrow(d)*15+80
  u <- 18
  png(filename='output/dueling.structure.png',width = w, height = h, units = "px", pointsize = u)
  
  par(mar = c(2,0,1,13))
  sharpshootR::dueling.dendrograms(as.phylo(as.hclust(t)), as.phylo(as.hclust(t.compound)))
  dev.off()
  
}
