library(stringr)
library(BiodiversityR)
library(cluster)
library(ape)
library(dendextend)
library(dplyr)
library(plyr)
library(dynamicTreeCut)
library(proxy)

substitutesymbols <- as.data.frame(cbind(SYMBOL2=c('2AG','TRBO2','SMLA3','CAREX'), form2=c('Nonvascular','Forb/Herb','Forb/Herb','Grass/grass-like (Graminoids)'), type2=c('Native','Native','Native','Native'), taxon=c('Chara','Lysimachia borealis','Smilax lasioneuron','Carex [ovales]')), stringsAsFactors = FALSE)


makeplot <- function(a,d,t,k){
  filename <- paste0('output/plot_',"_",a,'.png')
  
  #make cuts and reformat dendrogram

  groups <- cutree(t, k = k)
  
  plot <- names(groups)
  cluster <- unname(groups)
  groupdf <- as.data.frame(cbind(plot, cluster))
  groupdf$cluster <- (as.numeric(as.character(groupdf$cluster)))
  maxcluster <- max(groupdf$cluster)
  numberzeros <- nrow(groupdf[(groupdf$cluster == 0),])
  whichrecords <- which(groupdf$cluster == 0)
  if (nrow(groupdf[groupdf$cluster == 0,]) != 0){
    for (i in 1:numberzeros){ #assign all zero clusters to unique cluster number.
      groupdf[whichrecords[i],]$cluster <- maxcluster+i}}
  
  newlabels <- t$order.lab
  newlabels <- as.data.frame(newlabels)
  newlabels$row <- row(newlabels)
  newlabels <- merge(newlabels, groupdf, by.x='newlabels', by.y ='plot')
  newlabels$newlabels <- paste(newlabels$cluster, newlabels$newlabels)
  newlabels <- newlabels[order(newlabels$row),1]
  newtree <- t
  newtree$order.lab <- newlabels
  
  dend1 <- color_branches(as.hclust(newtree), k = k)
  dend1 <- color_labels(dend1, k = k)
  
  #output file
  
  w <- 800
  h <- nrow(d)*12+80
  u <- 12
  png(filename=filename,width = w, height = h, units = "px", pointsize = u)
  
  par(mar = c(2,0,1,13))
  plot(dend1, horiz = TRUE, main=paste('floristic simularity', a), font=1, cex=0.84)
  dev.off()
  
}
roundF<-function(p){
  p<-ifelse(p<0.5, floor(p/0.1+0.5)*0.1,ifelse(p<2, floor(p/0.5+0.5)*0.5, ifelse(p<5, floor(p+0.5),floor(p/5+0.5)*5)))
}
BA.fact10usc.metric<-function(p){
  p*10000/43560*10
}
BA.metric.usc.round<-function(p){
  round(p*43560/10000,1)
}
cm.in<-function(p){
  round(p/2.54,1)
}
m.ft.round<-function(p){
  round(p/0.3048,3)
}

#----
plots <- c(
  'GRR.GJS.2016.21',
  'GRR.GJS.2016.59',
  'GRR.GJS.2016.30',
  'GRR.GJS.2016.32',
  'GRR.2011.GJS.12',
  'GRR.GJS.2015.27',
  'GRR.GJS.2015.26',
  'GRR.GJS.2015.28',
  'GRR.GJS.2015.20',
  'GRR.GJS.2015.22',
  'GRR.GJS.2015.21',
  'GRR.GJS.2015.25',
  'GRR.GJS.2015.29',
  'GRR.GJS.2015.30',
  'GRR.GJS.2017.8',
  'GRR.GJS.2017.9',
  'GRR.GJS.2017.19',
  'GRR.GJS.2017.23',
  'GRR.GJS.2018.13',
  'GRR.GJS.2018.14',
  'GRR.GJS.2018.3',
  'GRR.GJS.2018.4',
  'GRR.GJS.2018.5',
  'GRR.GJS.2018.17',
  'GRR.GJS.2018.21',
  'GRR.GJS.2018.22',
  'GRR.GJS.2018.24',
  'GRR.GJS.2018.28',
  'GRR.GJS.2018.29',
  'GRR.GJS.2015.24',
  'GRR.GJS.2015.23',
  'GRR.GJS.2012.23',
  'GRR.GJS.2012.24',
  'GRR.GJS.2012.25',
  'GRR.GJS.2012.31',
  'GRR.GJS.2012.70',
  'GRR.GJS.2012.34',
  'GRR.GJS.2012.35',
  'GRR.GJS.2012.36',
  'GRR.GJS.2012.40',
  'GRR.GJS.2019.22',
  'GRR.GJS.2019.11',
  'GRR.GJS.2019.12',
  'GRR.GJS.2019.13',
  'GRR.GJS.2019.17',
  'GRR.GJS.2014.39',
  'GRR.GJS.2014.40',
  'GRR.GJS.2020.13',
  'GRR.GJS.2020.14'
)
#----
List_Habits <- read.delim("data/List_Habits.txt", na.strings = '', stringsAsFactors = FALSE)
Plant_Heights <- read.delim("data/Plant_Heights.txt", na.strings = '', stringsAsFactors = FALSE, encoding = 'UTF-8')
List_Habits[List_Habits$ESIS.Group %in% 'Grass/grass-like',]$ESIS.Group <- 'Grass/grass-like (Graminoids)'
handpicked  <- read.delim("data/handpicked.txt", na.strings = '', stringsAsFactors = FALSE)
listspp <- read.delim("data/List_Species2011.txt", encoding = 'UTF-8', na.strings = '', stringsAsFactors = FALSE)
#listspp <- readRDS("data/listspp.RDS")
listspp$Nativity <- ifelse(listspp$Eastern.North.America %in% 'N','Native',ifelse(listspp$Eastern.North.America %in% 'X','Introduced','Unknown'))
obs <- read.delim("data/Sites.txt")
obsspp <- read.delim("data/Observed_Species.txt", encoding = 'UTF-8', na.strings = '')
obsspp <- merge(obs[,c('Observation_ID','Observation_Label')],obsspp, by='Observation_ID')
obsspp[obsspp$AcSpecies == 'Phalaris' & !is.na(obsspp$AcTaxon) ,]$AcTaxon <- 'Phalaris arundinacea'
obsspp <- obsspp[!grepl("\\?", obsspp$AcTaxon) | is.na(obsspp$AcTaxon) ,]
#obsspp <- subset(obsspp, !is.na(specific_epithet) | AcTaxon %in% c('Sphagnum', 'Chara', 'Cladonia'))
obsspp <- merge(obsspp, List_Habits[,c('Form','Simple')], by.x = 'Habit', by.y = 'Form', all.x = TRUE)
obsspp <- subset(obsspp, Observation_ID %in% plots)

obsspp[obsspp$BA==0,]$BA <- NA
obsspp$BA <- BA.fact10usc.metric(obsspp$BA)
#----
#observed species means by plot

Com.Sp.sum<-aggregate(obsspp[,c('Field', 'Shrub', 'Subcanopy', 'Tree', 'BA')], by=list(obsspp$Observation_ID, obsspp$Observation_Label, obsspp$AcTaxon, obsspp$Simple), FUN=sum) #sum within plot
colnames(Com.Sp.sum)<-c('Observation_ID', 'Observation_Label', 'Species', 'Simple', 'Field', 'Shrub', 'Subcanopy', 'Tree', 'BA') #restore column names
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
Com.Sp.hts<-aggregate(obsspp[,c('Fmin', 'Fmax', 'Smin', 'Smax', 'SCmin', 'SCmax','Tmin', 'Tmax', 'Dmin', 'Dmax')], by=list(obsspp$Observation_Label, obsspp$AcTaxon, obsspp$Simple), FUN=mean, na.action = na.omit) #sum within plot
colnames(Com.Sp.hts)<-c('Observation_Label', 'Species', 'Simple', 'Fmin', 'Fmax', 'Smin', 'Smax', 'SCmin', 'SCmax','Tmin', 'Tmax', 'Dmin', 'Dmax') #restore column names

#Com.Sp.hts[,c('Fmin', 'Fmax', 'Smin', 'Smax', 'SCmin', 'SCmax','Tmin', 'Tmax', 'Dmin', 'Dmax')]<-   lapply(Com.Sp.hts[,c('Fmin', 'Fmax', 'Smin', 'Smax', 'SCmin', 'SCmax','Tmin', 'Tmax', 'Dmin', 'Dmax')], FUN = roundF)







Com.Sp.mean <- merge(Com.Sp.mean,Com.Sp.hts, by=c('Observation_Label', 'Species', 'Simple'), all.x = T)

#ensure not to exceed 100%
Com.Sp.mean$Field <- ifelse(Com.Sp.mean$Field > 100,100,Com.Sp.mean$Field)
Com.Sp.mean$Shrub <- ifelse(Com.Sp.mean$Shrub > 100,100,Com.Sp.mean$Shrub)
Com.Sp.mean$Subcanopy <- ifelse(Com.Sp.mean$Subcanopy > 100,100,Com.Sp.mean$Subcanopy)
Com.Sp.mean$Tree <- ifelse(Com.Sp.mean$Tree > 100,100,Com.Sp.mean$Tree)
#average overstory and understory
Com.Sp.mean$Total <- 100*(1-10^(apply(log10(1-(Com.Sp.mean[,c('Field', 'Shrub', 'Subcanopy', 'Tree')]/100.001)), MARGIN = 1, FUN='sum')))
Com.Sp.mean <-subset(Com.Sp.mean, !substr(Species,1,1) %in% '-'& !Species %in% '')

#----
#hand picked clusters
Com.Sp.groups <- merge(unique(handpicked[,c('Observation_ID', 'phase')]), Com.Sp.mean, by='Observation_ID')
#----
Com.Sp.groups <- merge(groupdf,  Com.Sp.mean, by='plot', all.x=TRUE, all.y = TRUE)

#----
#average spp by cluster

Com.Sp.groups.sum <- aggregate(Com.Sp.groups[,c('Field', 'Shrub', 'Subcanopy','Tree', 'Total')],
                               by=list(Com.Sp.groups$phase, Com.Sp.groups$Species), FUN='sum')
colnames(Com.Sp.groups.sum) <- c('phase', 'Species', 'Field', 'Shrub', 'Subcanopy','Tree', 'Total')
Com.Sp.groups.count <- aggregate(unique(Com.Sp.groups[c('phase', 'plot')])$plot, 
                                 by=list(unique(Com.Sp.groups[c('phase', 'plot')])$cluster), FUN='length')
colnames(Com.Sp.groups.count) <- c('phase', 'count')
Com.Sp.groups.mean <- merge(Com.Sp.groups.sum, Com.Sp.groups.count, by = 'phase')
Com.Sp.groups.mean[,c('Field', 'Shrub', 'Subcanopy','Tree', 'Total')] <- Com.Sp.groups.mean[,c('Field', 'Shrub', 'Subcanopy','Tree', 'Total')]/Com.Sp.groups.mean$count
rm(Com.Sp.groups.sum, Com.Sp.groups.count)

#----
#average heights by cluster

Com.Sp.groups.hts <- aggregate(Com.Sp.groups[,c('Fmin', 'Fmax', 'Smin', 'Smax', 'SCmin', 'SCmax', 'Tmin', 'Tmax', 'Dmin', 'Dmax')],
                               by=list(Com.Sp.groups$cluster, Com.Sp.groups$Species), FUN='mean')
colnames(Com.Sp.groups.hts) <- c('cluster', 'Species', 'Fmin', 'Fmax', 'Smin', 'Smax', 'SCmin', 'SCmax','Tmin', 'Tmax', 'Dmin', 'Dmax')

#----
#frequency spp by cluster
Com.Sp.prefreq <- Com.Sp.groups
Com.Sp.prefreq$Total <- ifelse(Com.Sp.prefreq$Total >0, 1,0)
Com.Sp.freq.sum <- aggregate(Com.Sp.prefreq$Total,
                             by=list(Com.Sp.prefreq$cluster, Com.Sp.prefreq$Species), FUN='sum')
colnames(Com.Sp.freq.sum) <- c('cluster', 'Species', 'freq')
Com.Sp.groups.count <- aggregate(unique(Com.Sp.prefreq[c('cluster', 'plot')])$plot, 
                                 by=list(unique(Com.Sp.prefreq[c('cluster', 'plot')])$cluster), FUN='length')
colnames(Com.Sp.groups.count) <- c('cluster', 'count')
Com.Sp.groups.freq <- merge(Com.Sp.freq.sum, Com.Sp.groups.count, by = 'cluster')
Com.Sp.groups.freq$freq <- Com.Sp.groups.freq$freq/Com.Sp.groups.freq$count*100
Com.Sp.groups.mean <- merge(Com.Sp.groups.mean, Com.Sp.groups.freq[,c('cluster', 'Species', 'freq')], by = c('cluster', 'Species'))
Com.Sp.groups.mean$freqcover <- (Com.Sp.groups.mean$Total+Com.Sp.groups.mean$freq*3)/4
rm(Com.Sp.freq.sum, Com.Sp.groups.count)
#percentile spp by cluster
Com.Sp.groups.pctl <- ddply(Com.Sp.groups, c('cluster', 'Species'), summarise,
                            f05 = quantile(Field, 0.05),
                            f75 = quantile(Field, 0.75),
                            s05 = quantile(Shrub, 0.05),
                            s75 = quantile(Shrub, 0.75),
                            sc05 = quantile(Subcanopy, 0.05),
                            sc75 = quantile(Subcanopy, 0.75),
                            t05 = quantile(Tree, 0.05),
                            t75 = quantile(Tree, 0.75),
                            b05 = quantile(BA, 0.05, na.rm = TRUE),
                            b95 = quantile(BA, 0.95, na.rm = TRUE)
                            )

Com.Sp.groups.pctl <- merge(Com.Sp.groups.pctl, Com.Sp.groups.freq, by.x=c('cluster', 'Species'), by.y=c('cluster', 'Species'), all.x = T)
Com.Sp.groups.pctl$f05 <- ifelse(Com.Sp.groups.pctl$freq < 75, 0,
                                 Com.Sp.groups.pctl$f05)
Com.Sp.groups.pctl$s05 <- ifelse(Com.Sp.groups.pctl$freq < 75, 0,
                                 Com.Sp.groups.pctl$s05)
Com.Sp.groups.pctl$sc05 <- ifelse(Com.Sp.groups.pctl$freq < 75, 0,
                                  Com.Sp.groups.pctl$sc05)
Com.Sp.groups.pctl$t05 <- ifelse(Com.Sp.groups.pctl$freq < 75, 0,
                                 Com.Sp.groups.pctl$t05)
Com.Sp.groups.pctl[!is.na(Com.Sp.groups.pctl$b05),]$b05 <- 
  ifelse(Com.Sp.groups.pctl[!is.na(Com.Sp.groups.pctl$b05),]$freq
         < 75,0,Com.Sp.groups.pctl[!is.na(Com.Sp.groups.pctl$b05),]$b05)


#Com.Sp.groups.pctl[,c('f05','s05','sc05','t05','f75','s75','sc75','t75','b05','b95')] <- lapply(Com.Sp.groups.pctl[,c('f05','s05','sc05','t05','f75','s75','sc75','t75','b05','b95')], FUN=roundF)

#----
#combine results
Group.Summary <- merge(Com.Sp.groups.pctl, Com.Sp.groups.hts, by=c('cluster','Species'), all.x = T)
Group.Summary <- merge(listspp[,c('Taxon', 'Symbol','Form','Nativity')], Group.Summary, by.x = 'Taxon', by.y = 'Species', all.y = TRUE)
Group.Summary <- merge(Group.Summary, substitutesymbols, by.x = 'Taxon', by.y='taxon', all.x = TRUE)
Group.Summary <- merge(List_Habits[,c('Form', 'ESIS.Group')], Group.Summary, by = 'Form', all.y = TRUE)
#----
#replacement values
Group.Summary[is.na(Group.Summary$Symbol) & !is.na(Group.Summary$SYMBOL2),]$Symbol <-
  Group.Summary[is.na(Group.Summary$Symbol) & !is.na(Group.Summary$SYMBOL2),]$SYMBOL2

Group.Summary[is.na(Group.Summary$ESIS.Group) & !is.na(Group.Summary$form2),]$ESIS.Group <-
  Group.Summary[is.na(Group.Summary$ESIS.Group) & !is.na(Group.Summary$form2),]$form2

Group.Summary[is.na(Group.Summary$Nativity) & !is.na(Group.Summary$type2),]$Nativity <-
  Group.Summary[is.na(Group.Summary$Nativity) & !is.na(Group.Summary$type2),]$type2



#----
Group.Summary <- Group.Summary[order(Group.Summary$cluster, Group.Summary$Taxon),c( 'cluster','Taxon', 'Symbol', 'ESIS.Group','Nativity', 'f05', 'f75', 's05', 's75', 'sc05', 'sc75', 't05', 't75', 'Fmin', 'Fmax', 'Smin', 'Smax', 'SCmin', 'SCmax', 'Tmin', 'Tmax','Dmin','Dmax','b05','b95')]

colnames(Group.Summary) <- c( 'cluster','Taxon', 'Plant.Symbol', 'Plant.Type','Nativity', 'f05', 'f75', 's05', 's75', 'sc05', 'sc75', 't05', 't75', 'Fmin', 'Fmax', 'Smin', 'Smax', 'SCmin', 'SCmax', 'Tmin', 'Tmax','Dmin','Dmax','b05','b95')

Group.Summary$maxHt <- pmax(Group.Summary$Fmax,Group.Summary$Smax,Group.Summary$SCmax,Group.Summary$Tmax, na.rm = TRUE)
Group.Summary <- merge(Group.Summary, Plant_Heights[,c('Scientific.Name', 'Ht_m')], by.x = 'Taxon', by.y='Scientific.Name', all.x=TRUE)
Group.Summary[is.na(Group.Summary$maxHt),]$maxHt <- Group.Summary[is.na(Group.Summary$maxHt),]$Ht_m
Group.Summary[Group.Summary$Plant.Type %in% 'Vine/Liana',]$maxHt <- NA


Forest.Overstory <- Group.Summary[Group.Summary$t75 >0,c( 'cluster','Taxon',  'Plant.Symbol', 'Plant.Type','Nativity', 't05', 't75','Tmin', 'Tmax','Dmin','Dmax','b05','b95', 'maxHt')]
colnames(Forest.Overstory) <- c('cluster','Taxon',  'Plant.Symbol', 'Plant.Type','Nativity', 'Cover.Low', 'Cover.High','Canopy.Bottom', 'Canopy.Top','Diam.Low','Diam.High','BA.Low','BA.High', 'maxHt')
Forest.Overstory[is.na(Forest.Overstory$Canopy.Bottom & Forest.Overstory$Plant.Type %in% c('Shrub/Subshrub','Vine/Liana')),]$
  Canopy.Bottom <- 5
Forest.Overstory[is.na(Forest.Overstory$Canopy.Bottom) & Forest.Overstory$Plant.Type %in% c('Tree'),]$
  Canopy.Bottom <- 10
Forest.Overstory[is.na(Forest.Overstory$Canopy.Top) & Forest.Overstory$Plant.Type %in% c('Shrub/Subshrub','Vine/Liana'),]$
  Canopy.Top <- 20
Forest.Overstory[is.na(Forest.Overstory$Canopy.Top) & Forest.Overstory$Plant.Type %in% c('Tree'),]$
  Canopy.Top <- 25
Forest.Overstory[is.na(Forest.Overstory$Canopy.Bottom),]$
  Canopy.Bottom <- 15

Forest.Overstory.sub<- Group.Summary[Group.Summary$sc75 >0,c( 'cluster','Taxon', 'Plant.Symbol', 'Plant.Type','Nativity','sc05', 'sc75','SCmin', 'SCmax', 'Dmin','Dmax','b05','b95', 'maxHt')]
Forest.Overstory.sub[,c('Dmin','Dmax','b05','b95')]<- NA
colnames(Forest.Overstory.sub) <- c('cluster','Taxon', 'Plant.Symbol', 'Plant.Type','Nativity', 'Cover.Low', 'Cover.High','Canopy.Bottom', 'Canopy.Top','Diam.Low','Diam.High','BA.Low','BA.High','maxHt')
Forest.Overstory.sub[is.na(Forest.Overstory.sub$Canopy.Bottom & Forest.Overstory.sub$Plant.Type %in% c('Shrub/Subshrub','Vine/Liana')),]$
  Canopy.Bottom <- 2
Forest.Overstory.sub[is.na(Forest.Overstory.sub$Canopy.Bottom & Forest.Overstory.sub$Plant.Type %in% c('Tree')),]$
  Canopy.Bottom <- 5
Forest.Overstory.sub[is.na(Forest.Overstory.sub$Canopy.Top) & Forest.Overstory.sub$Plant.Type %in% c('Shrub/Subshrub'),]$
  Canopy.Top <- 10
Forest.Overstory.sub[is.na(Forest.Overstory.sub$Canopy.Top) & Forest.Overstory.sub$Plant.Type %in% c('Tree','Vine/Liana'),]$
  Canopy.Top <- 15
Forest.Overstory.sub[is.na(Forest.Overstory.sub$Canopy.Bottom),]$
  Canopy.Bottom <- 5

Forest.Overstory <- rbind(Forest.Overstory, Forest.Overstory.sub)

Forest.Overstory$Canopy.Top <- ifelse(!is.na(Forest.Overstory$maxHt),
                                      pmin(Forest.Overstory$Canopy.Top,
                                           Forest.Overstory$maxHt),
                                      Forest.Overstory$Canopy.Top)
Forest.Overstory$Canopy.Bottom <- ifelse(!is.na(Forest.Overstory$maxHt),
                                      pmin(Forest.Overstory$Canopy.Bottom,
                                           Forest.Overstory$maxHt),
                                      Forest.Overstory$Canopy.Bottom)


Forest.Understory <- Group.Summary[Group.Summary$s75 >0,c( 'cluster','Taxon', 'Plant.Symbol', 'Plant.Type','Nativity', 's05', 's75','Smin', 'Smax', 'maxHt')]
colnames(Forest.Understory) <- c('cluster','Taxon','Plant.Symbol', 'Plant.Type','Nativity', 'Cover.Low', 'Cover.High','Canopy.Bottom', 'Canopy.Top', 'maxHt')
Forest.Understory[is.na(Forest.Understory$Canopy.Bottom) & Forest.Understory$Plant.Type %in% c('Shrub/Subshrub'),]$
  Canopy.Bottom <- 0.5
Forest.Understory[is.na(Forest.Understory$Canopy.Bottom) & Forest.Understory$Plant.Type %in% c('Tree','Vine/Liana'),]$
  Canopy.Bottom <- 1
Forest.Understory[is.na(Forest.Understory$Canopy.Top) & Forest.Understory$Plant.Type %in% c('Shrub/Subshrub'),]$
  Canopy.Top <- 2
Forest.Understory[is.na(Forest.Understory$Canopy.Top) & Forest.Understory$Plant.Type %in% c('Tree','Vine/Liana'),]$
  Canopy.Top <- 5
Forest.Understory[is.na(Forest.Understory$Canopy.Bottom),]$
  Canopy.Bottom <- 0.5


Forest.Understory.sub <- Group.Summary[Group.Summary$f75 >0,c( 'cluster','Taxon','Plant.Symbol', 'Plant.Type','Nativity', 'f05', 'f75','Fmin', 'Fmax', 'maxHt')]
colnames(Forest.Understory.sub) <- c('cluster','Taxon','Plant.Symbol', 'Plant.Type','Nativity', 'Cover.Low', 'Cover.High','Canopy.Bottom', 'Canopy.Top', 'maxHt')
Forest.Understory.sub[is.na(Forest.Understory.sub$Canopy.Bottom) & Forest.Understory.sub$Plant.Type %in% c('Nonvascular'),]$Canopy.Bottom <- 0.0
Forest.Understory.sub[is.na(Forest.Understory.sub$Canopy.Bottom) & Forest.Understory.sub$Plant.Type %in% c('Tree','Vine/Liana','Grass/grass-like (Graminoids)','Fern/fern ally','Forb/Herb', 'Nonvascular'),]$Canopy.Bottom <- 0.1
Forest.Understory.sub[is.na(Forest.Understory.sub$Canopy.Top) & Forest.Understory.sub$Plant.Type %in% c('Nonvascular'),]$Canopy.Top <- 0.05
Forest.Understory.sub[is.na(Forest.Understory.sub$Canopy.Top) & Forest.Understory.sub$Plant.Type %in% c('Grass/grass-like (Graminoids)','Fern/fern ally','Forb/Herb'),]$Canopy.Top <- 0.6
Forest.Understory.sub[is.na(Forest.Understory.sub$Canopy.Top) & Forest.Understory.sub$Plant.Type %in% c('Shrub/Subshrub'),]$Canopy.Top <- 0.3
Forest.Understory.sub[is.na(Forest.Understory.sub$Canopy.Top) & Forest.Understory.sub$Plant.Type %in% c('Tree','Vine/Liana'),]$Canopy.Top <- 0.5
Forest.Understory.sub[is.na(Forest.Understory.sub$Canopy.Bottom),]$
  Canopy.Bottom <- 0

Forest.Understory <- rbind(Forest.Understory, Forest.Understory.sub)

Forest.Understory$Canopy.Top <- ifelse(!is.na(Forest.Understory$maxHt),
                                      pmin(Forest.Understory$Canopy.Top,
                                           Forest.Understory$maxHt),
                                      Forest.Understory$Canopy.Top)
Forest.Understory$Canopy.Bottom <- ifelse(!is.na(Forest.Understory$maxHt),
                                         pmin(Forest.Understory$Canopy.Bottom,
                                              Forest.Understory$maxHt),
                                         Forest.Understory$Canopy.Bottom)

Forest.Overstory[,c('Cover.Low', 'Cover.High', 'Canopy.Bottom', 'Canopy.Top', 'Diam.Low', 'Diam.High', 'BA.Low', 'BA.High')] <- lapply(Forest.Overstory[,c('Cover.Low', 'Cover.High', 'Canopy.Bottom', 'Canopy.Top', 'Diam.Low', 'Diam.High', 'BA.Low', 'BA.High')], FUN=roundF)

Forest.Understory[,c('Cover.Low', 'Cover.High', 'Canopy.Bottom', 'Canopy.Top')] <- lapply(Forest.Understory[,c('Cover.Low', 'Cover.High', 'Canopy.Bottom', 'Canopy.Top')], FUN=roundF)

Sort_Habits <- read.delim("data/Sort_Habits.txt", na.strings = '', stringsAsFactors = FALSE)

Forest.Overstory1 <- merge(Forest.Overstory, Sort_Habits, by='Plant.Type', all.x = TRUE)
Forest.Understory1 <- merge(Forest.Understory, Sort_Habits, by='Plant.Type', all.x = TRUE)
Forest.Overstory1 <- Forest.Overstory1[order(
  Forest.Overstory1$cluster,
  Forest.Overstory1$order,
  -Forest.Overstory1$Cover.Low,
  -Forest.Overstory1$Cover.High,
  -Forest.Overstory1$Cover.High,
  Forest.Overstory1$Taxon
),c('cluster',	'Taxon', 'order', 'Plant.Symbol',	'Plant.Type',	'Nativity',	'Cover.Low',	'Cover.High',	'Canopy.Bottom',	'Canopy.Top',	'Diam.Low',	'Diam.High',	'BA.Low',	'BA.High')]

Forest.Understory1 <- Forest.Understory1[order(
  Forest.Understory1$cluster,
  Forest.Understory1$order,
  -Forest.Understory1$Cover.Low,
  -Forest.Understory1$Cover.High,
  -Forest.Understory1$Cover.High,
  Forest.Understory1$Taxon
),c('cluster',	'Taxon', 'order', 'Plant.Symbol',	'Plant.Type',	'Nativity',	'Cover.Low',	'Cover.High',	'Canopy.Bottom',	'Canopy.Top')]


write.csv(Forest.Overstory1, 'output/Forest.Overstory.csv', row.names = FALSE, na = "")
write.csv(Forest.Understory1, 'output/Forest.Understory.csv', row.names = FALSE, na = "")


