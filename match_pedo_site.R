library(stringr)
library(BiodiversityR)
library(cluster)
library(ape)
library(dendextend)
library(dplyr)
#----
NASISPEDONS <- read.delim("data/NASISPEDONS.txt")
List_Habits <- read.delim("data/List_Habits.txt", na.strings = '')
#Habit_Symbols <- read.delim("data/Habit_Symbols.txt", encoding = 'UTF-8', na.strings = '')
#FloraNorthAmerica <- read.delim("data/FloraNorthAmerica.txt", encoding = 'UTF-8', na.strings = '')
#colnames(FloraNorthAmerica)[1] <- c('Taxon') #This fixes the name of the first field, because Notepad adds "byte order mark" to begining of file if saving as UTF-8, which import as junk characters.
#write.table(FloraNorthAmerica, file = 'data/FloraNorthAmerica2.txt',  sep = "\t", quote = FALSE) 
obsspp <- read.delim("data/Observed_Species.txt", encoding = 'UTF-8', na.strings = '')
obsspp <- subset(obsspp, !substr(AcTaxon,1,1) %in% '-'& !AcTaxon %in% '' & !is.na(Habit))
obsspp <- merge(obsspp, List_Habits[,c('Form','Simple')], by.x = 'Habit', by.y = 'Form', all.x = TRUE)
obs <- read.delim("data/Sites.txt")
obsspp <- merge(obs[,c('Observation_ID','Observation_Label')],obsspp, by='Observation_ID')
obsspp <- subset(obsspp, Field+Shrub+Subcanopy+Tree > 0)

VEGOBS <- read.delim("data/VEGOBS.txt")
VEGOBS$pedon <- ""
VEGOBS$taxonname <- ""
VEGOBS$taxonclass <- ""
VEGOBS$pedondist <- NA
n <- nrow(VEGOBS)
for (i in 1:n){
NASISPEDONS$distance <- (((VEGOBS[i,]$Latitude - NASISPEDONS$Std.Latitude)/360*40041.47*1000)^2 +
  ((VEGOBS[i,]$Longitude - NASISPEDONS$Std.Longitude)/360*40041.47*1000*cos(VEGOBS[i,]$Latitude/2/360*2*3.141592))^2)^0.5
mindist <- min(NASISPEDONS$distance, na.rm = TRUE)

VEGOBS[i,]$pedondist <- round(mindist,1)
VEGOBS[i,]$pedon <- as.character(subset(NASISPEDONS, NASISPEDONS$distance == mindist)[1,]$User.Site.ID)
VEGOBS[i,]$taxonname <- as.character(subset(NASISPEDONS, NASISPEDONS$distance == mindist)[1,]$Current.Taxon.Name)
VEGOBS[i,]$taxonclass <- as.character(subset(NASISPEDONS, NASISPEDONS$distance == mindist)[1,]$Current.Taxonomic.Class)

}

VEGOBS[VEGOBS$pedondist > 50,]$taxonname <- ""
VEGOBS[VEGOBS$pedondist > 50,]$taxonclass <- ""
VEGOBS[VEGOBS$pedondist > 50,]$pedon <- ""

mu <- readRDS(file='data/mu.RDS')

VEGOBS_mukeys <- read.delim("data/VEGOBS_mukeys.txt")
VEGOBS_soilnames <- merge(VEGOBS_mukeys[,c('Observatio','RASTERVALU')], mu[,c('lmapunitiid', 'muname')], by.x='RASTERVALU', by.y= 'lmapunitiid')
VEGOBS <- merge(VEGOBS, VEGOBS_soilnames[,c('Observatio', 'muname')], by.x='Observation_Label', by.y= 'Observatio')
VEGOBS$Soil <- VEGOBS$taxonname
VEGOBS[VEGOBS$Soil %in% '',]$Soil <- str_split_fixed(VEGOBS[VEGOBS$Soil %in% '',]$muname, " ",2)[,1]

#----
#observed species

Com.Sp.sum<-aggregate(obsspp[,c('Field', 'Shrub', 'Subcanopy', 'Tree')], by=list(obsspp$Observation_Label, obsspp$AcTaxon, obsspp$Simple), FUN=sum) #sum within plot
colnames(Com.Sp.sum)<-c('Observation_Label', 'Species', 'Simple', 'Field', 'Shrub', 'Subcanopy', 'Tree') #restore column names

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
#cluster analysis
Com.Sp.mean$sqrttotal <- sqrt(Com.Sp.mean$Total)

plotinputs1 <- makecommunitydataset(Com.Sp.mean, row = 'soilplot', column = 'Species', value = 'sqrttotal', drop = TRUE)

jacdist1 <- as.data.frame(as.matrix(vegdist(plotinputs1,method='jaccard', binary=FALSE, na.rm=T)))

jactree <- agnes(jacdist1, method='average')

w <- 800
h <- 3000
u <- 12
ngroups <- 16
groups <- cutree(jactree, k = ngroups)

dend1 <- color_branches(as.hclust(jactree), k = ngroups)
dend1 <- color_labels(dend1, k = ngroups)

png(filename="output/jactree_flora.png",width = w, height = h, units = "px", pointsize = u)
par(mar = c(2,0,1,13))
plot(dend1, horiz = TRUE, main='floristic simularity - jaccard metric', font=1, cex=0.85)
rect.dendrogram(dend1, k = ngroups, horiz = TRUE)
dev.off()

#----
Com.Sp.Agg$sqrttotal <- Com.Sp.Agg$Total^0.5
plotinputs2 <- makecommunitydataset(Com.Sp.Agg, row = 'soilplot', column = 'Simple', value = 'Total', drop = TRUE)

jacdist2 <- as.data.frame(as.matrix(vegdist(plotinputs2,method='jaccard', binary=FALSE, na.rm=T)))

jactree <- agnes(jacdist2, method='average')

ngroups <- 8
groups <- cutree(jactree, k = ngroups)

dend1 <- color_branches(as.hclust(jactree), k = ngroups)
dend1 <- color_labels(dend1, k = ngroups)

w <- 800
h <- 3000
u <- 12
png(filename="output/jactree_simple.png",width = w, height = h, units = "px", pointsize = u)

par(mar = c(2,0,1,13))
plot(dend1, horiz = TRUE, main='floristic simularity - jaccard metric', font=1, cex=0.85)
rect.dendrogram(dend1, k = ngroups, horiz = TRUE)
dev.off()
#----
#jacdist <- (jacdist1*2+jacdist2*1)/3
plotinputs <- cbind(plotinputs1, plotinputs2)
jacdist <- as.data.frame(as.matrix(vegdist(plotinputs1, method='bray', binary=FALSE, na.rm=T)))
maxdist <- max(na.omit(jacdist))
jactree <- agnes(jacdist, method='average')
#jactree <- diana(jacdist)

ngroups <- 16
groups <- cutree(jactree, k = ngroups)
soilplot <- names(groups)
clust <- unname(groups)
groupdf <- as.data.frame(cbind(soilplot, clust))
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

w <- 800
h <- 3000
u <- 12
png(filename="output/jactree_comb.png",width = w, height = h, units = "px", pointsize = u)

par(mar = c(2,0,1,13))
plot(dend1, horiz = TRUE, main='floristic simularity - bray metric', font=1, cex=0.85)
rect.dendrogram(dend1, k = ngroups, horiz = TRUE)
dev.off()
#----
#group dominant and indicator species

soilplot <- names(groups)
clust <- unname(groups)
groupdf <- as.data.frame(cbind(soilplot, clust))


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
                       (((freqranks <= 3 | affranks <= 1) & stratum %in% c(1,2,3))|(subranks <= 1 & stratum == 3))
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
Com.Structure[order(as.numeric(as.character(Com.Structure$cluster))),c("cluster", "association", "Structure")]

nrow(unique(Com.Sp.groups[,c('Soil','clust')]))/length(unique(Com.Sp.groups$Soil))
#comparison of 16 clusters different matrix, distance, and tree building
clusterpersoilbrayflorahabitmatrixcomb <- 1.834862
clusterpersoileuclidflorahabitmatrixcomb <- 1.853211
clusterpersoilbraybinaryflorahabitmatrixcomb <- 1.862385
clusterpersoiljacflorahabitmatrixcomb <- 1.917431
clusterpersoiljacfloraonly <- 1.990826
clusterpersoiljacflorahabitdist2to1 <- 2.009174
clusterpersoiljachabitonly <- 2.036697
clusterpersoiljacdiana <- 2.110092
clusterpersoilbraynotsqrtflorahabitmatrixcomb <- 2.119266
clusterpersoiljacward <-2.275229
clusterpersoilbraynotsqrtfloraonly <- 2.073394
clusterpersoiljacnotsqrtfloraonly <- 2.018349
clusterpersoilbraysqrtbothflorahabitcomb <- 1.954128
clusterpersoilbrayfloraonly <- 1.963303
#optimum is with sqrt sp matrix with non-sqrt habits combined matrix then bray distance. But flora only matrix may have more cohesive community composition, but this defailts to mostly forests.
write.table(VEGOBS, 'output/VEGOBS-export.txt', row.names = FALSE, sep = "\t")
write.dbf(VEGOBS[,c(1,3:ncol(VEGOBS))], 'output/VEGOBS.dbf')

