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
Com.Sp.mean$total <- ((Com.Sp.mean$Field + Com.Sp.mean$Shrub + Com.Sp.mean$Subcanopy + Com.Sp.mean$Tree)/2)^0.5 #average overstory and understory and apply sqroot transform
Com.Sp.mean <- merge(Com.Sp.mean, VEGOBS[,c('Observation_Label', 'Soil')], by='Observation_Label')
Com.Sp.mean <-subset(Com.Sp.mean, !substr(Species,1,1) %in% '-'& !Species %in% '')
Com.Sp.mean$soilplot <- paste(Com.Sp.mean$Soil , Com.Sp.mean$Observation_Label)

#----
#need formula for aggregating within strata... 
Com.Sp.preagg <- Com.Sp.mean

Com.Sp.preagg$Field <- ifelse(Com.Sp.preagg$Field > 100,100,Com.Sp.preagg$Field)
Com.Sp.preagg$Shrub <- ifelse(Com.Sp.preagg$Shrub > 100,100,Com.Sp.preagg$Shrub)
Com.Sp.preagg$Subcanopy <- ifelse(Com.Sp.preagg$Subcanopy > 100,100,Com.Sp.preagg$Subcanopy)
Com.Sp.preagg$Tree <- ifelse(Com.Sp.preagg$Tree > 100,100,Com.Sp.preagg$Tree)
Com.Sp.preagg[,c('Field', 'Shrub', 'Subcanopy', 'Tree')] <- log10(1-(Com.Sp.preagg[,c('Field', 'Shrub', 'Subcanopy', 'Tree')]/100.001))
Com.Sp.Agg <- aggregate(Com.Sp.preagg[,c('Field', 'Shrub', 'Subcanopy', 'Tree')], by=list(Com.Sp.preagg$soilplot, Com.Sp.preagg$Simple), FUN = 'sum')
colnames(Com.Sp.Agg) <- c('soilplot', 'Simple', 'Field', 'Shrub', 'Subcanopy', 'Tree')
Com.Sp.Agg$Total <- Com.Sp.Agg$Field + Com.Sp.Agg$Shrub
Com.Sp.AggCanopy <- subset(Com.Sp.Agg,  Simple %in% c('Deciduous','Evergreen'))
Com.Sp.AggCanopy$Simple <- paste(Com.Sp.AggCanopy$Simple, 'Tree')
Com.Sp.AggCanopy$Total <- Com.Sp.AggCanopy$Subcanopy + Com.Sp.AggCanopy$Tree
Com.Sp.Agg$Simple <-as.character(Com.Sp.Agg$Simple)
Com.Sp.Agg[Com.Sp.Agg$Simple %in% c('Evergreen', 'Deciduous'),]$Simple <- 
  paste(Com.Sp.Agg[Com.Sp.Agg$Simple %in% c('Evergreen', 'Deciduous'),]$Simple, 'Shrub')
Com.Sp.Agg <- rbind(Com.Sp.Agg, Com.Sp.AggCanopy)
rm(Com.Sp.preagg, Com.Sp.AggCanopy)
Com.Sp.Agg$Total <- (10^(Com.Sp.Agg$Total)*-1+1)*100
#----
#cluster analysis


plotinputs1 <- makecommunitydataset(Com.Sp.mean, row = 'soilplot', column = 'Species', value = 'total', drop = TRUE)

jacdist1 <- as.data.frame(as.matrix(vegdist(plotinputs1,method='jaccard', binary=FALSE, na.rm=T)))

jactree <- agnes(jacdist1, method='average')

w <- 800
h <- 3000
u <- 12

groups <- cutree(jactree, k = 4)

dend1 <- color_branches(as.hclust(jactree), k = 4)
dend1 <- color_labels(dend1, k = 4)

png(filename="output/jactree_flora.png",width = w, height = h, units = "px", pointsize = u)
par(mar = c(2,0,1,13))
plot(dend1, horiz = TRUE, main='floristic simularity - jaccard metric', font=1, cex=0.85)
rect.dendrogram(dend1, k =4, horiz = TRUE)
dev.off()

#----
plotinputs2 <- makecommunitydataset(Com.Sp.Agg, row = 'soilplot', column = 'Simple', value = 'Total', drop = TRUE)

jacdist2 <- as.data.frame(as.matrix(vegdist(plotinputs2,method='jaccard', binary=FALSE, na.rm=T)))

jactree <- agnes(jacdist2, method='average')

groups <- cutree(jactree, k = 4)

dend1 <- color_branches(as.hclust(jactree), k = 4)
dend1 <- color_labels(dend1, k = 4)

w <- 800
h <- 3000
u <- 12
png(filename="output/jactree_simple.png",width = w, height = h, units = "px", pointsize = u)

par(mar = c(2,0,1,13))
plot(dend1, horiz = TRUE, main='floristic simularity - jaccard metric', font=1, cex=0.85)
rect.dendrogram(dend1, k =4, horiz = TRUE)
dev.off()
#----
jacdist <- (jacdist1+jacdist2)/2
maxdist <- max(na.omit(jacdist))
jactree <- agnes(jacdist, method='average')

groups <- cutree(jactree, k = 4)

dend1 <- color_branches(as.hclust(jactree), k = 4)
dend1 <- color_labels(dend1, k = 4)

w <- 800
h <- 3000
u <- 12
png(filename="output/jactree_comb.png",width = w, height = h, units = "px", pointsize = u)

par(mar = c(2,0,1,13))
plot(dend1, horiz = TRUE, main='floristic simularity - jaccard metric', font=1, cex=0.85)
rect.dendrogram(dend1, k =4, horiz = TRUE)
dev.off()
#----
#group dominant and indicator species

soilplot <- names(groups)
clust <- unname(groups)
groupdf <- as.data.frame(cbind(soilplot, clust))
Com.Sp.Agg$soilplot <- str_replace_all(Com.Sp.Agg$soilplot, ' ', '.')
Com.Sp.Agg$soilplot <- str_replace_all(Com.Sp.Agg$soilplot, '-', '.')
Com.Sp.Agg$soilplot <- str_replace_all(Com.Sp.Agg$soilplot, ',', '.')
Com.Sp.Agg$soilplot <- str_replace_all(Com.Sp.Agg$soilplot, ':', '.')
Com.Sp.Agg$soilplot <- str_replace_all(Com.Sp.Agg$soilplot, ';', '.')
Com.Sp.Agg$soilplot <- str_replace_all(Com.Sp.Agg$soilplot, '&', '.')

Com.Sp.groups <- merge(groupdf,  Com.Sp.Agg, by='soilplot', all.x=TRUE, all.y = TRUE)
#need to add placeholders for all zeros

unique(Com.Sp.mean$Species)

#----
Com.Sp.grouprank <- 
  Com.Sp.groups %>%
  group_by(clust) %>%
  mutate(my_ranks = order(order(Total, decreasing=TRUE)))

write.table(VEGOBS, 'output/VEGOBS-export.txt', row.names = FALSE, sep = "\t")
write.dbf(VEGOBS[,c(1,3:ncol(VEGOBS))], 'output/VEGOBS.dbf')

