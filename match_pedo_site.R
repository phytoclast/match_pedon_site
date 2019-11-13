library(stringr)
library(BiodiversityR)
library(cluster)
library(ape)
NASISPEDONS <- read.delim("data/NASISPEDONS.txt")
obsspp <- read.delim("data/Observed_Species.txt")
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

Com.Sp.sum<-aggregate(obsspp[,c('Field', 'Shrub', 'Subcanopy', 'Tree')], by=list(obsspp$Observation_Label, obsspp$AcTaxon), FUN=sum) #sum within plot
colnames(Com.Sp.sum)<-c('Observation_Label', 'Species','Field', 'Shrub', 'Subcanopy', 'Tree') #restore column names

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
#cluster analysis


plotinputs <- makecommunitydataset(Com.Sp.mean, row = 'soilplot', column = 'Species', value = 'total', drop = TRUE)

jacdist <- as.data.frame(as.matrix(vegdist(plotinputs,method='jaccard', binary=FALSE, na.rm=T)))

jactree <- agnes(jacdist, method='average')

w <- 800
h <- 3000
u <- 12
png(filename="output/jactree.png",width = w, height = h, units = "px", pointsize = u)
plot(as.phylo(as.hclust(jactree)), main='floristic simularity - jaccard metric',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()

#----


write.table(VEGOBS, 'output/VEGOBS-export.txt', row.names = FALSE, sep = "\t")
write.dbf(VEGOBS[,c(1,3:ncol(VEGOBS))], 'output/VEGOBS.dbf')

