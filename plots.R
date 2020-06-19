#spatial packages
require(sp)
require(raster)
require(rgdal)
require(rms)
require(RODBC)
#cluster analysis packages
library(vegan)
library(cluster)
library(ape)
#color pallette
library(RColorBrewer)
library(colorspace)
# Create a scatter plot matrix
library(ggplot2)
#calculate percentiles
library(plyr)
#classification tree packages
library(randomForest)
library(rpart)
library(rpart.plot)
#Discriminate Analysis
library(MASS)
#Composition in the form of boxes
library(treemap)


#0 import data-------------------------
obsspp <- read.delim("Observed_Species.txt")
synonymy <- read.delim("List_Species2011.txt")
growthhabit <- read.delim("List_Habits.txt")
clusters <- read.delim("Community_Clusters.txt")
obsspp$AcTaxon<- NULL #removing this column because we will create it as a demo
obsspp<- obsspp[obsspp$Module!='BA',] #removing this column because we will create it as a demo
growthhabit<-growthhabit[!(is.na(growthhabit$Form) | growthhabit$Form==""),] #remove blanks
#1 Synonymize species list------------------------------------------
#merge taxon parts into single field

obsspp$taxon<-ifelse(is.na(obsspp$var)|obsspp$var=="",
                     ifelse(is.na(obsspp$subsp)|obsspp$subsp=="",
                            paste(obsspp$Genus," ",obsspp$specific_epithet, sep=""),
                            paste(obsspp$Genus," ",obsspp$specific_epithet," ssp. ", obsspp$subsp, sep="")),
                     paste(obsspp$Genus," ",obsspp$specific_epithet," var. ", obsspp$var, sep=""))
str(obsspp)
#merge with synonym database, we will use the accepted binomial field for consistency sake (not everyone will use infrataxa)
#In this case we will keep non matches in the first table
#but first we got to fix case, by capitalizing only first letter
obsspp$taxon<-tolower(obsspp$taxon)
obsspp$taxon<-paste(toupper(substr(obsspp$taxon,1,1)), substr(obsspp$taxon, 2,nchar(obsspp$taxon)), sep="")

obsspp<-merge(obsspp,synonymy[,c("Taxon","AcTaxon","AcBinomial")], by.x="taxon", by.y="Taxon", all.x = T)
#move unknown taxa to the binomial field
obsspp$Species<-ifelse(is.na(obsspp$AcBinomial), obsspp$taxon,as.character(obsspp$AcBinomial)) #had to use as.character since the field populated as.factor level numbers
                    
#2 Add ESIS growth habits------------------------
obsspp<-merge(obsspp,growthhabit[,c("Form","Habitname","ESIS.Group")], by.x="Habit", by.y="Form")

#3 Add community Classification based on a the Sandy Lakeplain ESD community phases--------------------
obsspp<-merge(obsspp,clusters[,c("Plot","Community")], by.x="Observation_ID", by.y="Plot")

#4 Get single line per plot by averaging subplots if any exist (a maximum of 4)
Com.Sp.sum<-aggregate(obsspp[,c('Field', 'Shrub', 'Subcanopy', 'Tree')], by=list(obsspp$Observation_ID, obsspp$Community, obsspp$ESIS.Group, obsspp$Species), FUN=sum) #sum within plot
colnames(Com.Sp.sum)<-c('Observation_ID', 'Community', 'ESIS.Group', 'Species','Field', 'Shrub', 'Subcanopy', 'Tree') #restore column names

Com.Sp.freq<-aggregate(obsspp[,c('Species')], by=list(obsspp$Observation_ID, obsspp$Community, obsspp$ESIS.Group, obsspp$Species), FUN=length) #frequency within plot
colnames(Com.Sp.freq)<- c('Observation_ID', 'Community', 'ESIS.Group', 'Species', 'freq')

Com.max.freq<-aggregate(Com.Sp.freq[,c('freq')], by=list(Com.Sp.freq$Observation_ID, Com.Sp.freq$Community), FUN=max) #freq within plot
colnames(Com.max.freq)<- c('Observation_ID', 'Community', 'mfreq')
Com.max.freq$mfreq<-ifelse(Com.max.freq$mfreq>4,4,Com.max.freq$mfreq)#effectively ensureing values do not exceed 4. Species listed 5 times might occur if surveyer was unaware of species already counted in subplots, but this only adds a trace amount.

Com.Sp.mean<-merge(Com.Sp.sum, Com.max.freq[,c("Observation_ID","mfreq")], by="Observation_ID")
Com.Sp.mean$Field<-Com.Sp.mean$Field/Com.Sp.mean$mfreq
Com.Sp.mean$Shrub<-Com.Sp.mean$Shrub/Com.Sp.mean$mfreq
Com.Sp.mean$Subcanopy<-Com.Sp.mean$Subcanopy/Com.Sp.mean$mfreq
Com.Sp.mean$Tree<-Com.Sp.mean$Tree/Com.Sp.mean$mfreq
Com.Sp.ht<-aggregate(obsspp[,c('Fmin', 'Fmax', 'Smin', 'Smax', 'SCmin', 'SCmax', 'Tmin', 'Tmax')], by=list(obsspp$Observation_ID, obsspp$Community, obsspp$ESIS.Group, obsspp$Species), FUN=mean)
colnames(Com.Sp.ht)<-c('Observation_ID', 'Community', 'ESIS.Group', 'Species','Fmin', 'Fmax', 'Smin', 'Smax', 'SCmin', 'SCmax', 'Tmin', 'Tmax')
Com.Sp.mean<-merge(Com.Sp.mean, Com.Sp.ht, by=c('Observation_ID', 'Community', 'ESIS.Group', 'Species'))

#5 Determine frequency of species occuring in a plot for each community type------------------
Com.Sp.freq2<-aggregate(Com.Sp.mean[,c('Species')], by=list(Com.Sp.mean$Community, Com.Sp.mean$ESIS.Group, Com.Sp.mean$Species), FUN=length) #freq within community
colnames(Com.Sp.freq2)<-c('Community', 'ESIS.Group', 'Species', 'freq')
Countplots1<-unique(Com.Sp.mean[,c('Observation_ID', 'Community')])
Countplots2<-aggregate(Countplots1[,c('Community')], by=list(Countplots1$Community), FUN=length) #plots per community
colnames(Countplots2)<-c('Community', 'total')
Com.Sp.freq2<-merge(Com.Sp.freq2,Countplots2,by="Community")
Com.Sp.freq2$pfreq<-Com.Sp.freq2$freq/Com.Sp.freq2$total #This is the percent frequency of species within a given community type.
Com.Sp.mean<-merge(Com.Sp.mean, Com.Sp.freq2[,c('Community', 'Species', 'pfreq')], by=c('Community','Species'))

#6 Calculate upper and lower values for cover by stratum using percentiles for cover and means for canopy heights-----------------
Com.Sp.ref<-ddply(Com.Sp.mean, c("Community","Species"), summarise, Field05 = quantile(Field, .05), Field95 = quantile(Field, .95), Shrub05 = quantile(Shrub, .05), Shrub95 = quantile(Shrub, .95), Subcanopy05 = quantile(Subcanopy, .05), Subcanopy95 = quantile(Subcanopy, .95), Tree05 = quantile(Tree, .05), Tree95 = quantile(Tree, .95), Fmin = mean(Fmin), Fmax = mean(Fmax), Smin = mean(Smin), Smax = mean(Smax), SCmin = mean(SCmin), SCmax = mean(SCmax), Tmin = mean(Tmin), Tmax = mean(Tmax), pfreq = mean(pfreq))
#zero out minimum cover for species that are less frequent than 95% frequent 
Com.Sp.ref$Field05<-ifelse(Com.Sp.ref$pfreq<0.95,0,Com.Sp.ref$Field05)
Com.Sp.ref$Shrub05<-ifelse(Com.Sp.ref$pfreq<0.95,0,Com.Sp.ref$Shrub05)
Com.Sp.ref$Subcanopy05<-ifelse(Com.Sp.ref$pfreq<0.95,0,Com.Sp.ref$Subcanopy05)
Com.Sp.ref$Tree05<-ifelse(Com.Sp.ref$pfreq<0.95,0,Com.Sp.ref$Tree05)

#7 Rounding to avoid unnecessary appearence of precision in reference values-----------------

roundF<-function(p){
  p<-ifelse(p<0.5, floor(p/0.1+0.5)*0.1,ifelse(p<2, floor(p/0.5+0.5)*0.5, ifelse(p<5, floor(p+0.5),floor(p/5+0.5)*5)))
}

Com.Sp.round<-Com.Sp.ref
Com.Sp.round$Field05<-roundF(Com.Sp.round$Field05)
Com.Sp.round$Field95<-roundF(Com.Sp.round$Field95)
Com.Sp.round$Shrub05<-roundF(Com.Sp.round$Shrub05)
Com.Sp.round$Shrub95<-roundF(Com.Sp.round$Shrub95)
Com.Sp.round$Subcanopy05<-roundF(Com.Sp.round$Subcanopy05)
Com.Sp.round$Subcanopy95<-roundF(Com.Sp.round$Subcanopy95)
Com.Sp.round$Tree05<-roundF(Com.Sp.round$Tree05)
Com.Sp.round$Tree95<-roundF(Com.Sp.round$Tree95)
Com.Sp.round$Fmin<-roundF(Com.Sp.round$Fmin)
Com.Sp.round$Fmax<-roundF(Com.Sp.round$Fmax)
Com.Sp.round$Smin<-roundF(Com.Sp.round$Smin)
Com.Sp.round$Smax<-roundF(Com.Sp.round$Smax)
Com.Sp.round$SCmin<-roundF(Com.Sp.round$SCmin)
Com.Sp.round$SCmax<-roundF(Com.Sp.round$SCmax)
Com.Sp.round$Tmin<-roundF(Com.Sp.round$Tmin)
Com.Sp.round$Tmax<-roundF(Com.Sp.round$Tmax)

#Last minute conversion from meters to feet depending default units required in ESD database
Com.Sp.round$Fminft<-round(Com.Sp.round$Fmin/0.3048,1)
Com.Sp.round$Fmaxft<-round(Com.Sp.round$Fmax/0.3048,1)
Com.Sp.round$Sminft<-round(Com.Sp.round$Smin/0.3048,1)
Com.Sp.round$Smaxft<-round(Com.Sp.round$Smax/0.3048,1)
Com.Sp.round$SCminft<-round(Com.Sp.round$SCmin/0.3048,1)
Com.Sp.round$SCmaxft<-round(Com.Sp.round$SCmax/0.3048,1)
Com.Sp.round$Tminft<-round(Com.Sp.round$Tmin/0.3048,1)
Com.Sp.round$Tmaxft<-round(Com.Sp.round$Tmax/0.3048,1)

#8 Treemap plot species by community---------------------
Com.Sp.mean2<-ddply(Com.Sp.mean, c("Community",'ESIS.Group',"Species"), summarise, Field = mean(Field), Shrub = mean(Shrub), Subcanopy = mean(Subcanopy), Tree = mean(Tree), pfreq = mean(pfreq))
#sum mean covers and multiply by percent frequency
Com.Sp.mean2$Total<-(Com.Sp.mean2$Field+Com.Sp.mean2$Shrub+Com.Sp.mean2$Subcanopy+Com.Sp.mean2$Tree)*Com.Sp.mean2$pfreq

treeMapPlot<-treemap(Com.Sp.mean2, #Your data frame object
                     index=c("Community", "Species"),  #A list of your categorical variables
                     vSize = "Total",  #This is your quantitative variable
                     vColor= "Species", #This is a categorical variable
                     type="categorical", #Type sets the organization and color scheme of your treemap
                     palette = c("#CCCCCC", "#006400", "#00EE00", "#B23AEE", "#00868B", "#CDCD00", "#40E0D0", "#228B22", "#DB7093", "#EE799F", "#4169E1", "#8470FF", "#CD661D", "#EE0000", "#FFA54F", "#FFB6C1", "#FFE7BA", "#98F5FF", "#32CD32", "#79CDCD", "#458B00", "#9ACD32", "#90EE90"),
                     title="Frequency of Species by Community", #Customize your title
                     fontsize.labels=c(10,8), # size of labels. Give the size per level of                                                       aggregation: size for group, size for subgroup, sub                                                -subgroups...
                     fontcolor.labels=c("black","black","black"),    # Color of labels
                     fontface.labels=c(2,1,1),   # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
                     bg.labels=c("#CCCCCCDC"),              # Background color of labels
                     align.labels=list(
                       c("left", "top"), 
                       c("right", "bottom"), 
                       c("center", "center")
                     ), # Where to place labels in the rectangle?
                     overlap.labels=0.5,     # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
                     inflate.labels=F      # If true, labels are bigger when rectangle is bigger.
                     
)

print(treeMapPlot)


#9 plot treemap by growth habit----------------
treeMapPlot2<-treemap(Com.Sp.mean2, #Your data frame object
                     index=c("Community", "ESIS.Group"),  #A list of your categorical variables
                     vSize = "Total",  #This is your quantitative variable
                     vColor= "ESIS.Group", #This is a categorical variable
                     type="categorical", #Type sets the organization and color scheme of your treemap
                     palette = c("#CCCCCC", "#006400",  "#B23AEE", "#CDCD00","#00868B", "#00EE00", "#228B22","#40E0D0",  "#DB7093", "#EE799F", "#4169E1", "#8470FF", "#CD661D", "#EE0000", "#FFA54F", "#FFB6C1", "#FFE7BA", "#98F5FF", "#32CD32", "#79CDCD", "#458B00", "#9ACD32", "#90EE90"),
                     title="Frequency of Growth Habit by Community", #Customize your title
                     fontsize.labels=c(10,8), # size of labels. Give the size per level of                                                       aggregation: size for group, size for subgroup, sub                                                -subgroups...
                     fontcolor.labels=c("black","black","black"),    # Color of labels
                     fontface.labels=c(2,1,1),   # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
                     bg.labels=c("#CCCCCCDC"),              # Background color of labels
                     align.labels=list(
                       c("left", "top"), 
                       c("right", "bottom"), 
                       c("center", "center")
                     ), # Where to place labels in the rectangle?
                     overlap.labels=0.5,     # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
                     inflate.labels=F      # If true, labels are bigger when rectangle is bigger.
                     
)
print(treeMapPlot2)
