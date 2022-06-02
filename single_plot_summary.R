# packages ----
library(stringr)
library(BiodiversityR)
library(cluster)
library(ape)
library(dendextend)
library(dplyr)
library(dynamicTreeCut)
library(proxy)
library(Hmisc)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# functions ----
substitutesymbols <- as.data.frame(cbind(SYMBOL2=c('2AG','TRBO2','SMLA3','CAREX'), form2=c('Nonvascular','Forb/Herb','Forb/Herb','Grass/grass-like (Graminoids)'), type2=c('Native','Native','Native','Native'), taxon=c('Chara','Lysimachia borealis','Smilax lasioneuron','Carex [ovales]')), stringsAsFactors = FALSE)



  

roundF<-function(p){
  p<-ifelse(p<0.5, floor(p/0.1+0.5)*0.1,ifelse(p<5, floor(p/0.5+0.5)*0.5, ifelse(p<15, floor(p+0.5),floor(p/5+0.5)*5)))
}
tofoliar <- function(c){
  f <- (0.41*c/100 + 0.41*(c/100)^2)*100
  return(f)
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


# data ----
plots <- c(
  'GRR.GJS.2022.1',
  'GRR.GJS.2022.2',
  'GRR.GJS.2022.3',
  'GRR.GJS.2022.4',
  'GRR.GJS.2022.5',
  'GRR.GJS.2022.6',
  'GRR.GJS.2022.7',
  'GRR.GJS.2022.8',
  'GRR.GJS.2022.9',
  'GRR.GJS.2022.10',
  'GRR.GJS.2022.11',
  'GRR.GJS.2022.12'
)
# process data ----
List_Habits <- read.delim("data/List_Habits.txt", na.strings = '', stringsAsFactors = FALSE)
Plant_Heights <- read.delim("data/Plant_Heights.txt", na.strings = '', stringsAsFactors = FALSE, encoding = 'UTF-8')
List_Habits[List_Habits$ESIS.Group %in% 'Grass/grass-like',]$ESIS.Group <- 'Grass/grass-like (Graminoids)'
listspp <- read.delim("data/List_Species2011.txt", encoding = 'UTF-8', na.strings = '', stringsAsFactors = FALSE)
#listspp <- readRDS("data/listspp.RDS")
#fix for trees and shrub mismatch
SBD2 <- c('Ilex verticillata', 'Salix discolor', 'Salix interior', 'Staphylea trifolia', 
          'Salix bebbiana', 'Salix eriocephala', 'Salix petiolaris', 'Rhamnus cathartica', 'Frangula alnus',
          'Salix exigua', 'Elaeagnus angustifolia','Ptelea trifoliata', 'Toxicodendron vernix')
listspp[listspp$AcTaxon %in% SBD2,]$Form <- 'SBD2'

listspp$Nativity <- ifelse(listspp$Eastern.North.America %in% 'N','Native',ifelse(listspp$Eastern.North.America %in% 'X','Introduced','Unknown'))
obs <- read.delim("data/Sites.txt")
obsspp <- read.delim("data/Observed_Species.txt", encoding = 'UTF-8', na.strings = '')

obsspp <- merge(obs[,c('Observation_ID','Observation_Label','Latitude','Longitude')],obsspp, by='Observation_ID')
obsspp[obsspp$AcSpecies == 'Phalaris' & !is.na(obsspp$AcTaxon) ,]$AcTaxon <- 'Phalaris arundinacea'
obsspp <- obsspp[!grepl("\\?", obsspp$AcTaxon) | is.na(obsspp$AcTaxon) ,]
#obsspp <- subset(obsspp, !is.na(specific_epithet) | AcTaxon %in% c('Sphagnum', 'Chara', 'Cladonia'))
obsspp <- merge(obsspp, List_Habits[,c('Form','Simple')], by.x = 'Habit', by.y = 'Form', all.x = TRUE)
obsspp <- subset(obsspp, Observation_ID %in% plots)

obsspp[obsspp$BA==0,]$BA <- NA
obsspp$BA <- BA.fact10usc.metric(obsspp$BA)


# observed species means by plot ----
Com.Sp.sum<-aggregate(obsspp[,c('Field', 'Shrub', 'Subcanopy', 'Tree', 'BA')], by=list(Latitude=obsspp$Latitude,Longitude=obsspp$Longitude,Observation_ID=obsspp$Observation_ID, Observation_Label=obsspp$Observation_Label, Species=obsspp$AcTaxon, Habit=obsspp$Habit, Simple=obsspp$Simple), FUN=sum) #sum within plot
# freq in subplot ----
Com.Sp.freq<-aggregate(list(freq= obsspp$AcTaxon), by=list(Observation_Label=obsspp$Observation_Label, Species=obsspp$AcTaxon), FUN=length) #frequency within plot

Com.max.freq<-aggregate(list(mfreq=Com.Sp.freq$freq), by=list(Observation_Label=Com.Sp.freq$Observation_Label), FUN=max) #freq within plot
Com.max.freq$mfreq<-ifelse(Com.max.freq$mfreq>4,4,Com.max.freq$mfreq)#effectively ensuring values do not exceed 4. Species listed 5 times might occur if surveyer was unaware of species already counted in subplots, but this only adds a trace amount.
Com.Sp.mean<-merge(Com.Sp.sum, Com.max.freq[,c("Observation_Label","mfreq")], by="Observation_Label")
Com.Sp.mean$Field<-Com.Sp.mean$Field/Com.Sp.mean$mfreq
Com.Sp.mean$Shrub<-Com.Sp.mean$Shrub/Com.Sp.mean$mfreq
Com.Sp.mean$Subcanopy<-Com.Sp.mean$Subcanopy/Com.Sp.mean$mfreq
Com.Sp.mean$Tree<-Com.Sp.mean$Tree/Com.Sp.mean$mfreq
rm(Com.max.freq)
# heights ----
#transfer canopy heights to appropriate stratum
obsspp[!is.na(obsspp$Hmax) & obsspp$Hmax > 15 & is.na(obsspp$Tmax) & !grepl('^H',obsspp$Habit),]$Tmax <- 
  obsspp[!is.na(obsspp$Hmax) & obsspp$Hmax > 15 & is.na(obsspp$Tmax) & !grepl('^H',obsspp$Habit),]$Hmax

obsspp[!is.na(obsspp$Hmax) & obsspp$Hmax > 15 & is.na(obsspp$Tmin) & !grepl('^H',obsspp$Habit),]$Tmin <- 
  obsspp[!is.na(obsspp$Hmax) & obsspp$Hmax > 15 & is.na(obsspp$Tmin) & !grepl('^H',obsspp$Habit),]$Hmin

obsspp[!is.na(obsspp$Hmax) & obsspp$Hmax > 5 & obsspp$Hmax <= 15 & is.na(obsspp$SCmax) & !grepl('^H',obsspp$Habit),]$SCmax <- 
  obsspp[!is.na(obsspp$Hmax) & obsspp$Hmax > 5 & obsspp$Hmax <= 15 & is.na(obsspp$SCmax) & !grepl('^H',obsspp$Habit),]$Hmax

obsspp[!is.na(obsspp$Hmax) & obsspp$Hmax > 5 & obsspp$Hmax <= 15 & is.na(obsspp$SCmin) & !grepl('^H',obsspp$Habit),]$SCmin <- 
  obsspp[!is.na(obsspp$Hmax) & obsspp$Hmax > 5 & obsspp$Hmax <= 15 & is.na(obsspp$SCmin) & !grepl('^H',obsspp$Habit),]$Hmin

obsspp[!is.na(obsspp$Hmax) & obsspp$Hmax > 0.5 & obsspp$Hmax <= 5 & is.na(obsspp$Smax) & !grepl('^H',obsspp$Habit),]$Smax <- 
  obsspp[!is.na(obsspp$Hmax) & obsspp$Hmax > 0.5 & obsspp$Hmax <= 5 & is.na(obsspp$Smax) & !grepl('^H',obsspp$Habit),]$Hmax

obsspp[!is.na(obsspp$Hmax) & obsspp$Hmax > 0.5 & obsspp$Hmax <= 5 & is.na(obsspp$Smin) & !grepl('^H',obsspp$Habit),]$Smin <- 
  obsspp[!is.na(obsspp$Hmax) & obsspp$Hmax > 0.5 & obsspp$Hmax <= 5 & is.na(obsspp$Smin) & !grepl('^H',obsspp$Habit),]$Hmin

obsspp[!is.na(obsspp$Hmax) & obsspp$Hmax > 0 & obsspp$Hmax <= 0.5 & is.na(obsspp$Fmax) | grepl('^H',obsspp$Habit),]$Fmax <- 
  obsspp[!is.na(obsspp$Hmax) & obsspp$Hmax > 0 & obsspp$Hmax <= 0.5 & is.na(obsspp$Fmax) | grepl('^H',obsspp$Habit),]$Hmax

obsspp[!is.na(obsspp$Hmax) & obsspp$Hmax > 0 & obsspp$Hmax <= 0.5 & is.na(obsspp$Fmin)  | grepl('^H',obsspp$Habit),]$Fmin <- 
  obsspp[!is.na(obsspp$Hmax) & obsspp$Hmax > 0 & obsspp$Hmax <= 0.5 & is.na(obsspp$Fmin)  | grepl('^H',obsspp$Habit),]$Hmin

Com.Sp.hts<-aggregate(obsspp[,c('Fmin', 'Fmax', 'Smin', 'Smax', 'SCmin', 'SCmax','Tmin', 'Tmax', 'Dmin', 'Dmax')], by=list(Observation_Label=obsspp$Observation_Label, Species=obsspp$AcTaxon, Habit=obsspp$Habit, Simple=obsspp$Simple), FUN=mean, na.action = na.omit) #sum within plot

#Com.Sp.hts[,c('Fmin', 'Fmax', 'Smin', 'Smax', 'SCmin', 'SCmax','Tmin', 'Tmax', 'Dmin', 'Dmax')]<-   lapply(Com.Sp.hts[,c('Fmin', 'Fmax', 'Smin', 'Smax', 'SCmin', 'SCmax','Tmin', 'Tmax', 'Dmin', 'Dmax')], FUN = roundF)

Com.Sp.mean <- merge(Com.Sp.mean,Com.Sp.hts, by=c('Observation_Label', 'Species', 'Habit', 'Simple'), all.x = T)

#ensure not to exceed 100%
Com.Sp.mean$Field <- ifelse(Com.Sp.mean$Field > 100,100,Com.Sp.mean$Field)
Com.Sp.mean$Shrub <- ifelse(Com.Sp.mean$Shrub > 100,100,Com.Sp.mean$Shrub)
Com.Sp.mean$Subcanopy <- ifelse(Com.Sp.mean$Subcanopy > 100,100,Com.Sp.mean$Subcanopy)
Com.Sp.mean$Tree <- ifelse(Com.Sp.mean$Tree > 100,100,Com.Sp.mean$Tree)
#average overstory and understory
Com.Sp.mean$Total <- 100*(1-10^(apply(log10(1-(Com.Sp.mean[,c('Field', 'Shrub', 'Subcanopy', 'Tree')]/100.001)), MARGIN = 1, FUN='sum')))
Com.Sp.mean <-subset(Com.Sp.mean, !substr(Species,1,1) %in% '-'& !Species %in% '')
Com.Sp.mean <- left_join(Com.Sp.mean, List_Habits[,c('Form', 'ESIS.Group')], by=c("Habit"="Form"))

Com.Sp.mean$maxHt <- pmax(Com.Sp.mean$Fmax,Com.Sp.mean$Smax,Com.Sp.mean$SCmax,Com.Sp.mean$Tmax, na.rm = TRUE)
Com.Sp.mean <- merge(Com.Sp.mean, Plant_Heights[,c('Scientific.Name', 'Ht_m')], by.x = 'Species', by.y='Scientific.Name', all.x=TRUE)
Com.Sp.mean[is.na(Com.Sp.mean$maxHt),]$maxHt <- Com.Sp.mean[is.na(Com.Sp.mean$maxHt),]$Ht_m
Com.Sp.mean[grepl('^L',Com.Sp.mean$Habit),]$maxHt <- NA


Com.Sp.mean[Com.Sp.mean$Tree > 0 & is.na(Com.Sp.mean$Tmin) & Com.Sp.mean$ESIS.Group %in% c('Shrub/Subshrub','Vine/Liana'),]$
  Tmin <- 5
Com.Sp.mean[Com.Sp.mean$Tree > 0 & is.na(Com.Sp.mean$Tmin) & Com.Sp.mean$ESIS.Group %in% c('Tree'),]$
  Tmin <- 10
Com.Sp.mean[Com.Sp.mean$Tree > 0 & is.na(Com.Sp.mean$Tmax) & Com.Sp.mean$ESIS.Group %in% c('Shrub/Subshrub'),]$
  Tmax <- 17
Com.Sp.mean[Com.Sp.mean$Tree > 0 & is.na(Com.Sp.mean$Tmax) & Com.Sp.mean$ESIS.Group %in% c('Vine/Liana'),]$
  Tmax <- 20
Com.Sp.mean[Com.Sp.mean$Tree > 0 & is.na(Com.Sp.mean$Tmax) & Com.Sp.mean$ESIS.Group %in% c('Tree'),]$
  Tmax <- 25
Com.Sp.mean[Com.Sp.mean$Tree > 0 & is.na(Com.Sp.mean$Tmin),]$
  Tmin <- 15

Com.Sp.mean$Tmax <- ifelse(Com.Sp.mean$Tree > 0 & !is.na(Com.Sp.mean$maxHt) & Com.Sp.mean$maxHt > 15,
                                      pmin(Com.Sp.mean$Tmax,
                                           Com.Sp.mean$maxHt),
                                      Com.Sp.mean$Tmax)
Com.Sp.mean$Tmin <- ifelse(Com.Sp.mean$Tree > 0 & !is.na(Com.Sp.mean$maxHt) & Com.Sp.mean$maxHt > 15,
                                      pmin(Com.Sp.mean$Tmin,
                                           Com.Sp.mean$maxHt),
                                      Com.Sp.mean$Tmin)


Com.Sp.mean[Com.Sp.mean$Subcanopy > 0 & is.na(Com.Sp.mean$SCmin) & Com.Sp.mean$ESIS.Group %in% c('Shrub/Subshrub','Vine/Liana'),]$
  SCmin <- 2
Com.Sp.mean[Com.Sp.mean$Subcanopy > 0 & is.na(Com.Sp.mean$SCmin) & Com.Sp.mean$ESIS.Group %in% c('Tree'),]$
  SCmin <- 5
Com.Sp.mean[Com.Sp.mean$Subcanopy > 0 & is.na(Com.Sp.mean$SCmax) & Com.Sp.mean$ESIS.Group %in% c('Shrub/Subshrub'),]$
  SCmax <- 6
Com.Sp.mean[Com.Sp.mean$Subcanopy > 0 & is.na(Com.Sp.mean$SCmax) & Com.Sp.mean$ESIS.Group %in% c('Tree','Vine/Liana'),]$
  SCmax <- 15
Com.Sp.mean[Com.Sp.mean$Subcanopy > 0 & is.na(Com.Sp.mean$SCmin),]$
  SCmin <- 5

Com.Sp.mean$SCmax <- ifelse(Com.Sp.mean$Subcanopy > 0 & !is.na(Com.Sp.mean$maxHt) & Com.Sp.mean$maxHt > 5,
                                      pmin(Com.Sp.mean$SCmax,
                                           Com.Sp.mean$maxHt),
                                      Com.Sp.mean$SCmax)
Com.Sp.mean$SCmin <- ifelse(!is.na(Com.Sp.mean$maxHt) & Com.Sp.mean$maxHt > 5,
                                         pmin(Com.Sp.mean$SCmin,
                                              Com.Sp.mean$maxHt),
                                         Com.Sp.mean$SCmin)





Com.Sp.mean[Com.Sp.mean$Shrub > 0 & is.na(Com.Sp.mean$Smin) & Com.Sp.mean$ESIS.Group %in% c('Shrub/Subshrub'),]$
  Smin <- 0.5
Com.Sp.mean[Com.Sp.mean$Shrub > 0 & is.na(Com.Sp.mean$Smin) & Com.Sp.mean$ESIS.Group %in% c('Tree','Vine/Liana'),]$
  Smin <- 1
Com.Sp.mean[Com.Sp.mean$Shrub > 0 & is.na(Com.Sp.mean$Smax) & Com.Sp.mean$ESIS.Group %in% c('Shrub/Subshrub'),]$
  Smax <- 2
Com.Sp.mean[Com.Sp.mean$Shrub > 0 & is.na(Com.Sp.mean$Smax) & Com.Sp.mean$ESIS.Group %in% c('Tree','Vine/Liana'),]$
  Smax <- 5
Com.Sp.mean[Com.Sp.mean$Shrub > 0 & is.na(Com.Sp.mean$Smin),]$
  Smin <- 0.5

Com.Sp.mean$Smax <- ifelse(Com.Sp.mean$Shrub > 0 & !is.na(Com.Sp.mean$maxHt) & Com.Sp.mean$maxHt > 0.5,
                                      pmin(Com.Sp.mean$Smax,
                                           Com.Sp.mean$maxHt),
                                      Com.Sp.mean$Smax)
Com.Sp.mean$Smin <- ifelse(Com.Sp.mean$Shrub > 0 & !is.na(Com.Sp.mean$maxHt) & Com.Sp.mean$maxHt > 0.5,
                                         pmin(Com.Sp.mean$Smin,
                                              Com.Sp.mean$maxHt),
                                         Com.Sp.mean$Smin)


Com.Sp.mean[Com.Sp.mean$Field > 0 & is.na(Com.Sp.mean$Fmin) & Com.Sp.mean$ESIS.Group %in% c('Nonvascular'),]$Fmin <- 0.0

Com.Sp.mean[Com.Sp.mean$Field > 0 & is.na(Com.Sp.mean$Fmax) & !is.na(Com.Sp.mean$maxHt) & Com.Sp.mean$ESIS.Group %in% c('Grass/grass-like (Graminoids)','Fern/fern ally','Forb/Herb'),]$Fmax <- 
  pmin(Com.Sp.mean[is.na(Com.Sp.mean$Fmax) & !is.na(Com.Sp.mean$maxHt) & Com.Sp.mean$ESIS.Group %in% c('Grass/grass-like (Graminoids)','Fern/fern ally','Forb/Herb'),]$maxHt, 1.5)

Com.Sp.mean[Com.Sp.mean$Field > 0 & is.na(Com.Sp.mean$Fmin) & !is.na(Com.Sp.mean$Fmax) & Com.Sp.mean$ESIS.Group %in% c('Grass/grass-like (Graminoids)','Fern/fern ally','Forb/Herb'),]$Fmin <- Com.Sp.mean[is.na(Com.Sp.mean$Fmin) & !is.na(Com.Sp.mean$Fmax) & Com.Sp.mean$ESIS.Group %in% c('Grass/grass-like (Graminoids)','Fern/fern ally','Forb/Herb'),]$Fmax/2#4+.05

Com.Sp.mean[Com.Sp.mean$Field > 0 & is.na(Com.Sp.mean$Fmin) & Com.Sp.mean$ESIS.Group %in% c('Tree','Vine/Liana','Grass/grass-like (Graminoids)','Fern/fern ally','Forb/Herb', 'Nonvascular'),]$Fmin <- 0.1

Com.Sp.mean[Com.Sp.mean$Field > 0 & is.na(Com.Sp.mean$Fmax) & Com.Sp.mean$ESIS.Group %in% c('Nonvascular'),]$Fmax <- 0.05
Com.Sp.mean[Com.Sp.mean$Field > 0 & is.na(Com.Sp.mean$Fmax) & Com.Sp.mean$ESIS.Group %in% c('Grass/grass-like (Graminoids)','Fern/fern ally','Forb/Herb'),]$Fmax <- 0.6
Com.Sp.mean[Com.Sp.mean$Field > 0 & is.na(Com.Sp.mean$Fmax) & Com.Sp.mean$ESIS.Group %in% c('Shrub/Subshrub'),]$Fmax <- 0.3
Com.Sp.mean[Com.Sp.mean$Field > 0 & is.na(Com.Sp.mean$Fmax) & Com.Sp.mean$ESIS.Group %in% c('Tree','Vine/Liana'),]$Fmax <- 0.5
Com.Sp.mean[Com.Sp.mean$Field > 0 & is.na(Com.Sp.mean$Fmin),]$
  Fmin <- 0

Com.Sp.mean$Fmax <- ifelse(Com.Sp.mean$Field > 0 & !is.na(Com.Sp.mean$maxHt),
                                       pmin(Com.Sp.mean$Fmax,
                                            Com.Sp.mean$maxHt),
                                       Com.Sp.mean$Fmax)
Com.Sp.mean$Fmin <- ifelse(Com.Sp.mean$Field > 0 & !is.na(Com.Sp.mean$maxHt),
                                          pmin(Com.Sp.mean$Fmin,
                                               Com.Sp.mean$maxHt),
                                          Com.Sp.mean$Fmin)

#sorting
#
Com.Sp.mean <- arrange(Com.Sp.mean, Observation_Label, desc(Tree), desc(Subcanopy), desc(Shrub), desc(Field), Species)
Com.Sp.mean <-  subset(Com.Sp.mean, select=c("Latitude","Longitude","Observation_ID","Observation_Label","Species","Habit","ESIS.Group","Simple","Field", "Shrub","Subcanopy","Tree", "Fmin" , "Fmax","Smin","Smax","SCmin","SCmax","Tmin", "Tmax","Dmin","Dmax","BA","Total"))
colnames(Com.Sp.mean)





write.csv(Com.Sp.mean, 'output/plotsummary.csv', na = "", row.names = F)



