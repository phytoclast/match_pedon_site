
#----

if(F){
  newpedons <- soilDB::fetchNASIS(from = "pedons", rmHzErrors = FALSE, SS = FALSE)
  newsite <- aqp::site(newpedons)
  saveRDS(newsite, "data/NASISPEDONS2.RDS", )
  saveRDS(newpedons, 'data/fp2.RDS')
}
NASISPEDONS <- readRDS("data/NASISPEDONS2.RDS")
fp <- readRDS('data/fp2.RDS')

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
List_Habits[List_Habits$ESIS.Group %in% 'Grass/grass-like',]$ESIS.Group <- 'Grass/grass-like (Graminoids)'
# obsspp <- read.delim("data/Observed_Species.txt", encoding = 'UTF-8', na.strings = '')
obsspp <- read.delim("data/Observed_Species.txt", encoding = 'latin-1', na.strings = '')

obsspp <- obsspp %>% mutate(AcTaxon = 
                              case_when(
                                AcTaxon %in% 'Phalaris' ~ 'Phalaris arundinacea',
                                TRUE ~ AcTaxon
                              ))

# obsspp <- obsspp %>% mutate(AcTaxon = ifelse(AcTaxon %in% 'Phalaris','Phalaris arundinacea',AcTaxon))
# 
# idx <- obsspp$AcTaxon %in% 'Phalaris' %>% which()
# obsspp[idx,]$AcTaxon <- 'Phalaris arundinacea'[length(idx)>0]

obsspp <- obsspp %>% mutate(AcTaxon= ifelse(grepl('Lonicera.*bella',AcTaxon),'Lonicera.bella', AcTaxon))
obsspp <- obsspp[!grepl("\\?", obsspp$AcTaxon) | is.na(obsspp$AcTaxon) ,]
obsspp <- subset(obsspp, !is.na(specific_epithet) | AcTaxon %in% c('Sphagnum', 'Chara'))
obsspp[obsspp$BA %in% 0,]$BA <- NA
obsspp$BA <- BA.fact10usc.metric(obsspp$BA)

listspp <- readRDS("data/listspp.RDS")
#saveRDS(listspp, 'data/listspp.RDS')
SBD2 <- c('Ilex verticillata', 'Salix discolor', 'Salix interior', 'Staphylea trifolia', 
          'Salix bebbiana', 'Salix eriocephala', 'Salix petiolaris', 'Rhamnus cathartica', 'Frangula alnus',
          'Salix exigua', 'Elaeagnus angustifolia','Ptelea trifoliata', 'Toxicodendron vernix')
listspp[listspp$AcTaxon %in% SBD2,]$Form <- 'SBD2'
listspp$Nativity <- ifelse(listspp$Eastern.North.America %in% 'N','Native',ifelse(listspp$Eastern.North.America %in% 'X','Introduced','Unknown'))

obsspp <- subset(obsspp, !substr(AcTaxon,1,1) %in% '-'& !AcTaxon %in% '' & !is.na(Habit))
obsspp <- merge(obsspp, List_Habits[,c('Form','Simple')], by.x = 'Habit', by.y = 'Form', all.x = TRUE)

obs <- read.delim("data/Sites.txt")
obs <- subset(obs,Observer_Code %in% 'GRR.GJS' & Year >=2011 & !Observation_Type %in% c('Bogus', 'Floristics')) #filter out bogus points
unique(obs$Observation_Type)
# Obs.shape <- sf::st_as_sf(x = obs, 
#                           coords = c('Longitude', 'Latitude'),
#                           crs = "+proj=longlat +datum=WGS84")
# sf::write_sf(Obs.shape, 'data/Obs.shape.shp')

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

VEGOBS_mukeys <- readRDS('output/ssurgo.RDS')
VEGOBS_soilnames <- merge(VEGOBS_mukeys[,c('obs.id','mukey')], mu[,c('lmapunitiid', 'muname')], by.x='mukey', by.y= 'lmapunitiid')
VEGOBS <- merge(VEGOBS, VEGOBS_soilnames[,c('obs.id', 'muname')], by.x='Observation_ID', by.y= 'obs.id')

VEGOBS <- VEGOBS %>% mutate(
  eval=case_when(
  pedondist < 50 ~ "dump",
  pedondist < 50 ~ "keep1",
  TRUE~''),
  eval=case_when(
  (str_split_fixed(muname, " ",2)[,1]) %in% taxonname & pedondist < 1000 & eval  %in%  'dump' ~ "keep2",
  pedondate %in% Year & pedondist < 100 & eval %in% 'dump'~ "keep3",
  TRUE~eval),

  taxonname=ifelse(eval %in% 'dump',"", taxonname),
  taxonclass=ifelse(eval %in% 'dump',"", taxonclass),
  pedon=ifelse(eval %in% 'dump',"", pedon)
  )

VEGOBS$Soil <- VEGOBS$taxonname
for (i in 1:nrow(pedonoverride)){ #replace known soils that disagree with map unit and don't have pedon records
  VEGOBS[VEGOBS$Observation_Label %in% pedonoverride$sitelabel[i],]$Soil <- as.character(pedonoverride$soil[i])}
VEGOBS[VEGOBS$Soil %in% '',]$Soil <- str_split_fixed(VEGOBS[VEGOBS$Soil %in% '',]$muname, " ",2)[,1]
write.csv(VEGOBS, 'output/VEGOBS.csv', row.names = FALSE, na = "")


#----
#narrow to soil series
# soilgroup <- 'all'
# ngroups <- 18
# if (T){
#   remove <- c('2016MI037002',
#               'S12060501', 
#               'S12060502',
#               'S12060503',
#               'S12071304',
#               's20190701.02'
#   )
#   add <- 'S12062503'
#   #sortsoils <- unique(subset(s, (T150_OM >= 10 | grepl('histic',taxsubgrp) |grepl('histosols',taxorder)) & !grepl('Dysic',taxclname), select = 'compname'))[,1]
#   sortsoils <- unique(subset(s, (T150_OM >= 10 | grepl('histic',taxsubgrp) |grepl('histosols',taxorder)) , select = 'compname'))[,1]
#   #VEGOBS <- subset(VEGOBS,Soil %in% sortsoils & !Observation_Label %in% remove |Observation_Label %in% add )
#   VEGOBS <- subset(VEGOBS,Soil %in% sortsoils |Observation_Label %in% add )
#   ngroups <- 8
#   soilgroup <- 'euic_mucks'}


#narrow to soil series
sortsoils <- unique(
  subset(s, 
         hydricrating %in% 'yes' &
           (T50_sand < 70 & T150_sand < 80) &
           !flood %in% 'flood'&
           !(T150_OM >= 10 | grepl('histic',taxsubgrp) |grepl('histosols',taxorder))&
           !compname %in% 'Kingsville', 
         select = 'compname'))[,1]
#VEGOBS <- subset(VEGOBS,Soil %in% sortsoils & !Observation_Label %in% remove |Observation_Label %in% add )
#
obsids <- c('s20210805.02','s20220720.001','s20220719.001','s20220719.002')
VEGOBS <- subset(VEGOBS,Soil %in% sortsoils & !Observation_Type %in% c('Bogus', 'Floristic','Site Index') | Observation_Label %in% obsids)
ngroups <- 8
soilgroup <- 'wetloamy'
#----
#observed species

Com.Sp.sum<-aggregate(obsspp[,c('Field', 'Shrub', 'Subcanopy', 'Tree', 'BA')], by=list(obsspp$Observation_ID,obsspp$Observation_Label, obsspp$AcTaxon, obsspp$Simple), FUN=sum) #sum within plot
colnames(Com.Sp.sum)<-c('Observation_ID', 'Observation_Label', 'Species', 'Simple', 'Field', 'Shrub', 'Subcanopy', 'Tree', 'BA') #restore column names

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
obs$soilplot <- paste(obs$Soil , obs$Observation_Label)
obs$soilplot <- str_replace_all(obs$soilplot, ' ', '.')
obs$soilplot <- str_replace_all(obs$soilplot, '-', '.')
obs$soilplot <- str_replace_all(obs$soilplot, ',', '.')
obs$soilplot <- str_replace_all(obs$soilplot, ':', '.')
obs$soilplot <- str_replace_all(obs$soilplot, ';', '.')
obs$soilplot <- str_replace_all(obs$soilplot, '&', '.')
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

#plotdata <- makecommunitydataset(Com.Sp.mean, row = 'soilplot', column = 'Species', value = 'sqrttotal', drop = TRUE)


Com.Sp.mean$logtotal1 <- (log10(100*(1-10^(apply(log10(1-(Com.Sp.mean[,c('Field', 'Shrub', 'Subcanopy', 'Tree')]/100.001)), MARGIN = 1, FUN='sum'))))+2.1)
Com.Sp.mean$logtotal2 <- ifelse(apply(Com.Sp.mean[,c('Shrub', 'Subcanopy', 'Tree')], MARGIN = 1, FUN=sum) == 0,0, (log10(100*(1-10^(apply(log10(1-(Com.Sp.mean[,c('Shrub', 'Subcanopy', 'Tree')]/100.001)), MARGIN = 1, FUN='sum'))))+2.1))
Com.Sp.mean$logtotal3 <- ifelse(apply(Com.Sp.mean[,c('Subcanopy', 'Tree')], MARGIN = 1, FUN=sum) == 0,0,(log10(100*(1-10^(apply(log10(1-(Com.Sp.mean[,c('Subcanopy', 'Tree')]/100.001)), MARGIN = 1, FUN='sum'))))+2.1))
Com.Sp.mean$logtotal4 <- ifelse(Com.Sp.mean[,c('Tree')] == 0,0,(log10(100*(1-10^(log10(1-(Com.Sp.mean[,c( 'Tree')]/100.001)))))+2.1))

Com.Sp.mean$logtotal <- (Com.Sp.mean$logtotal1+Com.Sp.mean$logtotal2+Com.Sp.mean$logtotal3+Com.Sp.mean$logtotal4)/4

plotdata <- makecommunitydataset(Com.Sp.mean, row = 'soilplot', column = 'Species', value = 'logtotal', drop = TRUE)


overstorycover <- Com.Sp.mean %>% mutate(overstorycover = (apply(log10(1-(Com.Sp.mean[,c('Subcanopy', 'Tree')]/100.001)), MARGIN = 1, FUN='sum')),
                                         shrubberstorycover = (apply(log10(1-(Com.Sp.mean[,c('Shrub','Subcanopy', 'Tree')]/100.001)), MARGIN = 1, FUN='sum')))
overstorycover <- overstorycover %>% group_by(soilplot) %>% summarise(overstorycover = 100*(1-10^(sum(overstorycover))), 
                                                                      shrubberstorycover = 100*(1-10^(sum(shrubberstorycover))),
                                                                      shrubcover = shrubberstorycover - overstorycover)

