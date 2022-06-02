library(sf)
library(soilDB) 
library(sp)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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
obs <- subset(obs,Observer_Code %in% 'GRR.GJS' & Year >=2011 & !Observation_Type %in% c('Bogus', 'Floristics')) #



obs.geom <- st_as_sf(obs[,c('Observation_ID', 'Longitude', 'Latitude')], coords = c('Longitude', 'Latitude'))
st_crs(obs.geom) <-  st_crs('epsg:4326')
obs.geom2 <-  as_Spatial(obs.geom)
ssurgo=NULL
ssurgo.i<-NULL
for (i in 1:nrow(obs.geom2)){#i=400 length(db)
  obs.geom.i <- obs.geom2[i,]
  ssurgo.i <- SDA_spatialQuery(
    obs.geom.i,
    what = "mukey",
    geomIntersection = FALSE,
    db = c("SSURGO")
  )
  if(!is.null(ssurgo.i)){
    ssurgo.i <-  cbind(list(obs.id = obs.geom.i$Observation_ID[1]), ssurgo.i)}
  else{
    ssurgo.i <-  list(obs.id = obs.geom.i$Observation_ID[1], mukey = NA, muname = NA)}
  if(is.null(ssurgo)){
    ssurgo <- ssurgo.i
  }else{
    ssurgo <- rbind(ssurgo,ssurgo.i)
  }
}
saveRDS(ssurgo, 'output/ssurgo.RDS')
