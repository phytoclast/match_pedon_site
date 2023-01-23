
library(soilDB)
library(aqp)#load before dplyr
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
library(optpart)
library(dendsort)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


allsites <- soilDB::get_site_data_from_NASIS_db(SS=F)
obs <- read.delim("data/Sites.txt")
obs <- subset(obs,Observer_Code %in% 'GRR.GJS' & Year >=2011 & !Observation_Type %in% c('Bogus', 'Floristics')) #
plotids <- subset(obs, select=c("Site_Type","Observation_ID","Observation_Label","User_Plot_ID","Observation_Type","Latitude","Longitude", "Year","Mon","Day","State","County", "User_Pedon_ID"))

plotids$oldid <- plotids$User_Plot_ID
plotids$dist <- -999 %>% as.numeric()
rownames(plotids) <- 1:nrow(plotids)
n <- nrow(plotids)
for (i in 1:n){#i=418
  
  allsites$distance <- (((plotids[i,]$Latitude - allsites$y_std)/360*40041.47*1000)^2 +
                             ((plotids[i,]$Longitude - allsites$x_std)/360*40041.47*1000*cos(plotids[i,]$Latitude/2/360*2*3.141592))^2)^0.5


mindist <- min(allsites$distance, na.rm = TRUE)
allsites0 <- allsites %>% subset(distance %in% mindist)

# plotids[i,]$User_Plot_ID <-  ifelse(is.na(plotids[i,]$User_Plot_ID)|nchar(plotids[i,]$User_Plot_ID)<1, plotids[i,]$User_Plot_ID, ifelse(mindist<50, allsites0[1,]$site_id, NA))
# 
# plotids[i,]$dist <-  ifelse(is.na(plotids[i,]$User_Plot_ID)|nchar(plotids[i,]$User_Plot_ID)<1, plotids[i,]$User_Plot_ID, ifelse(mindist<50, mindist, NA) %>% as.numeric)

plotids[i,]$User_Plot_ID <-  ifelse(mindist<50, allsites0[1,]$site_id, NA)

plotids[i,]$dist <-  ifelse(mindist<50, mindist, NA) %>% as.numeric

}
plotids$dist <- as.numeric(plotids$dist)
write.csv(plotids, 'output/plotids.csv', row.names = F, na = "")

flagmismatch <- subset(plotids, !User_Plot_ID %in% oldid )
