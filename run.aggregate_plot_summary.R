substitutesymbols <- as.data.frame(cbind(SYMBOL2=c('2AG','TRBO2','SMLA3','CAREX'), form2=c('Nonvascular','Forb/Herb','Forb/Herb','Grass/grass-like (Graminoids)'), type2=c('Native','Native','Native','Native'), taxon=c('Chara','Lysimachia borealis','Smilax lasioneuron','Carex [ovales]')), stringsAsFactors = FALSE)





# data ----
Plant_Heights <- read.delim("data/Plant_Heights.txt", na.strings = '', stringsAsFactors = FALSE, encoding = 'UTF-8')




# heights ----
#transfer canopy heights to appropriate stratum
#
obsspp <- obsspp %>% mutate(
  Tmax = ifelse(
    !is.na(Hmax) & Hmax > 15 & is.na(Tmax) & !grepl('^H',Habit), Hmax, Tmax),
  Tmin = ifelse(
    !is.na(Hmax) & Hmax > 15 & is.na(Tmin) & !grepl('^H',Habit), Hmin, Tmin),
  SCmax = ifelse(
    !is.na(Hmax) & Hmax > 5 & Hmax <= 15 & is.na(SCmax) & !grepl('^H',Habit), Hmax, SCmax),
  SCmin = ifelse(
    !is.na(Hmax) & Hmax > 5 & Hmax <= 15 & is.na(SCmin) & !grepl('^H',Habit), Hmin, SCmin),
  Smax = ifelse(
    !is.na(Hmax) & Hmax > 0.5 & Hmax <= 5 & is.na(Smax) & !grepl('^H',Habit), Hmax, Smax),
  Smin = ifelse(
    !is.na(Hmax) & Hmax > 0.5 & Hmax <= 5 & is.na(Smin) & !grepl('^H',Habit), Hmin, Smin),
  Fmax = ifelse(
    !is.na(Hmax) & Hmax > 0 & Hmax <= 0.5 & is.na(Fmax) | grepl('^H',Habit), Hmax, Fmax),
  Fmin = ifelse(
    !is.na(Hmax) & Hmax > 0 & Hmax <= 0.5 & is.na(Fmin) | grepl('^H',Habit), Hmin, Fmin))


Com.Sp.hts <- obsspp %>% group_by(Observation_Label, AcTaxon, Simple) %>%
  summarise(
    Fmin=mean(Fmin, na.rm=T),
    Fmax=mean(Fmax, na.rm=T),
    Smin=mean(Smin, na.rm=T),
    Smax=mean(Smax, na.rm=T),
    SCmin=mean(SCmin, na.rm=T),
    SCmax=mean(SCmax, na.rm=T),
    Tmin=mean(Tmin, na.rm=T),
    Tmax=mean(Tmax, na.rm=T),
    Dmin=mean(Dmin, na.rm=T),
    Dmax=mean(Dmax, na.rm=T)
  )
  

Com.Sp.mean.ht <- Com.Sp.mean %>% left_join(Com.Sp.hts, by=c('Observation_Label'='Observation_Label', 'Species'='AcTaxon', 'Simple'='Simple'))

#ensure not to exceed 100%
Com.Sp.mean.ht$Field <- ifelse(Com.Sp.mean.ht$Field > 100,100,Com.Sp.mean.ht$Field)
Com.Sp.mean.ht$Shrub <- ifelse(Com.Sp.mean.ht$Shrub > 100,100,Com.Sp.mean.ht$Shrub)
Com.Sp.mean.ht$Subcanopy <- ifelse(Com.Sp.mean.ht$Subcanopy > 100,100,Com.Sp.mean.ht$Subcanopy)
Com.Sp.mean.ht$Tree <- ifelse(Com.Sp.mean.ht$Tree > 100,100,Com.Sp.mean.ht$Tree)
#average overstory and understory
Com.Sp.mean.ht$Total <- 100*(1-10^(apply(log10(1-(Com.Sp.mean.ht[,c('Field', 'Shrub', 'Subcanopy', 'Tree')]/100.001)), MARGIN = 1, FUN='sum')))
Com.Sp.mean.ht <-subset(Com.Sp.mean.ht, !substr(Species,1,1) %in% '-'& !Species %in% '')

##cluster derived phases ----


groups1 <- as.data.frame(cbind(soilplot=names(groups), phase=groups))
Com.Sp.groups <- merge(groups1, Com.Sp.mean.ht, by='soilplot')

##average spp by phase ----


Com.Sp.groups.sum <- aggregate(Com.Sp.groups[,c('Field', 'Shrub', 'Subcanopy','Tree', 'Total')],
                               by=list(Com.Sp.groups$phase, Com.Sp.groups$Species), FUN='sum')
colnames(Com.Sp.groups.sum) <- c('phase', 'Species', 'Field', 'Shrub', 'Subcanopy','Tree', 'Total')
Com.Sp.groups.count <- aggregate(unique(Com.Sp.groups[c('phase', 'Observation_ID')])$Observation_ID, 
                                 by=list(unique(Com.Sp.groups[c('phase', 'Observation_ID')])$phase), FUN='length')
colnames(Com.Sp.groups.count) <- c('phase', 'count')
Com.Sp.groups.mean <- merge(Com.Sp.groups.sum, Com.Sp.groups.count, by = 'phase')
Com.Sp.groups.mean[,c('Field', 'Shrub', 'Subcanopy','Tree', 'Total')] <- Com.Sp.groups.mean[,c('Field', 'Shrub', 'Subcanopy','Tree', 'Total')]/Com.Sp.groups.mean$count
rm(Com.Sp.groups.sum, Com.Sp.groups.count)
# aggregated total overstory and understory for species ranking
Com.Sp.groups.mean$over <- 100*(1-10^(apply(log10(1-(Com.Sp.groups.mean[,c('Subcanopy', 'Tree')]/100.001)), MARGIN = 1, FUN='sum')))
Com.Sp.groups.mean$under <- 100*(1-10^(apply(log10(1-(Com.Sp.groups.mean[,c('Field', 'Shrub')]/100.001)), MARGIN = 1, FUN='sum')))
##average heights by phase ----
Com.Sp.groups.hts <- Com.Sp.groups %>% group_by(phase, Species, Simple) %>%
  summarise(
    Fmin=mean(Fmin, na.rm=T),
    Fmax=mean(Fmax, na.rm=T),
    Smin=mean(Smin, na.rm=T),
    Smax=mean(Smax, na.rm=T),
    SCmin=mean(SCmin, na.rm=T),
    SCmax=mean(SCmax, na.rm=T),
    Tmin=mean(Tmin, na.rm=T),
    Tmax=mean(Tmax, na.rm=T),
    Dmin=mean(Dmin, na.rm=T),
    Dmax=mean(Dmax, na.rm=T)
  )
Com.Sp.groups.hts <-  as.data.frame(Com.Sp.groups.hts)

for (i in 1:ncol(Com.Sp.groups.hts)){#i=3
Com.Sp.groups.hts[is.nan(Com.Sp.groups.hts[,i]),i] <- NA
  }#clean up NaNs


#percentile spp by phase ----
#
countbyphase <- subset(Com.Sp.groups, select= c(soilplot, phase)) %>% unique() %>% group_by(phase) %>% summarise(total=length(soilplot))
Com.Sp.groups.pctl <- Com.Sp.groups %>% group_by(phase, Species) %>% summarise(
  f25 = quantile(Field, 0.15),
  f75 = quantile(Field, 0.85),
  s25 = quantile(Shrub, 0.15),
  s75 = quantile(Shrub, 0.85),
  sc25 = quantile(Subcanopy, 0.15),
  sc75 = quantile(Subcanopy, 0.85),
  t25 = quantile(Tree, 0.15),
  t75 = quantile(Tree, 0.85),
  b05 = quantile(BA, 0.05, na.rm = TRUE),
  b95 = quantile(BA, 0.95, na.rm = TRUE),
  freq = sum(Field+Shrub+Subcanopy+Tree>0)
) %>% left_join(countbyphase)
Com.Sp.groups.pctl$freq <- Com.Sp.groups.pctl$freq/Com.Sp.groups.pctl$total*100
# adjust max min values to frequency
Com.Sp.groups.pctl$f75 <- pmin(1,Com.Sp.groups.pctl$freq*3/200)*Com.Sp.groups.pctl$f75
Com.Sp.groups.pctl$s75 <- pmin(1,Com.Sp.groups.pctl$freq*3/200)*Com.Sp.groups.pctl$s75
Com.Sp.groups.pctl$sc75 <- pmin(1,Com.Sp.groups.pctl$freq*3/200)*Com.Sp.groups.pctl$sc75
Com.Sp.groups.pctl$t75 <- pmin(1,Com.Sp.groups.pctl$freq*3/200)*Com.Sp.groups.pctl$t75
Com.Sp.groups.pctl$b95 <- pmin(1,Com.Sp.groups.pctl$freq*3/200)*Com.Sp.groups.pctl$b95

Com.Sp.groups.pctl$f25 <- pmax(0,Com.Sp.groups.pctl$freq*3/200-1/2)*Com.Sp.groups.pctl$f25
Com.Sp.groups.pctl$s25 <- pmax(0,Com.Sp.groups.pctl$freq*3/200-1/2)*Com.Sp.groups.pctl$s25
Com.Sp.groups.pctl$sc25 <- pmax(0,Com.Sp.groups.pctl$freq*3/200-1/2)*Com.Sp.groups.pctl$sc25
Com.Sp.groups.pctl$t25 <- pmax(0,Com.Sp.groups.pctl$freq*3/200-1/2)*Com.Sp.groups.pctl$t25
Com.Sp.groups.pctl$b05 <- pmax(0,Com.Sp.groups.pctl$freq*3/200-1/2)*Com.Sp.groups.pctl$b05




#Com.Sp.groups.pctl[,c('f05','s05','sc05','t05','f75','s75','sc75','t75','b05','b95')] <- lapply(Com.Sp.groups.pctl[,c('f05','s05','sc05','t05','f75','s75','sc75','t75','b05','b95')], FUN=roundF)

##combine results ----

Group.Summary <- merge(Com.Sp.groups.pctl, Com.Sp.groups.hts, by=c('phase','Species'), all.x = T)
Group.Summary <- merge(listspp[,c('Taxon', 'Symbol','Form','Nativity')], Group.Summary, by.x = 'Taxon', by.y = 'Species', all.y = TRUE)
Group.Summary <- merge(Group.Summary, substitutesymbols, by.x = 'Taxon', by.y='taxon', all.x = TRUE)
Group.Summary <- merge(List_Habits[,c('Form', 'ESIS.Group')], Group.Summary, by = 'Form', all.y = TRUE)
Group.Summary <- merge(Group.Summary, Com.Sp.groups.mean[,c('Species','phase','over','under')], 
                       by.x = c('Taxon','phase'), by.y=c('Species','phase'), all.x = TRUE)
##replacement values ----

Group.Summary[is.na(Group.Summary$Symbol) & !is.na(Group.Summary$SYMBOL2),]$Symbol <-
  Group.Summary[is.na(Group.Summary$Symbol) & !is.na(Group.Summary$SYMBOL2),]$SYMBOL2

Group.Summary[is.na(Group.Summary$ESIS.Group) & !is.na(Group.Summary$form2),]$ESIS.Group <-
  Group.Summary[is.na(Group.Summary$ESIS.Group) & !is.na(Group.Summary$form2),]$form2

Group.Summary[is.na(Group.Summary$Nativity) & !is.na(Group.Summary$type2),]$Nativity <-
  Group.Summary[is.na(Group.Summary$Nativity) & !is.na(Group.Summary$type2),]$type2



# group summary export and fix missing heights ----
Group.Summary <- Group.Summary[order(Group.Summary$phase, Group.Summary$Taxon),c( 'phase','Taxon', 'Symbol', 'ESIS.Group','Nativity', 'f25', 'f75', 's25', 's75', 'sc25', 'sc75', 't25', 't75', 'Fmin', 'Fmax', 'Smin', 'Smax', 'SCmin', 'SCmax', 'Tmin', 'Tmax','Dmin','Dmax','b05','b95', 'over', 'under')]

colnames(Group.Summary) <- c( 'phase','Taxon', 'Plant.Symbol', 'Plant.Type','Nativity', 'f25', 'f75', 's25', 's75', 'sc25', 'sc75', 't25', 't75', 'Fmin', 'Fmax', 'Smin', 'Smax', 'SCmin', 'SCmax', 'Tmin', 'Tmax','Dmin','Dmax','b05','b95', 'over', 'under')

Group.Summary$maxHt <- pmax(Group.Summary$Fmax,Group.Summary$Smax,Group.Summary$SCmax,Group.Summary$Tmax, na.rm = TRUE)
Group.Summary <- merge(Group.Summary, Plant_Heights[,c('Scientific.Name', 'Ht_m')], by.x = 'Taxon', by.y='Scientific.Name', all.x=TRUE)
Group.Summary[is.na(Group.Summary$maxHt),]$maxHt <- Group.Summary[is.na(Group.Summary$maxHt),]$Ht_m
Group.Summary[Group.Summary$Plant.Type %in% 'Vine/Liana',]$maxHt <- NA


Forest.Overstory <- Group.Summary[Group.Summary$t75 >0,c( 'phase','Taxon',  'Plant.Symbol', 'Plant.Type','Nativity', 't25', 't75','Tmin', 'Tmax','Dmin','Dmax','b05','b95', 'maxHt', 'over')]
colnames(Forest.Overstory) <- c('phase','Taxon',  'Plant.Symbol', 'Plant.Type','Nativity', 'Cover.Low', 'Cover.High','Canopy.Bottom', 'Canopy.Top','Diam.Low','Diam.High','BA.Low','BA.High', 'maxHt', 'over')


Forest.Overstory <- Forest.Overstory %>% mutate(
  Canopy.Bottom = ifelse(
    is.na(Canopy.Bottom) & Plant.Type %in% c('Shrub/Subshrub','Vine/Liana'), 5, Canopy.Bottom),
  Canopy.Bottom = ifelse(
    is.na(Canopy.Bottom) & Plant.Type %in% c('Tree'), 10, Canopy.Bottom),
  Canopy.Top = ifelse(
    is.na(Canopy.Top) & Plant.Type %in% c('Shrub/Subshrub'), 17, Canopy.Top),
  Canopy.Top = ifelse(
    is.na(Canopy.Top) & Plant.Type %in% c('Vine/Liana'), 20, Canopy.Top),
  Canopy.Top = ifelse(
    is.na(Canopy.Top) & Plant.Type %in% c('Tree'), 25, Canopy.Top),
  Canopy.Bottom = ifelse(
    is.na(Canopy.Bottom), 15, Canopy.Bottom),
  Canopy.Top = ifelse(
    !is.na(maxHt) & maxHt > 15, pmin(Canopy.Top, maxHt),Canopy.Top),
  Canopy.Bottom = ifelse(
    !is.na(maxHt) & maxHt > 15, pmin(Canopy.Bottom, maxHt), Canopy.Bottom))


Forest.Overstory.sub<- Group.Summary[Group.Summary$sc75 >0,c( 'phase','Taxon', 'Plant.Symbol', 'Plant.Type','Nativity','sc25', 'sc75','SCmin', 'SCmax', 'Dmin','Dmax','b05','b95', 'maxHt', 'over')]
Forest.Overstory.sub[,c('Dmin','Dmax','b05','b95')]<- NA
colnames(Forest.Overstory.sub) <- c('phase','Taxon', 'Plant.Symbol', 'Plant.Type','Nativity', 'Cover.Low', 'Cover.High','Canopy.Bottom', 'Canopy.Top','Diam.Low','Diam.High','BA.Low','BA.High','maxHt', 'over')

Forest.Overstory.sub <- Forest.Overstory.sub %>% mutate(
  Canopy.Bottom = ifelse(
    is.na(Canopy.Bottom) & Plant.Type %in% c('Shrub/Subshrub','Vine/Liana'), 2, Canopy.Bottom),
  Canopy.Bottom = ifelse(
    is.na(Canopy.Bottom) & Plant.Type %in% c('Tree'), 5, Canopy.Bottom),
  Canopy.Top = ifelse(
    is.na(Canopy.Top) & Plant.Type %in% c('Shrub/Subshrub'), 6, Canopy.Top),
  Canopy.Top = ifelse(
    is.na(Canopy.Top) & Plant.Type %in% c('Vine/Liana','Tree'), 15, Canopy.Top),
  Canopy.Bottom = ifelse(
    is.na(Canopy.Bottom), 5, Canopy.Bottom),
  Canopy.Top = ifelse(
    !is.na(maxHt) & maxHt > 5, pmin(Canopy.Top, maxHt),Canopy.Top),
  Canopy.Bottom = ifelse(
    !is.na(maxHt) & maxHt > 5, pmin(Canopy.Bottom, maxHt), Canopy.Bottom))

Forest.Overstory <- rbind(Forest.Overstory, Forest.Overstory.sub)



Forest.Understory <- Group.Summary[Group.Summary$s75 >0,c( 'phase','Taxon', 'Plant.Symbol', 'Plant.Type','Nativity', 's25', 's75','Smin', 'Smax', 'maxHt', 'under')]
colnames(Forest.Understory) <- c('phase','Taxon','Plant.Symbol', 'Plant.Type','Nativity', 'Cover.Low', 'Cover.High','Canopy.Bottom', 'Canopy.Top', 'maxHt', 'under')

Forest.Understory <- Forest.Understory %>% mutate(
  Canopy.Bottom = ifelse(
    is.na(Canopy.Bottom) & Plant.Type %in% c('Shrub/Subshrub'), 0.5, Canopy.Bottom),
  Canopy.Bottom = ifelse(
    is.na(Canopy.Bottom) & Plant.Type %in% c('Tree','Vine/Liana'), 1, Canopy.Bottom),
  Canopy.Top = ifelse(
    is.na(Canopy.Top) & Plant.Type %in% c('Shrub/Subshrub'), 2, Canopy.Top),
  Canopy.Top = ifelse(
    is.na(Canopy.Top) & Plant.Type %in% c('Tree','Vine/Liana'), 5, Canopy.Top),
  Canopy.Bottom = ifelse(
    is.na(Canopy.Bottom), 0.5, Canopy.Bottom),
  Canopy.Top = ifelse(
    !is.na(maxHt) & maxHt > 0.5, pmin(Canopy.Top, maxHt),Canopy.Top),
  Canopy.Bottom = ifelse(
    !is.na(maxHt) & maxHt > 0.5, pmin(Canopy.Bottom, maxHt), Canopy.Bottom))


Forest.Understory.sub <- Group.Summary[Group.Summary$f75 >0,c( 'phase','Taxon','Plant.Symbol', 'Plant.Type','Nativity', 'f25', 'f75','Fmin', 'Fmax', 'maxHt', 'under')]
colnames(Forest.Understory.sub) <- c('phase','Taxon','Plant.Symbol', 'Plant.Type','Nativity', 'Cover.Low', 'Cover.High','Canopy.Bottom', 'Canopy.Top', 'maxHt', 'under')


Forest.Understory.sub <- Forest.Understory.sub %>% mutate(
  Canopy.Bottom = ifelse(
    is.na(Forest.Understory.sub$Canopy.Bottom) & Forest.Understory.sub$Plant.Type %in% c('Nonvascular'), 0, Canopy.Bottom),
  Canopy.Top = ifelse(
    is.na(Canopy.Top) & !is.na(maxHt) & Plant.Type %in% c('Grass/grass-like (Graminoids)','Fern/fern ally','Forb/Herb'), 1.5, Canopy.Top),
  Canopy.Bottom = ifelse(
    is.na(Canopy.Bottom) & !is.na(Canopy.Top) & Plant.Type %in% c('Grass/grass-like (Graminoids)','Fern/fern ally','Forb/Herb'), Canopy.Top/2, Canopy.Bottom),
  Canopy.Bottom = ifelse(
    is.na(Canopy.Bottom) & Plant.Type %in% c('Tree','Vine/Liana','Grass/grass-like (Graminoids)','Fern/fern ally','Forb/Herb', 'Nonvascular'), 0, Canopy.Bottom),
  Canopy.Top = ifelse(
    !is.na(maxHt) & maxHt > 0.5, pmin(Canopy.Top, maxHt),Canopy.Top),
  Canopy.Bottom = ifelse(
    !is.na(maxHt) & maxHt > 0.5, pmin(Canopy.Bottom, maxHt), Canopy.Bottom))



Forest.Understory <- rbind(Forest.Understory, Forest.Understory.sub)



Forest.Overstory[,c('Cover.Low', 'Cover.High', 'Canopy.Bottom', 'Canopy.Top', 'Diam.Low', 'Diam.High', 'BA.Low', 'BA.High')] <- lapply(Forest.Overstory[,c('Cover.Low', 'Cover.High', 'Canopy.Bottom', 'Canopy.Top', 'Diam.Low', 'Diam.High', 'BA.Low', 'BA.High')], FUN=roundF)

Forest.Understory[,c('Cover.Low', 'Cover.High', 'Canopy.Bottom', 'Canopy.Top')] <- lapply(Forest.Understory[,c('Cover.Low', 'Cover.High', 'Canopy.Bottom', 'Canopy.Top')], FUN=roundF)

Sort_Habits <- read.delim("data/Sort_Habits.txt", na.strings = '', stringsAsFactors = FALSE)

Forest.Overstory1 <- merge(Forest.Overstory, Sort_Habits, by='Plant.Type', all.x = TRUE)
Forest.Understory1 <- merge(Forest.Understory, Sort_Habits, by='Plant.Type', all.x = TRUE)
Forest.Overstory1 <- Forest.Overstory1[order(
  Forest.Overstory1$phase,
  Forest.Overstory1$order,
  -Forest.Overstory1$over,
  Forest.Overstory1$Taxon,
  -Forest.Overstory1$Canopy.Top
),c('phase',	'Taxon', 'order', 'over', 'Plant.Symbol',	'Plant.Type',	'Nativity',	'Cover.Low',	'Cover.High',	'Canopy.Bottom',	'Canopy.Top',	'Diam.Low',	'Diam.High',	'BA.Low',	'BA.High')]

Forest.Understory1 <- Forest.Understory1[order(
  Forest.Understory1$phase,
  Forest.Understory1$order,
  -Forest.Understory1$under,
  Forest.Understory1$Taxon,
  -Forest.Understory1$Canopy.Top
),c('phase',	'Taxon', 'order', 'under', 'Plant.Symbol',	'Plant.Type',	'Nativity',	'Cover.Low',	'Cover.High',	'Canopy.Bottom',	'Canopy.Top')]


write.csv(Forest.Overstory1, 'output/Forest.Overstory.csv', row.names = FALSE, na = "")
write.csv(Forest.Understory1, 'output/Forest.Understory.csv', row.names = FALSE, na = "")
#ESIS profile ----

ESIS <- rbind(Forest.Understory[,c('phase','Taxon','Plant.Symbol','Plant.Type','Nativity','Cover.Low','Cover.High','Canopy.Bottom','Canopy.Top')], Forest.Overstory[,c('phase','Taxon','Plant.Symbol','Plant.Type','Nativity','Cover.Low','Cover.High','Canopy.Bottom','Canopy.Top')])
ESIS <- ESIS %>% mutate(Plant.Type2 = case_when(
  Plant.Type %in% c('Grass/grass-like (Graminoids)') ~ 'Grass/Grasslike',
  Plant.Type %in% c('Vine/Liana', 'Shrub/Subshrub') ~ 'Shrub/Vine',
  Plant.Type %in% c('Tree', 'Tree Fern') ~ 'Tree',
  Plant.Type %in% c('Nonvascular') ~ 'Nonvascular',
  TRUE ~ 'Forb'))
                                                
                                      
ESIS.01 <- subset(ESIS, Canopy.Top <= 0.5| Plant.Type %in% c('Grass/grass-like (Graminoids)', 'Fern/fern ally', 'Forb/Herb'))
ESIS.01 <- merge(ESIS.01, Com.Sp.groups[,c('Observation_ID','phase','Observation_Label','Species','Simple','Field')], by.x = c('phase','Taxon'),  by.y = c('phase','Species') )
names(ESIS.01)[names(ESIS.01) == 'Field'] <- 'Cover'
ESIS.02 <- subset(ESIS, Canopy.Top > 0.5 & Canopy.Top <= 5 & !Plant.Type %in% c('Grass/grass-like (Graminoids)', 'Fern/fern ally', 'Forb/Herb'))
ESIS.02 <- merge(ESIS.02, Com.Sp.groups[,c('Observation_ID','phase','Observation_Label','Species','Simple','Shrub')], by.x = c('phase','Taxon'),  by.y = c('phase','Species') )
names(ESIS.02)[names(ESIS.02) == 'Shrub'] <- 'Cover'
ESIS.03 <- subset(ESIS, Canopy.Top > 5 & Canopy.Top <= 15 & !Plant.Type %in% c('Grass/grass-like (Graminoids)', 'Fern/fern ally', 'Forb/Herb'))
ESIS.03 <- merge(ESIS.03, Com.Sp.groups[,c('Observation_ID','phase','Observation_Label','Species','Simple','Subcanopy')], by.x = c('phase','Taxon'),by.y = c('phase','Species'))
names(ESIS.03)[names(ESIS.03) == 'Subcanopy'] <- 'Cover'
ESIS.04 <- subset(ESIS, Canopy.Top > 15 & !Plant.Type %in% c('Grass/grass-like (Graminoids)', 'Fern/fern ally', 'Forb/Herb'))
ESIS.04 <- merge(ESIS.04, Com.Sp.groups[,c('Observation_ID','phase','Observation_Label','Species','Simple','Tree')], by.x = c('phase','Taxon'),  by.y = c('phase','Species') )
names(ESIS.04)[names(ESIS.04) == 'Tree'] <- 'Cover'
ESIS <- rbind(ESIS.01, ESIS.02, ESIS.03, ESIS.04)


#set the stratum breaks
ESIS.Stratum <- c('f1', 'f2', 'f3', 's1', 's2', 't1', 't2', 't3', 't4')
Macrostratum <- c('under', 'under', 'under', 'under', 'under', 'over', 'over', 'over', 'over')
#ESIS.min <- c(0,.5,1,2,4.5,13,40,80,120)*.3048
#ESIS.max <- c(.5,1,2,4.5,13,40,80,120,379)*.3048
ESIS.min <- c(0,0.1,0.25, 0.5, 1.5, 5, 15, 25, 35)
ESIS.max <- c(0.1,0.25, 0.5, 1.5, 5, 15, 25, 35, 115)


ESIS.strata <- data.frame(cbind(as.data.frame(Macrostratum, stringsAsFactors = FALSE), as.data.frame(ESIS.Stratum, stringsAsFactors = FALSE), as.data.frame(cbind(ESIS.min, ESIS.max))))
#create blank first row
phase = 'x'
Observation_ID = 'x'
Plant.Type2 = 'x'
Cover = 0 
Macrostratum = 'x'
Stratum = 'x'
ESIS.z <- cbind(as.data.frame(cbind(phase, Observation_ID, Plant.Type2, Macrostratum, Stratum),stringsAsFactors = FALSE),as.data.frame(Cover))
#loop to sum the cover by stratum and plant type
i <- 1                        
for (i in 1:9){
ESIS.x <- subset(ESIS, Canopy.Top > ESIS.strata[i,]$ESIS.min & Canopy.Bottom <= ESIS.strata[i,]$ESIS.max )
if (nrow(ESIS.x) > 0){
ESIS.y <- aggregate(list(Cover = log10(1-(ESIS.x$Cover/100.001))), by=list(phase = ESIS.x$phase, Observation_ID = ESIS.x$Observation_ID, Plant.Type2 = ESIS.x$Plant.Type2),  FUN='sum')
ESIS.y$Cover <- 100*(1-10^(ESIS.y$Cover))
ESIS.y$Macrostratum <- ESIS.strata[i,]$Macrostratum
ESIS.y$Stratum <- ESIS.strata[i,]$ESIS.Stratum
ESIS.y <- subset(ESIS.y, select=c('phase', 'Observation_ID', 'Plant.Type2', 'Macrostratum', 'Stratum',  'Cover'))
ESIS.z <- rbind(ESIS.z, ESIS.y)
}}
#insert zeroes for missing values
ESIS.z <- ESIS.z[-1,]
ESIS.zero1 <- unique(ESIS.z[,c('phase', 'Observation_ID')])
ESIS.zero2 <- unique(ESIS.z[,c('Plant.Type2', 'Macrostratum', 'Stratum')])
ESIS.zero <- merge(ESIS.zero1, ESIS.zero2)
ESIS.zero$Cover <- 0
ESIS.z <- rbind(ESIS.z, ESIS.zero)
ESIS.z <- ESIS.z %>% group_by(phase, Observation_ID, Plant.Type2, Macrostratum, Stratum) %>% summarise(Cover = max(Cover))
#get percentiles
ESIS.quantile <- ESIS.z %>% group_by(phase,Plant.Type2,Macrostratum,Stratum) %>% summarise(c15 = quantile(Cover, 0.15),c85 = quantile(Cover, 0.85))

ESIS.quantile[,c('c15','c85')] <- lapply(ESIS.quantile[,c('c15','c85')], FUN = roundF) #rounding

#rearrange quantiles
ESIS.table <- unique(ESIS.quantile[,c('phase', 'Stratum')])
ESIS.table.Forb <- subset(ESIS.quantile, Plant.Type2 %in% 'Forb', select = c(phase, Stratum, c15, c85))
colnames(ESIS.table.Forb)[3:4] <- c('Forb.min', 'Forb.max')
ESIS.table.Grass <- subset(ESIS.quantile, Plant.Type2 %in% 'Grass/Grasslike', select = c(phase, Stratum, c15, c85))
colnames(ESIS.table.Grass)[3:4] <- c('Grass.min', 'Grass.max')
ESIS.table.Shrub <- subset(ESIS.quantile, Plant.Type2 %in% 'Shrub/Vine', select = c(phase, Stratum, c15, c85))
colnames(ESIS.table.Shrub)[3:4] <- c('Shrub.min', 'Shrub.max')
ESIS.table.Tree <- subset(ESIS.quantile, Plant.Type2 %in% 'Tree', select = c(phase, Stratum, c15, c85))
colnames(ESIS.table.Tree)[3:4] <- c('Tree.min', 'Tree.max')
ESIS.table <- merge(ESIS.table,ESIS.table.Forb, by=c('phase', 'Stratum'), all.x = TRUE)
ESIS.table <- merge(ESIS.table,ESIS.table.Grass, by=c('phase', 'Stratum'), all.x = TRUE)
ESIS.table <- merge(ESIS.table,ESIS.table.Shrub, by=c('phase', 'Stratum'), all.x = TRUE)
ESIS.table <- merge(ESIS.table,ESIS.table.Tree, by=c('phase', 'Stratum'), all.x = TRUE)
ESIS.table <- ESIS.table[,c('phase', 'Stratum', 'Tree.min', 'Tree.max', 'Shrub.min', 'Shrub.max', 'Grass.min', 'Grass.max', 'Forb.min', 'Forb.max')]
write.csv(ESIS.table, 'output/ESIS.table.csv', na = "", row.names = F)

#Stratiles ----
ESIS.stratiles.sum <- aggregate(list(cover.sum = ESIS.z$Cover), by=list(phase = ESIS.z$phase, Plant.Type2 = ESIS.z$Plant.Type2, Macrostratum = ESIS.z$Macrostratum, Stratum = ESIS.z$Stratum), FUN='sum') 
ESIS.stratiles.max <- aggregate(list(cover.max = ESIS.stratiles.sum$cover.sum), by=list(phase = ESIS.stratiles.sum$phase, Plant.Type2 = ESIS.stratiles.sum$Plant.Type2, Macrostratum = ESIS.stratiles.sum$Macrostratum), FUN='max')
ESIS.stratiles <- merge(ESIS.stratiles.sum, ESIS.stratiles.max, by=c('phase', 'Plant.Type2', 'Macrostratum'))
ESIS.stratiles$RCover <- ESIS.stratiles$cover.sum/(ESIS.stratiles$cover.max+0.000001)
ESIS.Macrostratum <- ESIS
ESIS.Macrostratum$Macrostratum <- ifelse(ESIS.Macrostratum$Canopy.Top >5, 'over','under')
ESIS.Macrostratum <- aggregate(list(Cover = log10(1-(ESIS.Macrostratum$Cover/100.001))), by=list(phase = ESIS.Macrostratum$phase, Observation_ID = ESIS.Macrostratum$Observation_ID, Plant.Type2 = ESIS.Macrostratum$Plant.Type2, Macrostratum = ESIS.Macrostratum$Macrostratum),  FUN='sum')
ESIS.Macrostratum$Cover <- 100*(1-10^(ESIS.Macrostratum$Cover))
ESIS.Macrostratum.quantile <- ESIS.Macrostratum %>% group_by(phase,Plant.Type2,Macrostratum) %>% 
  summarise(c15 = quantile(Cover, 0.15),c85 = quantile(Cover, 0.85)
)
ESIS.stratiles <- merge(ESIS.Macrostratum.quantile, ESIS.stratiles, by=c('phase', 'Plant.Type2', 'Macrostratum'))
ESIS.stratiles[,c('c15','c85')] <- ESIS.stratiles[,c('c15','c85')]*ESIS.stratiles$RCover
ESIS.stratiles[,c('c15','c85')] <- lapply(ESIS.stratiles[,c('c15','c85')], FUN = roundF) #rounding

#rearrange stratiles
ESIS.table2 <- unique(ESIS.stratiles[,c('phase', 'Stratum')])
ESIS.table.Forb2 <- subset(ESIS.stratiles, Plant.Type2 %in% 'Forb', select = c(phase, Stratum, c15, c85))
colnames(ESIS.table.Forb2)[3:4] <- c('Forb.min', 'Forb.max')
ESIS.table.Grass2 <- subset(ESIS.stratiles, Plant.Type2 %in% 'Grass/Grasslike', select = c(phase, Stratum, c15, c85))
colnames(ESIS.table.Grass2)[3:4] <- c('Grass.min', 'Grass.max')
ESIS.table.Shrub2 <- subset(ESIS.stratiles, Plant.Type2 %in% 'Shrub/Vine', select = c(phase, Stratum, c15, c85))
colnames(ESIS.table.Shrub2)[3:4] <- c('Shrub.min', 'Shrub.max')
ESIS.table.Tree2 <- subset(ESIS.stratiles, Plant.Type2 %in% 'Tree', select = c(phase, Stratum, c15, c85))
colnames(ESIS.table.Tree2)[3:4] <- c('Tree.min', 'Tree.max')
ESIS.table2 <- merge(ESIS.table2,ESIS.table.Forb2, by=c('phase', 'Stratum'), all.x = TRUE)
ESIS.table2 <- merge(ESIS.table2,ESIS.table.Grass2, by=c('phase', 'Stratum'), all.x = TRUE)
ESIS.table2 <- merge(ESIS.table2,ESIS.table.Shrub2, by=c('phase', 'Stratum'), all.x = TRUE)
ESIS.table2 <- merge(ESIS.table2,ESIS.table.Tree2, by=c('phase', 'Stratum'), all.x = TRUE)
ESIS.table2 <- ESIS.table2[,c('phase', 'Stratum', 'Tree.min', 'Tree.max', 'Shrub.min', 'Shrub.max', 'Grass.min', 'Grass.max', 'Forb.min', 'Forb.max')]
write.csv(ESIS.table2, 'output/ESIS.table2.csv', na = "", row.names = F)

# 
# # DWD and Litter ----
# obs.litter <- left_join(obs, groups1)
# 
# obs.litter.pctl <- obs.litter %>% group_by(phase) %>% summarise(
#                             litter.min = quantile(Litter_Cover, 0.15, na.rm = TRUE),
#                             litter.max = quantile(Litter_Cover, 0.85, na.rm = TRUE),
#                             logs.min = quantile(Log_Cover, 0.15, na.rm = TRUE),
#                             logs.max = quantile(Log_Cover, 0.85, na.rm = TRUE),
#                             moss.min = quantile(Moss_Cover, 0.15, na.rm = TRUE),
#                             moss.max = quantile(Moss_Cover, 0.85, na.rm = TRUE),
#                             lichen.min = quantile(Lichen_Cover, 0.15, na.rm = TRUE),
#                             lichen.max = quantile(Lichen_Cover, 0.85, na.rm = TRUE),
#                          water.min = quantile(Water_Cover, 0.15, na.rm = TRUE),
#                          water.max = quantile(Water_Cover, 0.85, na.rm = TRUE),
#                          hour1.min = quantile(DWD_Hits1*0.5/Transect_Length , 0.15, na.rm = TRUE),
#                          hour1.max = quantile(DWD_Hits1*0.5/Transect_Length , 0.85, na.rm = TRUE),
#                          hour10.min = quantile(DWD_Hits2*1.75/Transect_Length , 0.15, na.rm = TRUE),
#                          hour10.max = quantile(DWD_Hits2*1.75/Transect_Length , 0.85, na.rm = TRUE),
#                          hour100.min = quantile(DWD_Hits3*6.25/Transect_Length , 0.15, na.rm = TRUE),
#                          hour100.max = quantile(DWD_Hits3*6.25/Transect_Length , 0.85, na.rm = TRUE),
#                          hour1000.min = quantile(DWD_Hits4*17.5/Transect_Length , 0.15, na.rm = TRUE),
#                          hour1000.max = quantile(DWD_Hits4*17.5/Transect_Length , 0.85, na.rm = TRUE),
#                          hour10000.min = quantile(DWD_Hits5*37.5/Transect_Length , 0.15, na.rm = TRUE),
#                          hour10000.max = quantile(DWD_Hits5*37.5/Transect_Length , 0.85, na.rm = TRUE)
# )
# nums <- unlist(lapply(obs.litter.pctl, is.numeric)) #identify all numeric columns
# obs.litter.pctl[,nums] <- lapply(obs.litter.pctl[,nums], FUN = roundF) #rounding
# 
# write.csv(obs.litter.pctl, 'output/obs.litter.pctl.csv', na = "", row.names = F)
# 
# # Total Crown Cover for ground cover ----
# Ground <- subset(ESIS, (Plant.Type2 %in% 'Tree' & Canopy.Top > 5)|(!Plant.Type2 %in% 'Tree') )
# 
# 
# Ground.cover <- aggregate(list(Cover = log10(1-(Ground$Cover/100.001))), by=list(phase = Ground$phase, Observation_ID = Ground$Observation_ID, Plant.Type2 = Ground$Plant.Type2),  FUN='sum')
# 
# Ground.cover$Cover <- 100*(1-10^(Ground.cover$Cover))
# 
# Ground.cover.pctl <- Ground.cover %>% group_by(phase, Plant.Type2) %>% summarise(
#                            cover.min = quantile(Cover, 0.15, na.rm = TRUE),
#                            cover.mid = quantile(Cover, 0.5, na.rm = TRUE),
#                          cover.max = quantile(Cover, 0.85, na.rm = TRUE))
# Ground.cover.pctl$foliar.min <- tofoliar(Ground.cover.pctl$cover.min)
# Ground.cover.pctl$foliar.max <- tofoliar(Ground.cover.pctl$cover.max)
# 
# Ground.cover.pctl[,c('cover.min', 'cover.mid', 'cover.max', 'foliar.min', 'foliar.max')] <- lapply(Ground.cover.pctl[,c('cover.min', 'cover.mid', 'cover.max', 'foliar.min', 'foliar.max')], FUN = roundF) #rounding
# 
# write.csv(Ground.cover.pctl, 'output/Ground.cover.pctl.csv', na = "", row.names = F)

# # snags ----
# Snags <- subset(obsspp, grepl('snag', Genus)|grepl('<snag', Genus), select = c('Observation_ID','soilplot', 'Genus', 'BA', 'Dmin', 'Dmax'))
# Snags <- left_join(Snags, obs[,c('Observation_ID', 'soilplot')])
# Snags$D <- ifelse(!is.na(Snags$Dmax),ifelse(!is.na(Snags$Dmin), (Snags$Dmin+Snags$Dmax)/2, Snags$Dmax), 20)
# Snags$D.area <- (Snags$D/200)^2*3.141592
# Snags$Snags.Ha <- Snags$BA/Snags$D.area
# 
# Snags.pctl <- ddply(Snags, c('phase'), summarise,
#                     Snags.min = quantile(Snags.Ha, 0.15, na.rm = TRUE),
#                     Snags.mid = quantile(Snags.Ha, 0.5, na.rm = TRUE),
#                     Snags.max = quantile(Snags.Ha, 0.85, na.rm = TRUE))
# Snags.pctl[,c('Snags.min', 'Snags.mid', 'Snags.max')] <- lapply(Snags.pctl[,c('Snags.min', 'Snags.mid', 'Snags.max')], FUN = roundF) #rounding
# 
# write.csv(Snags.pctl, 'output/Snags.pctl.csv', na = "", row.names = F)
# 
