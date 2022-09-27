simplecover <- Com.Sp.mean %>% group_by(soilplot, Simple) %>% summarise(cover = 100*(1-10^(sum(log10(1-(Total/100.001))))))

exotics <- c(as.character(subset(listspp, Eastern.North.America %in% "X")$AcTaxon),  'Phalaris arundinacea') %>% unique()
exoticover <- Com.Sp.mean %>% mutate(exotic = ifelse(Species %in% exotics, 'Exotic','Native'))
exoticover <- exoticover %>% group_by(soilplot, exotic) %>% summarise(cover = 100*(1-10^(sum(log10(1-(Total/100.001))))))
exoticover <- subset(exoticover, exotic %in% 'Exotic') %>% mutate(Exotic=cover)
wetness <- Com.Sp.mean %>% left_join(unique(listspp[,c('AcTaxon','Wetness')]), by=c('Species'='AcTaxon'))
wetness <- wetness %>% group_by(soilplot) %>% summarise(wetness= weighted.mean(Wetness, w=Total, na.rm=TRUE))

functionalcover <- makecommunitydataset(as.data.frame(simplecover), row='soilplot', column = 'Simple', value =  'cover')
functionalcover <- functionalcover %>% mutate(soilplot = rownames(functionalcover))
functionalcover <- functionalcover %>% left_join(exoticover[c('soilplot', 'Exotic')]) %>% 
  left_join(overstorycover[c('soilplot', 'overstorycover')]) %>% 
  left_join(wetness)
functionalcover <- functionalcover %>% mutate(Exotic = ifelse(is.na(Exotic), 0, Exotic))

rownames(functionalcover) <- functionalcover$soilplot
functionalcover <- subset(functionalcover, select= -soilplot)
