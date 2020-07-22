library(soilDB)
library(aqp)
library(plyr)
site <- get_site_data_from_NASIS_db()
plot <- get_vegplot_from_NASIS_db()
SI <- get_vegplot_tree_si_summary_from_NASIS_db()

SI <- merge(SI, plot[,c('vegplotiid', 'legacysoilcompname')], by='vegplotiid')
SI <- subset(SI, !is.na(siteindexplotave))
SI.quantile <- ddply(SI, c('plantsym', 'plantsciname', 'legacysoilcompname'), summarise,
                     si.min = min(siteindexplotave),
                     si.15 = quantile(siteindexplotave, 0.15),
                     si.mean = mean(siteindexplotave),
                     si.85 = quantile(siteindexplotave, 0.85),
                     si.max = max(siteindexplotave),
                     si.count = length(siteindexplotave)
)

SI.quantile.species <- ddply(SI, c('plantsym', 'plantsciname'), summarise,
                     si.min = min(siteindexplotave),
                     si.15 = quantile(siteindexplotave, 0.15),
                     si.mean = mean(siteindexplotave),
                     si.85 = quantile(siteindexplotave, 0.85),
                     si.max = max(siteindexplotave),
                     si.count = length(siteindexplotave)
)
SI.quantile.soil <- ddply(SI, c('legacysoilcompname'), summarise,
                     si.min = min(siteindexplotave),
                     si.15 = quantile(siteindexplotave, 0.15),
                     si.mean = mean(siteindexplotave),
                     si.85 = quantile(siteindexplotave, 0.85),
                     si.max = max(siteindexplotave),
                     si.count = length(siteindexplotave)
)

write.csv(SI.quantile, 'output/SI.quantile.csv')
write.csv(SI.quantile.species, 'output/SI.quantile.species.csv')
write.csv(SI.quantile.soil, 'output/SI.quantile.soil.csv')