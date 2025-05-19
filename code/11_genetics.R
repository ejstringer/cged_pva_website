

# genetics ----------------------------------------------------------------

library(dartRverse)

ged2 <- readRDS('../em-pva/output/ged_genlight_filtered.rds')
# wild  -------------------------------------------------------------------

wild <- gl.keep.pop(ged2, pop.list = c('Cookanalla','Jerrabomberra East',
                                       'Jerrabomberra West','Majura Training Area (MTA)'))

wild2016 <- gl.keep.pop(wild, as.pop = 'year', pop.list = 2016)

gl.report.diversity(wild2016)

table(pop(wild2016))

sapply(pops, function(x) mean(gl.He(x)))

pops <- seppop(wild2016)

table(factor(wild@other$ind.metrics$pop),
      wild@other$ind.metrics$year)
