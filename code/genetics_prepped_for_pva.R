# Dispersal
library(tidyverse)
library(dartRverse)
library(otuSummary)
# genetics 

ged2 <- readRDS('./data/ged_genlight_filtered.rds')

# wild  -------------------------------------------------------------------

wild <- gl.keep.pop(ged2, pop.list = c('Cookanalla','Jerrabomberra East',
                                       'Jerrabomberra West','Majura Training Area (MTA)'))
wild@pop %>% table
wild@other$ind.metrics$grid %>% factor %>% table

keep_index <- grepl('_', wild@other$ind.metrics$grid)
wild <- gl.keep.ind(wild, ind.list = wild@ind.names[keep_index])

wild@other$ind.metrics$grid %>% factor %>% table

pop(wild) <- wild@other$ind.metrics$grid

## duplicates ---------------------------------------------------------

pop(wild) <- wild@other$ind.metrics$year
wildyear <- seppop(wild)

gl.download.binary(software = 'emibd9', out.dir = getwd())
rel <- lapply(wildyear, gl.run.EMIBD9)

rel[[1]]$rel
reldf <- lapply(rel, function(x) matrixConvert(x$rel, c('id1', 'id2', 'r'))) %>% 
  do.call('rbind', .) %>% 
  mutate(year = substr(rownames(.), 1,4)) %>% 
  filter(r > 0.4)

reldf$id1 %>% duplicated %>% table
reldf$id2 %>% duplicated %>% table
reldf$id2[duplicated(reldf$id2)]

reldf$id1 %in% reldf$id2

wild_unique <- gl.drop.ind(wild, ind.list = reldf$id2) 
saveRDS(wild_unique,'./output/wild_unique_genlight_tympos.rds')
