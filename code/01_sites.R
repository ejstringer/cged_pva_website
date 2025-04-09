
# load --------------------------------------------------------------------
source('./code/00_libraries.R')
grid1 <- vect('./data/grasslands/one_grid.shp')
grasslands <- vect("./data/grasslands/grasslands_no_extras.shp") 
ntg <- vect('./data/grasslands/NTG_5properties.shp')
exotic <- vect('./data/grasslands/exotic_map.shp')


# buffer around ntg (metres)
buffer_metres <- 100

# site size ---------------------------------------------------------------
## fix layers
exotic <- aggregate(exotic, by = 'COMMUNITY')
exotic <- buffer(exotic, 1)

grass<- grasslands
grass$land <- 'grassland'
grass <- aggregate(grass, by = 'land')
is.valid(grass)

## buffer ---------
ntg_buffer <- buffer(ntg, width=buffer_metres)

## crop ----------
ntg_buffer <- crop(ntg_buffer, grass)
ntg_buffer <- st_difference(st_as_sf(ntg_buffer), st_as_sf(exotic))

my_buffer <- vect(ntg_buffer)


# site map ----------------------------------------------------------------
plot_buffer <- sf::st_cast(st_as_sf(my_buffer), "MULTIPOLYGON")


sitemap <- ggplot()+
  geom_sf(data = st_as_sf(grasslands), colour = 'grey90')+
  geom_sf(data = plot_buffer, aes(fill = Proprty), alpha = 0.7)+
  geom_sf(data = st_as_sf(ntg),aes(fill = Proprty),
          colour = 'grey50')+
  scale_fill_manual(values = adegenet::virid(7)[2:6],
                    name = 'Sites',labels = c('Bonshaw','Cookanalla',
                                              'Jerra East', 'Jerra West',
                                              'Majura'))+
  ggtitle(paste('Dragon area use: NTG with', buffer_metres, 'buffer'))+
  theme_bw()+
  theme(legend.position = 'inside',
        legend.position.inside = c(0.2,0.85),
        legend.background = element_rect(colour = 'grey'))

sitemap

# area --------------------------------------------------------------------
propArea <-  c(expanse(my_buffer,unit = 'ha'),
               expanse(ntg, unit = 'ha'),
               expanse(grid1, unit = 'ha'))

dfarea <- data.frame(site = c(my_buffer$Proprty, ntg$Proprty, 'Grid'), 
                     area_ha = propArea, 
                     buffer = c(rep(T, 5), rep(F,6))) %>% 
  mutate(site = ifelse(site == 'Cookanela', 'Cookanalla', site),
         site = ifelse(site == 'MTA', 'Majura', site))

fxtb<- dfarea %>% mutate(area_of = ifelse(buffer, 
                                          paste('NTG buffer with',
                                                buffer_metres, 'm (ha)'),
                                          'NTG (ha)'),
                  area_of = ifelse(site == 'Grid', 
                                   'monitoring grid (ha)', area_of),
                  area_ha = round(area_ha, 2)) %>% 
  dplyr::select(-buffer) %>% 
  pivot_wider(values_from = area_ha, names_from = area_of) %>% 
  #filter(property != 'Bonshaw') %>% 
  flextable %>% autofit() %>% 
  theme_alafoli() %>% 
  bold(part  = 'header') %>% 
  italic(i = c(1)) %>% 
  hline(i = c(5),
        border = fp_border_default(width = 1.5, 
                                   style = 'dashed')) %>% 
  hline_bottom(border = fp_border_default(width = 2)) %>% 
  hline_top(border = fp_border_default(width = 2))

fxtb
area_used <- dfarea %>% filter(site != 'Bonshaw',
                               # property != 'grid',
                               buffer)

# area_used$area_ha[4] <- area_used$area_ha[4]*0.5 # adjusting MA ntg

write.csv(area_used, './output/area_used.csv', row.names = F)
