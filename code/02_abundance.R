
# load --------------------------------------------------------------------

source('./code/00_libraries.R')
caps <- read.csv('./data/N_estimates.csv')

# adult N -----------------------------------------------------------------

n_start <- caps %>% filter(year > 2012, site %in% c('Majura', 'JerraWest', 
                                                   'JerraEast', 'Cookanalla'),
                          !grepl('L', grid.id))


ngrids <- table(n_start$site[!duplicated(n_start$grid)])
n_sites <- n_start %>% group_by(year, site) %>% 
  summarise(N = (sum(N)),
            lcl = (sum(lcl)),
            ucl = (sum(ucl)),
            se_sum = sum(se),
            se_sum2 = sqrt(sum(se^2)),
            n = n())


grid.labs <- c("Grid one", "Grid two","Grid three", "Grid four")
names(grid.labs) <- c(1,2,3,4)

site.labs <- c("Aerial", "Airport"   ,    "Bonshaw"  ,     "Cookanalla"  ,  "Jerra East" ,    "Jerra West"  ,   "Majura" )
names(site.labs) <- c("AerialPaddock" ,"Airport"   ,    "Bonshaw"  ,     "Cookanalla"  ,  "JerraEast" ,    "JerraWest"  ,   "Majura" )

abundanceFig <- ggplot(n_start, aes(x = year, y = N, group = site))+
  geom_errorbar(aes(ymin = lcl, ymax = ucl), colour = "grey", width = 0) +
  geom_line()+
  geom_point(aes(colour = site), size = 2)+
  geom_point(shape = 1, size = 2, colour = "black")+
  facet_grid(site ~ grid.n, scales = "free",labeller=labeller(grid.n = grid.labs, site = site.labs))+
  scale_x_continuous(breaks = seq(2002, 2022, 2), 
                     labels = c(2002, '', 2006, '',2010, '', 2014, '', 2018, '', 2022 )) +
  labs(x="Year", y="Abundance (N)") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size=8, face = "bold" ),
        strip.text.y = element_text(size=9, face = "bold"),
        strip.background = element_blank(),
        axis.line.y = element_line(colour = 'black', linewidth = 1),
        legend.position = "none");abundanceFig


abundanceFig2 <- ggplot(n_sites, aes(x = year, y = N, group = site))+
  geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = "grey") +
  geom_line(show.legend = F)+
  geom_point(data = filter(n_sites, year == 2013), 
               aes(fill = factor(year)), pch = 21, colour = 'yellow',
             alpha = 0.8, size = 5)+
  geom_point(aes(colour = site), size = 2, show.legend = F)+
  geom_point(shape = 1, size = 2, colour = "black", show.legend = F)+
  scale_fill_manual(values = c('yellow'), 
                      name = 'Simulation Start')+
  facet_grid(~site, scales = "free",labeller=labeller(site = site.labs))+
  scale_x_continuous(breaks = seq(2002, 2022, 2)) +
  labs(x="Year", y="Abundance (N)") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'inside',
        legend.position.inside = c(0.12,0.9),
        legend.background = element_rect(colour = 'grey'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size=12, face = "bold" ),
        strip.background = element_blank(),
        axis.line.y = element_line(colour = 'black', linewidth = 1));abundanceFig2


write.csv(n_sites, './output/site_abundance.csv', row.names = F) 

