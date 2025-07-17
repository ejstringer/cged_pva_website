satelite <- read.csv('../em-pva/data/ged_grids_satellite2000-2020.csv')
satelite$LST_scales <- scale(satelite$LST)
satelite$NDVI_scales <- scale(satelite$NDVI)

phiASA  <-  param_selecteddf$mean[grepl('survival_ASA', param_selecteddf$name)]
phiJ  <-  param_selecteddf$mean[grepl('survival_J', param_selecteddf$name)]


mean_sat <- satelite %>% group_by(new_site) %>% 
  filter(new_site %in% c('Cookanalla', 'Jerra East', 'Jerra West', 'Majura'),
         grid_id != 'CA-03_L', year > 2012) %>% 
  summarise(ndvi = mean(NDVI_scales, na.rm = T), 
            lst = mean(LST_scales, na.rm = T),
            both = mean(NDVI_scales+LST_scales, na.rm = T),
            se = sd(NDVI_scales+LST_scales, na.rm = T)/sqrt(n()),
            n = n()) %>% 
  mutate(sur = phiASA,
         surJ = phiJ,
         logitsur = logit(surJ),
         lwr = both - se*1.96,
         upr = both + se*1.96) 
  ggplot(mean_sat,aes(ndvi, lst, colour = new_site))+
  geom_point(aes(size = sur))+
    theme_classic()

  
  ggplot(mean_sat,aes(lst+ndvi, logitsur, colour = new_site))+
    geom_point(size = 5)+
    geom_point(aes(y = sur))+
    theme_classic()+
    #geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0)+
    geom_smooth(method = 'lm', colour = 'black')+
    xlab('NDVI + LST (scaled)')+
    ylab('Adult Survival')+
    labs(colour = 'Grassland')+
    theme(legend.position = 'inside',
          legend.position.inside = c(0.2, 0.2),
          legend.background = element_rect(colour = 'grey'))

  cor(mean_sat$logitsur, mean_sat$lst)
  
  cor(mean_sat$logitsur, mean_sat$both)

  
  cor(mean_sat$logitsur, mean_sat$ndvi)
  lm(logitsur ~ both, data = mean_sat) %>% summary
  lm(logitsur ~ ndvi + lst, data = mean_sat) %>% summary
  
  x <- seq(-0.5,0.5, 0.001)
  y <- predict(m, newdata = data.frame(both = x))
  
  survival_habitat<-  data.frame(ndvi = x/2, lst = x/2, both = x, logitsur = y, sur = unlogit(y),
             lt = 'dashed') %>%
    mutate(lt = ifelse(both > min(mean_sat$both),
                      'dashed', lt),
           lt = ifelse(both > max(mean_sat$both), 'dotted', lt)) %>% 
    ggplot(aes(lst+ndvi, sur))+
    geom_line(aes(linetype = lt), show.legend = F)+
    geom_point(data = mean_sat, aes(colour = new_site), size = 2)+
    scale_linetype_manual(values = c(2,2,1))+
    theme_bw()+
    xlab('NDVI + LST (scaled)')+
    ylab('Adult Survival')+
    geom_hline(yintercept = c(0,1), colour = 'grey', linewidth = 0.25)+
    labs(colour = 'Grassland')+
    theme(legend.position = 'inside',
          legend.position.inside = c(0.85, 0.7),
          legend.background = element_rect(colour = 'grey'),
          panel.grid = element_blank());survival_habitat
  
  ggsave('./figures/habitat_survival_relationship2.png', plot = survival_habitat,
         width = 5, height = 3)


# extra -------------------------------------------------------------------

    
  satelite$both_scaled <- as.numeric(satelite$LST_scales + satelite$NDVI_scales)
  sat_dat <- satelite %>% 
    filter(new_site %in% c('Cookanalla', 'Jerra East', 'Jerra West', 'Majura'),
           grid_id != 'CA-03_L', year > 2012) %>% 
    left_join(data.frame(new_site = c('Cookanalla', 'Jerra East', 'Jerra West', 'Majura'),
                         sur = phiASA, surJ = phiJ)) %>%
    mutate(logitsur = logit(sur)) 
   
   
    ggplot(sat_dat,aes(both_scaled, logitsur, colour= new_site))+
    geom_point(alpha = 0.5)+
    geom_point(data = mean_sat, aes(x = both, fill = new_site), size = 4,
               colour = 'black', pch = 21)+
    geom_smooth(method = 'lm', colour = 'black')+
    theme_classic()+
    xlab('NDVI + LST (scaled)')+
    ylab('Adult Survival')+
    labs(colour = 'Grassland')
    
  sapply(1:5, function(x) plot(mm, which = x))

  mm <- lm(logitsur ~ LST_scales + NDVI_scales, data = sat_dat)
  summary(mm)
  
  
  lmer(both_scaled~ logitsur + (1|new_site), data = sat_dat) %>% summary
  mmboth <- lm(logitsur ~ both_scaled, data = sat_dat) 
x <- seq(-4,4, 0.001)
y <- predict(mmboth, newdata = data.frame(both_scaled = x), interval = "confidence")
yCI <- predict(mmboth, newdata = data.frame(both_scaled = x), )


data.frame(ndvi = x/2, lst = x/2, both_scaled = x, logitsur = y[,1], 
           lwr = y[,2],
           upr = y[,3],
           lt = 'dashed') %>%
  mutate(sur = unlogit(logitsur),
         lwru = unlogit(lwr),
         upru = unlogit(upr)) %>% 
  # mutate(lt = ifelse(both > min(mean_sat$both),
  #                    'solid', lt),
  #        lt = ifelse(both > max(mean_sat$both), 'dotted', lt)) %>% 
  ggplot(aes(both_scaled, sur))+
  geom_line()+
  geom_ribbon(aes(ymin = lwru, ymax = upru), alpha = 0.5, fill = 'grey80',
              colour = NA)+
  geom_point(data = sat_dat, aes(colour = new_site), size = 2)+
  scale_linetype_manual(values = c(2,2,1))+
  theme_bw()+
  xlab('NDVI + LST (scaled)')+
  ylab('Adult Survival')+
  labs(colour = 'Grassland')+
  theme(legend.position = 'inside',
        legend.position.inside = c(0.8, 0.8),
        legend.background = element_rect(colour = 'grey'),
        panel.grid = element_blank())+
  geom_hline(yintercept = c(0,1), colour = 'grey', linewidth = 0.25)





satelite %>% 
  filter(new_site %in% c('Cookanalla', 'Jerra East', 'Jerra West', 'Majura'),
         grid_id != 'CA-03_L', year > 2012) %>% 
  group_by(new_site, year) %>% 
  summarise(ndvi = mean(NDVI_scales, na.rm = T), 
            lst = mean(LST_scales, na.rm = T)) %>% 
  ggplot(aes(ndvi,  lst, colour = new_site))+
  geom_point()+
  facet_wrap(~new_site)+
  geom_point(data = mean_sat, aes(y = lst, x = ndvi, fill = new_site,
                                  size = sur),
             shape = 21, colour = 'black')+
  theme_classic()+
  ggforce::geom_mark_hull(
    aes(
      fill = new_site,
      label = new_site
    ),
    concavity = 2.5
    # if you get error, install "concaveman" package
  )
  ggforce::geom_mark_ellipse(aes(fill = factor(new_site)))

mean_sat %>% ggplot(aes(ndvi+lst, sur, colour = new_site))+
  geom_point(size = 3)+
  theme_classic()

mean_sat %>% ggplot(aes(lst, sur, colour = new_site))+
  geom_point(size = 3)+
  theme_classic()

mean_sat %>% ggplot(aes(ndvi, lst, colour = new_site))+
  geom_point(size = 3)+
  theme_classic()


mean_satyear <- satelite %>% group_by(new_site, year) %>% 
  filter(new_site %in% c('Cookanalla', 'Jerra East', 'Jerra West', 'Majura'),
         grid_id != 'CA-03_L', year > 2012) %>% 
  summarise(ndvi = sum(NDVI_scales < 0, na.rm = T)/n(), 
            lst = sum(LST_scales < 0, na.rm = T)/n(),
            both = sum(LST_scales < 0 & NDVI_scales < 0, 
                       na.rm = T)/n()) %>% 
  left_join(data.frame(new_site = c('Cookanalla', 'Jerra East', 'Jerra West', 'Majura'),
                       sur = phiASA, surJ = phiJ))
  ggplot(mean_satyear,aes(ndvi, sur, colour = new_site))+
  geom_point(size = 3)+
  theme_classic()
mean_sat %>% ggplot(aes(ndvi, sur, colour = new_site))+
  geom_point(size = 3)+
  theme_classic()

lm(mean_sat$sur ~ mean_sat$lst) %>% summary

satelite %>% group_by(new_site, year) %>% 
  filter(new_site %in% c('Cookanalla', 'Jerra East', 'Jerra West', 'Majura'),
         grid_id != 'CA-03_L', year > 2012) %>% 
  summarise(ndvi = mean(NDVI_scales, na.rm = T), 
            lst = mean(LST_scales, na.rm = T)) %>% 
  left_join(data.frame(new_site = c('Cookanalla', 'Jerra East', 'Jerra West', 'Majura'),
                       sur = phiASA, surJ = phiJ)) %>% 
  ggplot(aes(year, ndvi, colour = new_site, group = new_site))+
  geom_point(size = 3)+
  geom_line()+
  theme_classic()

cor(sat_dat$logitsur[complete.cases(sat_dat$both_scaled)],
    sat_dat$both_scaled[complete.cases(sat_dat$both_scaled)])
