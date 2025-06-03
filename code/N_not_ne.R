ne <-read.csv('./data/Ne_estimates_Ryan.csv') %>%
  separate(Cohort, into = c('site', 'year'), sep = '-') %>% 
  rename(Ne = Ne..LDNe.0.05.) %>% 
  mutate(site = factor(site),
         year = as.numeric(year),
         year = year-0) %>% 
  select(site, year, Ne,Sample.Size) %>% 
  as_tibble()

levels(ne$site) <- levels(factor(n_sites$site))


n_sites$site <- factor(n_sites$site)

n_both <- left_join(n_sites, ne)

n_both %>% #filter(site != 'JerraWest') %>% 
ggplot(aes(N, Ne, colour = site))+
  geom_point()+
  theme_classic()+
  geom_smooth(aes(group = 'no group'), se = F, method = 'lm')
