

# load --------------------------------------------------------------------

source('./code/00_libraries.R')
area_used <- read.csv('./output/area_used.csv')
n_sites <- read.csv('./output/site_abundance.csv')

# starting densities ------------------------------------------------------

n_sites %>% 
  filter(year == 2013) %>% 
  left_join(area_used) %>%
  mutate(N_ha = N/n,
         lcl_ha = lcl/n,
         ucl_ha = ucl/n,
         year = factor(year)) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  dplyr::select(-buffer, ) %>% 
  flextable() %>% 
  theme_alafoli() %>% 
  hline_top(border = fp_border_default(width = 2))
  
N_abundances <- n_sites %>% 
  filter(year == 2013) %>% 
  left_join(area_used)%>%
  dplyr::select(-buffer) %>% 
  rename(n_grids = n) 

fxtb_N <- N_abundances %>% 
  pivot_longer(cols = N:se_sum2) %>%
  mutate(value_ha = value/n_grids,
         year = as.character(year)) %>% 
  mutate(N_all = ifelse(name %in% c('N', 'lcl', 'ucl'),
                        round(value_ha*area_ha), NA)) %>%
  filter(name != 'se_sum2') %>% 
  mutate_if(is.numeric, round, 2)  %>%
  mutate(n_grids = ifelse(duplicated(site), '', n_grids)) %>% 
  mutate(across(c(year, site, area_ha), ~ifelse( duplicated(.x), '', .x))) %>%
  relocate(area_ha, .before = N_all) %>% 
  flextable() %>% 
  theme_alafoli() %>% 
  hline_top(border = fp_border_default(width = 2)) %>% 
  hline(i = seq(4,13,4),
        j = 2:8,
        border = fp_border_default(width = 1.4, style = 'dashed')) %>% 
  hline_bottom(border = fp_border_default(width = 2));fxtb_N


n_real <- n_sites %>% 
  left_join(area_used)%>%
  dplyr::select(-buffer) %>% 
  rename(n_grids = n) %>%
  mutate(N_all = (N/n_grids)*area_ha,
         lcl_all = (lcl/n_grids)*area_ha,
         ucl_all = (ucl/n_grids)*area_ha,
         site = factor(site))  
levels(n_real$site) <- c('CA', 'JE', 'JW', 'MA')
write.csv(n_real, './output/real_abundance_for_comparison.csv', row.names = F)
