
caps <- read.csv('../em-pva/data/ged_monitoring_grids2002-2023.csv')

caps %>% filter(year == 2013) %>% 
  mutate(age2 = ifelse(age == 'J', 'J', 'A')) %>% 
  group_by(age2) %>% 
  summarise(n = n())

37/(37+69)
0.35
