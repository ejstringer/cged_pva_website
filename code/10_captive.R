
tid <- read.csv('./data/tidbinbilla_colony.csv')

tid$Age

tid %>% 
  separate(col = Age, into = c('year', 'month', 'day'), sep = ',', remove = F) %>% 
  mutate(day = substr(trimws(day), 1,3),
         day = as.numeric(sub('D', '', day)),
         month = as.numeric(sub('M', '', month)),
         year = as.numeric(sub('Y', '', year)),
         days.alive = day+(month*30)+(year*365),
         status = case_when(
           grepl('Release', Status) ~ 'released',
           grepl('Dead', Status) ~ 'dead',
           grepl('Lost', Status) ~ 'lost',
           .default = Status
         )) %>% 
  filter(status != 'lost',
         status != 'released') %>% 
  select(days.alive, status) %>% 
  filter(days.alive < (365*6)) %>% 
  mutate(survival = ifelse(status == 'dead', 0, 1)) -> tid_colony


m<-glm(survival ~ days.alive, data = tid_colony, family = 'binomial')

x_days <- 0:2655

pred.m <- predict(m, newdata=data.frame(days.alive = x_days), type = 'response')

pred.surv <- data.frame(days.alive = x_days, 
           survival = pred.m) %>% 
  mutate(year = round(days.alive/365, 3),
         years.alive = days.alive/365,
         months.alive = years.alive*12)

pred.surv %>% 
  filter(year %in% c(1:7)) %>% 
  mutate(stage = c('J', paste0('A', 1:6)))

tid_colony %>% 
  mutate(years.alive = days.alive/365,
         months.alive = years.alive*12) %>% 
ggplot(aes(y = survival, x = months.alive))+
  geom_point()+
  theme_classic()+
  scale_x_continuous(breaks = seq(0,100, 12))+
  geom_line(data = pred.surv, colour = 'grey', lwd = 1)
