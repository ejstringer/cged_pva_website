library(survival)
# zims data -----------------------------------------------
tid <- read.csv('./data/tidbinbilla_colony.csv')
mel <- read.csv('./data/melbourne_colony.csv')

tid$Age
mel$Status %>% table

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

table(tid_colony$status[(tid_colony$days.alive/365)<1])

table(round(tid_colony$time), tid_colony$status) %>% 
  as.matrix %>% 
  as.data.frame() %>%
  pivot_wider(names_from = Var2, values_from = Freq) %>% 
  mutate(Var1 = as.numeric(as.character(Var1)),
         Var1 = ifelse(Var1 > 9, 9, Var1)) %>% 
  group_by(Var1) %>%
  summarise(Alive = sum(Alive),
            Dead = sum(dead)) %>% 
  mutate(All = Alive + Dead,
         survival = Alive/All,
         Age = paste0(Var1, '-', Var1 +1, 'yrs'))

pred.surv %>% 
  filter(year %in% c(1:7)) %>% 
  mutate(stage = c('J', paste0('A', 1:6)),
         survival_calc = 0.89^c(1,2,3,4,5,6,7),
         surv_diff = abs(survival - survival_calc))

tid_colony %>% 
  mutate(years.alive = days.alive/365,
         months.alive = years.alive*12) %>% 
ggplot(aes(y = survival, x = months.alive))+
  geom_point()+
  theme_classic()+
  scale_x_continuous(breaks = seq(0,100, 12))+
  geom_line(data = pred.surv, colour = 'grey', lwd = 1)



# survival analysis -------------------------------------------------------
#https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/
tid_colony$time <- tid_colony$days.alive/365
tid_colony$event <- ifelse(tid_colony$survival == 0, 1,0)
tid_colony$age <- ifelse(tid_colony$days.alive > (365*3), 'old', 'young')

km.curve <- with(tid_colony, Surv(time, event))
head(km.curve,80)


km_fit <- survfit(Surv(time, event) ~ 1, data=tid_colony)
summary(km_fit, times = c(1:7))


plot(km_fit, xlab="Years", main = 'Kaplan Meyer Plot')

km_AG_fit <- survfit(Surv(time, event) ~ age, data=tid_colony)
summary(km_AG_fit, times = c(1:5), age = 'old')
plot(km_AG_fit, xlab="Years", main = 'Kaplan Meyer Plot', col = 1:2)

# data from ryan ----------------------------------------------------------

library(tidyverse)
library(lubridate)
ged <- read.csv('./data/GED_Pedigree_for_Emily.csv') 
ged[grep('or 22', ged$sire.ID),]

ged$Notes[grep('or 22', ged$sire.ID)] <- paste(ged$Notes[grep('or 22', ged$sire.ID)], '; or sire.id = ged.id 22')
ged$sire.ID[grep('or 22', ged$sire.ID)] <- '105'

## fix dates ---------------------------------------------------------------
em.date.fix <- function(xdate){
  
  if(grepl('-', xdate)){
    xmonth <- which(sapply(month.abb, function(x) grepl(x, xdate)))
    xdatetidy <- gsub('-', '/', sub(month.abb[xmonth], xmonth, xdate))
    xdateformat<- as.Date(xdatetidy, tryFormats = c('%d/%m/%y'))
  }else{
    xdateformat <- as.Date(xdate, tryFormats = '%d/%m/%Y')
  }
  return(xdateformat)
}




ged$lay.date <- do.call('c', lapply(ged$lay.date, em.date.fix))

ged$hatch.date <- do.call('c', lapply(ged$hatch.date, em.date.fix))

## fix missing clutch ids --------------------------------------------------

ged %>% filter(is.na(clutch.ID) & complete.cases(sire.ID))
ged_miss <- ged[is.na(ged$clutch.ID) & complete.cases(ged$sire.ID),]

new_clutch_ids <- paste(ged_miss$colony, ged_miss$dam.ID, ged_miss$year, 
                        'x',
                        sep = '-')
new_clutch_ids[12:17] <- paste0(new_clutch_ids[12:17], 'x')

ged$clutch.ID[is.na(ged$clutch.ID) & complete.cases(ged$sire.ID)] <- new_clutch_ids

ged[169:188,]
ged$clutch.ID[183:188]
# 212 likely from same clutch as 204 - parents 177 and 188
# 215 likely from same clutch as 213 214 and 216 - parents 130 and 195
## explore -----------------------------------------------------------------

names(ged)
head(ged)
index <- paste(ged$clutch.ID, ged$egg.no) %>% duplicated
ged[index,]

ged %>% filter(sex == 'F', 
               complete.cases(clutch.ID)) %>%
  group_by(dam.ID) %>% 
  summarise(clutches = n()) %>% arrange(clutches) 

table(ged$Female.Clutch.Count, ged$year)
ged$egg.no %>% boxplot

filter(ged, dam.ID == 'B90783')


# data tbls ---------------------------------------------------------------


## dragon tbl --------------------------------------------------------------

dragons <- ged %>% 
  select(GED.ID, colony, sex, year, clutch.ID, hatch.date, Origin, Notes) %>% 
  rename(dragon.notes = Notes) %>% 
  filter(complete.cases(GED.ID)) %>% 
  rename_all(tolower)


## clutch tbl --------------------------------------------------------------
clutches <- ged %>% 
  filter(complete.cases(clutch.ID)) %>% 
  group_by(clutch.ID, dam.ID, sire.ID, lay.date, clutch.size) %>% 
  summarise(number.hatched = sum(hatch.success))%>% 
  rename_all(tolower)

clutches[grep(195, clutches$sire.id),]
## egg tbl -----------------------------------------------------------------
eggs <- ged %>% 
  select(clutch.ID, egg.no, hatch.success, GED.ID, Notes) %>% 
  rename(egg.notes = Notes) %>% 
  filter(complete.cases(egg.no))%>% 
  rename_all(tolower)


## AA tbl ------------------------------------------------------------------
AAnumber <- ged %>% select(GED.ID, AA_Number) %>% 
  filter(complete.cases(AA_Number))%>% 
  rename_all(tolower)




# save csvs ------------
#dir.create('./data/colony')

write.csv(dragons, './data/colony/dragons.csv', row.names = FALSE)
write.csv(clutches, './data/colony/clutches.csv', row.names = FALSE)
write.csv(eggs, './data/colony/eggs.csv', row.names = FALSE)
write.csv(AAnumber, './data/colony/AAsamples.csv', row.names = FALSE)


# assess births -----------------------------------------------------------

births <- dragons %>% filter(complete.cases(hatch.date)) %>% 
  left_join(clutches)
head(births)
births$hatch.date

momsbirth<-clutches %>% 
  left_join(dragons[,c('ged.id', 'hatch.date')], 
            by = c('dam.id' = 'ged.id')) %>% 
  #filter(complete.cases(hatch.date), complete.cases(lay.date)) %>% 
  filter(!grepl('UC', clutch.id)) %>% 
  mutate(age = lay.date - hatch.date,
         years = as.numeric(age/365),
         year = year(lay.date),
         months = years*12,
         wholeyear = round(years),
         not_hatched = clutch.size-number.hatched,
         per_hatched = number.hatched/clutch.size) 
momsbirth$years %>% max(., na.rm = T)
plot(momsbirth$years, momsbirth$number.hatched)
momsbirth$hatch.date %>% table
table(round(momsbirth$years,2))
boxplot(momsbirth$per_hatched~momsbirth$wholeyear)
hatched_age_model <- aov(number.hatched ~ factor(wholeyear), data = momsbirth)
summary(hatched_age_model)
TukeyHSD(hatched_age_model)

momsbirth %>%
  group_by(wholeyear) %>% 
  summarise(hatched = mean(number.hatched),
            sd = sd(number.hatched)) %>% 
  mutate(lcl = hatched - (sd*2/n()),
         ucl = hatched + (sd*2)/n())%>% 
ggplot(aes(wholeyear, hatched))+
  #geom_bar(stat = 'identity')
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.1)+
  geom_point(size = 2)+
  ylim(0,6)

barplot(table(momsbirth$number.hatched))
glm(number.hatched ~ (wholeyear)+I(wholeyear^2), 
    data = momsbirth, family ='poisson') %>% summary

#install.packages('pscl')
library(pscl)

summary(m1 <- zeroinfl(number.hatched ~ factor(wholeyear), data = momsbirth))

summary(m1 <- zeroinfl(number.hatched ~ wholeyear + I(wholeyear^2), data = momsbirth))

mnull <- update(m1, . ~ 1)

pchisq(2 * (logLik(m1) - logLik(mnull)), df = 6, lower.tail = FALSE)
AIC(m1, mnull)

mean(momsbirth$per_hatched)

sum(momsbirth$number.hatched)/sum(momsbirth$clutch.size)
mean(momsbirth$clutch.size*0.68)


(table(momsbirth$number.hatched[momsbirth$clutch.size %in% (3:8)],
       momsbirth$wholeyear[momsbirth$clutch.size %in% (3:8)]))
table(momsbirth$number.hatched)
barplot(table(momsbirth$number.hatched)[-1])
hist(momsbirth$not_hatched)
mean(momsbirth$not_hatched)
table(momsbirth$not_hatched>0)
barplot(table(momsbirth$clutch.size)[2:7])
table(momsbirth$clutch.size, momsbirth$number.hatched)

clutch_size_prop <- table(momsbirth$clutch.size)[2:7]
barplot(clutch_size_prop/sum(clutch_size_prop), ylim = c(0, 0.4))

momsbirth$clutch.size %>% table
shapiro.test(momsbirth$years)


tb_clutches <- table(momsbirth$dam.id[grepl('TNR', momsbirth$clutch.id)], 
      momsbirth$year[grepl('TNR', momsbirth$clutch.id)])

tb_clutches <- table(momsbirth$dam.id, 
                     momsbirth$year)
round(table(tb_clutches)[-c(1,5)]/sum(table(tb_clutches)[-c(1,5)]),4) 
 barplot(table(tb_clutches)[-1]/sum(table(tb_clutches)[-1]))
apply(tb_clutches, 2, function(x) table(x[x>0]))
px <- seq(0,100, 1)
plot(px/10,dpois(px, lambda = 15), type = 'l')

barplot(table(rpois(1000, lambda = 2))[-1])
