
library(tidyverse)
library(lubridate)
library(pscl)
# load --------------------------------------------------------------------

# zims data 
tid <- read.csv('./data/tidbinbilla_colony.csv')
mel <- read.csv('./data/melbourne_colony.csv')

# breeding 
mel2024 <- read.csv('./data/melbourn_zoo_tympo_breeding_2024.csv') %>% 
  filter(complete.cases(Number.of.clutches))

# colony data
ged <- read.csv('./data/GED_Pedigree_for_Emily_tidied.csv') %>% 
  rename_all(tolower)

# K ---------------------

Ktid <- table(tid$Status)[1]

# survival analysis -------------------------------------------------------
#https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/
tid_colony<- tid %>% 
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
  mutate(event = ifelse(status == 'dead', 1, 0),
         time = days.alive/365)

table(round(tid_colony$time))
km.curve <- with(tid_colony, Surv(time, event))
head(km.curve,80)


km_fit <- survfit(Surv(time, event) ~ 1, data=tid_colony)
survival_m <- summary(km_fit, times = c(1:7))


plot(km_fit, xlab="Years", main = 'Kaplan Meyer Plot')


# reproduction ----------------------------------------------------------------


# perc reproducing --------------------------------------------------------

tb_mel <- table(mel2024$Number.of.clutches)

perc_repro <- 1 - (tb_mel[1]/sum(tb_mel))

dups <- mel2024$Female[duplicated(mel2024$Female)]
mel2024[mel2024$Female %in% dups,]

## clutch distribution -----------------------------------------------------
clutches_mums <- ged %>% 
  filter(complete.cases(clutch.id)) %>% 
  group_by(clutch.id, dam.id, sire.id, lay.date, clutch.size) %>% 
  summarise(number.hatched = sum(hatch.success))%>% 
  rename_all(tolower)%>% 
  left_join(ged[,c('ged.id', 'hatch.date')], 
            by = c('dam.id' = 'ged.id')) %>% 
  filter(!grepl('UC', clutch.id)) %>% 
  ungroup() %>% 
  mutate(year = year(lay.date),
         age_days = ymd(lay.date) - ymd(hatch.date),
         age_years = as.numeric(age_days/365),
         wholeyear = round(age_years)) %>% 
  arrange(hatch.date) %>% 
  filter(complete.cases(hatch.date))

clutches_mums %>% head

### age diff ----------------------------------------------------------------
clutches_mums %>% 
  group_by(wholeyear) %>% 
  summarise(size = mean(clutch.size),
            sd = sd(clutch.size)) %>% 
  mutate(lcl = size - (sd*2/n()),
         ucl = size + (sd*2)/n())%>% 
  ggplot(aes(wholeyear, size))+
  # geom_jitter(data = clutches_mums, aes(y = clutch.size),
  #             colour = 'grey85', width = 0.1)+
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.1)+
  geom_point(size = 3)+
  theme_classic()

barplot(table(clutches_mums$clutch.size))
summary(m1 <- glm(number.hatched ~ factor(wholeyear),
                       data = clutches_mums), family = 'poisson')
anova(m1) # marginal but to keep it simple assume no difference
mnull <- update(m1, . ~ 1)

AIC(m1, mnull) # do not reject null model 

# assume no age difference

### clutch size --------------------

table(clutches_mums$clutch.size, clutches_mums$number.hatched)
# clutch size 2 and 9 never had hatched individuals - will remove them as viable clutch size
# clutch size of 8 only produced 7 so lets maintain clutch size max of 7 (for modelling)

c_size <- clutches_mums[clutches_mums$clutch.size %in% (3:7),]

barplot(table(c_size$clutch.size))

perc_size <- table(c_size$clutch.size)/sum(table(c_size$clutch.size))

clutch_beta <- dbeta(seq(0.1,0.9, length.out = length(perc_size)),2,3)

perc_beta <- clutch_beta/sum(clutch_beta)
sum(abs(perc_size-perc_beta))
names(perc_beta) <- names(perc_size)

round(perc_size,2)
round(perc_beta,2)
    


# number of clutches ------------------------------------------------------


tb_clutches <- table(clutches_mums$dam.id, 
                     clutches_mums$year)
table(tb_clutches)
# only one mother had 4 so will consider this an outlier and not model
n_clutches <- table(tb_clutches)[-c(1,5)]
perc_n <- n_clutches/sum(n_clutches)
round(perc_n,2) 

barplot(perc_n)


tb_clutches


# Juvenile survival -------------------------------------------------------

## age difference ----------------------------------------------------------
clutches_mums %>% 
  group_by(wholeyear) %>% 
  summarise(hatched = mean(number.hatched),
            sd = sd(number.hatched)) %>% 
  mutate(lcl = hatched - (sd*2/n()),
         ucl = hatched + (sd*2)/n())%>% 
  ggplot(aes(wholeyear, hatched))+
  geom_jitter(data = clutches_mums, aes(y = number.hatched),
              colour = 'grey85', width = 0.1)+
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.1)+
  geom_point(size = 3)+
  ylim(0,6)+
  theme_classic()

barplot(table(clutches_mums$number.hatched))
summary(m1 <- zeroinfl(number.hatched ~ factor(wholeyear),
                       data = clutches_mums))

mnull <- update(m1, . ~ 1)

pchisq(2 * (logLik(m1) - logLik(mnull)), df = 6, lower.tail = FALSE)
AIC(m1, mnull) # do not reject null model 

# No age difference


## survival ----------------------------------------------------------------



perc_hatch <- sum(momsbirth$number.hatched)/sum(momsbirth$clutch.size)

surv0_1 <- survival_m$surv[1]

juv_survival <- perc_hatch*surv0_1

juv_survival

adult_survival  <- survival_m$surv[2]
adult_survival


# Parameters --------------------------------------------------------------
Ktid # make 100
juv_survival
adult_survival

perc_n # probability of number of clutches

perc_beta # prob of clutch size

perc_repro # percent of females reproducing

