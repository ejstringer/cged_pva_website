
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

femaleHatches <- mel2024 %>% 
  group_by(Female) %>% 
  summarise(n_clutches = sum(total.offspring.successfully.hatched))

tb_mel <- table(femaleHatches$n_clutches)

females2024 <- mel2024 %>% 
  group_by(Female) %>% 
  summarise(n_clutches = sum(Number.of.clutches))

tb_mel <- table(females2024$n_clutches)


perc_repro <- 1 - (tb_mel[1]/sum(tb_mel))


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

table(clutches_mums$number.hatched==0)
15/(15+68)
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

perc_repro # percent of females reproducing

perc_n # probability of number of clutches

perc_beta # prob of clutch size


#### _________________ ####
# model -------------------------------------------------------------------
for(sub in 1:2){
reps <- 1000
seq_sup <- c(5,10,15)

supplement <- sub == 1

N_sim <- em.pva_simulator(populations = c('CA', 'JE', 'JW', 'MA','TR'),
                          stages = c('J', 'SA','A1','A2','A3', 'A4'),
                          stage_distribution = stage_distribution, 
                          initial_ab = c(52,34,233,146, 93), 
                          survival = c(0.45, 0.39, 0.27, 0.39, 0.88), 
                          survival_J = c(0.58, 0.52, 0.39, 0.52, 0.67),
                          survival_logit_sd = NULL,
                          site_adjust = NULL,
                          env_stoch = c(0.3, 1.9, 1.7, 2.1, 0.001), 
                          transition_mat = transition_mat, 
                          f_reproducing = c(rep(0.35,4), 0.6), 
                          clutch_sizes = 3:7,   
                          K = c(100,440, 740, 1000, 100),
                         # density_stage = c(2:5),
                          time_steps = 20, # time
                          replicates = reps,
                          supp = supplement,
                          n.supp = c(rep(25, 4),0), 
                          when.supp = seq_sup)


N_sim %>% length


# join --------------------------------------------------------------------

system.time(N_simulated <- mclapply(1:reps, em.extract_N, N_sim,
                                    mc.cores = 5) %>% 
              do.call('rbind', .) %>% 
              mutate_if(is.character, factor))

nlevels(N_simulated$stage)
nlevels(N_simulated$pop)
levels(N_simulated$stage) <- c('J', 'SA', 'A1', 'A2', 'A3', 'A4')
levels(N_simulated$pop) <- c('CA', 'JE', 'JW', 'MA', 'TR')


sim_N_sum <- N_simulated %>% 
  filter(stage != 'J') %>% 
  group_by(rep, tstep, pop) %>% 
  summarise(N = sum(N)) %>% 
  mutate(extinct = N < 2)
sum(sim_N_sum$extinct)

if(sub == 1) sim_N_sum_Sup <- sim_N_sum
if(sub == 2) sim_N_sum_base <- sim_N_sum


}
sim_N_sum_base$supplementation <- 'no'
sim_N_sum_Sup$supplementation <- 'yes'
# summarised plot ---------------------------------------------------------

sim_N_sum_base %>%
  bind_rows(sim_N_sum_Sup) %>% 
  group_by(tstep, pop, supplementation) %>% 
  summarise(N = mean(N)) %>% 
  mutate(year = tstep +2012) %>% 
  ggplot(aes(year, N, colour = supplementation))+
  geom_hline(yintercept = 0, colour = 'red', lty = 2, linewidth = 0.3)+
  geom_line(alpha = 1)+
  facet_wrap(~pop, scale = 'free')+
  theme_bw()+
  geom_vline(xintercept = seq_sup+2012, colour = 'grey60',
             lty = 2)+
  theme(panel.grid = element_blank(),
        legend.position = 'bottom',
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold'))


extinction_prob <- sim_N_sum_base %>%
  bind_rows(sim_N_sum_Sup) %>% ungroup() %>% 
  group_by(tstep,pop,supplementation) %>% 
  summarise(sims = n(),
            sim0 = sum(extinct),
            ext_p = sum(extinct)/n())
tail(extinction_prob)

extinction_prob %>% 
  mutate(year = tstep +2012,
         Sites = pop) %>%
  ggplot(aes(year, ext_p, colour = Sites))+
  geom_vline(xintercept = seq_sup+2012, lty = 3, colour = 'grey')+
  geom_line(aes(linetype = supplementation))+
  geom_point(aes(shape = supplementation))+
  theme_classic()+
  ylab('Extinction probability')+
  xlab('Year')+
  scale_x_continuous(breaks = c(2013, 2025,seq(2010,2060,10)))+
  scale_y_continuous(breaks = seq(0,1,0.2))
 
# delcine threshold ------------------------------------------------------
#https://kevintshoemaker.github.io/NRES-470/LECTURE12.html#Step_1:_conceptualize_the_life_history
final_N <- sim_N_sum_base %>%
  bind_rows(sim_N_sum_Sup) %>% 
  ungroup() %>% 
  filter(tstep == max(tstep),
         pop == "CA")


# plot probabilities of different severities of decline
Init_N <- 52
declines <- seq(0,100,by=2)
declineprob <- numeric(length(declines))

for(s in 1:length(declines)){
  declineprob[s] <- length(which(final_N$N <(Init_N-(declines[s]/100)*Init_N)))/length(final_N$N)
}

plot(declines[-101],declineprob[-101],type="l",lwd=2,xlab="Decline threshold (percent)",ylab="Probability of falling below threshold")

abline(v=25,col="red",lwd=2)
Init_N*0.25

# check juv ---------------------------------------------------------------


N_simulated %>% 
  mutate(age = ifelse(stage == 'J', 'juv', 'adult')) %>% 
  group_by(rep, tstep, pop, age) %>% 
  summarise(N = sum(N)) %>% 
  pivot_wider(names_from = age, values_from = N) %>% 
  filter(adult == 1) %>% summary

