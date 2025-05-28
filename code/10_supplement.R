figprefix <- ''
source('./code/00_libraries.R')
source('./code/05_model.R')
# load --------------------------------------------------------------------

base_params <- readRDS('./output/base_parameters.rds')
param_selecteddf <-readRDS('./output/selected_25models_parameters_df.RDS')
param_envdf <- readRDS('./output/selected_25models_parameters_df_environment.RDS')
area_ntg <- read.csv('./output/area_used.csv')

# parameters --------------------------------------------------------------

stages <- names(readRDS('./output/stage_distribution.rds'))
transition_mat <- readRDS('./output/transition_matrix.rds')
stage_distribution <- readRDS('./output/stage_distribution_base.rds')
grasslands <- c('CA', 'JE', 'JW', 'MA')
initial_N <-round(param_selecteddf$mean[grepl('N_init', param_selecteddf$name)])
phiASA  <-  param_selecteddf$mean[grepl('survival_ASA', param_selecteddf$name)]
phiJ  <-  param_selecteddf$mean[grepl('survival_J', param_selecteddf$name)]
fecundity <- rep(param_selecteddf$mean[grepl('F_', param_selecteddf$name)],4)
K <- base_params$K
clutch_sizes <- base_params$clutch
environmental <- param_envdf$mean


area_ntg$site <- grasslands 
# supplementation ---------------------------------------------------------
supplement <- T
number_suplementation <- rep(25,4)
seq_suplementation <- c(5,10,15)

# model -------------------------------------------------------------------
reps <- 1000
N_sim <- em.pva_simulator(populations = grasslands,
                          stages = stages,
                          stage_distribution = stage_distribution, 
                          initial_ab = initial_N, 
                          survival = phiASA, 
                          survival_J = phiJ,
                          survival_logit_sd = NULL,
                          site_adjust = NULL,
                          env_stoch = environmental, 
                          transition_mat = transition_mat, 
                          f_reproducing = fecundity, 
                          clutch_sizes = clutch_sizes,   
                          K = K,
                          time_steps = 50, # time
                          replicates = reps,
                          supp = supplement,
                          n.supp = number_suplementation, 
                          when.supp = seq_suplementation)

N_sim %>% length
system.time(N_simulated <- mclapply(1:reps, em.extract_N, N_sim,
                                    mc.cores = 5) %>% 
              do.call('rbind', .) %>% 
              mutate_if(is.character, factor))
nlevels(N_simulated$stage)
nlevels(N_simulated$pop)
levels(N_simulated$stage) <- stages
levels(N_simulated$pop) <- grasslands

sim_N_sum <- N_simulated %>% 
  filter(stage != 'J') %>% 
  group_by(rep, tstep, pop) %>% 
  summarise(N = sum(N)) %>% 
  mutate(extinct = N < 2,
         supplementation = supplement)
sum(sim_N_sum$extinct)

# abundance ---------------------------------------------------------------


sim_N_sum %>% 
  group_by(tstep, pop, supplementation) %>% 
  summarise(N = mean(N)) %>% 
  mutate(year = tstep +2012) %>% 
  ggplot(aes(year, N, colour = supplementation))+
  geom_hline(yintercept = 0, colour = 'red', lty = 2, linewidth = 0.3)+
  geom_line(alpha = 1)+
  facet_wrap(~pop, scale = 'free')+
  theme_bw()+
  geom_vline(xintercept = seq_suplementation+2012, colour = 'grey60',
             lty = 2)+
  theme(panel.grid = element_blank(),
        legend.position = 'bottom',
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold'))


# extinction --------------------------------------------------------------

extinction_prob <- sim_N_sum %>% ungroup() %>% 
  group_by(tstep,pop,supplementation) %>% 
  summarise(sims = n(),
            sim0 = sum(extinct),
            ext_p = sum(extinct)/n())
tail(extinction_prob)

extinction_prob %>% 
  mutate(year = tstep +2012,
         Sites = pop) %>%
  ggplot(aes(year, ext_p, colour = Sites))+
  geom_vline(xintercept = seq_suplementation+2012, lty = 3, colour = 'grey')+
  geom_line(aes(linetype = supplementation))+
  geom_point(aes(shape = supplementation))+
  theme_classic()+
  ylab('Extinction probability')+
  xlab('Year')+
  scale_x_continuous(breaks = c(2013, 2025,seq(2010,2060,10)))+
  scale_y_continuous(breaks = seq(0,1,0.2))

# delcine threshold ------------------------------------------------------
#https://kevintshoemaker.github.io/NRES-470/LECTURE12.html#Step_1:_conceptualize_the_life_history
final_N <- sim_N_sum %>%
  left_join(area_ntg, by = c('pop' = 'site')) %>% 
  ungroup() %>% 
  mutate(N_ha = N/area_ha) %>% 
  filter(tstep == 20) %>% 
  split(., .$pop)

# plot probabilities of different severities of decline

area_ntg$site
Init_N <- initial_N/area_ntg$area_ha
declines <- seq(0,1,by=0.02)
declineprob <- numeric(length(declines))

declineprob_list <- list()

for (pop in 1:length(grasslands)) {
  declineprob <- numeric(length(declines))

for(s in 1:length(declines)){
  N_pop <-final_N[[pop]]$N_ha
  N_under <- sum(N_pop < (Init_N-(declines[s])*Init_N[pop]))
  declineprob[s] <- N_under/length(N_pop)
}

  declineprob_list[[pop]] <- data.frame(pop = grasslands[pop],
                                        declines = declines,
                                        prob = declineprob)
  
}


declineprob_list %>%
  do.call('rbind', .) %>% 
  ggplot(aes(declines*100, prob, colour = pop))+
  geom_line()+
  theme_classic()+
  xlab("Decline threshold (percent)")+
  ylab("Probability of falling below threshold")+
  geom_vline(xintercept = 25, colour = 'black', alpha = 0.5)
