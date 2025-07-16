

# load --------------------------------------------------------------------
figprefix <- ''
source('./code/00_libraries.R')
source('./code/05_model.R')
stages <- names(readRDS('./output/stage_distribution.rds'))
transition_mat <- readRDS('./output/transition_matrix.rds')
base_params <- readRDS('./output/base_parameters.rds')
fper_repro <- base_params$fecundity
max_age <- base_params$age
K <- base_params$K
clutch_sizes <- base_params$clutch
transition <- base_params$transition

param_selected <- readRDS('./output/selected_25models_parameters.RDS')
param_selecteddf <-readRDS('./output/selected_25models_parameters_df.RDS')
param_envdf <- readRDS('./output/selected_25models_parameters_df_environment.RDS')
# stage distribution ------------------------------------------------------

prob_clutch <- dbeta(seq(0.1,0.9, length.out = length(clutch_sizes)), 2,3)
prob_clutch/sum(prob_clutch)

mean_clutch <- round(mean(sample(clutch_sizes, 100000, T, prob_clutch)),1)


fecundity <- fper_repro*0.5*mean_clutch
ASAsur <- mean(param_selecteddf$mean[grep('survival_A', param_selecteddf$name)])
Jt <- mean(param_selecteddf$mean[grep('survival_J', param_selecteddf$name)])

stage.mat <- matrix(0, 
                    nrow = length(stages),
                    ncol = length(stages),
                    byrow = TRUE,
                    dimnames = list(stages,
                                    stages))

# fecundity
fecund_ages <- which(stages %in% paste0('A', 1:(max_age)))
stage.mat[1,fecund_ages] <- fecundity
# juvenile survival + growth
stage.mat[2:3,which(stages == 'J')] <- Jt*c(transition, 1-transition)

# adult survival + growth
for (i in 2:(max_age+2)) {
  stg <- stages[i]
  if (stg == 'SA') stage.mat[i+2,i] <- ASAsur
  if (stg != 'SA' & i < (max_age+2)) stage.mat[i+1, i] <- ASAsur
}

stage.mat
eigen_analysis(stage.mat)

stage_distribution <- eigen_analysis(stage.mat)[[2]]
names(stage_distribution)<- stages
saveRDS(stage_distribution, './output/stage_distribution_base.rds')


as.data.frame(stage.mat) %>% 
  mutate_if(is.numeric, round, 3) %>% 
  mutate(stages = rownames(.)) %>% 
  rbind(c(round(stage_distribution,2),'stable:')) %>% 
  rbind(c(round(eigen_analysis(stage.mat)[[1]],3), rep(NA,max_age+1), 'lambda =')) %>% 
  relocate(stages) %>% 
  flextable() %>% 
  autofit() %>% 
 # border_remove() %>%
  vline(j = 1, i = 1:nrow(stage.mat), part = 'body', border = fp_border_default(width = 1.5)) %>% 
  vline(j = 1, part = 'header', border = fp_border_default(width = 1.5)) %>% 
  hline(i = nrow(stage.mat), border = fp_border_default(width = 1.5)) %>% 
  hline(i = nrow(stage.mat)+1, j=1:2, border = fp_border_default(width = 1.5)) %>% 
  vline(i = nrow(stage.mat)+2, j = c(2),border = fp_border_default(width = 1.5), part = 'body') %>% 
  align(j = 1,  align = 'right') %>% 
  align(part = 'header', j  = 1,align = 'right') %>% 
  hline(i = nrow(stage.mat)+1, border = fp_border_default(width = 1.5)) %>% 
  hline_bottom(j = 3:(nrow(stage.mat)+1),border = fp_border_default(width = 0)) %>%  
  vline_left(i = nrow(stage.mat)+2, part = 'body', border = fp_border_default(width = 1.5)) %>% 
  bold(part = 'header') %>% 
  bold(j = 2, i = nrow(stage.mat)+2, part = 'body') %>% 
  bold(j= 1) %>% 
  font(part = 'all', fontname='calibri') %>% 
  saveRDS('./output/new_stage_distribution.rds')


# model inputs ------------------------------------------------------------
top25 <- param_selecteddf
paramlist <- list(populations = c('CA', 'JE', 'JW', 'MA'),
initial_ab = round(top25$mean[grep('N_', top25$name)]), # adult abundance for populations
survival = top25$mean[grep('_ASA_', top25$name)], # survival of adults and SA at sites
survival_J = top25$mean[grep('_J_', top25$name)], # juvenile survival
env_stoch = param_envdf$mean, # sd on survival
f_reproducing = fper_repro, # proportion of females reproducing
clutch_sizes = clutch_sizes, # clutch size range   
K = K # carrying capcity applied to adults and SA
)




paramlist[which(names(paramlist) != 'clutch_sizes')] %>% 
  do.call('cbind', .) %>% 
  as.data.frame %>% #dplyr::select(-clutch_sizes) %>% 
  pivot_longer(initial_ab:K) %>% 
  mutate(value = as.numeric(value),
         value = round(value, 2),
         name = ifelse(name =='survival', 'survival_A_SA', name),
         name = ifelse(name == 'initial_ab', 'N_initial', name),
         name = gsub('_', ' ', name)) %>% 
  pivot_wider(names_from = populations, values_from = value) %>% 
  mutate(name = ifelse(name == 'env stoch', 'env (logit)', name)) %>% 
  rename(Parameter = name) %>% 
  mutate_if(is.double, as.character) %>% 
  flextable() %>% 
  autofit() %>%
  theme_zebra() %>%  
  hline_bottom(border = fp_border_default(width = 2)) %>% 
  hline_top(border = fp_border_default(width = 1)) %>%  
  bold(part = 'header') %>% 
  bold(j = 1) %>% 
  font(part = 'all', fontname='calibri') %>% 
  saveRDS('./output/base_parameter_inputs.rds')


# base model --------------------------------------------------------------
reps_base <- 1000
#top25 <- apply(param_selected, 2, summary)[4,]

N_sim <- em.pva_simulator(populations = c('CA', 'JE', 'JW', 'MA'),
                          stages = stages,
                          stage_distribution = stage_distribution, 
                          initial_ab = round(paramlist$initial_ab), # adult abundance for populations
                          survival = paramlist$survival, # survival of adults and SA at sites
                          survival_J = paramlist$survival_J, # juvenile survival
                          survival_logit_sd = NULL,
                          site_adjust = NULL,
                          env_stoch = paramlist$env_stoch, # sd on survival
                          transition_mat = transition_mat, # transition prob to SA
                          f_reproducing = paramlist$f_reproducing, # proportion of females reproducing
                          clutch_sizes = clutch_sizes, # clutch size range   
                          K = K, # carrying capcity applied to adults and SA
                          time_steps = 50, # time
                          replicates = reps_base)


# extract N ---------------------------------------------------------------


N_sim %>% length


system.time(N_simulated <- mclapply(1:reps_base, em.extract_N, N_sim,
                                    mc.cores = 10) %>% 
  do.call('rbind', .) %>% 
  mutate_if(is.character, factor))

nlevels(N_simulated$stage)
nlevels(N_simulated$pop)
levels(N_simulated$stage) <- stages
levels(N_simulated$pop) <- c('CA', 'JE', 'JW', 'MA')


ggplot(N_simulated, aes(tstep, N+1, group = rep, colour= rep))+
  geom_line(alpha = 0.1)+
  facet_grid(stage~pop, scale = 'free')+
  theme_bw()+
  #scale_y_log10()+
  theme(strip.text.y = element_text(angle = 0),
        panel.grid = element_blank(),
        legend.position = 'none')

ggsave(paste0('./figures/',figprefix, 'base_model_Pops_Stages.png'),dpi = 300,
       height = 7, width = 7, units = 'in')  



# extinction --------------------------------------------------------------

sim_N_sum <- N_simulated %>% 
  group_by(rep, tstep, pop) %>% 
  summarise(N = sum(N)) %>% 
  mutate(extinct = N < 2)
sum(sim_N_sum$extinct)
extinction_prob <- sim_N_sum %>% ungroup() %>% 
  group_by(tstep,pop) %>% 
  summarise(sims = n(),
            sim0 = sum(extinct),
            ext_p = sum(extinct)/n())
extinction_prob
sim_N_sum %>% 
  mutate(year = tstep +2012) %>% 
  ggplot(aes(year, N, group = rep))+
 # geom_hline(yintercept = 100, colour = 'red', lty = 2)+
  geom_line(alpha = 0.1, colour = 'grey50')+
  facet_wrap(~pop, scale = 'free')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold'))
 #scale_y_log10()
tail(extinction_prob)

ggsave(paste0('./figures/',figprefix, 'base_model_Pops.png'),dpi = 300,
       height = 4, width = 7, units = 'in')  

extinction_prob %>% 
  mutate(year = tstep +2012,
         Sites = pop) %>%
  ggplot(aes(year, ext_p, colour = Sites))+
  geom_vline(xintercept = c(2023, 2025), lty = c(3,2), colour = 'grey')+
  geom_line()+
  geom_point()+
  theme_classic()+
  ylab('Extinction probability')+
  xlab('Year')+
  scale_x_continuous(breaks = c(2013, 2025,seq(2010,2060,10)))+
  scale_y_continuous(breaks = seq(0,1,0.2))

ggsave(paste0('./figures/',figprefix, 'base_model_extinction_noerror.png'),dpi = 300,
       height = 3.8, width = 7, units = 'in') 


extinction_prob_sample <- sim_N_sum %>% ungroup() %>% 
  mutate(subsample = rep %% 10) %>% 
  group_by(tstep,pop, subsample) %>% 
  summarise(sims = n(),
            sim0 = sum(extinct),
            ext_p = sum(extinct)/n())
table(extinction_prob_sample$subsample)
extinction_prob_sample %>%
  ungroup() %>% 
  group_by(tstep, pop) %>%
  summarise(p = mean(ext_p),
            sd = sd(ext_p),
            lcl = quantile(ext_p, 0.025),
            ucl = quantile(ext_p, 0.975),
            #lcl = p - sd*15,
            #ucl = p + sd*15
            ) %>% 
  mutate(year = tstep +2012,
         Sites = pop) %>%
  ggplot(aes(year, p, colour = Sites))+
  geom_hline(yintercept = c(0,1), linewidth = 0.5, colour = 'grey50')+
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0, alpha = 0.5)+
  geom_vline(xintercept = c(2023, 2025), lty = c(3,2), colour = 'grey')+
  geom_line()+
  geom_point()+
  theme_classic()+
  ylab('Extinction probability')+
  xlab('Year')+
  scale_x_continuous(breaks = c(2013, 2025,seq(2010,2060,10)))+
  scale_y_continuous(breaks = seq(0,1,0.2))
ggsave(paste0('./figures/',figprefix, 'base_model_extinction.png'),dpi = 300,
       height = 3.8, width = 7, units = 'in') 



# summary stats -----------------------------------------------------------

## growth rate --------------
base_summary <- N_simulated %>% 
  left_join(extinction_prob) %>% 
  group_by(rep, tstep, pop) %>% 
  summarise(N = sum(N),
            ext_p = mean(ext_p, na.rm = TRUE)) %>% 
  mutate(ext100 = ext_p == 1) %>% 
  ungroup() %>% 
  group_by(pop, rep) %>% 
  mutate(diff_year = tstep - lag(tstep),
         diff_growth = N - lag(N),
         rate_percent = (diff_growth /diff_year)/ lag(N) * 100,
         time2ext = min(tstep[ext100 | pop == 'CA']),
         time2ext = ifelse(time2ext ==1, NA, time2ext)) %>%
  ungroup() %>% 
  group_by(pop) %>% 
  summarise(base_growth_rate = mean(rate_percent, na.rm = T),
            extinction_base = mean(ext_p),
            time_to_extinction = mean(time2ext)) %>% 
  mutate(year_of_extinction = time_to_extinction + 2012)
base_summary
saveRDS(base_summary, './output/base_summary_r_ext_p.rds')
