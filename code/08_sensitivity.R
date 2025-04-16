original_paramlist <- paramlist
sim_summary
paramlist 

names(paramlist)
sim_summary <- N_simulated %>% 
  group_by(rep, tstep, pop) %>% 
  summarise(N = sum(N)) %>% 
  ungroup() %>% 
  group_by(pop, rep) %>% 
  mutate(diff_year = tstep - lag(tstep),
         diff_growth = N - lag(N),
         rate_percent = (diff_growth /diff_year)/ lag(N) * 100) %>%
  ungroup() %>% 
  group_by(pop) %>% 
  summarise(growth_rate_base = mean(rate_percent, na.rm = T)) 
sim_summary <- NULL
seqencedf <- data.frame(parameter = rep(names(paramlist)[2:6], each = 2), ten = c(0.9, 1.1))

# base model --------------------------------------------------------------
reps_base <- 1000
for (i in 1:nrow(seqencedf)) {
  
#top25 <- apply(param_selected, 2, summary)[4,]
param_adjust <- seqencedf$parameter[i]
adjust <- seqencedf$ten[i]
paramlist10 <- paramlist
paramlist10[[param_adjust]] <- paramlist10[[param_adjust]]*adjust  


N_sim <- em.pva_simulator(populations = c('CA', 'JE', 'JW', 'MA'),
                          stages = c('J', 'SA','A1','A2','A3', 'A4'),
                          stage_distribution = stage_distribution, 
                          initial_ab = round(paramlist10$initial_ab), # adult abundance for populations
                          survival = paramlist10$survival, # survival of adults and SA at sites
                          survival_J = paramlist10$survival_J, # juvenile survival
                          survival_logit_sd = NULL,
                          site_adjust = NULL,
                          env_stoch = paramlist10$env_stoch, # sd on survival
                          transition_mat = transition_mat, # transition prob to SA
                          f_reproducing = paramlist10$f_reproducing, # proportion of females reproducing
                          clutch_sizes = clutch_sizes, # clutch size range   
                          K = K, # carrying capcity applied to adults and SA
                          time_steps = 50, # time
                          replicates = reps_base)


# extract N ---------------------------------------------------------------


N_sim %>% length


system.time(N_simulated <- mclapply(1:reps_base, 
                                    em.extract_N, N_sim, mc.cores = 10) %>% 
              do.call('rbind', .) %>% 
              mutate_if(is.character, factor))
nlevels(N_simulated$stage)
nlevels(N_simulated$pop)
levels(N_simulated$stage) <- c('J', 'SA', 'A1', 'A2', 'A3', 'A4')
levels(N_simulated$pop) <- c('CA', 'JE', 'JW', 'MA')

sim_adjusted_lambda <- N_simulated %>% 
  group_by(rep, tstep, pop) %>% 
  summarise(N = sum(N)) %>% 
  ungroup() %>% 
  group_by(pop, rep) %>% 
  mutate(diff_year = tstep - lag(tstep),
         diff_growth = N - lag(N),
         rate_percent = (diff_growth /diff_year)/ lag(N) * 100) %>%
  ungroup() %>% 
  group_by(pop) %>% 
  summarise(growth_rate = mean(rate_percent, na.rm = T)) %>% 
  left_join(base_summary)

sim_N_sum <- N_simulated %>% 
  group_by(rep, tstep, pop) %>% 
  summarise(N = sum(N)) %>% 
  mutate(extinct = N < 2)
sum(sim_N_sum$extinct)
sim_adjusted_lambda <- sim_N_sum %>% ungroup() %>% 
  group_by(tstep,pop) %>% 
  summarise(sims = n(),
            sim0 = sum(extinct),
            ext_p = sum(extinct)/n())%>% 
  group_by(pop) %>% 
  summarise(extinction = mean(ext_p)) %>% 
  left_join(base_summary)

sim_adjusted_lambda$parameter <- param_adjust
sim_adjusted_lambda$adjusted <- adjust

sim_summary <- bind_rows(sim_summary, sim_adjusted_lambda)
}
sim_summary %>% 
  pivot_wider(names_from = adjusted, values_from = growth_rate) %>% 
  mutate(S = (`1.1`-`0.9`)/base_growth_rate,
         param2 = rep(c('N', 'Sur', 'Juv', 'Env', 'fem'), each = 4)) %>% 
  ggplot(aes(param2, S, fill = pop)) +
  geom_bar(stat = 'identity')+
  facet_wrap(~pop, scale = 'free')

sim_summary$param2 <- factor(sim_summary$parameter) 

levels(sim_summary$param2)<- c( 'Env','F','N',  'Sur', 'Juv')
sim_summary %>%  
  #mutate(S = growth_rate-base_growth_rate) %>%
  pivot_wider(names_from = adjusted, 
              values_from = growth_rate) %>% 
  mutate(S = (`1.1`-`0.9`)/(0.1*base_growth_rate)) %>% 
  ggplot(aes(param2, S, fill = param2)) +
  geom_bar(stat = 'identity')+
  facet_wrap(~pop, scale = 'free_x')+
  theme_bw()

sim_summary %>% 
  pivot_wider(names_from = adjusted, 
              values_from = growth_rate) %>% 
  group_by(parameter, pop) %>% 
  summarise(lower = mean(`0.9`),
            baseline = mean(base_growth_rate),
            upper = mean(`1.1`)) %>% 
  pivot_longer(lower:upper) %>% 
  ggplot(aes(value,parameter))+
  geom_line(linewidth = 1)+
  geom_point(aes(size = name))+
  theme_bw()+
  facet_grid(~pop, scale = 'free')+
  scale_size_manual(values = c(3,0,0))

sim_summary %>% 
  pivot_wider(names_from = adjusted, 
              values_from = extinction) %>% 
  group_by(parameter, pop) %>% 
  summarise(lower = mean(`0.9`),
            baseline = mean(extinction_base),
            upper = mean(`1.1`)) %>% 
  # mutate(S = (upper-lower)/baseline) %>% 
  # ggplot(aes(parameter, S, fill = pop))+
  # geom_bar(stat = 'identity', position = position_dodge())+
  # facet_wrap(~pop, scale = 'free')
  pivot_longer(lower:upper) %>% 
  ggplot(aes(value,parameter))+
  geom_line(linewidth = 0.7)+
  geom_point(aes(size = name))+
  theme_bw()+
  facet_grid(~pop, scale = 'free')+
  scale_size_manual(values = c(3,0,0))



sim_summary %>% 
  pivot_wider(names_from = adjusted, 
              values_from = extinction) %>% 
  group_by(param2, pop) %>% 
  summarise(lower = mean(`0.9`),
            baseline = mean(extinction_base),
            upper = mean(`1.1`)) %>% 
  mutate(S = (upper-lower)/baseline) %>%
  group_by(param2) %>% 
  summarise(S = mean(S)) %>% 
  arrange(S) %>% 
  mutate(n = row_number()) %>% 
  ggplot(aes(n, S, fill = param2))+
  geom_bar(stat = 'identity')+
  theme_bw()
 # filter(param2 != 'N', param2 != 'Env') %>% 
  ggplot(aes(param2, S, fill = pop))+
  geom_bar(stat = 'identity', position = position_dodge())+
  facet_wrap(~pop, scale = 'free')
