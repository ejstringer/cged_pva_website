

# environmental stochasticity ---------------------------------------------

abundance <- readRDS('./output/top_25_best_fitting_models.rds') %>% 
  filter(year >= 2013)
abundance$tstep <- abundance$year-2012
abundance$log_Nsim <- log(abundance$Nsim+1)
abundance$log_Nall <- log(abundance$N_all+1)
head(abundance)

ggplot(abundance, aes(tstep, Nsim, colour = site, group = run))+
  geom_line()+
  facet_grid(~site)+
  scale_y_log10()+
  geom_smooth(method = 'lm', aes(group = site), colour = 'red')+
  theme_bw()

sim_models <- lapply(unique(abundance$site),
                     function(x) (lm(log_Nsim ~ tstep, 
                                       data = abundance[abundance$site == x,])))

names(sim_models) <- unique(abundance$site)
lapply(sim_models, summary)
lapply(sim_models, coef)
sim_models
sapply(lapply(sim_models, resid), hist)

sim_resid <- do.call(c, lapply(sim_models, resid))
abundance$sim_resid <- sim_resid
abundance_real <- abundance %>% 
  filter(run =='AFAG5151',
         complete.cases(N_all)) %>% 
  select(log_Nall, N_all, tstep, site)

ggplot(abundance_real, aes(tstep, N_all+1, colour = site))+
  geom_point()+
  geom_line()+
  facet_grid(~site)+
  geom_hline(yintercept = 1)+
  scale_y_log10()+
  geom_smooth(method = 'lm', aes(group = site), colour = 'red')+
  theme_bw()


real_models <- lapply(unique(abundance_real$site),
                     function(x) (lm(log_Nall ~ tstep, 
                                     data = abundance_real[abundance_real$site == x,])))

names(real_models) <- unique(abundance$site)
lapply(real_models, summary)
lapply(real_models, coef)

real_resid <- do.call(c, lapply(real_models, resid))

abundance_real$real_resid <- real_resid

abundance_resid <- abundance %>% 
  left_join(abundance_real[,colnames(abundance_real)!='run'])

ggplot(abundance_resid, aes(exp(sim_resid)))+
  geom_histogram(alpha = 1, fill = 'pink')+
  geom_histogram(alpha = 0.5, fill = 'forestgreen',
                 aes(x = exp(real_resid)))+
  facet_wrap(~site)+
  scale_x_log10()+
  theme_bw()

abundance_real2 <- abundance_real %>% mutate(run = 'real')
ggplot(abundance, aes(tstep, Nsim+1, colour = site, group = run))+
  geom_line(linewidth = 0.2)+
  geom_line(aes(y = N_all+1), data = abundance_real2, colour = 'black')+
  facet_grid(~site)+
  scale_y_log10()+
  geom_hline(yintercept = 1)+
  theme_bw()



abundance %>% 
  mutate(xdiff = N_all-Nsim) %>% 
  ggplot(aes(abs(xdiff), fill = site))+
  geom_histogram(colour = 'grey', linewidth = 0.25)+
  facet_wrap(~site)+
  theme_bw()


abundance %>% 
  group_by(site, year) %>% 
  summarise(N_all = mean(N_all),
            Nsim = mean(Nsim)) %>% 
  mutate(diff = N_all - Nsim) %>%  
  pivot_longer(N_all:Nsim) %>% 
  filter(complete.cases(value)) %>% 
  ggplot(aes(year, value, colour = name))+
  geom_line()+
  facet_wrap(~site, scale = 'free')+
  theme_bw()+
  scale_y_log10()

abundance %>% 
  group_by(site, year) %>% 
  summarise(N_all = mean(N_all),
            Nsim = mean(Nsim)) %>% 
  mutate(diff = log(N_all+1) - log(Nsim+1)) %>% 
  filter(complete.cases(diff)) %>% 
  ungroup() %>% 
  ggplot(aes(diff)) + 
  geom_histogram()+
  facet_wrap(~site)+
  theme_bw()

ggplot(abundance, aes(log(N_all)))+
  geom_histogram()


abundance$N_all %>% log %>% sd

# time to extinction ------------------------------------------------------

abundance



