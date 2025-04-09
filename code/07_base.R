
em.extract_N <- function(repx, sim){
  tstep <- dim(sim)[3]
  nstage <- dim(sim)[1]
  npop <- dim(sim)[2]
  
  df_N_overtime <- NULL
  for (ts in 1:tstep) {
    mat <- sim[,,ts,repx]
    rownames(mat) <- paste0('stage', 1:nstage)
    colnames(mat) <- paste0('pop', 1:npop)
    
    dfmat <- as.data.frame(mat) %>% 
      mutate(stage = rownames(mat)) %>% 
      pivot_longer(cols = colnames(mat),names_to = 'pop', values_to = 'N')
    dfmat$tstep <- ts
    dfmat$rep <- repx
    df_N_overtime <- rbind(df_N_overtime, dfmat)
    
  }
  return(df_N_overtime)
}

# base model --------------------------------------------------------------
reps_base <- 1000
top25 <- apply(param_selected, 2, summary)[4,]

N_sim <- em.pva_simulator(populations = c('CA', 'JE', 'JW', 'MA'),
                          stages = c('J', 'SA','A1','A2','A3', 'A4'),
                          stage_distribution = stage_distribution, 
                          initial_ab = round(top25[1:4]), # adult abundance for populations
                          survival = top25[6:9], # survival of adults and SA at sites
                          survival_J = top25[5], # juvenile survival
                          env_stoch = top25[11:14], # sd on survival
                          transition_mat = transition_mat, # transition prob to SA
                          f_reproducing = top25[10], # proportion of females reproducing
                          clutch_sizes = clutch_sizes, # clutch size range   
                          K = K, # carrying capcity applied to adults and SA
                          time_steps = 50, # time
                          replicates = reps_base)


# extract N ---------------------------------------------------------------


N_sim %>% length


N_simulated <- lapply(1:reps_base, em.extract_N, N_sim) %>% 
  do.call('rbind', .) %>% 
  mutate_if(is.character, factor)
nlevels(N_simulated$stage)
nlevels(N_simulated$pop)
levels(N_simulated$stage) <- c('J', 'SA', 'A1', 'A2', 'A3', 'A4')
levels(N_simulated$pop) <- c('CA', 'JE', 'JW', 'MA')


ggplot(N_simulated, aes(tstep, N+1, group = rep, colour= rep))+
  geom_line(alpha = 0.02)+
  facet_grid(stage~pop, scale = 'free')+
  theme_bw()+
  #scale_y_log10()+
  theme(strip.text.y = element_text(angle = 0),
        panel.grid = element_blank(),
        legend.position = 'none')

ggsave(paste0('./figures/',figprefix, '_base_model_Pops_Stages.png'),dpi = 300,
       height = 4, width = 7, units = 'in')  

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
  geom_hline(yintercept = 100, colour = 'red', lty = 2)+
  geom_line(alpha = 0.1, colour = 'grey50')+
  facet_wrap(~pop, scale = 'free')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold'))
 #scale_y_log10()
tail(extinction_prob)

ggsave(paste0('./figures/',figprefix, '_base_model_Pops.png'),dpi = 300,
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

ggsave(paste0('./figures/',figprefix, '_base_model_extinction.png'),dpi = 300,
       height = 3.8, width = 7, units = 'in') 
