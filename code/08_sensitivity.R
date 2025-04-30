
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

eigen_analysis <- function(A) {
  ev <- eigen(A)
  lmax <- which(Re(ev$values) == max(Re(ev$values)))
  lambda <- Re(ev$values[lmax])
  W <- ev$vectors
  w <- abs(Re(W[, lmax]))
  stable_age <- w/sum(w)
  return(list(lambda, stable_age))
}



# load --------------------------------------------------------------------
figprefix <- 'ntg100m'
source('./code/00_libraries.R')
source('./code/05_model.R')
param_starting_estimates <- readRDS('./output/starting_values.rds')
stage_distribution <- readRDS('./output/stage_distribution.rds')
transition_mat <- readRDS('./output/transition_matrix.rds')
K <- readRDS('./output/carrying_capcity.rds')
clutch_sizes <- 4:7
param_selected <- readRDS('./output/selected_25models_parameters.RDS')
param_selecteddf <-readRDS('./output/selected_25models_parameters_df.RDS')
# stage distribution ------------------------------------------------------



fecundity <- (mean(param_selecteddf$mean[param_selecteddf$name == 'F_reproduction'])*0.5)*5.5 
ASAsur <- mean(param_selecteddf$mean[grep('survival_A', param_selecteddf$name)])
Jt <- mean(param_selecteddf$mean[grep('survival_J', param_selecteddf$name)])

stage.mat <- matrix(c(0, 0, rep(fecundity,3),0,
                      Jt*0.194, 0, 0, 0,0,0,
                      Jt*0.806, 0,0,0,0,0,
                      0, ASAsur,ASAsur, 0,0,0,
                      0, 0, 0,ASAsur,0,0,
                      0, 0, 0,0,ASAsur,0), nrow = 6, ncol = 6, byrow = TRUE,
                    dimnames = list(c('J', 'SA','A1','A2','A3', 'A4'),
                                    c('J', 'SA','A1','A2','A3', 'A4')))

stage_distribution <- eigen_analysis(stage.mat)[[2]]
names(stage_distribution)<- c('J', 'SA','A1','A2','A3', 'A4')

# model inputs ------------------------------------------------------------
top25 <- param_selecteddf
paramlist <- list(populations = c('CA', 'JE', 'JW', 'MA'),
                  initial_ab = round(top25$mean[2:5]), # adult abundance for populations
                  survival = top25$mean[10:13], # survival of adults and SA at sites
                  survival_J = top25$mean[14:17], # juvenile survival
                  env_stoch = top25$mean[6:9], # sd on survival
                  f_reproducing = top25$mean[1], # proportion of females reproducing
                  clutch_sizes = clutch_sizes, # clutch size range   
                  K = K # carrying capcity applied to adults and SA
)



paramlist

# sensitivity -------------------------------------------------------------

original_paramlist <- paramlist

paramlist 

names(paramlist)

sim_summary <- NULL
seqencedf <- data.frame(parameter = rep(names(paramlist)[2:6], each = 2), 
                        ten = c(0.9, 1.1))
base_summary <- readRDS('./output/base_summary_r_ext_p.rds')
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

sim_N_sum <- N_simulated %>% 
  group_by(rep, tstep, pop) %>% 
  summarise(N = sum(N)) %>% 
  mutate(extinct = N < 2)

sum(sim_N_sum$extinct)

sim_adjusted_lambda <- sim_N_sum %>% 
  ungroup() %>% 
  group_by(pop, rep) %>% 
  mutate(diff_year = tstep - lag(tstep),
         diff_growth = N - lag(N),
         rate_percent = (diff_growth /diff_year)/ lag(N) * 100) %>%
  ungroup() %>% 
  group_by(tstep,pop) %>%
  summarise(sims = n(),
            sim0 = sum(extinct),
            ext_p = sum(extinct)/n(),
            rate = mean(rate_percent, na.rm=T))%>%
  ungroup() %>%
  group_by(pop) %>% 
  summarise(growth_rate = mean(rate, na.rm = T),
            extinction = mean(ext_p)) %>% 
  left_join(base_summary)

sim_adjusted_lambda$parameter <- param_adjust
sim_adjusted_lambda$adjusted <- adjust

sim_summary <- bind_rows(sim_summary, sim_adjusted_lambda)
}
sim_summary %>% 
  dplyr::select(-extinction, -extinction_base) %>% 
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
