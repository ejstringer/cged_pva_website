
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




# load --------------------------------------------------------------------
figprefix <- 'ntg100m'
source('./code/05_model.R')
param_starting_estimates <- readRDS('./output/starting_values.rds')
stage_distribution <- readRDS('./output/stage_distribution.rds')
transition_mat <- readRDS('./output/transition_matrix.rds')
K <- readRDS('./output/carrying_capcity.rds')
clutch_sizes <- 4:7
param_selected <- readRDS('./output/selected_25models_parameters.RDS')

# stage distribution ------------------------------------------------------

eigen_analysis <- function(A) {
  ev <- eigen(A)
  lmax <- which(Re(ev$values) == max(Re(ev$values)))
  lambda <- Re(ev$values[lmax])
  W <- ev$vectors
  w <- abs(Re(W[, lmax]))
  stable_age <- w/sum(w)
  return(list(lambda, stable_age))
}

fecundity <- (mean(param_selected[,'F_reproduction'])*0.5)*5.5 
ASAsur <- mean(param_selected[,grep('survival_A', colnames(param_selected))])
Jt <- mean(param_selected[,'survival_juv'])
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

as.data.frame(stage.mat) %>% 
  mutate_if(is.numeric, round, 3) %>% 
  mutate(stages = rownames(.)) %>% 
  rbind(c(round(stage_distribution,2),'stable:')) %>% 
  rbind(c(round(eigen_analysis(stage.mat)[[1]],3), rep(NA,5), 'lambda =')) %>% 
  relocate(stages) %>% 
  flextable() %>% 
  autofit() %>% 
 # border_remove() %>%
  vline(j = 1, i = 1:6, part = 'body', border = fp_border_default(width = 1.5)) %>% 
  vline(j = 1, part = 'header', border = fp_border_default(width = 1.5)) %>% 
  hline(i = 6, border = fp_border_default(width = 1.5)) %>% 
  hline(i = 7, j=1:2, border = fp_border_default(width = 1.5)) %>% 
  vline(i = 8, j = c(2),border = fp_border_default(width = 1.5), part = 'body') %>% 
  align(j = 1,  align = 'right') %>% 
  align(part = 'header', j  = 1,align = 'right') %>% 
  hline(i = 7, border = fp_border_default(width = 1.5)) %>% 
  hline_bottom(j = 3:7,border = fp_border_default(width = 0)) %>%  
  vline_left(i = 8, part = 'body', border = fp_border_default(width = 1.5)) %>% 
  bold(part = 'header') %>% 
  bold(j = 2, i = 8, part = 'body') %>% 
  bold(j= 1) %>% 
  font(part = 'all', fontname='calibri') %>% 
  saveRDS('./output/new_stage_distribution.rds')


# model inputs ------------------------------------------------------------

paramlist <- list(populations = c('CA', 'JE', 'JW', 'MA'),
initial_ab = round(top25[1:4]), # adult abundance for populations
survival = top25[6:9], # survival of adults and SA at sites
survival_J = top25[5], # juvenile survival
env_stoch = top25[11:14], # sd on survival
f_reproducing = top25[10], # proportion of females reproducing
clutch_sizes = clutch_sizes, # clutch size range   
K = K # carrying capcity applied to adults and SA
)



paramlist %>% 
  do.call('cbind', .) %>% 
  as.data.frame %>% dplyr::select(-clutch_sizes) %>% 
  pivot_longer(initial_ab:K) %>% 
  mutate(value = as.numeric(value),
         value = round(value, 2),
         name = ifelse(name =='survival', 'survival_A_SA', name),
         name = ifelse(name == 'initial_ab', 'N_initial', name),
         name = gsub('_', ' ', name)) %>% 
  pivot_wider(names_from = populations, values_from = value) %>% 
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
