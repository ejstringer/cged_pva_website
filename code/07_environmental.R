
# load --------------------------------------------------------------------
figprefix <- ''
source('./code/00_libraries.R')
source('./code/05_model.R')
param_starting_estimates <- readRDS('./output/starting_values.rds')
stages <- names(readRDS('./output/stage_distribution.rds'))
transition_mat <- readRDS('./output/transition_matrix.rds')
base_params <- readRDS('./output/base_parameters.rds')
max_age <- base_params$age
K <- base_params$K
clutch_sizes <- base_params$clutch
transition <- base_params$transition
# survival_unlogit <- base_params$survival
# survival_logit_sd <- base_params$surv_sd

param_selected <- readRDS('./output/selected_25models_parameters.RDS')
param_selecteddf <-readRDS('./output/selected_25models_parameters_df.RDS')

n_real <- read.csv('./output/real_abundance_for_comparison.csv')
n_real$ucl_all[n_real$N_all < 1 &n_real$ucl_all==0] <- n_real$N[n_real$N_all < 1 &n_real$ucl_all==0]*n_real$area_ha[n_real$N_all < 1 &n_real$ucl_all==0]

# stage distribution ------------------------------------------------------

fecundity <- (mean(param_selecteddf$mean[param_selecteddf$name == 'F_reproduction'])*0.5)*5.5 
ASAsur <- mean(param_selecteddf$mean[grep('survival_A', param_selecteddf$name)])
Jt <- mean(param_selecteddf$mean[grep('survival_J', param_selecteddf$name)])

stage.mat <- matrix(0, 
                    nrow = length(stages),
                    ncol = length(stages),
                    byrow = TRUE,
                    dimnames = list((stages),
                                    (stages)))

# fecundity
fecund_ages <- which(stages %in% paste0('A', 1:(max_age-1)))
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



# model inputs ------------------------------------------------------------
top25 <- param_selecteddf
paramlist <- list(populations = c('CA', 'JE', 'JW', 'MA'),
                  initial_ab = round(top25$mean[2:5]), # adult abundance for populations
                  survival = top25$mean[6:9], # survival of adults and SA at sites
                  survival_J = top25$mean[10:13], # juvenile survival
                  #env_stoch = top25$mean[6:9], # sd on survival
                  f_reproducing = top25$mean[1], # proportion of females reproducing
                  clutch_sizes = clutch_sizes, # clutch size range   
                  K = K # carrying capcity applied to adults and SA
)






# fitting -----------------------------------------------------------------


# env prior ---------------------------------------------------------------
reps = 50000
# sd applied to survival yearly  (logit)
stoch <- seq(0,3,0.1)
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(x = NULL, y = NULL, xlim = c(0,3), ylim = c(0,1), ylab = 'survival', xlab = 'sd')
abline(h = 0.975, col = 'grey', lty = 2)
abline(h = 0.025, col = 'grey', lty = 2)
abline(h = 0, col = 'grey')
abline(h = 1, col = 'grey')
for (i in 1:4) {
site <- i
m2 <- sapply(stoch, function(x) unlogit(qnorm(c(0.025), mean = logit(paramlist$survival[site]), sd = x)))
m3 <- sapply(stoch, function(x) unlogit(qnorm(c(0.975), mean = logit(paramlist$survival[site]), sd = x)))
lines(c(stoch), c(m2), col = i)
lines(c(stoch), c(m3), col = i)
}


plot(x = NULL, y = NULL, xlim = c(0,3), ylim = c(-8,8), ylab = 'survival', xlab = 'sd')
abline(h = logit(0.975), col = 'grey', lty = 2)
abline(h = logit(0.025), col = 'grey', lty = 2)

abline(h = logit(0.995), col = 'pink', lty = 3)
abline(h = logit(0.005), col = 'pink', lty = 3)

for (i in 1:4) {
  site <- i
  m2 <- sapply(stoch, function(x) (qnorm(c(0.025), mean = logit(paramlist$survival[site]), sd = x)))
  m3 <- sapply(stoch, function(x) (qnorm(c(0.975), mean = logit(paramlist$survival[site]), sd = x)))
  lines(c(stoch), c(m2), col = i)
  lines(c(stoch), c(m3), col = i)
}



env_lower <- 0
env_upper <- 2 # logit scale


sample_data <- randomLHS(reps, length(paramlist$populations))
colnames(sample_data) <-paste0('env_stoch_', paramlist$populations)
env_stochasticity <- (env_lower + sample_data * (env_upper - env_lower))

param_est <- env_stochasticity
set.seed(54614)
rownames(param_est) <- paste0(
  stri_rand_strings(nrow(sample_data), 4, "[A-Z]"),
  stri_rand_strings(nrow(sample_data), 4, "[0-9]"))

head(param_est)
param_starting_estimates <- param_est

## setup for run -----------------------------------------------------------

parameter_estimates <- list()

parameter_estimates[[1]] <- param_starting_estimates
no.runs <- c(1,1,1,2,1)*10000
topModels <- 100
nrow(parameter_estimates[[1]])
head(parameter_estimates[[1]])
n_cores <- 10

# run model ---------------------------------------------------------------

system.time({for (r in 1:length(no.runs)) {
  
  
  # define parameters
  
  param_est <- parameter_estimates[[r]]
  
  
  ## run multiple runs -------------------------------------------------------
  
  
  param_dist <- split(as.data.frame(param_est), rownames(param_est))
  sapply(lapply(param_dist, function(param) param$env_stoch), class)
  N_sim <- mclapply(param_dist, 
                  function(param) em.pva_simulator(populations = c('CA', 'JE', 'JW', 'MA'),
                                                   stages = stages,
                                                   stage_distribution = stage_distribution, 
                                                   initial_ab = round(paramlist$initial_ab), # adult abundance for populations
                                                   survival = paramlist$survival, # survival of adults and SA at sites
                                                   survival_J = paramlist$survival_J, # juvenile survival
                                                   survival_logit_sd = NULL,
                                                   site_adjust = NULL,
                                                   env_stoch = unlist(param), # sd on survival
                                                   transition_mat = transition_mat, # transition prob to SA
                                                   f_reproducing = paramlist$f_reproducing, # proportion of females reproducing
                                                   clutch_sizes = clutch_sizes, # clutch size range   
                                                   K = K, # carrying capcity applied to adults and SA
                                                   time_steps = 11, # time
                                                   replicates = 1), mc.cores = n_cores
  )
  print(paste('model', r, 'run'))
  ## model df ----------------------------------------------------------------
  sim_adults <- mclapply(N_sim, em.extract_ASA, mc.cores = n_cores)
  
  
  simdf_list <- list()
  for (i in 1:length(sim_adults)) {
    sim <- sim_adults[[i]]
    Nsim <- do.call('c', lapply(1:nrow(sim), function(x) sim[x,]))
    simdf_list[[i]] <-  data.frame(Nsim, 
                                   run = names(sim_adults)[i],
                                   year = 2013:2023, 
                                   site = rep(rownames(sim),
                                              each = length(2013:2023)))
  }
  
  
  sim_sd <- do.call('rbind', simdf_list) %>%
    left_join(n_real, by = c('year', 'site')) %>% 
    #filter(Nsim > 0) %>% 
    group_by(run, site) %>% 
    # summarise(sd_real = sd(log(N_all+1), na.rm = T),
    #        sd_sim = sd(log(Nsim+1), na.rm = T)) %>% 
    summarise(sd_real = sd(log(N_all+1), na.rm = T),
              sd_sim = sd(log(Nsim+1), na.rm = T)) %>% 
    mutate(sd_diff = abs(sd_sim - sd_real))
  
  simdf <- do.call('rbind', simdf_list) %>%
    left_join(n_real, by = c('year', 'site')) %>% 
    left_join(sim_sd)
  
  ## best models -------------------------------------------------------------
  
  
  
  best_model_fits <- simdf %>% filter(complete.cases(sd_diff)) %>% 
    group_by(run) %>% 
    summarise(diff = mean(sd_diff)) %>% 
    mutate(best = diff == min(diff))
  
  m5perc <- quantile(best_model_fits$diff, topModels/length(simdf_list))
  
  simdf_best <- left_join(simdf, best_model_fits, by = c('run')) %>% 
    mutate(bestx = diff < m5perc)
  
  simdf_best %>% 
    ggplot(aes(y = diff))+
    geom_boxplot(fill = 'grey')+
    geom_boxplot(data = simdf_best[simdf_best$bestx,],
                 colour= 'blue', size = 1)+
    theme_classic()
  
  ## new params from best models ---------------------------------------------
  models_selected <- simdf_best$run[which(simdf_best$bestx)] %>% unique()
  param_selected <- param_est[which(rownames(param_est) %in% models_selected),]
  
  param_selected %>% nrow
  
  apply(param_selected, 2, summary)
  apply(param_selected, 2, summary)[4,]
  
  
  
  ## estimate new distributions ----------------------------------------------
  
  param_prior <- lapply(1:ncol(param_selected), function(x) param_selected[,x])
  names(param_prior) <- colnames(param_selected)
  param_est2 <- lapply(param_prior, 
                       em.sample.distribution, samplesize = no.runs[r]) %>% 
    do.call('cbind', .)

  
  set.seed(54616)
  rownames(param_est2) <- paste0(
    stri_rand_strings(nrow(param_est2), 4, "[A-Z]"),
    stri_rand_strings(nrow(param_est2), 4, "[0-9]"))
  param_est2 
  
  parameter_estimates[[r+1]] <- param_est2
}# end of r
}) #system time


# plot best models --------------------------------------------------------

m5perc <- quantile(best_model_fits$diff, 25/length(simdf_list))


simdf_best <- left_join(simdf, best_model_fits) %>% 
  mutate(bestx = diff < m5perc)

sim_sample <- simdf_best$run %in% unique(simdf_best$run)[sample(1:no.runs[4],10000)]

fig_val <- ggplot(n_real, aes(year, N_all, fill = site))+
  geom_ribbon(aes(ymax = ucl_all, ymin = lcl_all), alpha = 0.7)+
  theme_light()+
  geom_line()+
  facet_wrap(~site, ncol = 2)+
  scale_x_continuous(breaks = seq(2013,2025,2))+
  geom_line(data = simdf_best,
            aes(y = Nsim, group = run, alpha =  +(bestx),
                colour = factor(bestx)))+
  scale_fill_manual(values = viridis::viridis(6)[1:4],
                    name = 'Real N estimates')+
  scale_colour_manual(values = c('grey20', 'red'),
                      name = 'Top simulated N')+
  guides(alpha = 'none')+
  scale_alpha(range = c(0.005,1))+
  theme(panel.grid = element_blank())

fig_val + scale_y_log10(labels = label_comma())
ggsave(paste0('./figures/',figprefix, 'environmental_models_runs.png'),
       dpi = 300,
       height = 6, width = 10, units = 'in')

simdf_best_save <- filter(simdf_best, bestx)
saveRDS(simdf_best_save, './output/top_25_best_fitting_models_environmental.rds')

models_selected <- simdf_best$run[which(simdf_best$bestx)] %>% unique()
param_selected <- param_est[which(rownames(param_est) %in% models_selected),]

param_selected %>% nrow

apply(param_selected, 2, summary)

saveRDS(param_selected, './output/selected_25models_parameters_env_stoch.RDS')
apply(param_selected, 2, summary)[4,]

m5perc
simdf_best %>% 
  ggplot(aes(y = diff))+
  geom_boxplot(fill = 'grey')+
  geom_boxplot(data = simdf_best[simdf_best$bestx,],
               colour= 'blue', size = 1)+
  theme_classic()+
  ylim(0, max(simdf_best$diff))+
  theme(axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  xlab('Difference in standard deviations')

ggsave(paste0('./figures/',figprefix, 'environmental_models_sd_diff.png'),dpi = 300,
       height = 4, width = 4, units = 'in')
# parameter optermisation ------------------------------------------------
parameter_estimates %>% length
parameter_estimates[[r+2]] <- param_selected
lapply(parameter_estimates, colnames)
sapply(parameter_estimates, nrow)
parameter_estimates2 <- parameter_estimates
for (j in 1:length(parameter_estimates)) parameter_estimates2[[j]] <- cbind(parameter_estimates[[j]], round = j-1)

em.parameters_to_df <- function(param_est){

  
  param_est %>% as.data.frame %>%
    relocate(round) %>% 
    pivot_longer(-round) %>% 
    mutate(round = as.character(round)) 
  
}

paramer_optermisation <- lapply(parameter_estimates2, em.parameters_to_df) %>% 
  do.call('rbind', .)


paramer_optermisation$name %>% table
n_p <- length(parameter_estimates)-1
paramer_optermisation %>% filter(round != n_p) %>%  
  mutate(round_reordered = (n_p-1) - as.numeric(round)) %>% 
  ggplot(aes(y=round_reordered, x=value, fill = round))+
  geom_boxplot()+
  facet_wrap(~name, scale = 'free_y', ncol = 2)+
  scale_fill_manual(values = adegenet::virid(n_p+3)[4:(n_p+3)],
                    name = 'Simulation round')+
  theme_bw()+
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.position.inside = c(0.75,0.06),
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold'),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank())

ggsave(paste0('./figures/',figprefix, 'environmental_parameter_estimation.png'),
       dpi = 300, height = 6, width = 6, units = 'in')

paramer_optermisation %>% filter(round == n_p) %>% 
  group_by(name) %>% summarise(mean = mean(value)) %>% 
  as.data.frame %>% 
  saveRDS('./output/selected_25models_parameters_df_environment.RDS')
apply(param_selected, 2, summary)[4,]



# survival plot -----------------------------------------------------------
paramlist$survival

surv_plot <- data.frame(survival = logit(paramlist$survival), 
           sd = apply(param_selected, 2, summary)[4,]) 

qp <- seq(0.001, 0.999, 0.001)

mat_surv <- sapply(1:4, function(x) dnorm(logit(qp), 
                              mean = surv_plot$survival[x],
                              sd = surv_plot$sd[x])) %>% 
  as.data.frame


colnames(mat_surv) <- rownames(surv_plot)

mat_surv %>% mutate(survival = qp,survival_logit= logit(qp)) %>% 
  pivot_longer(env_stoch_CA:env_stoch_MA) %>% 
  pivot_longer(survival:survival_logit, names_to = 'x', values_to = 'surv') %>% 
  mutate(site = sub('env_stoch_', '', name)) %>% 
  ggplot(aes(x = surv, y = value, colour = site))+
  geom_line()+
  facet_grid(site~x, scale = 'free')+
  theme_bw()+
  theme(strip.background.y = element_rect(fill = 'grey50', linewidth = 0),
        strip.background.x =  element_blank(),
        strip.text.y.right = element_text(angle = 0, face = 'bold', colour = 'white'),
        strip.text.x.top = element_text(face = 'bold', size = 12),
        panel.grid = element_blank(),
        legend.position = 'none')+
  xlab('Survival')+
  ylab('Density')

ggsave(paste0('./figures/',figprefix, 
              'environmental_stoch_survival_distribution.png'),
       dpi = 300, height = 6, width = 6, units = 'in')



