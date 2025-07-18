figprefix <- ''

# load --------------------------------------------------------------------
source('./code/00_libraries.R')
source('./code/05_model.R')
param_starting_estimates <- readRDS('./output/starting_values.rds')
param_starting_estimates <- param_starting_estimates[,-ncol(param_starting_estimates)]
stage_distribution <- readRDS('./output/stage_distribution.rds')
transition_mat <- readRDS('./output/transition_matrix.rds')
base_params <- readRDS('./output/base_parameters.rds')
fper_repro <- base_params$fecundity
K <- base_params$K
clutch_sizes <- base_params$clutch
survival_unlogit <- base_params$survival
survival_logit_sd <- base_params$surv_sd

n_real <- read.csv('./output/real_abundance_for_comparison.csv')
n_real$ucl_all[n_real$N_all < 1 &n_real$ucl_all==0] <- n_real$N[n_real$N_all < 1 &n_real$ucl_all==0]*n_real$area_ha[n_real$N_all < 1 &n_real$ucl_all==0]
# setup for run -----------------------------------------------------------

parameter_estimates <- list()

parameter_estimates[[1]] <- param_starting_estimates
no.runs <- c(1,1,1,2,1)*10000
topModels <- 100
nrow(parameter_estimates[[1]])
head(param_starting_estimates)

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
                                                 stages = names(stage_distribution),
                                                 stage_distribution = stage_distribution, 
                                                 initial_ab = unlist(param[1,grep('init', colnames(param))]), # adult abundance for populations
                                                 survival = survival_unlogit, # survival of adults and SA at sites
                                                 survival_J = param$survival_juv,# juvenile survival
                                                 survival_logit_sd = survival_logit_sd,
                                                 site_adjust = unlist(param[1,grep('adjust', colnames(param))]),
                                                 env_stoch = NULL, # sd on survival
                                                 transition_mat = transition_mat, # transition prob to SA
                                                 f_reproducing = fper_repro,#param$F_reproduction, # proportion of females reproducing
                                                 clutch_sizes = clutch_sizes, # clutch size range   
                                                 K = K, # carrying capcity applied to adults and SA
                                                 time_steps = 11, # time
                                                 replicates = 1), mc.cores = n_cores 
)
print(paste('model', r, 'run'))
## model df ----------------------------------------------------------------
sim_adults <- lapply(N_sim, em.extract_ASA)


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


simdf <- do.call('rbind', simdf_list) %>%
  left_join(n_real, by = c('year', 'site')) %>% 
  mutate(squared_diff = sqrt((N_all - (Nsim))^2))

## best models -------------------------------------------------------------



best_model_fits <- simdf %>% filter(complete.cases(squared_diff)) %>% 
  group_by(run) %>% 
  summarise(diff = mean(squared_diff)) %>% 
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
param_est2[,1:4] <- round(param_est2[,1:4])

set.seed(54614)
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
ggsave(paste0('./figures/',figprefix, 'validation_models_runs.png'),dpi = 300,
       height = 6, width = 10, units = 'in')

simdf_best_save <- filter(simdf_best, bestx)
saveRDS(simdf_best_save, './output/top_25_best_fitting_models.rds')

models_selected <- simdf_best$run[which(simdf_best$bestx)] %>% unique()
param_selected <- param_est[which(rownames(param_est) %in% models_selected),]

param_selected %>% nrow

apply(param_selected, 2, summary)

saveRDS(param_selected, './output/selected_25models_parameters.RDS')
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
  xlab('Mean squared difference')

ggsave(paste0('./figures/',figprefix, 'validation_models_mean_squred_diff.png'),dpi = 300,
       height = 4, width = 4, units = 'in')
# parameter optermisation ------------------------------------------------
parameter_estimates %>% length
parameter_estimates[[r+2]] <- param_selected
lapply(parameter_estimates, colnames)
sapply(parameter_estimates, nrow)
parameter_estimates2 <- parameter_estimates
for (j in 1:length(parameter_estimates)) parameter_estimates2[[j]] <- cbind(parameter_estimates[[j]], round = j-1)

em.parameters_to_df <- function(param_est){
  adjusted_param <- param_est[,grep('adjust', colnames(param_est))]
  back_transform_survival <- unlogit(qnorm(adjusted_param, mean = logit(survival_unlogit), sd = survival_logit_sd))
  colnames(back_transform_survival) <- sub('adjust','ASA', colnames(back_transform_survival))
  
  juv_surv <- param_est[,grep('_juv', colnames(param_est))]
  
  juv_site_logit <- apply(adjusted_param, 2, function(x) qnorm(x, mean = logit(juv_surv), sd = survival_logit_sd))
  juv_site_survival <- unlogit(juv_site_logit)
  colnames(juv_site_survival) <- sub('adjust','J', colnames(juv_site_survival))
  param_est2 <- cbind(param_est, back_transform_survival)
  param_est3 <- cbind(param_est2, juv_site_survival)
  param_est3 %>% as.data.frame %>%
    relocate(round) %>% 
    pivot_longer(-round) %>% 
    mutate(type = ifelse(grepl('init', name), 'inititial N', 'value'),
           type = ifelse(grepl('quality', name), 'site quality', type),
           type = ifelse(name %in% c('survival_juv', 'F_reproduction'),
                         'x','site'),
           round = as.character(round)) 
  
}

paramer_optermisation <- lapply(parameter_estimates2, em.parameters_to_df) %>% 
  do.call('rbind', .)


paramer_optermisation$name %>% table
n_p <- length(parameter_estimates)-1
paramer_optermisation %>% filter(round != n_p, 
                                 !grepl('adjust', name),
                                 name != 'survival_ju', 
                                 name != 'F_reproductio') %>%  
  mutate(round_reordered = 5 - as.numeric(round),
         name = ifelse(name == 'F_reproduction',
                       'xF_reproduction', name),
         name = ifelse(name == 'survival_juv',
                       'xSurvival_juv', name)) %>% 
  ggplot(aes(y=round_reordered, x=value, fill = round))+
  geom_boxplot()+
  facet_wrap(~name, scale = 'free', ncol = 4)+
  scale_fill_manual(values = adegenet::virid(n_p+3)[4:(n_p+3)],
                    name = 'Simulation round')+
  theme_bw()+
  theme(legend.position = 'inside',
        legend.direction = 'horizontal',
        legend.position.inside = c(0.75,0.09),
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold'),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank())

ggsave(paste0('./figures/',figprefix, 
              'validation_parameter_estimation.png'),dpi = 300,
       height = 7, width = 8, units = 'in')

paramer_optermisation %>% filter(round == '6') %>% 
  group_by(name) %>% summarise(mean = mean(value)) %>% 
  as.data.frame %>% 
  saveRDS('./output/selected_25models_parameters_df.RDS')
apply(param_selected, 2, summary)[4,]



as.dist(cor(param_selected)) %>% 
  otuSummary::matrixConvert(colname = c('parameter 1', 'parameter 2', 'r')) %>% 
  filter(abs(r)>0.4) %>% 
  mutate(r = round(r, 2)) %>% 
  arrange(abs(r)) %>% 
  flextable() %>% 
  autofit() %>% 
  theme_alafoli() %>% 
  bg(bg = 'honeydew') %>% 
  bold(part = 'header') %>% 
  font(fontname = 'Calibri', part = 'all') %>% 
  saveRDS('./output/top_correlated_params.RDS')
fx_r <- readRDS('./output/top_correlated_params.RDS')
fx_r %>% save_as_image(paste0('./figures/',figprefix, 
                              'validation_parameter_cor_top.png'),
                       res=300)

data.frame(dist = as.dist(cor(param_selected))) %>% 
  ggplot(aes(dist))+
  geom_histogram(bins = 30, colour = 'black', fill = 'lightgreen')+
  scale_x_continuous(limits = c(-0.8,0.8),
                     breaks = seq(-0.8, 0.8, 0.2))+
  xlab('r')+
  theme_classic()

ggsave(paste0('./figures/',figprefix, 'validation_parameter_correlation.png'),
       dpi = 300,
       height = 6,
       width = 10,
       units = 'in')


