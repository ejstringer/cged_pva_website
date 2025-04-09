

# load --------------------------------------------------------------------
em.extract_ASA <- function(sim){
  tstep <- dim(sim)[3]
  
  asa <- sapply(1:tstep, function(x) colSums(sim[2:6,,x,]))
  colnames(asa) <- 2013:(2012+tstep)
  rownames(asa) <- c('CA', 'JE', 'JW', 'MA')
  asa
}

em.sample.distribution <- function(x, samplesize = 10000){
  xhist=hist(x,freq=FALSE, col=rgb(0,0,1,1/4))
  #sample from it
  bins=with(xhist,sample(length(mids),samplesize,p=density,replace=TRUE)) # choose a bin
  result=runif(length(bins),xhist$breaks[bins],xhist$breaks[bins+1]) # sample a uniform in it
  hist(result,freq=FALSE,add=TRUE,bord=1, col = rgb(1,0,0,1/4))
  
  return(result)
}


source('./code/05_model.R')
param_starting_estimates <- readRDS('./output/starting_values.rds')
stage_distribution <- readRDS('./output/stage_distribution.rds')
transition_mat <- readRDS('./output/transition_matrix.rds')
K <- readRDS('./output/carrying_capcity.rds')
clutch_sizes <- 4:7

n_real <- read.csv('./output/real_abundance_for_comparison.csv')
n_real$ucl_all[n_real$N_all < 1 &n_real$ucl_all==0] <- n_real$N[n_real$N_all < 1 &n_real$ucl_all==0]*n_real$area_ha[n_real$N_all < 1 &n_real$ucl_all==0]
# setup for run -----------------------------------------------------------

parameter_estimates <- list()

parameter_estimates[[1]] <- param_starting_estimates
no.runs <- c(1,1,1,1,1)*10000
topModels <- 100
nrow(param_starting_estimates)
head(param_starting_estimates)

system.time({for (r in 1:length(no.runs)) {
  
  
  # define parameters
  
  param_est <- parameter_estimates[[r]]


# run multiple runs -------------------------------------------------------


param_dist <- split(as.data.frame(param_est), rownames(param_est))
sapply(lapply(param_dist, function(param) param$survival_juv), class)
N_sim <- lapply(param_dist, 
                function(param) em.pva_simulator(populations = c('CA', 'JE', 'JW', 'MA'),
                                                 stages = c('J', 'SA','A1','A2','A3', 'A4'),
                                                 stage_distribution = stage_distribution, 
                                                 initial_ab = unlist(param[1,grep('init', colnames(param))]), # adult abundance for populations
                                                 survival = unlist(param[1,grep('ASA', colnames(param))]), # survival of adults and SA at sites
                                                 survival_J = param$survival_juv, # juvenile survival
                                                 env_stoch = unlist(param[1,grep('env_stoch', colnames(param_est))]), # sd on survival
                                                 transition_mat = transition_mat, # transition prob to SA
                                                 f_reproducing = param$F_reproduction, # proportion of females reproducing
                                                 clutch_sizes = clutch_sizes, # clutch size range   
                                                 K = K, # carrying capcity applied to adults and SA
                                                 time_steps = 11, # time
                                                 replicates = 1)
)
print(paste('model', r, 'run'))
# model df ----------------------------------------------------------------
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

# best models -------------------------------------------------------------



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

# new params from best models ---------------------------------------------
models_selected <- simdf_best$run[which(simdf_best$bestx)] %>% unique()
param_selected <- param_est[which(rownames(param_est) %in% models_selected),]

param_selected %>% nrow

apply(param_selected, 2, summary)
apply(param_selected, 2, summary)[4,]



# estimate new distributions ----------------------------------------------

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
  guides(alpha = FALSE)+
  scale_alpha(range = c(0.005,1))+
  theme(panel.grid = element_blank())

fig_val + scale_y_log10(labels = label_comma())
ggsave(paste0('./figures/',figprefix, '_validation_models_runs.png'),dpi = 300,
       height = 6, width = 10, units = 'in')


models_selected <- simdf_best$run[which(simdf_best$bestx)] %>% unique()
param_selected <- param_est[which(rownames(param_est) %in% models_selected),]

param_selected %>% nrow

apply(param_selected, 2, summary)
apply(param_selected, 2, summary)[4,]

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

ggsave(paste0('./figures/',figprefix, '_validation_models_mean_squred_diff.png'),dpi = 300,
       height = 4, width = 4, units = 'in')
# parameter optermisation ------------------------------------------------
parameter_estimates %>% length
parameter_estimates[[r+2]] <- param_selected
lapply(parameter_estimates, colnames)
parameter_estimates2 <- parameter_estimates
for (j in 1:length(parameter_estimates)) parameter_estimates2[[j]] <- cbind(parameter_estimates[[j]], round = j-1)

em.parameters_to_df <- function(param_est){
  
  param_est %>% as.data.frame %>%
    pivot_longer(cols = N_init_CA:env_stoch_MA) %>% 
    mutate(type = ifelse(grepl('init', name), 'inititial N', 'value'),
           type = ifelse(grepl('quality', name), 'site quality', type),
           type = ifelse(name %in% c('survival_juv', 'F_reproduction'),
                         'x','site'),
           round = as.character(round)) 
  
}

paramer_optermisation <- lapply(parameter_estimates2, em.parameters_to_df) %>% 
  do.call('rbind', .)


paramer_optermisation %>% 
  ggplot(aes(y=round, x=value, fill = round))+
  geom_boxplot()+
  facet_wrap(type~name, scale = 'free')+
  scale_fill_manual(values = virid(9)[3:9],
                    name = 'Simulation round')+
  theme_bw()+
  theme(legend.position = 'inside',
        legend.direction = 'horizontal',
        legend.position.inside = c(0.75,0.1),
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold'),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank())

ggsave(paste0('./figures/',figprefix, '_validation_parameter_estimation.png'),dpi = 300,
       height = 6, width = 10, units = 'in')

apply(param_selected, 2, summary)

corparams <- data.frame(cor(param_selected))
corparams[which(abs(corparams) > 0.5)]

c5 <- apply(corparams, 2, function(x) ifelse(length(which(abs(x)>0.5 & x<1))==0,0,
                                       which(abs(x)>0.5 & x<1))) %>% 
  .[. >0]
corparams[c5, names(c5)]

data.frame(dist = as.dist(cor(param_selected))) %>% 
  ggplot(aes(dist))+
  geom_histogram(bins = 30, colour = 'black', fill = 'lightgreen')+
  scale_x_continuous(limits = c(-0.8,0.8),
                     breaks = seq(-0.8, 0.8, 0.2))+
  theme_classic()


corrplot(param_selected)

