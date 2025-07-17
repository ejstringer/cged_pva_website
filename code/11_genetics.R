figprefix <- ''
source('./code/00_libraries.R')
source('./code/05_model.R')
library(dartRverse)
# genetics ----------------------------------------------------------------

snps <- readRDS('./data/wild_unique_genlight_tympos.rds') %>% 
  gl.filter.monomorphs() %>% 
  gl.filter.maf(threshold = 0.2) %>% 
  gl.filter.callrate(threshold = 1)

table(factor(snps@other$ind.metrics$pop),snps@pop)

snps <- gl.keep.pop(snps, pop.list = c(2009:2016))
snps <- gl.keep.loc(snps, loc.list = sample(snps@loc.names, 100))


pop(snps) <- str_sub(snps@other$ind.metrics$grid, 1,2)
snps@pop %>% table

# parameters --------------------------------------------------------------

stage_distribution <- readRDS('./output/stage_distribution_base.rds')

stages <- names(readRDS('./output/stage_distribution.rds'))
transition_mat <- readRDS('./output/transition_matrix.rds')
base_params <- readRDS('./output/base_parameters.rds')
fper_repro <- base_params$fecundity
max_age <- base_params$age
K <- base_params$K
clutch_sizes <- base_params$clutch
transition <- base_params$transition

param_selecteddf <-readRDS('./output/selected_25models_parameters_df.RDS')
param_envdf <- readRDS('./output/selected_25models_parameters_df_environment.RDS')

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

t_steps <- 50
supp_every <- 1
supplement <- rep(c(F, T),10)
seq_suplementation <- seq(14,t_steps, supp_every)
n_supplement <- c(2, 10, 16, 10)

n_cores <- 10

genetics <- T
# PVA ---------------------------------------------------------------------
reps_base <- 1

 
system.time(N_sim_genetics <- mclapply(supplement, 
                          function(suppTF) em.pva_simulator(populations = paramlist$populations,
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
                          time_steps = t_steps, # time
                          replicates = reps_base,
                          supp = suppTF,
                          stage.supp = 1,
                          n.supp = n_supplement, 
                          when.supp = seq_suplementation,
                          GENETICS = genetics,
                          snps = snps), mc.cores = n_cores))
sim_genetics <- lapply(N_sim_genetics, function(x) x$G)

sim_genetics

em.extract_G <- function(sim_G, run){

  

    sim_G %>% as.data.frame() %>% 
  mutate(pop = paramlist$populations) %>% 
  pivot_longer(-pop, names_to = 'generation',
               values_to = 'diversity') %>% 
  mutate(generation = as.numeric(sub('V', '', generation)),
         beforesupp = ifelse(generation < min(seq_suplementation),
                             'yes', 'no')) %>%
  group_by(pop) %>% 
  mutate(percent = diversity/diversity[generation == 1],
         run = run) 

}


# summarise ---------------------------------------------------------------
sim_geneticslist <- lapply(1:length(sim_genetics), function(x) em.extract_G(sim_genetics[[x]], run = x))
sim_geneticsdf <- do.call('rbind', sim_geneticslist) %>% 
  mutate(supplemented = (run %% 2)== 0)%>% 
    mutate(supplemented = ifelse(supplemented, 'Supplemented', 'Not supplemented'),
           diversity2 = ifelse(is.na(diversity), 0, diversity),
           #beforesupp = ifelse(is.na(diversity), 'N = 0', beforesupp),
           `Before supplementation` = beforesupp)
  table(sim_geneticsdf$run)/50/4
  
  

# plot --------------------------------------------------------------------

  
  Heplot <-  sim_geneticsdf %>%    
    ggplot(aes(generation, diversity2, colour = `Before supplementation`))+
    geom_line(aes(group = run), linewidth =  0.2, alpha = 0.5)+
    geom_point(size = 0.5, alpha = 0.5)+
    geom_smooth(data = filter(sim_geneticsdf, complete.cases(diversity)),
                method = 'lm', colour = 'black', 
                aes(linetype = `Before supplementation`), show.legend = F)+
    facet_grid(pop~supplemented)+
    theme_bw()+
    ylab('Expected Heterozygosity')+
    theme(legend.position = 'bottom',
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.background.y = element_rect(fill = 'grey90', linewidth = 0),
          strip.text.y = element_text(angle = 0),
          strip.text = element_text(face = 'bold')); Heplot 

  ggsave(paste0('./figures/',figprefix, 'genetics_heterozygosity_supplementation.png'),
         plot = Heplot, dpi = 300,
         height = 15, width = 13, units = 'cm')  
  