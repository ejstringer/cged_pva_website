
# load --------------------------------------------------------------------
source('./code/00_libraries.R')

logit <- function(p) log(p/(1-p))
unlogit <- function(x) exp(x) / (1 + exp(x))

area_used <- read.csv('./output/area_used.csv')
n_sites <- read.csv('./output/site_abundance.csv')

N_start <- n_sites %>% 
  filter(year == 2013) %>% 
  left_join(area_used)%>%
  dplyr::select(-buffer) %>% 
  rename(n_grids = n) %>% 
  mutate(site = c('CA', 'JE', 'JW', 'MA'))

# breeding
mel2024 <- read.csv('./data/melbourn_zoo_tympo_breeding_2024.csv') %>% 
  filter(complete.cases(Number.of.clutches))

# K 
caps <- read.csv('./data/N_estimates.csv')

# replicates --------------------------------------------------------------

reps <- 50000

# initial N ---------------------------------------------------------------
N_start

sample_q <- lapply(1:nrow(N_start), function(x)
                  qnorm(runif(50000,
                  pnorm(N_start$lcl[x], mean = N_start$N[x], sd = N_start$se_sum[x]),
                  pnorm(N_start$ucl[x], mean = N_start$N[x], sd = N_start$se_sum[x]))))
range(sample_q)
sample_N <- lapply(1:nrow(N_start), 
                   function(x) (sample_q[[x]]*N_start$se_sum[x])+N_start$N[x])

sample_N_all <- lapply(1:nrow(N_start), 
                       function(x) round((sample_N[[x]]/N_start$n_grids[x])*N_start$area_ha[x]))

initial_N <- sample_N_all %>% do.call('cbind', .)
colnames(initial_N) <- paste('N_init', N_start$site, sep = '_')

initial_N %>% head

# survival ----------------------------------------------------------------

## Juvenile ----------
juv0_lower <- 0.15 # wendy estimates 
juv0_upper <- 0.6 # wendy winter survival?

sample_data <- randomLHS(reps, 1)
colnames(sample_data) <- 'survival_juv'
juv0 <- (juv0_lower + sample_data * (juv0_upper - juv0_lower))

## site adjust ------------------

survival_logit <- (-0.5738307)
survival_unlogit <- unlogit(-0.5738307)
survival_logit_sd <- 0.2515

sample_data <- randomLHS(reps, nrow(N_start))

site_survival <- sapply(1:nrow(N_start), 
                          function(x) unlogit(qnorm(sample_data[,x], 
                                            mean = survival_logit, 
                                            sd = survival_logit_sd)))

head(site_survival)
site_survival <- sample_data
class(site_survival)
colnames(site_survival) <- paste('survival_adjust', N_start$site, sep = '_')

site_survival %>% head

#par(mfrow = c(2,2)); apply(site_survival, 2, hist);par(mfrow = c(1,1))

# transition --------------------------------------------------------------

transition  = 0.194
#transition_sd = 0.071 # fixed parameter


# clutch size -------------------------------------------------------------

clutch_sizes <- 3:7

prob_clutch <- dbeta(seq(0.1,0.9, length.out = length(clutch_sizes)), 2,3)
prob_clutch/sum(prob_clutch)

mean_clutch <- mean(sample(clutch_sizes, 100000, T, prob_clutch))

# reproductive females ----------------------------------------------------

repro_fixed <- T # 0.65 from captive estimates

# fixed
femaleHatches <- mel2024 %>% 
  group_by(Female) %>% 
  summarise(n_clutches = sum(total.offspring.successfully.hatched))

tb_mel <- table(femaleHatches$n_clutches)

females2024 <- mel2024 %>% 
  group_by(Female) %>% 
  summarise(n_clutches = sum(Number.of.clutches))

tb_mel <- table(females2024$n_clutches)


mel_perc_repro <- 1 - (tb_mel[1]/sum(tb_mel))

# not fixed 

repro_lower <- 0
repro_upper <- 1

sample_data <- randomLHS(reps, 1)
colnames(sample_data) <- 'F_reproduction'
reproduction <- (repro_lower + sample_data * (repro_upper - repro_lower))


# carrying capacity -------------------------------------------------------

# max density 
max_density <- caps %>% 
  group_by(site, year) %>% 
  summarise(N = sum(N), n = n(), Ng = N/n()) %>%
  arrange(desc(N)) %>% 
  ungroup() %>% 
  filter(Ng == max(Ng))

K <- round(max_density$Ng*N_start$area_ha)
names(K) <-  paste('K', N_start$site, sep = '_')


# environmental stoch -----------------------------------------------------
# sd applied to survival yearly  (logit)


# new env stoch after other parameters
env_lower <- 0
env_upper <- 2.5 # logit scale


sample_data <- randomLHS(reps, length(N_start$site))
colnames(sample_data) <-paste0('env_stoch_', N_start$site)
env_stochasticity <- (env_lower + sample_data * (env_upper - env_lower))
# variable params ---------------------------------------------------------


param_est <- cbind(initial_N, 
                   juv0,
                   site_survival,
                   reproduction#, 
                   #env_stochasticity
                   )

if(repro_fixed) param_est[,'F_reproduction'] <- mel_perc_repro
set.seed(54612)
rownames(param_est) <- paste0(
  stri_rand_strings(nrow(sample_data), 4, "[A-Z]"),
  stri_rand_strings(nrow(sample_data), 4, "[0-9]"))

head(param_est)
param_starting_estimates <- param_est

saveRDS(param_starting_estimates, './output/starting_values.rds')

# parameter table ---------------------------------------------------------

parameter_table <- data.frame(site = N_start$site,
                          N = round(apply(initial_N, 2, mean)),
                          N_sd = round(apply(initial_N, 2, sd)),
                          N_lower = round(apply(initial_N, 2, min)),
                          N_upper = round(apply(initial_N, 2, max)),
                          survivalJ_lower = juv0_lower,
                          survivalJ_upper= juv0_upper,
                          survivalJ_sd = (survival_logit_sd),
                          survivalSA = unlogit(survival_logit),
                          survivalSA_sd = (survival_logit_sd),
                          survivalA = unlogit(survival_logit),
                          survivalA_sd = (survival_logit_sd),
                          transition,
                          reproduction_lower = repro_lower,
                          reproduction_upper = repro_upper,
                          clutches = round(mean_clutch,1),
                          clutches_lower = min(clutch_sizes),
                          clutches_upper = max(clutch_sizes),
                          K = K,
                          environmental_lower = env_lower,
                          environmental_upper = env_upper) %>% 
  pivot_longer(N:environmental_upper,
               names_to = 'Parameter',values_to = 'estimate') %>% 
  separate(Parameter, into = c('Parameter', 'type'), sep = '_') %>% 
  mutate(type = ifelse(is.na(type), 'value', type)) %>%
  pivot_wider(names_from = type, values_from = estimate) %>% 
  arrange(Parameter, site) %>% 
  mutate(site = ifelse(Parameter %in% c('K', 'N'), site, 'All'),
         type = ifelse(is.na(value), 'variable', 'fixed'),
         type = ifelse(Parameter %in% c('survivalJ', 'survivalSA','survivalA'),
                       'site variable', type)) %>% 
  rename_all(tolower) %>% 
  mutate_if(is.numeric, round,3) %>% 
  mutate(across(everything(), ~ replace(.x, is.na(.x), ""))) %>% 
  unique() %>% 
  mutate(parameter = factor(parameter, 
                            levels = c('N', 'survivalJ', 'survivalSA', 'survivalA',
                                       'transition', 'reproduction','clutches',
                                       'environmental','K')),
         type = ifelse(parameter == 'N', 'variable', type)) %>% 
  arrange(parameter)%>% 
  mutate(distribution = c(rep('truncated normal ', 4), 
                          'uniform/logit normal', rep('logit normal',2),
                          '', 'uniform', '~beta(2,3)', 'logit uniform',
                          rep('',4))) 
if (repro_fixed) {
  r_index <- parameter_table$parameter== 'reproduction'
  parameter_table$value[r_index] <- round(mel_perc_repro, 3)
  parameter_table$type[r_index]  <- 'fixed'
  col_index <- which(colnames(parameter_table) %in% c('distribution',
                                                      'lower','upper'))
  parameter_table[r_index, col_index] <- ''

}
paramter_fxtb <- parameter_table %>% 
  flextable() %>% 
autofit() %>%
  theme_zebra() %>% 
  color(i = which(grepl('variable', parameter_table$type)),j = 3:8,color = 'grey50') %>% 
  bold(i = which(grepl('variable', parameter_table$type)),j = 3:8) %>% 
  hline_bottom(border = fp_border_default(width = 2)) %>% 
  hline_top(border = fp_border_default(width = 1)) %>% 
  hline(i = c(4,11),border = fp_border_default(width = 0.9)) %>% 
  bold(part = 'header') %>% 
  font(part = 'all', fontname='calibri');paramter_fxtb

paramter_fxtb %>% save_as_image('./figures/parameter_priors_table.png',
                                res = 300)
# variable distributions ----------------------------
param.labs <- c('initial N', 'environmental stoch', 'survival', 'fecundity')
names(param.labs) <- c('N', 'e','s', 'F')



distribution_fig <- param_est %>% cbind(env_stochasticity) %>%
  as.data.frame %>% 
  rename(survival_J= survival_juv) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(type = factor(substr(name, 1,1),levels = c('N', 'e', 's','F')),
         variable = sub('N_init_', '', name),
         variable = sub('survival_', '', variable),
         variable = sub('env_stoch_', '', variable),
         variable = sub('F_reproduction', '%', variable)) %>% 
  ggplot(aes(x = value, fill = name))+
  geom_histogram(colour = 'black', bins = 20)+
  facet_wrap(type~variable, scale = 'free', 
             labeller = labeller(type = param.labs))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold'));distribution_fig

ggsave('./figures/parameter_priors.png',
       plot = distribution_fig, dpi = 300,
       height = 7, width = 10, units = 'in')



# stages ------------------------------------------------------------------
max_age <- 3 

stages <- c('J', 'SA',paste0('A', 1:max_age))

# stage matrix ------------------------------------------------------------


stage.mat <- matrix(0, 
                    nrow = max_age+2, ncol = max_age+2, 
                    byrow = TRUE,
                    dimnames = list(stages,stages))

# fecundity
fecund_ages <- which(stages %in% paste0('A', 1:(max_age-1)))
stage.mat[1,fecund_ages] <- 2.38
# juvenile survival + growth
stage.mat[2:3,which(stages == 'J')] <- round(survival_unlogit*c(transition, 1-transition),2)

# adult survival + growth
for (i in 2:(max_age+2)) {
  stg <- stages[i]
  if (stg == 'SA') stage.mat[i+2,i] <- round(survival_unlogit,2)
  if (stg != 'SA' & i < (max_age+2)) stage.mat[i+1, i] <- round(survival_unlogit,2)
}

stage.mat
stagemat_lambda <- eigen_analysis(stage.mat)[[1]]
stagemat_lambda
stage_distribution <- eigen_analysis(stage.mat)[[2]]

names(stage_distribution)<- stages

saveRDS(stage_distribution, './output/stage_distribution.rds')

fx_stages <- as.data.frame(stage.mat) %>% 
  mutate(stages = rownames(.)) %>% 
  rbind(c(round(stage_distribution,2),'stable:')) %>%
  relocate(stages) %>% 
flextable() %>% 
  autofit() %>% 
  vline(j = 1, i = 1:nrow(stage.mat), part = 'body', border = fp_border_default(width = 1.5)) %>% 
  vline(j = 1, part = 'header', border = fp_border_default(width = 1.5)) %>% 
  hline(i = nrow(stage.mat), border = fp_border_default(width = 1.5)) %>% 
  align(j = 1,  align = 'right') %>% 
  align(part = 'header', j  = 1,align = 'right') %>% 
  bold(part = 'header') %>% 
  bold(j= 1) %>% 
  font(part = 'all', fontname='calibri') %>% 
  add_footer_lines(paste('lambda:', round(stagemat_lambda,3))) %>%
  align(part = 'footer', align = 'right') %>% 
  font(part = "footer", fontname = "Consolas") %>% 
  bold(bold = T, part = "footer") %>%  
  set_caption(caption = paste('CGED Stage Matrix'));fx_stages

fx_stages %>% save_as_image('./figures/parameter_stage_matrix.png', res = 300)


transition_mat <- matrix(0, 
                    nrow = max_age+2, ncol = max_age+2, 
                    byrow = TRUE,
                    dimnames = list(stages,stages))

# juvenile transition
transition_mat[2:3,which(stages == 'J')] <- c(transition, 1-transition)

# adult survival + growth
for (i in 2:(max_age+2)) {
  stg <- stages[i]
  if (stg == 'SA') transition_mat[i+2,i] <- 1
  if (stg != 'SA' & i < (max_age+2)) transition_mat[i+1, i] <- 1
}

transition_mat


saveRDS(transition_mat, './output/transition_matrix.rds')

fx_transition<- transition_mat %>% 
  as.data.frame() %>% 
  mutate(stages = rownames(.)) %>% 
  relocate(stages) %>% 
  flextable() %>% 
  autofit() %>% 
  vline(j = 1, i = 1:nrow(transition_mat), part = 'body', border = fp_border_default(width = 1.5)) %>% 
  vline(j = 1, part = 'header', border = fp_border_default(width = 1.5)) %>% 
  align(j = 1,  align = 'right') %>% 
  align(part = 'header', j  = 1,align = 'right') %>% 
  bold(part = 'header') %>% 
  bold(j= 1) %>% 
  font(part = 'all', fontname='calibri')  %>% 
  set_caption(caption = 'CGED Transition Matrix');fx_transition

fx_transition %>% 
  save_as_image('./figures/parameter_transitions.png', res = 300)


# save parameters ---------------------------------------------------------

saveRDS(list(age = max_age, K = K, survival = survival_unlogit, 
             surv_sd = survival_logit_sd,
             transition = transition,
             clutch = clutch_sizes, beta = c(2,3),
             fecundity = mel_perc_repro), 
        './output/base_parameters.rds')


