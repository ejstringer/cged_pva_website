
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

## site adult ------------------

survival_logit <- -0.5738307
survival_logit_sd <- 0.2515

sample_data <- randomLHS(reps, nrow(N_start))

site_survival <- sapply(1:nrow(N_start), 
                          function(x) unlogit(qnorm(sample_data[,x], 
                                            mean = survival_logit, 
                                            sd = survival_logit_sd)))
head(site_survival)
class(site_survival)
colnames(site_survival) <- paste('survival_ASA', N_start$site, sep = '_')

site_survival %>% head

#par(mfrow = c(2,2)); apply(site_survival, 2, hist);par(mfrow = c(1,1))

# transition --------------------------------------------------------------

transition  = 0.194
#transition_sd = 0.071 # fixed parameter


# clutch size -------------------------------------------------------------

clutch_sizes <- 4:7

# reproductive females ----------------------------------------------------

repro_lower <- 0.98
repro_upper <- 1

sample_data <- randomLHS(reps, 1)
colnames(sample_data) <- 'F_reproduction'
reproduction <- (repro_lower + sample_data * (repro_upper - repro_lower))


# carrying capacity -------------------------------------------------------

K <- 50*N_start$area_ha %>% round
names(K) <-  paste('K', N_start$site, sep = '_')
saveRDS(K, './output/carrying_capcity.rds')

# environmental stoch -----------------------------------------------------
# sd applied to survival yearly  (logit)

# approx 10% variation in survival 
var10 <- c(unlogit(survival_logit)-0.1,unlogit(survival_logit)+0.1)
survival_logit -logit(var10)
logit(var10)

env_lower <- 0
env_upper <- 0.24 # logit scale

m <- unlogit(qnorm(c(0.025, 0.975), mean = survival_logit, env_upper))
var10;m;var10-m
#===================================__
# approx 20% variation in survival 
var10 <- c(unlogit(survival_logit)-0.2,unlogit(survival_logit)+0.2)
survival_logit -logit(var10)
logit(var10)

env_lower <- 0
env_upper <- 0.43 # logit scale

m <- unlogit(qnorm(c(0.025, 0.975), mean = survival_logit, env_upper))
var10;m;var10-m

# sample

sample_data <- randomLHS(reps, nrow(N_start))
colnames(sample_data) <-paste0('env_stoch_', N_start$site)
env_stochasticity <- (env_lower + sample_data * (env_upper - env_lower))


# variable params ---------------------------------------------------------


param_est <- cbind(initial_N, 
                   juv0,
                   site_survival,
                   reproduction, 
                   env_stochasticity)

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
                          survivalSA = unlogit(survival_logit),
                          survivalSA_sd = unlogit(survival_logit_sd),
                          survivalA = unlogit(survival_logit),
                          survivalA_sd = unlogit(survival_logit_sd),
                          transition,
                          reproduction_lower = repro_lower,
                          reproduction_upper = repro_upper,
                          clutches = mean(clutch_sizes),
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
         type = ifelse(Parameter %in% c('survivalSA','survivalA'),
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
  arrange(parameter)

paramter_fxtb <- parameter_table %>% 
  mutate(distribution = c(rep('truncated normal ', 4), 
                          'uniform', rep('logit normal',2),
                          '', 'uniform', '~beta(2,2)', 'uniform',
                          rep('',4))) %>% 
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


# variable distributions ----------------------------
param.labs <- c('initial N', 'environmental stoch', 'survival', 'fecundity')
names(param.labs) <- c('N', 'e','s', 'F')



distribution_fig <- param_est %>% as.data.frame %>% 
  rename(survival_J= survival_juv) %>% 
  pivot_longer(cols = N_init_CA:env_stoch_MA) %>% 
  mutate(type = factor(substr(name, 1,1),levels = c('N', 'e','s', 'F')),
         variable = sub('N_init_', '', name),
         variable = sub('survival_', '', variable),
         variable = sub('env_stoch_', '', variable),
         variable = sub('F_reproduction', '', variable)) %>% 
  ggplot(aes(x = value, fill = name))+
  geom_histogram(colour = 'black', bins = 30)+
  facet_wrap(type~variable, scale = 'free', 
             labeller = labeller(type = param.labs))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold'));distribution_fig


# stage matrix ------------------------------------------------------------
eigen_analysis <- function(A) {
  ev <- eigen(A)
  lmax <- which(Re(ev$values) == max(Re(ev$values)))
  lambda <- Re(ev$values[lmax])
  W <- ev$vectors
  w <- abs(Re(W[, lmax]))
  stable_age <- w/sum(w)
  return(list(lambda, stable_age))
}

stage.mat <- matrix(c(0, 0, 2.15,2.15,2.15,0,
                      0.072, 0, 0, 0,0,0,
                      0.288, 0,0,0,0,0,
                      0, 0.36,0.36, 0,0,0,
                      0, 0, 0,0.36,0,0,
                      0, 0, 0,0,0.36,0), nrow = 6, ncol = 6, byrow = TRUE,
                    dimnames = list(c('J', 'SA','A1','A2','A3', 'A4'),
                                    c('J', 'SA','A1','A2','A3', 'A4')))

stage_distribution <- eigen_analysis(stage.mat)[[2]]
names(stage_distribution)<- c('J', 'SA','A1','A2','A3', 'A4')

saveRDS(stage_distribution, './output/stage_distribution.rds')

fx_stages <- as.data.frame(stage.mat) %>% 
  mutate(stages = rownames(.)) %>% 
  rbind(c(round(stage_distribution,2),'stable:')) %>% 
  relocate(stages) %>% 
flextable() %>% 
  autofit() %>% 
  vline(j = 1, i = 1:6, part = 'body', border = fp_border_default(width = 1.5)) %>% 
  vline(j = 1, part = 'header', border = fp_border_default(width = 1.5)) %>% 
  hline(i = 6, border = fp_border_default(width = 1.5)) %>% 
  align(j = 1,  align = 'right') %>% 
  align(part = 'header', j  = 1,align = 'right') %>% 
  bold(part = 'header') %>% 
  bold(j= 1) %>% 
  font(part = 'all', fontname='calibri') %>% 
  set_caption(caption = paste('CGED Stage Matrix ( lamba =',
                              round(eigen_analysis(stage.mat)[[1]],3), ')'));fx_stages



transition_mat <- matrix(c(0, 0, 0,0,0,0,
                      0.194, 0, 0, 0,0,0,
                      0.806, 0,0,0,0,0,
                      0, 1,1, 0,0,0,
                      0, 0, 0,1,0,0,
                      0, 0, 0,0,1,0), nrow = 6, ncol = 6, byrow = TRUE,
                    dimnames = list(c('J', 'SA','A1','A2','A3', 'A4'),
                                    c('J', 'SA','A1','A2','A3', 'A4')))

saveRDS(transition_mat, './output/transition_matrix.rds')

fx_transition<- transition_mat %>% 
  as.data.frame() %>% 
  mutate(stages = rownames(.)) %>% 
  relocate(stages) %>% 
  flextable() %>% 
  autofit() %>% 
  vline(j = 1, i = 1:6, part = 'body', border = fp_border_default(width = 1.5)) %>% 
  vline(j = 1, part = 'header', border = fp_border_default(width = 1.5)) %>% 
  align(j = 1,  align = 'right') %>% 
  align(part = 'header', j  = 1,align = 'right') %>% 
  bold(part = 'header') %>% 
  bold(j= 1) %>% 
  font(part = 'all', fontname='calibri')  %>% 
  set_caption(caption = 'CGED Transition Matrix');fx_transition

