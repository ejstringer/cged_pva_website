figprefix <- ''
source('./code/00_libraries.R')
source('./code/05_model.R')
# load --------------------------------------------------------------------

base_params <- readRDS('./output/base_parameters.rds')
param_selecteddf <-readRDS('./output/selected_25models_parameters_df.RDS')
param_envdf <- readRDS('./output/selected_25models_parameters_df_environment.RDS')
area_ntg <- read.csv('./output/area_used.csv')

# parameters --------------------------------------------------------------

stages <- names(readRDS('./output/stage_distribution.rds'))
transition_mat <- readRDS('./output/transition_matrix.rds')
stage_distribution <- readRDS('./output/stage_distribution_base.rds')
grasslands <- c('CA', 'JE', 'JW', 'MA')
initial_N <-round(param_selecteddf$mean[grepl('N_init', param_selecteddf$name)])
phiASA  <-  param_selecteddf$mean[grepl('survival_ASA', param_selecteddf$name)]
phiJ  <-  param_selecteddf$mean[grepl('survival_J', param_selecteddf$name)]
fecundity <- rep(param_selecteddf$mean[grepl('F_', param_selecteddf$name)],4)
K <- base_params$K
clutch_sizes <- base_params$clutch
environmental <- param_envdf$mean


area_ntg$site <- grasslands 
# # supplementation ---------------------------------------------------------
# supplement <- T
# total_suplemented <- 100
# number_suplementation <- rep(total_suplemented/4,4)
# seq_suplementation <- c(5,10,15)
# 
# em.sup.distribution <- function(N_supp, n.pops = 4, precision = 1){
#  
#   N_dist <- numeric(n.pops)
#   for (i in 1:(n.pops-1)){
#     N_dist[i] <- round(rnorm(1, N_supp/n.pops, precision))
#     N_supp <- N_supp-N_dist[i]
#   } 
#   N_dist[n.pops] <- N_supp
#   N_dist[sample(1:n.pops, n.pops)]
# }
# em.sup.distribution(100, 4,10)
# par(mfrow = c(2,2))
# lapply(1:4, function(x) hist(sapply(rep(100,1000), em.sup.distribution, 
#                                     precision = 10)[x,],
#                              breaks = 10))
# reps <- 10000
# 50*1000
# total_suplemented <- (seq(12,500, 4))
# length(total_suplemented)
# 
# supplement_distribution<- t(sapply(rep(4,reps), em.sup.distribution))
# colnames(supplement_distribution) <- paste('sup', grasslands, sep = '_')
# supplement_sample <- list()
# for (i in 1:length(total_suplemented)) {
#   supplement_distribution*total_suplemented[i]
# supplement_distribution$N_supp <- total_suplemented[i]  
# supplement_sample[[i]] <- supplement_distribution
# }
# 
# supplement_sample_df <- do.call('rbind', supplement_sample)
# head(supplement_sample_df) 
# nrow(supplement_sample_df)

# supp take two -----------------------------------------------------------
# 
# em.sup_across_pops <- function(n.supp = 100, n.pops = 4){
#   x <- seq(0,n.supp, length.out = 6)
#   x_round <- round(x)
#   
#   if (sum(!(x == x_round))>0) print('rounded x') 
#   if (sum(!(x == x_round))>0) x_round <- unique(x_round)
#   
#   xlist <- list()
#   for (i in 1:4) xlist[[i]] <- x_round
#   all_combs <- expand.grid(xlist)
#   
#   all_combs1 <- all_combs[rowSums(all_combs)==n.supp,]
#   
#   all_combs1
# }
# seq(0,100,2)
# em.sup_across_pops(20,4) %>% nrow()
# supplement_sample <- lapply(c(20,40,60,80,100,200,300,400, 500, 600), em.sup_across_pops)
# sum(sapply(supplement_sample, nrow))
# 
# 
# 
# supplement_sample_df <- do.call('rbind', supplement_sample) %>% 
#   data.frame(n = 1:(sum(sapply(supplement_sample, nrow))*100), .,
#            row.names = NULL) %>% 
#   select(-n)
# 
# nrow(supplement_sample_df)
# set.seed(13457)
# rownames(supplement_sample_df) <- paste0(
#   stri_rand_strings(nrow(supplement_sample_df), 4, "[A-Z]"),
#   stri_rand_strings(nrow(supplement_sample_df), 4, "[0-9]"))
# 
# colnames(supplement_sample_df) <- grasslands
# supplement_sample_df

# take three ----------------
n_rep<- 30
df_supp <- data.frame(rep(1:200, n_rep), 
                      rep(1:200, n_rep),
                      rep(1:200, n_rep),
                      rep(1:200, n_rep))
df_supp
nrow(df_supp)

set.seed(13457)
rownames(df_supp) <- paste0(
  stri_rand_strings(nrow(df_supp), 4, "[A-Z]"),
  stri_rand_strings(nrow(df_supp), 4, "[0-9]"))

colnames(df_supp) <- grasslands
df_supp

# model setup -------------------------------------------------------------
t_steps <- 50
supp_every <- 1
supplement <- T
seq_suplementation <- seq(14,t_steps, supp_every)

supplement_sample_df <- df_supp

# model -------------------------------------------------------------------
reps <- nrow(supplement_sample_df)
sup_dist <- split(as.data.frame(supplement_sample_df),
                  rownames(supplement_sample_df))
sup_dist[[1]]
N_sim <- mclapply(sup_dist, 
                function(param_sup) em.pva_simulator(populations = grasslands,
                          stages = stages,
                          stage_distribution = stage_distribution, 
                          initial_ab = initial_N, 
                          survival = phiASA, 
                          survival_J = phiJ,
                          survival_logit_sd = NULL,
                          site_adjust = NULL,
                          env_stoch = environmental, 
                          transition_mat = transition_mat, 
                          f_reproducing = fecundity, 
                          clutch_sizes = clutch_sizes,   
                          K = K,
                          time_steps = t_steps, # time
                          replicates = 1,
                          supp = supplement,
                          stage.supp = 1,
                          n.supp = unname(unlist(param_sup)), 
                          when.supp = seq_suplementation), mc.cores = 5)
length(N_sim)
system.time(sim_N <- mclapply(N_sim, function(x) em.extract_N(1,x), mc.cores = 5))
for (i in 1:length(N_sim)) sim_N[[i]]$run <- names(sim_N)[i]

N_simulated <- do.call('rbind', sim_N)%>% 
  mutate_if(is.character, factor)

nlevels(N_simulated$stage)
nlevels(N_simulated$pop)
levels(N_simulated$stage) <- stages
levels(N_simulated$pop) <- grasslands


sim_N_sum <- N_simulated %>% 
 # filter(stage != 'J') %>% 
  mutate(age = ifelse(stage == 'J', 'N_Juv', 'N_adult')) %>% 
  group_by(run, tstep, pop, age) %>% 
  summarise(N = sum(N)) %>% 
  mutate(year = tstep +2012) %>% 
  pivot_wider(names_from = age, values_from = N) %>% 
  mutate(N = N_Juv + N_adult)


sim_N_sum_supp <- supplement_sample_df %>% 
  as.data.frame() %>% 
  mutate(run = rownames(.),
         supp_dist = paste(CA,JE, JW, MA)) %>% 
  pivot_longer(-c(run, supp_dist), names_to = 'pop', values_to = 'n_supp') %>% 
  right_join(sim_N_sum)  %>% 
  mutate(supp_group = cut(n_supp, #https://www.r-bloggers.com/2020/09/how-to-convert-continuous-variables-into-categorical-by-creating-bins/
                          breaks = unique(quantile(n_supp, 
                                                   probs=seq.int(0,1, by=1/40)), 
                                          include.lowest=TRUE)),
         supp_group = ifelse(is.na(supp_group), '(-1,0]',
                             as.character(supp_group)),
         supp_group = factor(supp_group),
         supp = sub('*\\(', '', supp_group),
         supp = sub(']', '', supp)) %>% 
  separate(supp, into = c('min_supp', 'max_supp'), sep = ',') %>% 
  mutate(max_supp = as.numeric(max_supp), 
         extinct = N < 2)

head(sim_N_sum_supp)


# N ~ tsteps  -------------------------------------------------------------
summary(1:100)
quantile(1:100, 0.25)
pg1 <- sim_N_sum_supp %>% 
  filter(n_supp %in% c(1,10,100),
         tstep > 8) %>% 
  group_by(tstep, year, pop, n_supp) %>% 
  summarise(N_mean = mean(N),
            N_sd = sd(N),
            q95 = quantile(N, 0.95),
            q05 = quantile(N, 0.05)) %>% 
  left_join(data.frame(pop = grasslands, K = K)) %>%
  mutate(K = ifelse(K < (q95), K, NA)) %>% 
  ggplot(aes(year, N_mean, colour = pop, fill = pop))+
  facet_grid(n_supp~pop, scale = 'free')+
  geom_ribbon(aes(ymin = q05, ymax=q95), alpha = 0.2, colour = NA)+
  geom_point()+
  geom_line()+
    theme_bw();pg1
  #geom_hline(aes(yintercept = K), colour = 'grey', lwd = 0.75, lty = 1)+
  #geom_hline(aes(yintercept = K), colour = 'black', lwd = 1, lty = 3)
  
ggsave('./figures/supplementation_440400_supp.png',
       plot = pg1, units = 'cm',
       height = 16, width = 20)

# N ~ n.supp --------------------------------------------------------------
sim_N_sum_supp_2050 <- sim_N_sum_supp %>% 
  filter(year == 2050) %>% 
  mutate(extinct = N < 2) %>% 
  group_by(n_supp, pop) %>% 
  summarise(N_mean = mean(N),
            N_sd = sd(N),
            q95 = quantile(N, 0.95),
            q05 = quantile(N, 0.05),
            N_adult_mean = mean(N_adult),
            N_adult_sd = sd(N_adult),
            q95_adult = quantile(N_adult, 0.95),
            q05_adult = quantile(N_adult, 0.05),
            ext_prob = sum(extinct)/n())

supp_seq <- sort(unique(sim_N_sum_supp$n_supp)) 
supp_seq[!supp_seq %in% c(64, 180, 320, 360, 480)]

pg2 <- sim_N_sum_supp_2050 %>% 
  filter(n_supp != 1) %>% 
  mutate(x_N = n_supp) %>% 
  ggplot(aes(n_supp, N_mean, colour = pop, fill = pop))+
  geom_line()+
  geom_point()+
  geom_line(aes(y = x_N), colour = 'red')+
  facet_wrap(~pop, scale = 'free')+
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.3,
              lwd = 0.2)+
#  scale_x_log10(breaks = supp_seq[!supp_seq %in% c(36, 100, 200, 500,64, 180, 320, 360, 480)])+
 # scale_x_continuous(n.breaks = 20)+
 # scale_x_sqrt()+
  theme_bw()+
  #geom_hline(yintercept = c(50), lty =2, colour = 'gold')+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = 'none')+
  geom_smooth(method = 'gam', colour = 'black'); pg2
plotly::ggplotly(pg2)
ggsave('./figures/supplementation_2050_supp.png',
       units = 'cm',
       height = 16, width = 20)

# Ext ~ n.supp
sim_N_sum_supp_2050 %>% 
  filter(n_supp > 0) %>% 
  mutate(prob_0 = ext_prob == 0) %>% 
  ggplot(aes(n_supp, ext_prob, colour = prob_0))+
  geom_line()+
  geom_point()+
  facet_grid(rows = vars(pop), scale = 'free_x')+
  theme_bw()+
  geom_hline(yintercept = 0.05, lty = 2)+
 # scale_x_log10(breaks = supp_seq[!supp_seq %in% c(64, 160,36, 40, 100, 240, 500,180, 320, 360, 480)])+
  scale_x_continuous(n.breaks = 20)+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  geom_vline(xintercept = 0)+
  geom_smooth(method = 'gam', colour = 'black')

ggsave('./figures/supplementation_2050_supp_extprob.png',
       units = 'cm',
       height = 16, width = 20)
