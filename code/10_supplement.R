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
# supplementation ---------------------------------------------------------
supplement <- T
total_suplemented <- 100
number_suplementation <- rep(total_suplemented/4,4)
seq_suplementation <- c(5,10,15)

em.sup.distribution <- function(N_supp, n.pops = 4, precision = 1){
 
  N_dist <- numeric(n.pops)
  for (i in 1:(n.pops-1)){
    N_dist[i] <- round(rnorm(1, N_supp/n.pops, precision))
    N_supp <- N_supp-N_dist[i]
  } 
  N_dist[n.pops] <- N_supp
  N_dist[sample(1:n.pops, n.pops)]
}
em.sup.distribution(100, 4,10)
par(mfrow = c(2,2))
lapply(1:4, function(x) hist(sapply(rep(100,1000), em.sup.distribution, 
                                    precision = 10)[x,],
                             breaks = 10))
reps <- 10000
50*1000
total_suplemented <- (seq(12,500, 4))
length(total_suplemented)

supplement_distribution<- t(sapply(rep(4,reps), em.sup.distribution))
#colnames(supplement_distribution) <- paste('sup', grasslands, sep = '_')
supplement_sample <- list()
for (i in 1:length(total_suplemented)) {
  supplement_distribution*total_suplemented[i]
supplement_distribution$N_supp <- total_suplemented[i]  
supplement_sample[[i]] <- supplement_distribution
}

supplement_sample_df <- do.call('rbind', supplement_sample)
head(supplement_sample_df) 
nrow(supplement_sample_df)

# supp take two -----------------------------------------------------------

em.sup_across_pops <- function(n.supp = 100, n.pops = 4){
  x <- seq(0,n.supp, length.out = 11)
  x_round <- round(x)
  
  if (sum(!(x == x_round))>0) print('rounded x') 
  if (sum(!(x == x_round))>0) x_round <- unique(x_round)
  
  xlist <- list()
  for (i in 1:4) xlist[[i]] <- x_round
  all_combs <- expand.grid(xlist)
  
  all_combs1 <- all_combs[rowSums(all_combs)==n.supp,]
  
  all_combs1
}
seq(0,100,2)
em.sup_across_pops(20,4) %>% nrow()
supplement_sample <- lapply(c(20,40,60,80,100,200,300,400, 500, 600), em.sup_across_pops)
sum(sapply(supplement_sample, nrow))



supplement_sample_df <- do.call('rbind', supplement_sample) %>% 
  data.frame(n = 1:(sum(sapply(supplement_sample, nrow))*10), .,
           row.names = NULL) %>% 
  select(-n)

nrow(supplement_sample_df)
set.seed(13456)
rownames(supplement_sample_df) <- paste0(
  stri_rand_strings(nrow(supplement_sample_df), 4, "[A-Z]"),
  stri_rand_strings(nrow(supplement_sample_df), 4, "[0-9]"))

supplement_sample_df
# model setup -------------------------------------------------------------
t_steps <- 50
supp_every <- 1
supplement <- T
seq_suplementation <- seq(14,t_steps, supp_every)



# model -------------------------------------------------------------------
reps <- nrow(supplement_sample_df)
sup_dist <- split(as.data.frame(supplement_sample_df),
                  rownames(supplement_sample_df))
sup_dist[[1]]
N_sim <- lapply(sup_dist, 
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
                          when.supp = seq_suplementation))

system.time(sim_N <- mclapply(N_sim, function(x) em.extract_N(1,x), mc.cores = 10))
for (i in 1:length(N_sim)) sim_N[[i]]$run <- names(sim_N)[i]

N_simulated <- do.call('rbind', sim_N)%>% 
  mutate_if(is.character, factor)
# system.time(N_simulated <- mclapply(1:reps, em.extract_N, N_sim,
#                                     mc.cores = 10) %>% 
#               do.call('rbind', .) %>% 
#               mutate_if(is.character, factor))
nlevels(N_simulated$stage)
nlevels(N_simulated$pop)
levels(N_simulated$stage) <- stages
levels(N_simulated$pop) <- grasslands


sim_N_sum <- N_simulated %>% 
  filter(stage != 'J', tstep == 20) %>% 
  group_by(run, tstep, pop) %>% 
  summarise(N = sum(N)) %>% 
  mutate(extinct = N < 2,
         year = tstep +2012)
sum(sim_N_sum$extinct)


# plots -------------------------------------------------------------------
pop.labs <- paste(grasslands, K, sep = ':')
names(pop.labs) <- grasslands
colnames(supplement_sample_df) <- grasslands

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
  mutate(max_supp = as.numeric(max_supp))

table(table(sim_N_sum_supp$run))
table(table(sim_N_sum_supp$supp_dist))

supplement_sample_df %>% 
  mutate(supp_dist = paste(CA,JE, JW, MA)) %>% 
  group_by(supp_dist) %>% 
  summarise_if(is.numeric, list(mean = mean,sd = sd))

table(sim_N_sum_supp$n_supp == 0)
table((sim_N_sum_supp$supp_group))

  ggplot(sim_N_sum_supp,aes(n_supp, N, colour = extinct))+
  geom_point(aes(size = extinct))+
  geom_smooth(method = 'gam',  se = F, colour = 'black')+
 # geom_hline(yintercept = K, colour = 'grey', lty = 2)+
  theme_light()+
  scale_size_manual(values = c(0.5, 2))+
  facet_grid(rows = vars(pop), scale = 'free', labeller = labeller(pop = pop.labs))+
  theme(strip.text.y.right = element_text(angle = 0))

  
extinction_prob <- sim_N_sum_supp %>% ungroup() %>%  
  group_by(tstep,pop, max_supp) %>% 
  summarise(sims = n(),
            sim0 = sum(extinct),
            ext_p = sum(extinct)/n()) %>% 
  mutate(extinct_prob = ext_p > 0)
tail(extinction_prob)

extinction_prob %>% 
  #filter(max_supp >0) %>% 
  ggplot(aes(max_supp, ext_p, colour = extinct_prob))+
  geom_point(alpha = 0.7)+
  geom_line()+
  facet_wrap(~pop)+
  scale_x_log10()+
  theme_light()+
  geom_hline(yintercept = 0.05) -> g;g

plotly::ggplotly(g)

###### abundance ---------------------------------------------------------------


sim_N_sum %>% 
  mutate(supplementation = T) %>% 
  group_by(tstep, pop, supplementation) %>% 
  summarise(N = mean(N)) %>% 
  mutate(year = tstep +2012,
         supp_year = ifelse(tstep %in% seq_suplementation, 'yes','no')) %>% 
  filter(year > 2012) %>% 
  ggplot(aes(year, N, colour = supplementation))+
  # geom_line(data = filter(sim_N_sum, rep != 1000),
  #           aes(group = rep), alpha = 0.4, colour = 'grey',
  #           linewidth = 0.25)+
  geom_point(aes(colour = supp_year))+
  geom_hline(yintercept = 0, colour = 'red', lty = 2, linewidth = 0.3)+
  geom_line(alpha = 1)+
  facet_wrap(~pop, scale = 'free')+
  theme_bw()+
  #geom_vline(xintercept = seq_suplementation+2012, colour = 'grey60',lty = 2)+
  scale_color_manual(values = c(NA, 'black', 'blue'))+
  theme(panel.grid = element_blank(),
        legend.position = 'bottom',
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold'))


## extinction --------------------------------------------------------------

extinction_prob <- sim_N_sum %>% ungroup() %>% 
  group_by(tstep,pop,supplementation) %>% 
  summarise(sims = n(),
            sim0 = sum(extinct),
            ext_p = sum(extinct)/n())
tail(extinction_prob)

extinction_prob %>% 
  mutate(year = tstep +2012,
         Sites = pop) %>%
  ggplot(aes(year, ext_p, colour = Sites))+
  geom_vline(xintercept = seq_suplementation+2012, lty = 3, colour = 'grey')+
  geom_line(aes(linetype = supplementation))+
  geom_point(aes(shape = supplementation))+
  theme_classic()+
  ylab('Extinction probability')+
  xlab('Year')+
  scale_x_continuous(breaks = c(2013, 2025,seq(2010,2060,10)))+
  scale_y_continuous(breaks = seq(0,1,0.2))

## extinction all --------------------------------------------------------------

extinction_prob <- sim_N_sum %>% ungroup() %>% 
  group_by(tstep,supplementation) %>% 
  summarise(sims = n(),
            sim0 = sum(extinct),
            ext_p = sum(extinct)/n())
tail(extinction_prob)

extinction_prob %>% 
  mutate(year = tstep +2012,
         Sites = pop) %>%
  ggplot(aes(year, ext_p))+
  geom_vline(xintercept = seq_suplementation+2012, lty = 3, colour = 'grey')+
  geom_line(aes(linetype = supplementation))+
  geom_point(aes(shape = supplementation))+
  theme_classic()+
  ylab('Extinction probability')+
  xlab('Year')+
  scale_x_continuous(breaks = c(2013, 2025,seq(2010,2060,10)))+
  scale_y_continuous(breaks = seq(0,1,0.2))

## delcine threshold ------------------------------------------------------
#https://kevintshoemaker.github.io/NRES-470/LECTURE12.html#Step_1:_conceptualize_the_life_history
final_N <- sim_N_sum %>%
  left_join(area_ntg, by = c('pop' = 'site')) %>% 
  ungroup() %>% 
  mutate(N_ha = N/area_ha) %>% 
  #filter(tstep == 20) %>% 
  split(., .$pop)

# plot probabilities of different severities of decline

gen_test <- 10:50
area_ntg$site
Init_N <- initial_N#/area_ntg$area_ha
declines <- seq(0,1,by=0.02)
declineprob <- numeric(length(declines))

declineprob_list <- list()

for (pop in 1:length(grasslands)) {
  declineprob <- matrix(data = NA, ncol = length(declines),
                        nrow = length(gen_test))
  colnames(declineprob) <- paste0('d', declines)
  rownames(declineprob) <- paste0('gen', gen_test)
for (gen in gen_test) {
for(s in 1:length(declines)){
  site <- final_N[[pop]]
  N_pop <- site$N[site$tstep == gen]
  N_under <- sum(N_pop < (Init_N[pop]-(declines[s])*Init_N[pop]))
  declineprob[gen-min(gen_test)+1,s] <- N_under/length(N_pop)
}
  declineprob_df <- declineprob %>% 
    as.data.frame %>% 
    mutate(gen = rownames(declineprob)) %>% 
    mutate_all(as.character()) %>% 
    pivot_longer(-gen) %>% 
    mutate(gen = sub('gen', '', gen),
           name = sub('d', '', name)) %>% 
    mutate_all(as.numeric) %>% 
    rename(declines = name, prob = value) %>% 
    mutate(pop = grasslands[pop])
  
  declineprob_list[[pop]] <- declineprob_df
} 
}
names(declineprob_list) <- grasslands


declineprob_list %>%
  do.call('rbind', .) %>%
  filter(gen == 20) %>% 
  ggplot(aes(declines*100, prob, colour = pop, group = paste(pop,gen)))+
  geom_line()+
  theme_classic()+
  xlab("Decline threshold (percent)")+
  ylab("Probability of falling below threshold")+
  geom_vline(xintercept = 25, colour = 'black', 
             alpha = 0.5, linetype = 'dashed')

declineprob_list %>%
  do.call('rbind', .) %>% 
  group_by(declines, gen) %>% 
  summarise(prob = mean(prob)) %>% 
  filter(declines != 1, gen %in% c(10,12,14,16,18,20, 22,24,26,28)) %>% 
  ggplot(aes(declines*100, prob, colour = gen, group = gen))+
  geom_line()+
  theme_classic()+
  xlab("Decline threshold (percent)")+
  ylab("Probability of falling below threshold")+
  geom_vline(xintercept = 25, colour = 'black', 
             alpha = 0.5, linetype = 'dashed')


# all ---------------------------------------------------------------------


sim_N_sum_supp %>% 
  group_by(run, supp_dist, tstep, year) %>% 
  summarise(N = sum(N),
            supp = sum(n_supp),
            extinct = sum(N < 2)) %>% 
  ungroup() %>% 
  group_by(supp_dist, tstep, year) %>% 
  summarise(supp = mean(supp),
            N_mean = mean(N),
            N_sd = sd(N),
            n = n(),
            ext_p = sum(extinct)/(n()*4)) %>% 
  mutate(lcl = N_mean - (N_sd/sqrt(n)),
         ucl = N_mean + (N_sd/sqrt(n))) %>% 
  ungroup() %>% 
  # group_by(supp) %>% 
  # filter(lcl == max(lcl)) %>% 
  # arrange(supp)
  ggplot(aes(supp, N_mean, group = factor(supp)))+
  geom_boxplot()+
  #scale_x_log10()+
  theme_light()+
  geom_hline(yintercept = 0.05)
