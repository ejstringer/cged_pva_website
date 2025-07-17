figprefix <- ''
source('./code/00_libraries.R')
library(mgcv)
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
fecundity <- rep(base_params$fecundity,4) #rep(param_selecteddf$mean[grepl('F_', param_selecteddf$name)],4)
K <- base_params$K
clutch_sizes <- base_params$clutch
environmental <- param_envdf$mean


area_ntg$site <- grasslands 
# # supplementation ---------------------------------------------------------
## take three ----------------
n_rep<- 100
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
n_cores <- 10
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
                          when.supp = seq_suplementation), mc.cores = n_cores)
length(N_sim)
system.time(sim_N <- mclapply(N_sim, function(x) em.extract_N(1,x), 
                              mc.cores = n_cores))
for (i in 1:length(N_sim)) sim_N[[i]]$run <- names(sim_N)[i]

N_simulated <- do.call('rbind', sim_N)%>% 
  mutate_if(is.character, factor)

nlevels(N_simulated$stage)
nlevels(N_simulated$pop)
levels(N_simulated$stage) <- stages
levels(N_simulated$pop) <- grasslands

## summarise ---------------
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


# plot --------------------------------------------------------------------


## N ~ tsteps  -------------------------------------------------------------
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

## N ~ n.supp --------------------------------------------------------------
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
#plotly::ggplotly(pg2)
ggsave('./figures/supplementation_2050_supp.png',
       units = 'cm',
       height = 16, width = 20)

# Ext ~ n.supp
pg3 <- sim_N_sum_supp_2050 %>% 
  filter(n_supp > 0, n_supp <101) %>% 
  mutate(prob_0 = ext_prob == 0) %>% 
  ggplot(aes(n_supp, ext_prob, colour = prob_0))+
  geom_line()+
  geom_point()+
  facet_grid(rows = vars(pop), scale = 'free_x')+
  theme_bw()+
  geom_hline(yintercept = 0.05, lty = 2)+
 # scale_x_log10(breaks = supp_seq[!supp_seq %in% c(64, 160,36, 40, 100, 240, 500,180, 320, 360, 480)])+
  scale_x_continuous(n.breaks = 10)+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  geom_vline(xintercept = 0)+
  geom_smooth(method = 'gam', colour = 'black')

ggsave('./figures/supplementation_2050_supp_extprob.png', plot = pg3,
       units = 'cm',
       height = 16, width = 20)


## extinction --------------------------------------------------------------
supp_run_code <- supplement_sample_df %>% 
  mutate(run = rownames(.),
         supp = CA) %>% 
  select(run, supp)

extinction_prob <- sim_N_sum %>% 
  mutate(extinct = N < 2,
         year = tstep +2012) %>% 
  left_join(supp_run_code) %>%
  ungroup() %>% 
  group_by(year, tstep,pop, supp) %>% 
  summarise(sims = n(),
            sim0 = sum(extinct),
            ext_p = sum(extinct)/n())
extinction_prob

pg4<-extinction_prob %>% 
  filter(year %in% seq(2025,2100, 10), supp > 10) %>% 
  ggplot(aes(supp, ext_p, colour = pop))+
  geom_point(alpha = 0.4)+
  geom_smooth(method = 'gam')+
  theme_bw()+
  facet_wrap(~year, scale = 'free')+
  ylab('Extinction probability')+
  xlab('Yearly Supplementation (>10)')+
  geom_hline(yintercept = 0.05, linetype = 'dashed');pg4

ggsave('./figures/supplementation_extinction_prob.png', plot = pg4,
       units = 'cm',
       height = 16, width = 20)



### gam model ---------------------------------------------------------------

pg5<-extinction_prob %>% 
  filter(year %in% seq(2035,2100, 10),
         ext_p >0, supp < 101) %>% 
  ggplot(aes(supp, ext_p, colour = pop))+
  geom_point(alpha = 0.4)+
  geom_smooth(method = 'gam', se = F, colour = 'grey30')+
  geom_smooth(method = 'lm', linetype = 'dashed', se = F, colour = 'red')+
  theme_bw()+
  scale_y_log10()+
  scale_x_log10()+
  facet_grid(pop~year, scale = 'free')+
  ylab('Extinction probability')+
  xlab('Yearly Supplementation')+
  geom_hline(yintercept = 0.05, linetype = 'dashed');pg5

plotly::ggplotly(pg4)


extinction_prob %>% 
  ungroup() %>% 
  filter(ext_p <= 0.05,
         year == 2050) %>% 
  arrange(pop,supp) %>% 
  group_by(year, pop, ext_p) %>% 
  summarise(max(supp))


extinction_prob %>% 
  filter(year %in% seq(2026,2100, 4),
         ext_p >0) %>% 
  ggplot(aes(supp, ext_p, colour = pop))+
  geom_point(alpha = 0.4)+
  geom_smooth(method = 'lm', se = F)+
  #geom_smooth(method = 'lm', aes(group = year))+
  theme_bw()+
  scale_y_log10()+
  scale_x_log10()+
  facet_wrap(~year, scale = 'free_x')+
  ylab('Extinction probability')+
  xlab('Yearly Supplementation')+
  geom_hline(yintercept = 0.05, linetype = 'dashed')

# years > 33 seem to converge

ext_data <- extinction_prob %>% 
  filter(year > 2035,
       #  ext_p >0, 
         supp < 100) %>% 
  mutate(logsupp = log(supp),
         logextp = log(ext_p),
         logitextp = logit(ext_p))

pg6 <-  ggplot(ext_data,aes(supp, ext_p))+
  geom_point(alpha = 0.1, colour = 'grey60')+
  # geom_smooth(method = 'lm', se = F, colour = 'grey30')+
  geom_smooth(method = 'gam', linetype = 'dashed',
              aes(group = year, colour = year), se = F, linewidth = 0.2)+
  theme_bw()+
  #scale_y_log10()+
  scale_x_log10()+
  facet_wrap(~pop, scale = 'free_x')+
  ylab('Extinction probability')+
  xlab('Yearly Supplementation')+
  geom_hline(yintercept = 0.05, linetype = 'dashed'); pg6


ggsave('./figures/supplementation_extinction_prob_year.png', plot = pg6,
       units = 'cm',
       height = 16, width = 20)


  lm(ext_p~logsupp*pop, data = ext_data[ext_data$year %in% c(2035:2040),]) %>% 
    summary
  
  
  gam_supp <- gam(ext_p ~ s(logsupp,
                            by = interaction(pop, year), k = 10) + pop + year, 
                  data = ext_data[ext_data$year %in% c(2035:2040),]) 
  
  summary(gam_supp)
plot(gam_supp$residuals)
  
gam_supp <- gam(ext_p ~ s(logsupp,
                          by = interaction(pop), k = 10) + pop, 
                data = ext_data) 

summary(gam_supp)
plot(gam_supp$residuals)

  ext_data$predicted <- predict(gam_supp)  

 pg7 <-  ggplot(ext_data,aes(supp, ext_p, colour = pop))+
    geom_point(alpha = 0.4)+
    #geom_smooth(method = 'lm', se = F, colour = 'grey30')+
    #geom_smooth(method = 'gam', colour = 'red', linetype = 'dashed')+
    theme_bw()+
    geom_line(aes(y = predicted), colour = 'black')+
    #scale_y_log10()+
    scale_x_log10()+
    facet_wrap(~pop, scale = 'free_x')+
    ylab('Extinction probability')+
    xlab('Yearly Supplementation')+
    geom_hline(yintercept = 0.05, linetype = 'dashed');pg7

  
 ggsave('./figures/supplementation_extinction_prob_gam.png', 
        plot = pg7,
        units = 'cm',
        height = 16, width = 20)
 
 xx <- seq(1, 60, 0.01)
  pred_extdata <- data.frame(logsupp = log(xx), supp = xx,
                                 pop = rep(levels(ext_data$pop), 
                                           each = length(xx)))
  pred_extdata$ext_p <- predict(gam_supp, newdata = pred_extdata)   
### flextable ---------
extinction_table <-  pred_extdata %>% 
    mutate(ext_p2 = round(ext_p, 2)) %>% 
    filter(ext_p2 <= 0.1) %>% 
    group_by(pop, ext_p2) %>% 
    summarise(sup = round(mean(supp))) %>% 
    pivot_wider(names_from = pop, values_from = sup) %>% 
  rename(`Extinction Probability` = ext_p2)

extinction_table %>% 
  mutate_all(~replace(., is.na(.), '-')) %>% 
  flextable() %>%
  border_remove() %>% 
  align(j = 1, align = 'center') %>% 
  hline_bottom(border = fp_border_default(width = 2)) %>%  
  hline(part = 'header', border = fp_border_default(width = 2)) %>% 
  vline(j = 1, part = 'body', border = fp_border_default(width = 1.5)) %>% 
  bold(part = 'header') %>% 
  font(part = 'all', fontname='calibri') %>% 
  bg(j = 1,bg = 'grey90') %>% 
  bg(part = 'header',bg = 'grey90') -> ext_fltb;ext_fltb

saveRDS(ext_fltb, './output/supplementation_extinct_prob_tbl.rds') 
  