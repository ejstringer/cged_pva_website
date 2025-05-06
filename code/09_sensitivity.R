sensitivity <- function(A) {
  d <- eigen(A)$values   # eigen values
  w <- eigen(A)$vectors  # right eigen vectors
  v <- Conj(solve(w))    # complex conjugate of left eigen vectors
  # output of eigenvalues is decreasingly sorted.
  v[1,] %*% t(w[,1])
}
eigen_analysis <- function(A) {
  ev <- eigen(A)
  lmax <- which(Re(ev$values) == max(Re(ev$values)))
  lambda <- Re(ev$values[lmax])
  W <- ev$vectors
  w <- abs(Re(W[, lmax]))
  stable_age <- w/sum(w)
  return(list(lambda, stable_age))
}
source('./code/00_libraries.R')
source('./code/05_model.R')
param_starting_estimates <- readRDS('./output/starting_values.rds')
stage_distribution <- readRDS('./output/stage_distribution.rds')
transition_mat <- readRDS('./output/transition_matrix.rds')
K <- readRDS('./output/carrying_capcity.rds')
clutch_sizes <- 4:7
param_selected <- readRDS('./output/selected_25models_parameters.RDS')
param_selecteddf <-readRDS('./output/selected_25models_parameters_df.RDS')
# stage distribution ------------------------------------------------------



# fecundity <- (mean(param_selected[,'F_reproduction'])*0.5)*5.5 
# ASAsur <- mean(param_selected[,grep('survival_A', colnames(param_selected))])
# Jt <- mean(param_selected[,'survival_juv'])

fecundity <- (mean(param_selecteddf$mean[param_selecteddf$name == 'F_reproduction'])*0.5)*5.5 
ASAsur <- mean(param_selecteddf$mean[grep('survival_A', param_selecteddf$name)])
Jt <- mean(param_selecteddf$mean[grep('survival_J', param_selecteddf$name)])

stage.mat <- matrix(c(0, 0, rep(fecundity,3),0,
                      Jt*0.194, 0, 0, 0,0,0,
                      Jt*0.806, 0,0,0,0,0,
                      0, ASAsur,ASAsur, 0,0,0,
                      0, 0, 0,ASAsur,0,0,
                      0, 0, 0,0,ASAsur,0.01), nrow = 6, ncol = 6, byrow = TRUE,
                    dimnames = list(c('J', 'SA','A1','A2','A3', 'A4'),
                                    c('J', 'SA','A1','A2','A3', 'A4')))
mat <- stage.mat


vals <- c(1.01, 1.05, 1.1)
perc_labs <- paste0('+', (vals-1)*100, '%')


results <- list()
for (x in 1:length(vals)) {
  testlam <- matrix(0, nrow=nrow(mat), ncol=nrow(mat))
  
  for (i in c(1:nrow(mat))) {
    for (j in c(1:ncol(mat))) {
      if (mat[i,j] == 0) {testlam[i,j] <- 0} else {
        tempmat <- mat
        new.val <- mat[i,j]*vals[x]
        if (i !=1 & new.val > 1) new.val <- 1 
        tempmat[i,j] <- mat[i,j]*vals[x]
        testlam[i,j] <- Re(eigen(tempmat)$values[1])
      }
    }
  }
  results[[x]] <- testlam
}

output <- matrix(NA, nrow=length(vals), ncol=length(which(as.vector(t(results[[1]])) >0)))
for (w in c(1: length(vals))) {
  output[w,] <- as.vector(t(results[[w]]))[which(as.vector(t(results[[w]])) >0)]
}

colnames(output) <- c('fecA1', 'fecA2', 'fecA3', 'JsurvSA', 'JsurvA', 
                      'SAsurvA1', 'A1surv','A2surv','A3surv','A4surv')

lower <- 0.8
lambda <- Re(eigen(mat)$values[1])

lambda

output/lambda

par(mar=c(5, 4, 2, 2))
barplot(output-lower, beside=TRUE, ylim=c(lower, 1.1), 
        ylab="Population growth rate",
        col=c("#007765","#6AA341","gray90"), offset = lower)
abline(h=lower, lwd=2)
abline(h=1.0, lwd=2, col="blue")
abline(h= (Re(eigen(mat)$values[1])), lwd=2, col="orange")
legend("topleft", inset=c(0.01,0.01), perc_labs,
       fill=c("#007765","#6AA341","gray90"), cex=0.8)
# text(x=c(3.7,7.7,11.65,15.8, 20), y=0.78,
#      labels=c("AdultFec","JsurvSA", "JsurvA", "SubSurv","AdultSurv"),
#      pos=2, xpd=TRUE, cex=0.9)

output %>% as.data.frame() %>% 
  mutate(diff = perc_labs) %>% 
  pivot_longer(-diff, names_to = 'stage', values_to = 'r') %>% 
  mutate(lambda = lambda,
         increase = r - lambda) %>% 
  pivot_longer(cols = c(r, increase)) %>%
  mutate(lambda = ifelse(name == 'r', lambda, 0),
         line = ifelse(name == 'r', 1, 0.05)) %>% 
  filter(name != 'r') %>% 
  ggplot(aes(stage, value, fill = diff))+
  geom_bar(stat = 'identity', position = position_dodge())+
  #facet_grid(~(name), scale = 'free')+
  scale_fill_manual(values = c("#007765","#6AA341","gray90"),
                    name = NULL)+
  theme_classic()+
  geom_hline(aes(yintercept = lambda), colour = 'orange')

emat <- popbio::elasticity(mat)
sumfec <- sum(emat[1,3:6])
sumgrow <- sum(emat[2:3,1])
sumgrow2 <- sum(emat[4,2])
sumstasis <- sum(emat[4:6,4:6])

par(mar=c(2, 2, 2, 2))
pie(c(sumfec,sumgrow,sumgrow2, sumstasis), col = c("gray50", "#6AA341",'#6AA990', "gray90"),
    labels=c("fecundity","growth+survival J","growth+survival SA", "growth+survival A"))




# junk --------------------------------------------------------------------


mat <- matrix(c(0,0,fecundity,
                Jt*0.194, 0,0,
                Jt*0.906, ASAsur,ASAsur), nrow = 3, ncol =3, byrow = T)

eigen_analysis(mat)
eigen(mat)
eigen_analysis(stage.mat)
eigen(stage.mat)

sensitivity(mat)

elasticity <- function() S1.spring() * P() / lambda(A.spring())



# code --------------------------------------------------------------------
#https://tomizonor.wordpress.com/2014/08/31/periodic-matrix-model/
#https://www.r-bloggers.com/2014/10/sensitivity-and-elasticity-of-seasonal-matrix-model/

#https://compadre-db.org/Education/article/sensitivity-and-elasticity-matrices
giraffe <- matrix(c(0, 0, 0.24, 0.57, 0, 0, 0, 0.79, 0.84), nrow=3, byrow=TRUE,
                  dimnames=list(c("calf","subadult","adult"),
                                c("calf","subadult","adult")))

mat <- giraffe

w <- eigen(mat)$vectors
v <- Conj(solve(w))

senmat <- Re(v[1,] %*% t(w[,1]))
emat <- (1/(Re(eigen(mat)$values[1]))) * senmat * mat

popbio::sensitivity(mat)
popbio::elasticity(mat)

sumfec <- sum(emat[1,2:3])
sumgrow <- sum(emat[2,1],emat[3,1],emat[3,2])
sumstasis <- sum(diag(emat),emat[2,3])



# big mat
sumfec <- sum(emat[1,2:6])
sumgrow <- sum(emat[2:3,1],emat[3,1],emat[3,2])
sumstasis <- sum(diag(emat),emat[2,3])

emat

par(mar=c(2, 2, 2, 2))
pie(c(sumgrow,sumfec,sumstasis), col = c("gray50", "#6AA341", "gray90"),
    labels=c("growth", "fecundity", "stasis"))

vals <- c(1.01, 1.1, 1.2)
perc_labs <- paste0('+', (vals-1)*100, '%')


results <- list()
for (x in 1:length(vals)) {
  testlam <- matrix(0, nrow=nrow(mat), ncol=nrow(mat))
  
  for (i in c(1:nrow(mat))) {
    for (j in c(1:ncol(mat))) {
      if (mat[i,j] == 0) {testlam[i,j] <- 0} else {
        tempmat <- mat
        tempmat[i,j] <- mat[i,j]*vals[x]
        testlam[i,j] <- Re(eigen(tempmat)$values[1])
      }
    }
  }
  results[[x]] <- testlam
}

output <- matrix(NA, nrow=length(vals), ncol=length(which(as.vector(t(results[[1]])) >0)))
for (w in c(1: length(vals))) {
  output[w,] <- as.vector(t(results[[w]]))[which(as.vector(t(results[[w]])) >0)]
}

colnames(output) <- c('fec', 'JsurvSA', 'JsurvA', 'SAsurv', 'Asurv')

colnames(output) <- c('fecA1', 'fecA2', 'fecA3', 'JsurvSA', 'JsurvA', 'SAsurvA1', 'A1surv','A2surv','A3surv','A4surv')

colnames(output) <- c('fec', 'calfsurv', 'subsurv', 'Asurv')
lower <- 0.8
Re(eigen(mat)$values[1])
par(mar=c(5, 4, 2, 2))
barplot(output-lower, beside=TRUE, ylim=c(lower, 1.12), ylab="Population growth rate",
        col=c("#007765","#6AA341","gray90"), offset= lower)
abline(h=lower, lwd=2)
abline(h=1.0, lwd=2, col="blue")
abline(h= (Re(eigen(mat)$values[1])), lwd=2, col="orange")
legend("topleft", inset=c(0.01,0.01), perc_labs,
       fill=c("#007765","#6AA341","gray90"), cex=0.8)
# text(x=c(3.7,7.7,11.65,15.8, 20), y=0.78,
#      labels=c("AdultFec","JsurvSA", "JsurvA", "SubSurv","AdultSurv"),
#      pos=2, xpd=TRUE, cex=0.9)



emat <- popbio::elasticity(mat)
sumfec <- sum(emat[1,3:6])
sumgrow <- sum(emat[2:3,1])
sumgrow2 <- sum(emat[4,2])
sumstasis <- sum(emat[4:6,4:6])

par(mar=c(2, 2, 2, 2))
pie(c(sumfec,sumgrow,sumgrow2, sumstasis), col = c("gray50", "#6AA341",'#6AA990', "gray90"),
    labels=c("fecundity","growth+survival J","growth+survival SA", "growth+survival A"))




# old ---------------------------------------------------------------------


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

eigen_analysis <- function(A) {
  ev <- eigen(A)
  lmax <- which(Re(ev$values) == max(Re(ev$values)))
  lambda <- Re(ev$values[lmax])
  W <- ev$vectors
  w <- abs(Re(W[, lmax]))
  stable_age <- w/sum(w)
  return(list(lambda, stable_age))
}



# load --------------------------------------------------------------------
figprefix <- 'ntg100m'
source('./code/00_libraries.R')
source('./code/05_model.R')
param_starting_estimates <- readRDS('./output/starting_values.rds')
stage_distribution <- readRDS('./output/stage_distribution.rds')
transition_mat <- readRDS('./output/transition_matrix.rds')
K <- readRDS('./output/carrying_capcity.rds')
clutch_sizes <- 4:7
param_selected <- readRDS('./output/selected_25models_parameters.RDS')
param_selecteddf <-readRDS('./output/selected_25models_parameters_df.RDS')
# stage distribution ------------------------------------------------------



fecundity <- (mean(param_selecteddf$mean[param_selecteddf$name == 'F_reproduction'])*0.5)*5.5 
ASAsur <- mean(param_selecteddf$mean[grep('survival_A', param_selecteddf$name)])
Jt <- mean(param_selecteddf$mean[grep('survival_J', param_selecteddf$name)])

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

# model inputs ------------------------------------------------------------
top25 <- param_selecteddf
paramlist <- list(populations = c('CA', 'JE', 'JW', 'MA'),
                  initial_ab = round(top25$mean[2:5]), # adult abundance for populations
                  survival = top25$mean[10:13], # survival of adults and SA at sites
                  survival_J = top25$mean[14:17], # juvenile survival
                  env_stoch = top25$mean[6:9], # sd on survival
                  f_reproducing = top25$mean[1], # proportion of females reproducing
                  clutch_sizes = clutch_sizes, # clutch size range   
                  K = K # carrying capcity applied to adults and SA
)



paramlist

# sensitivity -------------------------------------------------------------

original_paramlist <- paramlist

paramlist 

names(paramlist)

sim_summary <- NULL
seqencedf <- data.frame(parameter = rep(names(paramlist)[2:6], each = 2), 
                        ten = c(0.9, 1.1))
base_summary <- readRDS('./output/base_summary_r_ext_p.rds')
# base model --------------------------------------------------------------
reps_base <- 1000
for (i in 1:nrow(seqencedf)) {
  
#top25 <- apply(param_selected, 2, summary)[4,]
param_adjust <- seqencedf$parameter[i]
adjust <- seqencedf$ten[i]
paramlist10 <- paramlist
paramlist10[[param_adjust]] <- paramlist10[[param_adjust]]*adjust  


N_sim <- em.pva_simulator(populations = c('CA', 'JE', 'JW', 'MA'),
                          stages = c('J', 'SA','A1','A2','A3', 'A4'),
                          stage_distribution = stage_distribution, 
                          initial_ab = round(paramlist10$initial_ab), # adult abundance for populations
                          survival = paramlist10$survival, # survival of adults and SA at sites
                          survival_J = paramlist10$survival_J, # juvenile survival
                          survival_logit_sd = NULL,
                          site_adjust = NULL,
                          env_stoch = paramlist10$env_stoch, # sd on survival
                          transition_mat = transition_mat, # transition prob to SA
                          f_reproducing = paramlist10$f_reproducing, # proportion of females reproducing
                          clutch_sizes = clutch_sizes, # clutch size range   
                          K = K, # carrying capcity applied to adults and SA
                          time_steps = 50, # time
                          replicates = reps_base)


# extract N ---------------------------------------------------------------


N_sim %>% length


system.time(N_simulated <- mclapply(1:reps_base, 
                                    em.extract_N, N_sim, mc.cores = 10) %>% 
              do.call('rbind', .) %>% 
              mutate_if(is.character, factor))
nlevels(N_simulated$stage)
nlevels(N_simulated$pop)
levels(N_simulated$stage) <- c('J', 'SA', 'A1', 'A2', 'A3', 'A4')
levels(N_simulated$pop) <- c('CA', 'JE', 'JW', 'MA')

sim_N_sum <- N_simulated %>% 
  group_by(rep, tstep, pop) %>% 
  summarise(N = sum(N)) %>% 
  mutate(extinct = N < 2)

sum(sim_N_sum$extinct)

sim_adjusted_lambda <- sim_N_sum %>% 
  ungroup() %>% 
  group_by(pop, rep) %>% 
  mutate(diff_year = tstep - lag(tstep),
         diff_growth = N - lag(N),
         rate_percent = (diff_growth /diff_year)/ lag(N) * 100) %>%
  ungroup() %>% 
  group_by(tstep,pop) %>%
  summarise(sims = n(),
            sim0 = sum(extinct),
            ext_p = sum(extinct)/n(),
            rate = mean(rate_percent, na.rm=T))%>%
  ungroup() %>%
  group_by(pop) %>% 
  summarise(growth_rate = mean(rate, na.rm = T),
            extinction = mean(ext_p)) %>% 
  left_join(base_summary)

sim_adjusted_lambda$parameter <- param_adjust
sim_adjusted_lambda$adjusted <- adjust

sim_summary <- bind_rows(sim_summary, sim_adjusted_lambda)
}
sim_summary %>% 
  dplyr::select(-extinction, -extinction_base) %>% 
  pivot_wider(names_from = adjusted, values_from = growth_rate) %>% 
  mutate(S = (`1.1`-`0.9`)/base_growth_rate,
         param2 = rep(c('N', 'Sur', 'Juv', 'Env', 'fem'), each = 4)) %>% 
  ggplot(aes(param2, S, fill = pop)) +
  geom_bar(stat = 'identity')+
  facet_wrap(~pop, scale = 'free')

sim_summary$param2 <- factor(sim_summary$parameter) 

levels(sim_summary$param2)<- c( 'Env','F','N',  'Sur', 'Juv')
sim_summary %>%  
  #mutate(S = growth_rate-base_growth_rate) %>%
  pivot_wider(names_from = adjusted, 
              values_from = growth_rate) %>% 
  mutate(S = (`1.1`-`0.9`)/(0.1*base_growth_rate)) %>% 
  ggplot(aes(param2, S, fill = param2)) +
  geom_bar(stat = 'identity')+
  facet_wrap(~pop, scale = 'free_x')+
  theme_bw()

sim_summary %>% 
  pivot_wider(names_from = adjusted, 
              values_from = growth_rate) %>% 
  group_by(parameter, pop) %>% 
  summarise(lower = mean(`0.9`),
            baseline = mean(base_growth_rate),
            upper = mean(`1.1`)) %>% 
  pivot_longer(lower:upper) %>% 
  ggplot(aes(value,parameter))+
  geom_line(linewidth = 1)+
  geom_point(aes(size = name))+
  theme_bw()+
  facet_grid(~pop, scale = 'free')+
  scale_size_manual(values = c(3,0,0))

sim_summary %>% 
  pivot_wider(names_from = adjusted, 
              values_from = extinction) %>% 
  group_by(parameter, pop) %>% 
  summarise(lower = mean(`0.9`),
            baseline = mean(extinction_base),
            upper = mean(`1.1`)) %>% 
  # mutate(S = (upper-lower)/baseline) %>% 
  # ggplot(aes(parameter, S, fill = pop))+
  # geom_bar(stat = 'identity', position = position_dodge())+
  # facet_wrap(~pop, scale = 'free')
  pivot_longer(lower:upper) %>% 
  ggplot(aes(value,parameter))+
  geom_line(linewidth = 0.7)+
  geom_point(aes(size = name))+
  theme_bw()+
  facet_grid(~pop, scale = 'free')+
  scale_size_manual(values = c(3,0,0))



sim_summary %>% 
  pivot_wider(names_from = adjusted, 
              values_from = extinction) %>% 
  group_by(param2, pop) %>% 
  summarise(lower = mean(`0.9`),
            baseline = mean(extinction_base),
            upper = mean(`1.1`)) %>% 
  mutate(S = (upper-lower)/baseline) %>%
  group_by(param2) %>% 
  summarise(S = mean(S)) %>% 
  arrange(S) %>% 
  mutate(n = row_number()) %>% 
  ggplot(aes(n, S, fill = param2))+
  geom_bar(stat = 'identity')+
  theme_bw()
 # filter(param2 != 'N', param2 != 'Env') %>% 
  ggplot(aes(param2, S, fill = pop))+
  geom_bar(stat = 'identity', position = position_dodge())+
  facet_wrap(~pop, scale = 'free')
