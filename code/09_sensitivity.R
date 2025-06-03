
#https://tomizonor.wordpress.com/2014/08/31/periodic-matrix-model/
#https://www.r-bloggers.com/2014/10/sensitivity-and-elasticity-of-seasonal-matrix-model/

#https://compadre-db.org/Education/article/sensitivity-and-elasticity-matrices

# load --------------------------------------------------------------------


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
stages <- names(readRDS('./output/stage_distribution.rds'))
transition_mat <- readRDS('./output/transition_matrix.rds')
# K <- readRDS('./output/carrying_capcity.rds')
# clutch_sizes <- 4:7
base_params <- readRDS('./output/base_parameters.rds')
max_age <- base_params$age
K <- base_params$K
clutch_sizes <- base_params$clutch
transition <- base_params$transition

param_selected <- readRDS('./output/selected_25models_parameters.RDS')
param_selecteddf <-readRDS('./output/selected_25models_parameters_df.RDS')
# stage distribution ------------------------------------------------------



# fecundity <- (mean(param_selected[,'F_reproduction'])*0.5)*5.5 
# ASAsur <- mean(param_selected[,grep('survival_A', colnames(param_selected))])
# Jt <- mean(param_selected[,'survival_juv'])

fecundity <- (mean(param_selecteddf$mean[param_selecteddf$name == 'F_reproduction'])*0.5)*mean(clutch_sizes)
ASAsur_site <- (param_selecteddf$mean[grep('survival_A', param_selecteddf$name)])
Jt_site <- (param_selecteddf$mean[grep('survival_J', param_selecteddf$name)])

stage.mat_list <- list()
for (pop in 1:4) {
  ASAsur <- ASAsur_site[pop]
  Jt <- Jt_site[pop]
  
  stage.mat <- matrix(0, 
                      nrow = length(stages),
                      ncol = length(stages),
                      byrow = TRUE,
                      dimnames = list(stages,
                                      stages))
  
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
  
  stage.mat[length(stages), length(stages)] <- 0.001
  stage.mat_list[[pop]] <- stage.mat
  
}
output_list <- list()
for (pop in 1:4) {
  
mat <- stage.mat_list[[pop]]


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

output

colnames(output) <- c(paste0('fec', stages[which(testlam[1,] >0)]),
                      paste0('Jsur', stages[which(testlam[,1] >0)]),
                      paste0('SAsur', stages[which(testlam[,2] >0)]),
                      paste0(stages[!stages %in% c('J', 'SA')], 'surv'))



lambda <- Re(eigen(mat)$values[1])
lower <- 0.70
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

output_df <- output %>% as.data.frame() %>% 
  mutate(perc_diff = factor(perc_labs, levels = perc_labs)) %>% 
  pivot_longer(-perc_diff, names_to = 'stage', values_to = 'r') %>% 
  mutate(lambda = Re(eigen(mat)$values[1]),
         increase = r - lambda,
         stage = factor(stage, levels = colnames(output))) %>% 
  pivot_longer(cols = c(r, increase)) %>%
  mutate(#lambda = ifelse(name == 'r', lambda, 0),
         line = ifelse(name == 'r', 1, 0.05)) 

output_df$pop <- c('CA', 'JE', 'JW', 'MA')[pop]

output_list[[pop]] <- output_df
}


output_r <- output_list %>%
  do.call('rbind', .) %>% 
  filter(name == 'r') %>% 
  mutate(lower = 0.7,#lambda - 0.1,
         upper = lambda + 0.03,
         hwid = 0.4,
         site = factor(pop, levels = c('CA', 'MA', 'JE', 'JW')))
output_r %>% 
  #filter(perc_diff == '+10%') %>% 
  #pivot_longer(cols = c(lambda, value), names_to = 'change', values_to = 'r') %>% 
  # mutate(perc_diff = ifelse(change == 'lambda', '0%',
  #                           as.character(perc_diff))) %>% 
ggplot(aes(stage, value, fill = perc_diff))+
  # geom_bar(stat = 'identity', position = position_dodge(),
  #          colour = 'black', lwd = 0.2)+
  geom_hline(aes(yintercept = upper), lty = 2, colour = 'white')+
  geom_hline(aes(yintercept = lambda))+
  geom_rect(aes(xmin = as.numeric(stage)-hwid, 
                xmax = as.numeric(stage)+hwid,  ymin = lambda, ymax = value),
            colour = 'black', lwd = 0.2,
            position = position_dodge())+
  facet_grid(rows = vars(pop),scale = 'free')+
  theme_bw()+
  scale_fill_manual(values = c("#007765","#6AA341", "gray90"),
                    name = NULL)+
  #geom_hline(aes(yintercept = c(lambda)), lty = 2)+
  theme(strip.background = element_rect(fill = '#31688EFF', linewidth = 0),
        strip.text.y.right = element_text(angle = 0, face = 'bold', colour = 'white'),
        panel.grid.major.x = element_blank(),
        legend.position = 'inside',
        legend.position.inside = c(0.93,0.91),
        legend.background = element_rect(colour = 'grey'))+
  ylab(expression(italic(r)))+
  xlab('Stage adjusted')

ggsave(paste0('./figures/', 'sensitivity_lambda_diff.png'), 
       dpi = 300,
       height = 6, width = 8, units = 'in')


# average across sites ----------------------------------------------------

fecundity <- (mean(param_selecteddf$mean[param_selecteddf$name == 'F_reproduction'])*0.5)*mean(clutch_sizes)
ASAsur<- mean(param_selecteddf$mean[grep('survival_A', param_selecteddf$name)])
Jt<- mean(param_selecteddf$mean[grep('survival_J', param_selecteddf$name)])


  
stage.mat_mean <- matrix(0, 
                    nrow = length(stages),
                    ncol = length(stages),
                    byrow = TRUE,
                    dimnames = list(stages,
                                    stages))

# fecundity
fecund_ages <- which(stages %in% paste0('A', 1:(max_age-1)))
stage.mat_mean[1,fecund_ages] <- fecundity
# juvenile survival + growth
stage.mat_mean[2:3,which(stages == 'J')] <- Jt*c(transition, 1-transition)

# adult survival + growth
for (i in 2:(max_age+2)) {
  stg <- stages[i]
  if (stg == 'SA') stage.mat_mean[i+2,i] <- ASAsur
  if (stg != 'SA' & i < (max_age+2)) stage.mat_mean[i+1, i] <- ASAsur
}

stage.mat_mean[length(stages), length(stages)] <- 0.001

stage.mat_mean
eigen_analysis(stage.mat_mean)

emat <- popbio::elasticity(stage.mat_mean)

sumfec <- sum(emat[1,!(stages %in% c('J', 'SA'))])
sumgrow <- sum(emat[2:3,1])
sumgrow2 <- sum(emat[4,2])
sumstasis <- sum(emat[3:length(stages),3:length(stages)])


par(mar=c(5.1, 4.1, 4.1, 2.1))
barplot(c(fecundity = sumfec, Juv = sumgrow, SA = sumgrow2, Adult = sumstasis), 
        col = c("gray50", "#6AA341",'#6AA990', "gray90"), ylim = c(0,0.4),
        xlab = 'Fecundity / growth + Survival', ylab = 'Elasticity')
par(mar=c(0, 2, 0, 4))
pie(c(sumfec,sumgrow,sumgrow2, sumstasis), col = c("gray50", "#6AA341",'#6AA990', "gray90"),
    labels=c("fecundity","growth+survival J","growth+survival SA", "growth+survival A"))

png(filename = './figures/sensitivity_pie.png',res = 300, width = 7, 
    units = 'in', height = 5)
par(mar=c(0, 2, 0, 4))
pie(c(sumfec,sumgrow,sumgrow2, sumstasis), col = c("gray50", "#6AA341",'#6AA990', "gray90"),
    labels=c("fecundity","growth+survival J","growth+survival SA", "growth+survival A"))

dev.off()
par(mar=c(5.1, 4.1, 4.1, 2.1))

ematround <- round(emat, 3)
ematround[mat == 0] <- ''
fx_elasticity <- ematround %>% as.data.frame %>% 
  mutate(Elasticity = rownames(.)) %>% relocate(Elasticity) %>% 
  flextable() %>% 
  autofit() %>% 
    theme_box() %>% 
  font(fontname = 'calibri', part = 'all') %>% 
  bold(j = 1) %>% 
  hline(part = 'header', border = fp_border_default(width = 2)) %>% 
  vline(j = 1, border = fp_border_default(width = 2)) %>% 
  border_outer(border = fp_border_default(width = 2)) %>% 
  bg(j = -1,bg = 'grey90') %>% 
  bg(j = 1, part = 'header',bg = 'grey90');fx_elasticity

saveRDS(fx_elasticity,'./output/sensitivity_elasticity_table.rds')

fx_mat <- round(stage.mat_mean,3) %>% as.data.frame %>% 
  mutate_all(~ifelse(. == 0, '0', .)) %>% 
  mutate(Stage = rownames(.)) %>% relocate(Stage) %>% 
  flextable() %>% 
  autofit() %>% 
  theme_box() %>% 
  font(fontname = 'calibri', part = 'all') %>% 
  bold(j = 1) %>% 
  hline(part = 'header', border = fp_border_default(width = 2)) %>% 
  vline(j = 1, border = fp_border_default(width = 2)) %>% 
  border_outer(border = fp_border_default(width = 2)) %>% 
  bg(j = -1,bg = 'grey90') %>% 
  bg(j = 1, part = 'header',bg = 'grey90');fx_mat

saveRDS(fx_mat,'./output/sensitivity_stage_mat.rds')
