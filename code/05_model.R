
# funtions -------------------------------------------------

logit <- function(p) log(p/(1-p))
unlogit <- function(x) exp(x) / (1 + exp(x))

em.initial.N <- function(init, stage_distribution, p){
  
  adult.stages <- stage_distribution[-1]
  juv.stage <- stage_distribution[1]
  
  initial_abundance <- round(rbind(J=(init/(1-juv.stage))*juv.stage,
              sapply(init, function(x) (adult.stages/sum(adult.stages))*x)))
  
   
  rownames(initial_abundance) <- names(stage_distribution)
  colnames(initial_abundance) <- p
  initial_abundance
  
}


em.juveniles <- function(mothers, clutch_sizes){
  
  n <- length(clutch_sizes)
  prob_seq <- round(seq(0.1,0.9, length.out = n),2)
  beta2 <- dbeta(prob_seq,
                 2,2)
  clutch_prob <-beta2/sum(beta2)
  if(mothers >0){
    sample(clutch_sizes,
           size = mothers,
           prob = clutch_prob,
           replace = T)
  }else{
    0
  }
}

# pva model ---------------------------------------------------------------


em.pva_simulator <- function(populations = c('CA', 'JE', 'JW', 'MA'),
                             stages = c('J', 'SA','A1','A2','A3', 'A4'),
                             stage_distribution, 
                             initial_ab, # adult abundance for populations
                             survival, # survival of adults and SA at sites
                             survival_J, # juvenile survival
                             env_stoch, # sd on survival
                             transition_mat, # transition prob to SA
                             f_reproducing, # proportion of females reproducing
                             clutch_sizes, # clutch size range   
                             K = K, # carrying capcity applied to adults and SA
                             time_steps = 8, # time
                             replicates = 1){
  
  initial_abundance <- em.initial.N(initial_ab, stage_distribution, populations)
  n_stages <- nrow(initial_abundance)
  n_populations <- ncol(initial_abundance)
  
  
  N <- array(dim = c(n_stages, n_populations, time_steps, replicates))
  N[, , 1, ] <- initial_abundance
  
 
  # for loop replicate (i)
  for (i in 1:replicates) {
    # for loop year (yr)
    for (yr in 2:time_steps) {
      
      ## survival ---------
      ### environmental stoch (adults and SA) ---------------------------
      survival_env <- unlogit(sapply(1:n_populations,
                                     function(x) rnorm(1, 
                                                       logit(survival[x]),
                                                       env_stoch[x])))
      
      p_environment_effect <-sapply(1:n_populations,
                                    function(x) pnorm(logit(survival_env)[x], 
                                                              logit(survival[x]),
                                                              env_stoch[x]))
      
      survival_J_env <- unlogit(sapply(1:n_populations, 
              function(x) qnorm(p_environment_effect[x],
                                logit(survival_J),
                                env_stoch[x])))
      
      
      survival_year <- cbind(J= survival_J_env,
                        sapply(1:(n_stages-2), function(x) survival_env),
                        A4 = 0)
      
      for (sites in 1:n_populations) {
        surv <- survival_year[sites,]
        prev.N <- N[,sites , yr-1, i]
        new.N <- sapply(1:length(prev.N), function(x) sum(rbinom(prev.N[x],1,surv[x])))
        N[, sites, yr, i] <- new.N
          }
      
      ## transition ---------
      
      
      tran.mat <- transition_mat
      tran.mat[tran.mat<1 & tran.mat>0] <- 0
      SAtran <- transition_mat[2,1]
      transitioning_J <- N[1, , yr, i]
      new.tran <-  (tran.mat %*% N[, , yr, i])
      new.SA <- sapply(transitioning_J, function(x) sum(rbinom(x,1, SAtran)))
      new.A <- transitioning_J-new.SA
      new.tran[2,] <- new.SA
      new.tran[3,] <- new.A
      N[, , yr, i] <- new.tran
      
      ## reproduction ---------
      females <- round(colSums(N[3:n_stages,,yr,i])/2)
      mothers <- sapply(females, function(x) sum(rbinom(x, 1, f_reproducing)))
      clutches <- lapply(mothers, em.juveniles, clutch_sizes)
      
      juvs <- sapply(clutches, sum)
      
      N[1, , yr, i]<- juvs
      
      
      ## carrying capcity ---------
      adults <- colSums(N[2:n_stages,,yr,i])
      
      if(sum(adults > K)>0){
      adults.dist <- stage_distribution[-1]/sum(stage_distribution[-1])
        K.stages <- round(sapply(K, function(x) adults.dist*x))
        for(age in 2:n_stages){
          current_N <-  N[age, , yr, i]
          N[age, , yr, i]<- ifelse(current_N > K.stages[age-1,],
                                   K.stages[age-1,], current_N) 
        }
        
      }
      
    }# }yr
  }# }i
  
  return(N)
   
  }



