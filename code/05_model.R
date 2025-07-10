
# funtions -------------------------------------------------

# logit and unlogit functions
logit <- function(p) log(p/(1-p))
unlogit <- function(x) exp(x) / (1 + exp(x))


# initialise starting population size across stages based on stable stage
# distributions derived from best guess stage matrix. initial population size
# is estimated for adults and sub-adults only, so juvenile abundance added.

em.initial.N <- function(init_A_SA,               
                         stage_distribution, 
                         population_names){                 
  
  # summaries adult stable stages and juvenile stage
  adult.stages <- stage_distribution[-1] 
  adult.stages.proportion <- adult.stages/sum(adult.stages)
  juv.stage <- stage_distribution[1] #0.35              
  
  # Calculate N
  total_N <- init_A_SA/(1-juv.stage)
  juv_N <- (total_N)*juv.stage
  adult_N <- sapply(init_A_SA, function(x) (adult.stages.proportion)*x)
  
  initial_abundance <- round(rbind(juv_N, adult_N)) 
  
  # adding names to matrix 
  rownames(initial_abundance) <- names(stage_distribution)
  colnames(initial_abundance) <- population_names
  initial_abundance
  
}

# Function to determine clutch sizes of mothers
em.juveniles <- function(mothers, clutch_sizes){
  
  n <- length(clutch_sizes)
  prob_seq <- seq(0.1,0.9, length.out = n)
  beta2 <- dbeta(prob_seq,2,3) 
  clutch_prob <-beta2/sum(beta2)
  
  if(mothers >0){# necessary because mothers = 0 returns nothing
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
                             survival_J, # juvenile survival per site
                             survival_logit_sd = NULL,
                             site_adjust = NULL,
                             env_stoch = NULL, # sd on survival
                             transition_mat, # transition prob to SA
                             f_reproducing, # proportion of females reproducing
                             clutch_sizes, # clutch size range   
                             K = K, # carrying capcity applied to adults and SA
                             time_steps = 8, # time
                             replicates = 1,
                             supp = FALSE,
                             n.supp = c(10,10,10,10),
                             stage.supp = 3,
                             when.supp = seq(2,5)){
  
  ## initialise model
  # initial population sizes across stages and sites
  initial_abundance <- em.initial.N(initial_ab, stage_distribution, populations)
  n_stages <- nrow(initial_abundance)
  n_populations <- ncol(initial_abundance)
  
  if (length(f_reproducing) == 1) f_reproducing <- rep(f_reproducing, n_populations)
  
  # abundance array
  N <- array(dim = c(n_stages, n_populations, time_steps, replicates))
  N[, , 1, ] <- initial_abundance
  
  
  # survival adjustment for sites (calibration only)
  
  if(!is.null(site_adjust)){
    
    survival <- unlogit(sapply(site_adjust, 
                   qnorm, 
                   mean = logit(survival), 
                   sd = survival_logit_sd))
    survival_J <- unlogit(sapply(site_adjust, 
                               qnorm, 
                               mean = logit(survival_J), 
                               sd = survival_logit_sd))
    
  }
 
  
  for (i in 1:replicates) {
    for (yr in 2:time_steps) {
      
      Nprevious <- N[, , yr-1, i]
      ## supplementation -----
      if (supp) {
        if(yr %in% when.supp){
          for (sites in 1:n_populations) {
            
            prev.N <- Nprevious[,sites]
            # demographic stochasitcity (binomial/bernoulli distribution)
            new.N <- prev.N[stage.supp] + n.supp[sites]
            Nprevious[stage.supp,sites] <- new.N
          }
        }
      }
      
      ## survival ---------
      ### environmental stoch (adults and SA) ---------------------------
      
      survival_env <- survival
      survival_J_env <- survival_J
      if(!is.null(env_stoch)){
      survival_env <- unlogit(sapply(1:n_populations,
                                     function(x) rnorm(1, 
                                                       logit(survival[x]),
                                                       env_stoch[x])))
      
      p_environment_effect <-sapply(1:n_populations,
                                    function(x) pnorm(logit(survival_env)[x], 
                                                              logit(survival[x]),
                                                              env_stoch[x]))
      ### environmnetal stoch J - same effect as adults
      survival_J_env <- unlogit(sapply(1:n_populations, 
              function(x) qnorm(p_environment_effect[x],
                                logit(survival_J[x]),
                                env_stoch[x])))
      }
      
      ### yearly environmental adjusted survival 
      survival_year <- cbind(J= survival_J_env,
                        sapply(1:(n_stages-2), function(x) survival_env),
                        oldest = 0)
      # surviving
      for (sites in 1:n_populations) {
        surv <- survival_year[sites,]
        prev.N <- Nprevious[,sites]
        # demographic stochasitcity (binomial/bernoulli distribution)
        new.N <- sapply(1:length(prev.N), function(x) sum(rbinom(prev.N[x],1,surv[x])))
        N[, sites, yr, i] <- new.N
          }
      
      ## transition ---------
      new.tran <-  (transition_mat %*% N[, , yr, i])
      
      # Juvenile stage transition 
          stage_juv <- transition_mat[,1]
          stage_Njuv <- N[1, ,yr,i]
          juv_tran <- stage_juv[stage_juv>0]
          # demographic stochasitcity (binomial/bernoulli distribution)
          stage1 <- sapply(stage_Njuv, function(x) sum(rbinom(x,1, juv_tran[1])))
          stage2 <- stage_Njuv-stage1
          new.tran[stage_juv>0,] <- rbind(stage1, stage2)
      
          # update N
      N[, , yr, i] <- new.tran
      
      ## reproduction ---------
      females <- round(colSums(N[3:n_stages,,yr,i])/2)
      mothers <- sapply(1:n_populations, 
                        function(x) sum(rbinom(females[x], 1, f_reproducing[x])))
      clutches <- lapply(mothers, em.juveniles, clutch_sizes)
      
      juvs <- sapply(clutches, sum)
      
      N[1, , yr, i]<- juvs
      
      
      ## carrying capacity ---------
      # adults <- colSums(N[2:n_stages,,yr,i])
      # 
      # if(sum(adults > K)>0){
      #   sites_K <- which(adults>K)
      #     
      #   adult.stages <- stage_distribution[-1] 
      #   adults.dist <- adult.stages/sum(adult.stages)
      #   K.stages <- round(sapply(K[sites_K], function(x) adults.dist*x))
      #   for(age in 2:n_stages){
      #     current_N <-  N[age, sites_K, yr, i]
      #     N[age, sites_K, yr, i]<- ifelse(current_N > K.stages[age-1,],
      #                              K.stages[age-1,], current_N) 
      #   }
      #   
      # }
      
      yearN <- N[, , yr, i]
      
      yearNstages <- sapply(1:ncol(yearN), function(x) rep(stages, yearN[,x]))
      
      total <- colSums(yearN)
   
      extra_N <- ifelse((total - K)>0, total - K, 0)
      
      
      K_correct <- sapply(1:n_populations, 
                          function(x) table(factor(sample(yearNstages[[x]], 
                                                      size = extra_N[x],
                                                      replace = F), 
                                               levels = stages)))
      N[, , yr, i] <- N[, , yr, i] - unname(K_correct)
    
    }
    
  }

  return(N)
}

  


