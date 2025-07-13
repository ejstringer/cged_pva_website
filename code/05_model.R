
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

em.initial.Ngl <- function(gl, initial_popN, popname){
  
  initial_adults <- initial_popN[-1]
  N_size <- sum(initial_adults)
  gl2 <- gl.sim.ind(gl, n = N_size)
  pop(gl2) <- rep(names(initial_adults), initial_adults)
  gl2@other$sex <- rep(c('female', 'male'), 
                       length.out = N_size)[sample(1:N_size, N_size)]
  
  
  gljuv <- gl.sim.offspring(mothers = gl2[gl2@other$sex == 'female',],
                            fathers = gl2[gl2@other$sex == 'male',],
                            noffpermother = mean(clutch_sizes))[1:initial_popN[1],]
  pop(gljuv) <- rep('J', nInd(gljuv))
  glall <- rbind(gl2, gljuv)
  glall@other$ind.metrics <- data.frame(sex = sample(c('female', 'male'), 
                                             nInd(glall), T),
                                stage = glall@pop,
                                pop = rep(popname, nInd(glall)),
                                status = 'alive')
  rel <- gl.grm(glall)

  for (i in 1:nrow(rel)) glall@other$ind.metrics$f[i] <- rel[i,i] - 1
 
  return(glall)
}

em.gl.reproduce <- function(gl, mums, babies){
  if (mums > 0) {
    
 
  meta <- gl@other$ind.metrics
  f_index <- which(meta$sex == 'female' & meta$stage != 'SA')
  m_index <- which(meta$sex == 'male' & meta$stage != 'SA')
  mum_index <- sample(f_index,mums)
  if(length(mum_index) != length(babies)) stop('fix mums and bubs!')
  glbabies <- list()
  for (i in 1:length(mum_index)) {
    dad <- sample(m_index, 1)
    glbabies[[i]] <- gl.sim.offspring(mothers = gl[mum_index[i],],
                     fathers = gl[dad,], 
                     noffpermother = babies[i])
  }
  
  babysexes <- do.call('c', lapply(glbabies, function(x) x@other$sex))
  glJuvs <- do.call('rbind', glbabies)
  glJuvs@other$ind.metrics <- data.frame(sex = babysexes, 
                                         stage = 'J',
                                         pop = gl@other$ind.metrics$pop[1],
                                         status = 'alive',
                                         f = NA)
  glJuvs
  glall <- rbind(gl, glJuvs)
  
  glall@other$ind.metrics <- rbind(gl@other$ind.metrics, glJuvs@other$ind.metrics)
  
  rel <- gl.grm(glall)
  
  for (i in 1:nrow(rel)) glall@other$ind.metrics$f[i] <- rel[i,i] - 1
  
  }
  if(mums == 0) glall <- gl
  pop(glall)<- glall@other$ind.metrics$stage
   return(glall)
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
                             when.supp = seq(2,5),
                             GENETICS = FALSE,
                             snps){
  
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
 
  
  if(GENETICS){
    snpspop <- seppop(snps)
    if(sum(!(names(snpspop) == populations))>0) stop('genlight pops do not match')
    
    glN <- lapply(1:n_populations,
                  function(x) em.initial.Ngl(snpspop[[x]],
                                             initial_abundance[,x], 
                                             popname = populations[x]))
    names(glN) <- populations
    glN

      }
  
  for (i in 1:replicates) {
    for (yr in 2:time_steps) {
      print('supp')
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
      print('sur')
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
      # genetics survival
      if(GENETICS){
      for (sites in 1:n_populations) {
        if(!is.null(glN[[sites]])){
      survivors <- N[, sites, yr, i]
      names(survivors) <- stages
      diers <- N[, sites, yr-1, i] - survivors
      glNpop <- glN[[sites]]
      for (dd in 1:n_stages) {
        
        sindex <- which(glNpop@other$ind.metrics$stage == stages[dd])
        
        sindexDead <- sindex[sample(1:length(sindex), diers[dd])]
        glNpop@other$ind.metrics$status[sindexDead] <- 'dead'
        
      }
      if('alive' %in% glNpop@other$ind.metrics$status){
      glN[[sites]] <- glNpop[which(glNpop@other$ind.metrics$status == 'alive'),]
      }else{
        glN[[sites]] <- NULL
      }
        }
      }
      ## transition ---------
      print('trans')
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
      
      # Genetics
      if(GENETICS){
        for (sites in 1:n_populations) {
          if(!is.null(glN[[sites]])){
          glNtran <- glN[[sites]]
          new_stages <- glNtran@other$ind.metrics$stage
          
          
          for (stran in stages) {
            
            stage_tran_mat <- transition_mat[,stran]
            newstage <- names(stage_tran_mat[stage_tran_mat>0])
            current_stage <- which(glNtran@other$ind.metrics$stage==stran)
            if(length(newstage)>1){
              n_new_stage <- new.tran[stage_tran_mat>0,sites]
              new_stages[current_stage] <- rep(newstage, n_new_stage)
            }
            
            if(length(newstage) == 1){
            
            new_stages[current_stage] <- newstage 
            }
            
            
            
          }
          glNtran@other$ind.metrics$stage <- new_stages
          glN[[sites]] <- glNtran
          
          }
        }
        
      }
      
      
      ## reproduction ---------
      print('repro')
      females <- round(colSums(N[3:n_stages,,yr,i])/2)
      if(GENETICS) {
        
       gfemales<- sapply(glN, function(x) table(x@other$ind.metrics$sex,
                                      x@other$ind.metrics$stage)[1,])
       females <- colSums(gfemales[rownames(gfemales) %in% stages[3:n_stages],])
      }
      mothers <- sapply(1:n_populations, 
                        function(x) sum(rbinom(females[x], 1, f_reproducing[x])))
      clutches <- lapply(mothers, em.juveniles, clutch_sizes)
      
      if(GENETICS){
        
           glNew <- lapply(1:n_populations,
                                     function(x) em.gl.reproduce(gl = glN[[x]],
                                                                 mums = mothers[x],
                                                                 babies = clutches[[x]]))
        
        glN <- glNew
      }
      
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
      
      # Genetics 
      if(GENETICS){
        for(sites in 1:n_populations){
          if(!is.null(glN[[sites]])){
        glNK <- glN[[sites]]
        Ksite <- K_correct[,sites]
        for (xstage in 1:n_stages) {
          sindex <- which(glNK@other$ind.metrics$stage == stages[xstage])
          
          sindexDead <- sindex[sample(1:length(sindex), Ksite[xstage])]
          glNK@other$ind.metrics$status[sindexDead] <- 'dead'
          
        }
        glN[[sites]] <- glNK[which(glNK@other$ind.metrics$status == 'alive'),]
          }
        }
      }
      
      
    print(yr)
    }
    
  }

  return(N)
}

  


