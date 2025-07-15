
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
                                new_stage = NA,
                                pop = rep(popname, nInd(glall)),
                                status = 'alive',
                                mother = FALSE,
                                eggs = NA,
                                f = NA,
                                year = 1)
 # rel <- gl.grm(glall, verbose = 0)

 #  for (i in 1:nrow(rel)) glall@other$ind.metrics$f[i] <- rel[i,i] - 1
 
  pop(glall) <- glall@other$ind.metrics$pop
  return(glall)
}

em.gl.reproduce <- function(gl){

  meta <- gl@other$ind.metrics
  popname <- unique(meta$pop)
  
  if(sum(meta$status == 'alive')>0){
  if(length(popname) > 1) stop('too many pop names')
  if(sum(meta$mother)>0){
  mum_index <- which(meta$mother)
  
  m_index <- which(meta$sex == 'male' & 
                     !(meta$stage %in%  c('J', 'SA')) & 
                     meta$status != 'dead')
 

  glbabies <- list()
  for (clutch in 1:length(mum_index)) {
    dad <- sample(m_index, 1)
    
    glbabies[[clutch]] <- gl.sim.offspring(mothers = gl[mum_index[clutch],],
                     fathers = gl[dad,], 
                     noffpermother = meta$eggs[mum_index[clutch]])
  }
  
  babysexes <- do.call('c', lapply(glbabies, function(x) x@other$sex))
  glJuvs <- do.call('rbind', glbabies)
  glJuvs@other$ind.metrics <- data.frame(sex = babysexes, 
                                         stage = 'J',
                                         new_stage = NA,
                                         pop = popname,
                                         status = 'alive',
                                         mother = FALSE,
                                         eggs = NA,
                                         f = NA,
                                         year = NA)
  glJuvs
  pop(glJuvs) <- glJuvs@other$ind.metrics$pop
  }

  glAdulst <- gl[meta$status == 'alive',]
  
  if(sum(meta$mother)>0){
  glAll <- rbind(glAdulst, glJuvs)
  glAll@other$ind.metrics <- rbind(glAdulst@other$ind.metrics,
                                   glJuvs@other$ind.metrics)
  }
  
  if(sum(meta$mother)== 0) glAll <- glAdulst
  glAll@other$ind.metrics$mother <- F
  glAll@other$ind.metrics$eggs <- NA
  glAll@other$ind.metrics$new_stage <- factor(glAll@other$ind.metrics$new_stage,
                                              levels = levels(glAll@other$ind.metrics$stage))
  
  }
  
  if(sum(meta$status == 'alive')==0) glAll <- NULL
   return(glAll)
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
  
  ### initialise model ----
  # initial population sizes across stages and sites
  initial_abundance <- em.initial.N(initial_ab, stage_distribution, populations)
  n_stages <- nrow(initial_abundance)
  n_populations <- ncol(initial_abundance)
  
  
  if (length(f_reproducing) == 1) f_reproducing <- rep(f_reproducing, n_populations)
  
  # abundance array
  N <- array(dim = c(n_stages, n_populations, time_steps, replicates))
  N[, , 1, ] <- initial_abundance
  
  # genetics 
  suppGENETICS <- GENETICS
  G <- array(dim = c(n_populations, time_steps, replicates))
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
    
    glN_init <- lapply(1:n_populations,
                  function(x) em.initial.Ngl(snpspop[[x]],
                                             initial_abundance[,x], 
                                             popname = populations[x]))
    names(glN_init) <- populations
    glN <- do.call('rbind', glN_init)
    glN_meta <- lapply(glN_init, function(x) x@other$ind.metrics)
    meta_init <- do.call('rbind', glN_meta)
    meta_init$new_stage <- factor(meta_init$new_stage, levels = stages)
    meta_init$pop <- factor(meta_init$pop, levels = populations)
    meta_init$sex <- factor(meta_init$sex) 
    glN@other$ind.metrics <- meta_init
      }
  
  for (i in 1:replicates) {
    if (GENETICS) {
      he_start <- sapply(glN_init, function(x) mean(gl.He(x)))
      G[, 1,i] <- he_start
      genetics_meta <- meta_init
    }
    
    for (yr in 2:time_steps) {
      print(yr)
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
          
          if(suppGENETICS){
            
            glsup <- gl.sim.ind(snps, n = sum(n.supp))
            
            glsup@other$ind.metrics <- data.frame(sex = sample(c('female', 'male'), 
                                                               nInd(glsup), T),
                                                  stage = stages[stage.supp],
                                                  new_stage = NA,
                                                  pop = rep(populations, n.supp),
                                                  status = 'alive',
                                                  mother = FALSE,
                                                  eggs = NA,
                                                  f = NA,
                                                  year = yr)
            pop(glsup) <- glsup@other$ind.metrics$pop
            
            if(GENETICS){
            meta_notsupp <- glN@other$ind.metrics
            glN <- rbind(glN, glsup)
            glN@other$ind.metrics <- rbind(meta_notsupp, glsup@other$ind.metrics)
            genetics_meta <- glN@other$ind.metrics
            }
            if(!GENETICS){
              glN <- glsup
              genetics_meta <- glN@other$ind.metrics
              GENETICS <- TRUE
            }
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
      
      # genetics survival
      if(GENETICS){
        print('surv')
      for (sitename_gl in populations) {
      xsites <- which(sitename_gl == populations)
      survivors <- N[, xsites, yr, i]
      names(survivors) <- stages
  
      
      diers <- Nprevious[,xsites] - survivors
      
      for (dd in 1:n_stages) {
        
        sindex <- which(genetics_meta$stage == stages[dd] & genetics_meta$pop == sitename_gl)
        
        sindexDead <- sindex[sample(1:length(sindex), diers[dd])]
        genetics_meta$status[sindexDead] <- 'dead'
        
      }

      }
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
      
      # Genetics
      if(GENETICS){
        print('transition')
       
        new_stages <- genetics_meta$new_stage
        
        for (sitename_gl in populations) {
            xsites <- which(sitename_gl == populations)
            
          for (stran in stages) {
            
            stage_tran_mat <- transition_mat[,stran]
            newstage <- names(stage_tran_mat[stage_tran_mat>0])
            current_stage <- which(genetics_meta$stage==stran & 
                                     genetics_meta$pop == sitename_gl &
                                     genetics_meta$status == 'alive')
            if(length(newstage)>1){
              n_new_stage <- new.tran[stage_tran_mat>0,xsites]
              new_stages[current_stage] <- rep(newstage, n_new_stage)
            }
            
            if(length(newstage) == 1){
            
            new_stages[current_stage] <- newstage 
            }
 
          }
          
          
        }
        genetics_meta$new_stage <- new_stages
      }
      
      
      ## reproduction ---------
     
      females <- round(colSums(N[3:n_stages,,yr,i])/2)
      
      if(GENETICS) {
        print('repo')
        females<- table(genetics_meta$pop[genetics_meta$sex == 'female' &
                        genetics_meta$status == 'alive' &
                        genetics_meta$new_stage %in% stages[3:n_stages]])
        
        males <- table(genetics_meta$pop[genetics_meta$sex == 'male' &
                         genetics_meta$status == 'alive' &
                         genetics_meta$new_stage %in% stages[3:n_stages]])
        females[males == 0] <- 0

      }
      
      mothers <- sapply(1:n_populations, 
                        function(x) sum(rbinom(females[x], 1, f_reproducing[x])))
      clutches <- lapply(mothers, em.juveniles, clutch_sizes)
      
      if(GENETICS){
        
        for (sitename_gl in populations) {
            xsites <- which(sitename_gl == populations)
            
            if(sitename_gl %in% populations[mothers>0]){
          breeders <- which(genetics_meta$new_stage != 'SA' &
                                  genetics_meta$pop == sitename_gl &
                                  genetics_meta$status == 'alive' &
                             genetics_meta$sex == 'female')
         
          successful_breeders <- sample(breeders, mothers[xsites])
          if(length(breeders)==1) successful_breeders <- breeders
          print(paste('length B', length(breeders)))
          genetics_meta$mother[successful_breeders] <- TRUE
          
          genetics_meta$eggs[successful_breeders] <- clutches[[xsites]]
         
        
         
         }
            }
           
      
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
        print('K')
    
        for(sitename_gl in populations){
    
            xsites <- which(sitename_gl == populations)
        Ksite <- K_correct[,xsites]
        for (xstage in 2:n_stages) {
          sindex <- which(genetics_meta$new_stage == stages[xstage] & 
                            genetics_meta$pop == sitename_gl & 
                            genetics_meta$status == 'alive')
          
          sindexDead <- sindex[sample(1:length(sindex), Ksite[xstage])]
          genetics_meta$status[sindexDead] <- 'deadK'
          
        }
        if (Ksite[1] > 0) {
          sindex <- which(complete.cases(genetics_meta$eggs) & 
                            genetics_meta$pop == sitename_gl)
          
          juvDead <- sample(rep(sindex, genetics_meta$eggs[sindex]), Ksite[1])
          
          juvDeadCount <- table(juvDead)
          juvDeadmum <- as.numeric(names(juvDeadCount))
          genetics_meta$eggs[juvDeadmum] <- genetics_meta$eggs[juvDeadmum]-juvDeadCount
          
          genetics_meta$mother[complete.cases(genetics_meta$eggs) & genetics_meta$eggs ==  0] <- F
          
        }
          
        }
        
      }
      
      
      
      if(GENETICS){
      # adults
      testA <-genetics_meta %>%
        filter(status == 'alive') %>%
        group_by(pop, new_stage) %>%
        summarise(n = n()) %>%
        pivot_wider(names_from = pop, values_from = n)
      
      # juves
      testJ <- genetics_meta %>%
        group_by(pop) %>%
        summarise(J = sum(eggs, na.rm=T)) %>%
        pivot_wider(names_from = pop, values_from = J)
      
      testMat <- N[, , yr, i]
      testMat[testMat >0] <- 0
      
      testMat[1,populations %in% names(testJ)] <- as.numeric(testJ)
      
      testAJ <- testA[,-1]
      testAJ[is.na(testAJ)]  <- 0
      
      testMat[stages %in% testA$new_stage,
              populations %in% names(testAJ)] <- as.matrix(testAJ[order(testA$new_stage),])
      
     
      ## checks ----
      if(sum(!(N[, , yr, i] == testMat))>0) stop('N not equal to nG')
      N[, , yr, i] == testMat
      N[, , yr, i]
      testMat
      
      }
      
      ## Actual genetics ----
      if(GENETICS){
        print('genetics!!')
        glNsave <- glN
        glN@other$ind.metrics <- genetics_meta
        glN@other$ind.metrics$stage <- genetics_meta$new_stage
        glN@other$ind.metrics$new_stage <- NA
        pop(glN) <- glN@other$ind.metrics$pop
        glN_pops <- seppop(glN)
        
        next_generation <- lapply(glN_pops, em.gl.reproduce)
        
        notnull <- !sapply(next_generation, is.null)
        next_generation_occ <- next_generation[notnull]
        if('genlight' %in% sapply(next_generation_occ, class)){
        he <- sapply(next_generation_occ, function(x) mean(gl.He(x)))
        hepops <- which(populations %in% names(he))
        G[hepops,yr, i] <- he
        
        new_meta <- lapply(next_generation_occ, function(x) x@other$ind.metrics)
        glN <- do.call('rbind', next_generation_occ)
        glN@other$ind.metrics <- do.call('rbind', new_meta)
        glN@other$ind.metrics$year <- yr
        genetics_meta <- glN@other$ind.metrics
        }else{
          cat('No more pops GENETICS set to false \n')
          GENETICS <- FALSE
        }
        
      }

    }
    
  }

  return(list(N = N, G = G))
}





