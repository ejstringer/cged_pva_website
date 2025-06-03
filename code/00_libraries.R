
# libraries ---------------------------------------------------------------
library(terra)
library(tidyverse)
library(sf)
library(lhs)
library(stringi)
library(parallel)
library(flextable)
library(scales)
library(lubridate)
library(pscl)
# load --------------------------------------------------------------------

# functions ---------------------------------------------------------------

eigen_analysis <- function(A) {
  ev <- eigen(A)
  lmax <- which(Re(ev$values) == max(Re(ev$values)))
  lambda <- Re(ev$values[lmax])
  W <- ev$vectors
  w <- abs(Re(W[, lmax]))
  stable_age <- w/sum(w)
  return(list(lambda, stable_age))
}

em.extract_ASA <- function(sim){
  tstep <- dim(sim)[3]
  nstg <- dim(sim)[1]
  asa <- sapply(1:tstep, function(x) colSums(sim[2:nstg,,x,]))
  colnames(asa) <- 2013:(2012+tstep)
  rownames(asa) <- c('CA', 'JE', 'JW', 'MA')
  asa
}

em.sample.distribution <- function(x, samplesize = 10000){
  #https://stats.stackexchange.com/questions/191725/sample-from-distribution-given-by-histogram
  xhist=hist(x,freq=FALSE, col=rgb(0,0,1,1/4))
  # sample from it
  bins=with(xhist,sample(length(mids),samplesize,p=density,replace=TRUE)) # choose a bin
  result=runif(length(bins),xhist$breaks[bins],xhist$breaks[bins+1]) # sample a uniform in it
  hist(result,freq=FALSE,add=TRUE,bord=1, col = rgb(1,0,0,1/4))
  
  return(result)
}

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
