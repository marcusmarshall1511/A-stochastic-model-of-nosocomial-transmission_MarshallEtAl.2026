library(ggplot2)

library(rstudioapi)




current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))



#Gillespie code for sanity check of results


################################################################################

#Rate function
Rates <- function(x){
  l <- x[[1]]
  EHQ <- x[[2]][1]
  IHQ <- x[[2]][2]
  EHG <- x[[3]][1]
  IHG <- x[[3]][2]
  EP <- x[[4]][1]
  IP <- x[[4]][2]
  
  mu_IP <<- IP * deltas$deltaP 
  mu_EP <<- EP * deltas$deltaP 
  
  delta_iso <<- deltas$Hyp[[1]][l+1]
  
  lambda_HQ <<- (pops$NHQ - EHQ - IHQ)*(betas$Beta_O[[1]][l+1] + IHQ*betas$Beta_HQ + IHG*betas$Beta_H + IP*betas$Beta_HQP)
  lambda_HG <<- (pops$NHG -EHG -IHG)*(IHG*betas$Beta_HG + IHQ*betas$Beta_H + IP*betas$Beta_HGP)
  lambda_P <<- (pops$NP -EP -IP)*(IHQ*betas$Beta_HQP + IHG*betas$Beta_HGP + IP*betas$Beta_P)
  
  eta_HQ <<- EHQ * zetas$zeta_HQ
  eta_HG <<- EHG * zetas$zeta_HG
  eta_P <<- EP * zetas$zeta_P
  
  gamma_rate <<- IHQ * gammas$gamma_HQ + IHG * gammas$gamma_HG + IP * gammas$gamma_P
  
  Theta <<- mu_IP + mu_EP + delta_iso + lambda_HQ + lambda_HG + lambda_P + eta_HQ + eta_HG + eta_P + gamma_rate
}


# Stochastic process information
{# Uses list for stoch. process, e.g., list(K, c(EHQ,IHQ),c(EHG,IHG),c(EP,IP))
# For:
# K - Removal stage of Isolated patient
# EX <- Number of exposed individuals in subpopulation X
# IX <- Number of infected individuals in subpopulation X
}


one_gillespie <- function(params=params, replicate=NULL){
  
  if (missing(replicate)){
    replicate <- 1
  }
  
  t <- 0 #Initialise the process at time = 0
  init <-  params[[6]]
  x <-init #x is how we measure the stochastic process
  EplusI <- sum(x[[1]],x[[2]],x[[3]],x[[4]])
  
  detected <- FALSE # Start undetected outbreak
  
  
  #Tracking, see legend above
  
  #1.Overall
  Y_times <- c()
  #2.-#4 Marginals
  YHQ_times <- c()
  YHG_times <- c()
  YP_times <- c()
  
  
  
  while (t < params[[5]] & !(detected) & EplusI > 0){ #While time limit is not reached, outbreak not detected, and pathogen present
    
    Rates(x) #Compute the transition rates
    
    t <- t + -log(runif(1,0,1))/Theta #time increment
    
    {
    i <- sample(1:10, 1, prob = c(mu_EP,mu_IP, delta_iso, lambda_HQ,lambda_HG,lambda_P, eta_HQ, eta_HG, eta_P,gamma_rate)/Theta)
    if (i==1){ x <- replace( x,4,list( c(x[[4]][1]-1, x[[4]][2]) ) ) } #removal of exposed patient
    else if (i==2){x <- replace( x,4,list( c(x[[4]][1], x[[4]][2]-1) ) ) } #removal of infectious patient
    else if (i==3){x <- replace( x,1,x[[1]]-1) } #progression of stage to removal of isolated patient
    else if (i==4){x <- replace( x,2,list( c(x[[2]][1]+1, x[[2]][2]) ) ) #exposure of isolation-care HCW, HQ
                                Y_times <- append(Y_times,t)
                                YHQ_times <- append(YHQ_times,t)
                                } 
    else if (i==5){x <- replace( x,3,list( c(x[[3]][1]+1, x[[3]][2]) ) ) #exposure of general-care HCW, HG
                                Y_times <- append(Y_times,t)
                                YHG_times <- append(YHG_times,t)
                                } 
    else if (i==6){x <- replace( x,4,list( c(x[[4]][1]+1, x[[4]][2]) ) ) #exposure of general patient, P
                                Y_times <- append(Y_times,t)
                                YP_times <- append(YP_times,t)
                                }
    else if (i==7){x <- replace( x,2,list( c(x[[2]][1]-1, x[[2]][2]+1) ) ) } #maturation of exposed HQ into infectious
    else if (i==8){x <- replace( x,3,list( c(x[[3]][1]-1, x[[3]][2]+1) ) ) } #maturation of exposed HG into infectious
    else if (i==9){x <- replace( x,4,list( c(x[[4]][1]-1, x[[4]][2]+1) ) ) } #maturation of exposed P into infectious
    else if (i==10){detected <- TRUE } #outbreak is detected
    } # Choosing event
    
    EplusI <- sum(x[[1]],x[[2]],x[[3]],x[[4]]) #Checking existence of pathogen
  }
  
  # print(Y_times)
  # print(YHQ_times)
  # print(YHG_times)
  # print(YP_times)
  # print(init[[1]])
  data.frame(Y = length(Y_times), YHQ = length(YHQ_times), YHG = length(YHG_times), YP = length(YP_times), T = t, replicate = replicate, init = init[[1]])
  
  
} #A single gillespie simulation

gillespie_trials <- function(trials){
  runs <- data.frame(replicate = NULL, Y = NULL, YHQ = NULL, YHG = NULL, YP = NULL, T = NULL)
  for (trial in 1:trials){
    runs <- rbind(runs, one_gillespie(params=params,replicate=trial))
    if ((100*trial/trials)%%10 == 0){
    cat("trials: ", trial," out of ",trials," (",100*trial/trials,"%) \n")
    }
  }
  runs
  
} #Numerous gillespie simulations

gillespie_outbreaks <- function(o){
  indx <- 1
  j <- 0
  runs <- data.frame(replicate = NULL, Y = NULL, YHQ = NULL, YHG = NULL, YP = NULL, T = NULL)
  while (j < o){
    run <- one_gillespie(params=params,replicate=indx)
    runs <- rbind(runs, run)
    if (run$Y > 0){
      j <- j+1
    if ((100*j/o)%%10 == 0){
      cat("Outbreaks occured: ", j," out of ",o," (",100*j/o,"%) \n")
    }
    }
    indx <- indx + 1
    
    
  }
  runs
  runs_outbreak <- runs[runs$Y>0,]
  runs_outbreak$replicate <- 1:length(runs_outbreak$Y)
  return(runs_outbreak)
} #Simulating until 'X = o' number outbreaks

################################################################################

{ #Populations
  { #If using cohorting configurations, uncomment following section
  ################################################################################
  # uncomment desired cohorting type
  Cohorting <- "S" #Strong cohorting
  # Cohorting <- "M" #Medium cohorting
  # Cohorting <- "N" #No cohorting
  
  if (Cohorting == "S"){ H <- c(2,11) } else if (Cohorting=="M"){ H <- c(6,7) } else if (Cohorting=="N"){ H <- c(13,0) }
  
  pops <- data.frame(
    NHQ = H[1], # Number of Isolation-care HCWs
    NHG = H[2], #Number of General-care HCWs
    NP = 12 #Number of General Patients
  )
  ################################################################################


    
  # For manual population setting uncomment following section
  ################################################################################
  # pops <- data.frame(
  #   NHQ = 5, # Number of Isolation-care HCWs
  #   NHG = 5, #Number of General-care HCWs
  #   NP = 12 #Number of General Patients
  # )
  ################################################################################
 }
  
  #Pathogen infectivity, using R=beta*12*7 (Infect R individuals out of 12, over 7 day period)
  { R <- 10
  beta <- R/(12*7)}
  
  #Cohorting protection
  {rI <- 0.95
  rC <- 0.9}
  
  #Inf. rates
  {
    # eps <- 1-rC #Careful
    eps <- 1 #NOT careful
    
    
    betas <- data.frame(
    Beta_O = I(list(c(0, (1-rI)*beta*2/pops$NHQ, (1-rI)*beta*2/pops$NHQ, (1-rI)*beta*2/pops$NHQ, (1-rI)*beta*2/pops$NHQ))),
    Beta_H = eps*2/3*beta, Beta_HQ = 2/3*beta, Beta_HG = 2/3*beta,
    Beta_HQP = (1-rC)*beta, Beta_HGP = beta,
    Beta_P = 0
  )}
  
  #Latent rates
  {zeta <- 1/3
  zetas <- data.frame(
    zeta_HQ = zeta, #Incubation rate of Isolation-care HCWs
    zeta_HG = zeta, #Incubation rate of General-care HCWs
    zeta_P = zeta #Incubation rate of General Patients
  )}
  
  #Removal rates
  {deltas <- data.frame(
    deltaP = 1/6, #General patients
    Hyp = I(list(c(0,0.56, 0.56, 0.56, 6.69))) #Hypoexponential for Iso patient
  )}
  
  #Detection rates
  { gamma <- 1/5
    gammas <- data.frame(
    gamma_HQ = gamma, #Detection rate of infectious Isolation-care HCWs
    gamma_HG = gamma, #Detection rate of infectious General-care HCWs
    gamma_P = gamma #Detection rate of General Patients
  )}} #PARAMETERS

{tmax <- 100000
init <- list( length(deltas$Hyp[[1]]) -1 , c(0,0), c(0,0), c(0,0) )
params <- list(betas, zetas, deltas, gammas, tmax, init)} #tmax, init, compiling params


#RUN TRIALS AND SIFT THROUGH OUTBREAKS
################################################################################
# o <- 1000
# runs_o <- gillespie_outbreaks(o)
# Y <- sum(runs_o$Y)/o
# YHQ <- sum(runs_o$YHQ)/o
# YP <- sum(runs_o$YP)/o
# 
# 
# if (pops$NHG >0){
#   YHG <- sum(runs_o$YHG)/o
# 
#   cat("On average, ", Y," exposures occured with, E(YHQ)=",YHQ,", E(YHG,)=",YHG,", E(YP,)=",YP)
# } else {
#   cat("On average, ", Y," exposures occured with, E(YHQ)=",YHQ," and E(YP,)=",YP)
# }


# CONDITION ON OUTBREAK ANALYTICALLY
################################################################################

trials <- 100000


#The following is for the isolated patient case, where no other individuals are exposed/infectous


# P(isolated patient in stage l when first HQ exposure occurs | outbreak occurs)
{iso_patient_prob <- function(params){ 
 
  K <- length(params[[3]]$Hyp[[1]]) - 1
  Hyp <- params[[3]]$Hyp[[1]][2:(K+1)]          
  Beta_O <- params[[1]]$Beta_O[[1]][2:(K+1)]      
  lambda <- pops$NHQ * Beta_O
  
  surv <- Hyp / (Hyp + lambda)       # P(progress past stage K without exposure)
  expo <- lambda / (Hyp + lambda)    # P(exposure at stage K | reached stage K)
  
  # cum_surv[l] = prod_{k=l+1}^{K} surv[k]  (empty product = 1 when l=K)
  cum_surv <- rev(cumprod(c(1, rev(surv)))[1:K])
  
  P_first_exposure_at_l <- cum_surv * expo   # indexed l = K, K-1, ..., 1
  names(P_first_exposure_at_l) <- K:1
  
  p_outbreak <- sum(P_first_exposure_at_l)
  probs <- P_first_exposure_at_l / p_outbreak
  
  return(probs)
} }
  
#Uses probability of iso patient being in stage l, samples init state, then simulates gillespie
gillespie_trials_outbreak_analytical <- function(trials, probs){
    runs <- data.frame(replicate = NULL, Y = NULL, YHQ = NULL, YHG = NULL, YP = NULL, T = NULL)
    for (trial in 1:trials){
      
      K <- length(params[[3]]$Hyp[[1]]) - 1
      l <- sample(K:1,1,prob=iso_patient_prob(params = params))
      init <- list(l, c(1,0), c(0,0), c(0,0) )
      
      params[[6]] <- init
    
      runs <- rbind(runs, one_gillespie(params=params,replicate=trial))
      if ((100*trial/trials)%%10 == 0){
        cat("trials: ", trial," out of ",trials," (",100*trial/trials,"%) \n")
      }
    }
    runs
    
} 

sourceCpp("Y_engine_truncated.cpp")

{runs_o <- gillespie_trials_outbreak_analytical(trials, params)

Y <- sum(runs_o$Y)/trials +1
YHQ <- sum(runs_o$YHQ)/trials +1
YP <- sum(runs_o$YP)/trials

if (pops$NHG >0){
  YHG <- sum(runs_o$YHG)/trials
  
  cat("On average, ", Y," exposures occured with, E(YHQ)=",YHQ,", E(YHG,)=",YHG,", E(YP)=",YP)
} else {
  cat("On average, ", Y," exposures occured with, E(YHQ)=",YHQ," and E(YP)=",YP)
}}


# Y_distribution

runsY <- runs_o$Y


Y_dist <- c()
nmaxmax <- 20
for (n in 0:nmaxmax){
  Y_dist <- c(Y_dist,print(length(runsY[runsY==n])/trials))
  YHQ_dist <- c(YHQ_dist,print(length(runsYHQ[runsYHQ==n])/trials))
}
sum(Y_dist)
# fast_trial_cpp(int nmaxmax, double gamma, double rC, String cohort, double beta, double zeta, double eps)

Y_analyt <- fast_trial_cpp(nmaxmax = 20, gamma = gamma, rC = rC, cohort = "S", beta = beta, zeta = zeta, eps = eps)
Y_analyt <- Y_analyt[2:length(Y_analyt)]/(1-Y_analyt[1])
Y_analyt
