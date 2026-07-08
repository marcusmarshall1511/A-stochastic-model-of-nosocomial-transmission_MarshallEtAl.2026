library(expm)
library(readxl)
library(dplyr)

# set.seed(100)

# output_file <- "IRISH_ICU_results_ordered_NEWNEWNEWNEW.txt"
# writeLines("", output_file)

data <- read_excel("COVID-19_ICU_cases.xlsx")
raw_data <- data$`Before april 2021` %>%
  na.omit %>%
  as.vector()

raw_data

# Hypoexponential density


dhypoexp <- function(x, par) {
  alpha <- c(1, numeric(length(par) - 1))
  Theta <- matrix(0, nrow = length(par), ncol = length(par))
  for (i in 1:length(par)) {
    Theta[i, i] = -par[i]
    if (i < length(par)) {
      Theta[i, i + 1] = par[i]
    }
  }
  ones <- rep(1, length(par))
  as.numeric(-alpha %*% expm(x * Theta) %*% Theta %*% ones)
}

############################################################
# Random generator and minimum-of-two
############################################################

rhypoexp <- function(n, par){
  samples <- numeric(n)
  for (i in 1:n){
    samples[i] <- sum(rexp(length(par), rate = par))
  }
  samples
}

minimum_of_two <- function(n, par){
  pmin(
    rhypoexp(n, par),
    rgamma(n, shape = 5, scale = 1.4)
  )
}


############################################################
# Negative log-likelihood 
############################################################
counter <- 0

nll_hypoexp <- function(par, data) {
  
  
  
  vals <- sapply(data, function(x) dhypoexp(x, par))
  
  
  vals[vals <= 0] <- 1e-10
  
  nll <- -sum(log(vals))

  
  nll
  
}

############################################################
# Basic optim fit (single run)
############################################################
est_pars <- numeric(0)
NLLs <- numeric(0)
AICs <- numeric(0)


df <- tibble(
  k = integer(),
  nll = numeric(),
  rates = list()
)


kmax<-15
for (k in 1:kmax){
mean_x <- mean(raw_data)

lower <- rep(k / mean_x * 0.01,k)
upper <- rep(k / mean_x * 50, k)

init  <- runif(k, lower, upper)           # random start





res <-numeric(0)
library(parallel)

cl <- makeCluster(10)

clusterExport(cl, c("nll_hypoexp", "raw_data", "lower", "upper", "k", "dhypoexp"))
clusterEvalQ(cl, library(expm))

run_single <- function(init) {

  fit <- optim(
    par     = init,
    fn      = nll_hypoexp,
    data    = raw_data,
    method  = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = list(
      trace = 0,
      maxit = 5000,
      pgtol = 0
    )
  )
  list(value = fit$value, par = fit$par)
}

starts <- replicate(50, runif(k, lower, upper), simplify = FALSE)

results <- parLapply(cl, starts, run_single)

stopCluster(cl)

best <- results[[ which.min(sapply(results, `[[`, "value")) ]]

best$value
best$par



  rates <- list(best$par)

df_res <- tibble(
  k = k,
  nll = best$value,
  rates = rates )   # ALWAYS wrap in list()


df <- bind_rows(df, df_res)


logfile <- "optim_log_free_params.txt"

cat("starting run\n", file=logfile, append=TRUE)

cat(
  "finished ", k, " out of ", kmax,
  ". nll: ", best$value,
  " rates: ", paste(unlist(rates), collapse=", "),
  "\n",
  file = logfile,
  append = TRUE
)
}


################################################################################

true_par <- best$par
  # c(0.0703227823411154, 8.3198836017103, 8.30400755726981, 8.64881497333471, 8.72619937561971, 8.7971877006993, 10.3122997921632, 9.51076894142005)

raw_data <- minimum_of_two(n = 1000, true_par)



nll_hypoexp <- function(par, data) {
  
  # enforce positivity (safety)
  if (any(par <= 0)) return(1e10)
  
  
  
  vals <- sapply(data, function(x) dhypoexp(x, par))
  vals[vals <= 0] <- 1e-10
  
  nll <- -sum(log(vals))
  
  
  nll
}

############################################################
# Basic optim fit
############################################################


df <- tibble(
  k = integer(),
  nll = numeric(),
  rates = list()
)


kmax<-10
for (k in 1:kmax){
  mean_x <- mean(raw_data)
  
  lower <- rep(k / mean_x * 0.01, k)
  upper <- rep(k / mean_x * 50,  k)
  
  init  <- runif(k, lower, upper)           # random start
  
  
  
  
  
  res <-numeric(0)
  library(parallel)
  
  cl <- makeCluster(10)
  
  clusterExport(cl, c("nll_hypoexp", "raw_data", "lower", "upper", "k", "dhypoexp"))
  clusterEvalQ(cl, library(expm))
  
  run_single <- function(init) {
    
    fit <- optim(
      par     = init,
      fn      = nll_hypoexp,
      data    = raw_data,
      method  = "L-BFGS-B",
      lower   = lower,
      upper   = upper,
      control = list(
        trace = 0,
        maxit = 5000,
        pgtol = 0
      )
    )
    list(value = fit$value, par = fit$par)
  }
  
  starts <- replicate(50, runif(k, lower, upper), simplify = FALSE)
  
  results <- parLapply(cl, starts, run_single)
  
  stopCluster(cl)
  
  best <- results[[ which.min(sapply(results, `[[`, "value")) ]]
  
  best$value
  best$par
  
  
  
  rates <- list(best$par)
  df_res <- tibble(
    k = k,
    nll = best$value,
    rates = rates )   # ALWAYS wrap in list()
  
  
  df <- bind_rows(df, df_res)
  
  
  logfile <- "optim_log_motwo_free_params.txt"

  cat("starting run\n", file=logfile, append=TRUE)

  cat(
    "finished ", k, " out of ", kmax,
    ". nll: ", best$value,
    " rates: ", paste(unlist(rates), collapse=", "),
    "\n",
    file = logfile,
    append = TRUE
  )
}
print(df)