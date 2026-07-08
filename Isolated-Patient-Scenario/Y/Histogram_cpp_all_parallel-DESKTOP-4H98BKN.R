library(parallel)
library(ggplot2)
library(scales)
library(dplyr)
library(Rcpp)
setwd("C:/Users/Marcus/OneDrive - University of Leeds/Project1_Compartmental_model/Rcode_NEWASOFAPRIL/Y_R_rC/c++/extendedp")
writeLines("", "logs/progress.txt")

HQ <- function(nmax, gamma, rC, beta, zeta, eps){
message("Compiling C++ simulation engine...")
sourceCpp("fixed_HQ.cpp")
message("Compilation successful. Starting simulations...")

dir.create("logs", showWarnings = FALSE)
progress_file <- "logs/progress.txt"
if (!file.exists(progress_file)) {
  writeLines("", progress_file) 
}
# 2. Wrapper function to run all three cohorts
run_cohorts <- function(n, gamma, rC, beta, zeta, eps) {
  S <- fast_trial_cpp(2, gamma, rC, "S", beta, zeta, eps)
  #M <- fast_trial_cpp(6, gamma, rC, "M", beta, zeta, eps)
  N <- fast_trial_cpp(13, gamma, rC, "N", beta, zeta, eps)
  S <- S[2:length(S)]/(1-S[1])
  #M <- M[2:length(M)]/M[1]
  N <- N[2:length(N)]/(1-N[1])
  df <- rbind(
    data.frame(x = 1:(length(S)), y = S, category = "S"),
    #data.frame(x = 1:(length(M)), y = M, category = "M"),
    data.frame(x = 1:(length(N)), y = N, category = "N")
  )
  return(df)
}


df <- run_cohorts(nmax, gamma, rC, beta, zeta, eps)



df$category <- factor(df$category, 
                       levels = c("S","M","N"))

return(df)
}

HG <- function(nmax, gamma, rC, beta, zeta, eps){
  message("Compiling C++ simulation engine...")
  sourceCpp("fixed_HG_conditioned_on_outbreak.cpp")
  message("Compilation successful. Starting simulations...")
  
  dir.create("logs", showWarnings = FALSE)
  progress_file <- "logs/progress.txt"
  if (!file.exists(progress_file)) {
    writeLines("", progress_file) 
  }
  # 2. Wrapper function to run all three cohorts
  run_cohorts <- function(n, gamma, rC, beta, zeta, eps) {
    S <- fast_trial_cpp(11, gamma, rC, "S", beta, zeta, eps)
    #M <- fast_trial_cpp(7, gamma, rC, "M", beta, zeta, eps)
    N <- NA
    
    df <- rbind(
      data.frame(x = 0:(length(S)-1), y = S, category = "S"),
      #data.frame(x = 0:(length(M)-1), y = M, category = "M"),
      data.frame(x = 0:(length(N)-1), y = N, category = "N")
    )
    return(df)
  }
  
  df <- run_cohorts(nmax, gamma, rC, beta, zeta, eps)
  
  df$category <- factor(df$category, 
                        levels = c("S","M","N"))
  return(df)
}

P <- function(nmax, gamma, rC, beta, zeta, eps){
  message("Compiling C++ simulation engine...")
  sourceCpp("fixed_P_conditioned_on_outbreak.cpp")
  message("Compilation successful. Starting simulations...")
  
  dir.create("logs", showWarnings = FALSE)
  progress_file <- "logs/progress.txt"
  if (!file.exists(progress_file)) {
    writeLines("", progress_file) 
  }
  # 2. Wrapper function to run all three cohorts
  run_cohorts <- function(n, gamma, rC, beta, zeta, eps) {
    S <- fast_trial_cpp(n, gamma, rC, "S", beta, zeta, eps)
    #M <- fast_trial_cpp(n, gamma, rC, "M", beta, zeta, eps)
    N <- fast_trial_cpp(n, gamma, rC, "N", beta, zeta, eps)
    
    df <- rbind(
      data.frame(x = 0:(length(S)-1), y = S, category = "S"),
      #data.frame(x = 0:(length(M)-1), y = M, category = "M"),
      data.frame(x = 0:(length(N)-1), y = N, category = "N")
    )
    return(df)
  }
  
  df <- run_cohorts(nmax, gamma, rC, beta, zeta, eps)
  
  df$category <- factor(df$category, 
                        levels = c("S","M","N"))
  return(df)
}



run_one <- function(R, rC, eps,gamma) {
  message(paste("Running simulations for R =", R, ", rC =", rC, ", eps =", eps,", gamma =", gamma))

  nmax <- 20
  #gamma <- 1/3
  beta <- R/84
  zeta <- 1/3

  HQ_df <- HQ(nmax, gamma, rC, beta, zeta, eps)
  HQ_df$compartment <- "HQ"
  write(paste(Sys.time(), "HQ" , "- DONE R =", R, ", rC =", rC, ", eps =", eps,", gamma =", gamma, ", max_y = ", max(HQ_df$y, na.rm = TRUE)), 
        file = "logs/progress.txt", append = TRUE)
  HG_df <- HG(nmax, gamma, rC, beta, zeta, eps)
  HG_df$compartment <- "HG"
  write(paste(Sys.time(), "HG" , "- DONE R =", R, ", rC =", rC, ", eps =", eps,", gamma =", gamma, ", max_y = ", max(HG_df$y, na.rm = TRUE)), 
        file = "logs/progress.txt", append = TRUE)

  P_df <- P(nmax, gamma, rC, beta, zeta, eps)
  P_df$compartment <- "P"
  write(paste(Sys.time(), "P" , "- DONE R =", R, ", rC =", rC, ", eps =", eps,", gamma =", gamma, ", max_y = ", max(P_df$y, na.rm = TRUE)), 
        file = "logs/progress.txt", append = TRUE)

  all_df <- rbind(HQ_df, HG_df, P_df)

  #ylimit <- max(all_df$y[!is.na((all_df$y))])+0.1
  ylimit <- 0.9
  # Save plots
  Plotter <- function(cohorttype) {
  df <- all_df[all_df$category==cohorttype,]
  df$compartment <- factor(df$compartment, 
                        levels = c("HQ","P","HG"))
  
  legend_breaks <- if (cohorttype == "N") c("HQ", "P") else c("HQ", "P", "HG")
  
  bar_width <- .7
  p <- ggplot(data=df, aes(x=x, y=y, fill=compartment))+
    theme_bw()+
    geom_col(
      width = bar_width, 
      position = position_dodge(width = bar_width, preserve = "single"), 
      color = "black"
    )+
    scale_fill_discrete(breaks = legend_breaks) +
    scale_x_continuous(breaks=seq(0,15,2),limits = c(NA, 13))+
    scale_y_continuous(labels = label_number(accuracy = 0.01))+
    labs(x=expression(paste("Number of exposures during outbreak")), y="Probability"
         # , title = bquote(R == .(R))
    )+
    theme(legend.position = c(.7,0.8), text = element_text(size = 20), legend.title=element_blank(),
          # Centered and 'plain' to match your axis titles
          plot.title = element_text(hjust = 0.5, face = "plain",size=25),
          # This helps expressions align properly in the legend box
          legend.text.align = 0)+
    # xlim(0,15)+
    ylim(0,ylimit)
  return(p)
  ggsave(paste0("Histogram_", cohorttype,"_rC",rC,"_R",R,"_eps",eps,"_gamma",gamma,".pdf"), plot=p, width= 5, height= 5, unit="in", dpi=300, path="Histograms/")
}

  # At the bottom of run_one...
plot_S <- Plotter("S")
plot_N <- Plotter("N")

# Return everything as a named list
return(list(
  data = all_df,
  plot_S = plot_S,
  plot_N = plot_N
))

  
}




params <- do.call(
  rbind,
  lapply(c(0.9), function(rC) {
    expand.grid(
      R     = c(10),
      rC    = rC,
      eps   = c(1 - rC, 1),
      gamma = c(1/1, 1/5)
    )
  })
)


library(pbapply)
library(parallel)

cl <- makeCluster(detectCores() - 1)

clusterEvalQ(cl, {
  library(Rcpp)
  library(ggplot2)
  library(dplyr)
  library(scales)
  setwd("C:/Users/Marcus/OneDrive - University of Leeds/Project1_Compartmental_model/Rcode_NEWASOFAPRIL/Y_R_rC/c++/extendedp")
})

clusterExport(cl, c("HQ", "HG", "P", "run_one", "params"))


results <- parLapply(
  cl,
  seq_len(nrow(params)),
  function(i) {
    p <- params[i, ]
    run_one(p$R, p$rC, p$eps, p$gamma)
  }
)


stopCluster(cl)

saveRDS(results, "simulation_plots_and_data.rds")








