library(parallel)
library(ggplot2)
library(scales)
library(dplyr)
library(Rcpp)
library(RcppArmadillo)

library(patchwork)
library(cowplot)
library(gridExtra)


library(rstudioapi)

current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))

writeLines("", "logs/progress.txt")

HQ <- function(nmax, gamma, rC, beta, zeta, eps){
message("Compiling C++ simulation engine...")
sourceCpp("YHQ_engine.cpp")
message("Compilation successful. Starting simulations...")

dir.create("logs", showWarnings = FALSE)
progress_file <- "logs/progress.txt"
if (!file.exists(progress_file)) {
  writeLines("", progress_file) 
}
# 2. Wrapper function to run all three cohorts
run_cohorts <- function(n, gamma, rC, beta, zeta, eps) {
  S <- fast_trial_cpp_HQ(2, gamma, rC, "S", beta, zeta, eps)
  #M <- fast_trial_cpp_HQ(6, gamma, rC, "M", beta, zeta, eps)
  N <- fast_trial_cpp_HQ(13, gamma, rC, "N", beta, zeta, eps)
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
  sourceCpp("YHG_engine.cpp")
  message("Compilation successful. Starting simulations...")
  
  dir.create("logs", showWarnings = FALSE)
  progress_file <- "logs/progress.txt"
  if (!file.exists(progress_file)) {
    writeLines("", progress_file) 
  }
  # 2. Wrapper function to run all three cohorts
  run_cohorts <- function(n, gamma, rC, beta, zeta, eps) {
    S <- fast_trial_cpp_HG(11, gamma, rC, "S", beta, zeta, eps)
    #M <- fast_trial_cpp_HG(7, gamma, rC, "M", beta, zeta, eps)
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
  sourceCpp("YP_engine.cpp")
  message("Compilation successful. Starting simulations...")
  
  dir.create("logs", showWarnings = FALSE)
  progress_file <- "logs/progress.txt"
  if (!file.exists(progress_file)) {
    writeLines("", progress_file) 
  }
  # 2. Wrapper function to run all three cohorts
  run_cohorts <- function(n, gamma, rC, beta, zeta, eps) {
    S <- fast_trial_cpp_P(n, gamma, rC, "S", beta, zeta, eps)
    #M <- fast_trial_cpp_P(n, gamma, rC, "M", beta, zeta, eps)
    N <- fast_trial_cpp_P(n, gamma, rC, "N", beta, zeta, eps)
    
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

  nmax <- 15
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
  # current_path = rstudioapi::getActiveDocumentContext()$path 
  # setwd(dirname(current_path ))
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



####################################################################################################################################




yshift_ <- 0.1


p1 <- results[[1]]$plot_N+
  theme(axis.title.x = element_blank()) +
  scale_x_continuous(labels=NULL, breaks = seq(0,20,2), limits = c(NA,13))+
  ylim(0,1)+
  theme(legend.position = c(0.85,0.85))
show(p1)

p2 <- results[[3]]$plot_N+
  theme(axis.title.x = element_blank()) +
  scale_x_continuous(labels=NULL, breaks = seq(0,20,2), limits = c(NA,13))+
  theme(axis.title.y = element_blank(), legend.position = c(0.85,0.85)) +
  scale_y_continuous(labels=NULL, breaks = seq(0,10,0.25),limits=c(0,1))
show(p2)



p3 <- results[[1]]$plot_S+
  theme(axis.title.x = element_blank(), legend.position = c(0.85,0.85)) +
  scale_x_continuous(labels=NULL, breaks = seq(0,13,2), limits = c(NA,13))+
  ylim(0,1)
show(p3)
xlim_zoom <- c(2.5, 13)
ylim_zoom <- c(0, 0.003)
zoom_plot <- p3 +
  xlim(xlim_zoom) +
  ylim(ylim_zoom) +
  theme_bw()+
  theme(legend.position ="none")+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),plot.margin = margin(b=5, 0, 0, 0),
        panel.spacing = unit(0, "lines"))+
  scale_x_continuous(breaks=seq(4,14,2), limits=xlim_zoom)
zoom_plot
yshift <- 0.05
p3 <- p3+
  annotation_custom(
    grob = ggplotGrob(zoom_plot), 
    xmin = 2.5, xmax = 13, 
    ymin = 0.08+ yshift, ymax = 0.6+ yshift
  ) +
  annotate("rect", xmin = xlim_zoom[1], xmax = xlim_zoom[2], ymin = ylim_zoom[1], ymax = ylim_zoom[2], 
           color = "black", fill = NA, linetype = "dashed")
show(p3)



p4 <- results[[3]]$plot_S+
  theme(axis.title.x = element_blank(), legend.position = c(0.85,0.85)) +
  scale_x_continuous(labels=NULL, breaks = seq(0,13,2), limits = c(NA,13))+
  theme(axis.title.y = element_blank()) +
  scale_y_continuous(labels=NULL, breaks = seq(0,1,0.25),limits=c(0,1))
show(p4)
xlim_zoom <- c(2.5, 13)
ylim_zoom <- c(0, 0.05)
zoom_plot <- p4 +
  xlim(xlim_zoom) +
  ylim(ylim_zoom) +
  theme_bw()+
  theme(legend.position ="none")+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),plot.margin = margin(b=5, 0, 0, 0),
        panel.spacing = unit(0, "lines"))+
  scale_x_continuous(breaks=seq(4,14,2), limits=xlim_zoom)
zoom_plot

p4 <- p4+
  annotation_custom(
    grob = ggplotGrob(zoom_plot), 
    xmin = 2.5, xmax = 13, 
    ymin = 0.08+ yshift, ymax = 0.6+ yshift
  ) +
  annotate("rect", xmin = xlim_zoom[1], xmax = xlim_zoom[2], ymin = ylim_zoom[1], ymax = ylim_zoom[2], 
           color = "black", fill = NA, linetype = "dashed")
show(p4)


p5 <- results[[2]]$plot_S+theme(axis.title.x = element_blank(), legend.position = c(0.85,0.85))+
  ylim(0,1)
show(p5)
xlim_zoom <- c(2.5, 13)
ylim_zoom <- c(0, 0.06)
zoom_plot <- p5 +
  xlim(xlim_zoom) +
  ylim(ylim_zoom) +
  theme_bw()+
  theme(legend.position ="none")+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),plot.margin = margin(b=5, 0, 0, 0),
        panel.spacing = unit(0, "lines"))+
  scale_x_continuous(breaks=seq(4,14,2), limits=xlim_zoom)
zoom_plot

p5 <- p5+
  annotation_custom(
    grob = ggplotGrob(zoom_plot), 
    xmin = 2.5, xmax = 13, 
    ymin = 0.08+ yshift, ymax = 0.6+ yshift
  ) +
  annotate("rect", xmin = xlim_zoom[1], xmax = xlim_zoom[2], ymin = ylim_zoom[1], ymax = ylim_zoom[2], 
           color = "black", fill = NA, linetype = "dashed")
show(p5)

p6 <- results[[4]]$plot_S+
  theme(axis.title.y = element_blank(), legend.position = c(0.85,0.85)) +
  scale_y_continuous(labels=NULL, breaks = seq(0,1,0.25),limits=c(0,1))+theme(axis.title.x = element_blank())
show(p6)
xlim_zoom <- c(2.5, 13)
ylim_zoom <- c(0, 0.13)
zoom_plot <- p6 +
  xlim(xlim_zoom) +
  ylim(ylim_zoom) +
  theme_bw()+
  theme(legend.position ="none")+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),plot.margin = margin(b=5, 0, 0, 0),
        panel.spacing = unit(0, "lines"))+
  scale_x_continuous(breaks=seq(4,14,2), limits=xlim_zoom)
zoom_plot


p6 <- p6+
  annotation_custom(
    grob = ggplotGrob(zoom_plot), 
    xmin = 2.5, xmax = 13, 
    ymin = 0.08+ yshift, ymax = 0.6+ yshift
  ) +
  annotate("rect", xmin = xlim_zoom[1], xmax = xlim_zoom[2], ymin = ylim_zoom[1], ymax = ylim_zoom[2], 
           color = "black", fill = NA, linetype = "dashed")
show(p6)


expected_value <- function(df) {
  # Identify the columns
  xcol <- df[[1]] 
  ycol <- df[[2]]          
  
  
  # Expected value
  sum(xcol * ycol)
}



d1HQ <- results[[1]]$plot_N$data[ results[[1]]$plot_N$data$compartment=="HQ",]
d1P <- results[[1]]$plot_N$data[ results[[1]]$plot_N$data$compartment=="P",]

d2HQ <- results[[3]]$plot_N$data[ results[[3]]$plot_N$data$compartment=="HQ",]
d2P<- results[[3]]$plot_N$data[ results[[3]]$plot_N$data$compartment=="P",]

d3HQ <- results[[1]]$plot_S$data[ results[[1]]$plot_S$data$compartment=="HQ",]
d3P<- results[[1]]$plot_S$data[ results[[1]]$plot_S$data$compartment=="P",]
d3HG <- results[[1]]$plot_S$data[ results[[1]]$plot_S$data$compartment=="HG",]

d4HQ <- results[[3]]$plot_S$data[ results[[3]]$plot_S$data$compartment=="HQ",]
d4P<- results[[3]]$plot_S$data[ results[[3]]$plot_S$data$compartment=="P",]
d4HG <- results[[3]]$plot_S$data[ results[[3]]$plot_S$data$compartment=="HG",]

d5HQ <- results[[2]]$plot_S$data[ results[[2]]$plot_S$data$compartment=="HQ",]
d5P<- results[[2]]$plot_S$data[ results[[2]]$plot_S$data$compartment=="P",]
d5HG <- results[[2]]$plot_S$data[ results[[2]]$plot_S$data$compartment=="HG",]

d6HQ <- results[[4]]$plot_S$data[ results[[4]]$plot_S$data$compartment=="HQ",]
d6P<- results[[4]]$plot_S$data[ results[[4]]$plot_S$data$compartment=="P",]
d6HG <- results[[4]]$plot_S$data[ results[[4]]$plot_S$data$compartment=="HG",]

m1HQ <- expected_value(d1HQ)
m1P <- expected_value(d1P)

m2HQ <- expected_value(d2HQ)
m2P <- expected_value(d2P)

m3HQ <- expected_value(d3HQ)
m3P <- expected_value(d3P)
m3HG <- expected_value(d3HG)

m4HQ <- expected_value(d4HQ)
m4P <- expected_value(d4P)
m4HG <- expected_value(d4HG)

m5HQ <- expected_value(d5HQ)
m5P <- expected_value(d5P)
m5HG <- expected_value(d5HG)

m6HQ <- expected_value(d6HQ)
m6P <- expected_value(d6P)
m6HG <- expected_value(d6HG)


xshift <- 3

p1 <- p1+ theme(legend.position = c(0.85-0.4,0.85))+
  geom_text(
    aes(
      x = 5-0.1+xshift, y = 0.825+yshift_,
      label = paste0("E*group('[', Y[H[Q]]*'|'*Y>0, ']') == ", sprintf("%.2f", m1HQ))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )+
  geom_text(
    aes(
      x = 5-0.1+xshift, y = 0.825-0.065+yshift_,
      label = paste0("E*group('[', Y[P]*phantom('o')*'|'*Y>0, ']')== ", sprintf("%.2f", m1P))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )+
  geom_segment(
    aes(x = 7.5, xend = 13, y = 0.825 - 0.097+yshift_, yend = 0.825 - 0.097+yshift_),
    inherit.aes = FALSE,
    linetype = "dashed",  # or "solid"
    colour = "grey40"
  ) +
  geom_text(
    aes(
      x = 5-0.1+xshift, y = 0.825-0.13-0.01+yshift_,
      label = paste0("E*group('[',Y*phantom('oii')*'|'*Y>0, ']')  == ", sprintf("%.2f", m1HQ+m1P))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )

p2 <- p2+ theme(legend.position = c(0.85-0.4,0.85))+
  geom_text(
    aes(
      x = 5-0.1+xshift, y = 0.825+yshift_,
      label = paste0("E*group('[', Y[H[Q]]*'|'*Y>0, ']') == ", sprintf("%.2f", m2HQ))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )+
  geom_text(
    aes(
      x = 5-0.1+xshift, y = 0.825-0.065+yshift_,
      label = paste0("E*group('[', Y[P]*phantom('o')*'|'*Y>0, ']')== ", sprintf("%.2f", m2P))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )+
  geom_segment(
    aes(x = 7.5, xend =13, y = 0.825 - 0.097+yshift_, yend = 0.825 - 0.097+yshift_),
    inherit.aes = FALSE,
    linetype = "dashed",  # or "solid"
    colour = "grey40"
  ) +
  geom_text(
    aes(
      x = 5-0.1+xshift, y = 0.825-0.13-0.01+yshift_,
      label = paste0("E*group('[',Y*phantom('oii')*'|'*Y>0, ']')  == ", sprintf("%.2f", m2HQ+m2P))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )



p3 <- p3+ theme(legend.position = c(0.85-0.4,0.85))+
  geom_text(
    aes(
      x = 4.5+xshift, y = 0.86+yshift_,
      label = paste0("E*group('[', Y[H[Q]]*'|'*Y>0, ']') == ", sprintf("%.2f", m3HQ))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )+
  geom_text(
    aes(
      x = 4.5+xshift, y = 0.86-0.065+yshift_,
      label = paste0("E*group('[', Y[P]*phantom('o')*'|'*Y>0, ']')== ", sprintf("%.2f", m3P))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )+
  geom_text(
    aes(
      x = 4.5+xshift, y = 0.86-0.13+yshift_,
      label = paste0("E*group('[', Y[HG]*'|'*Y>0, ']')== ", sprintf("%.2f", m3HG))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )+
  geom_segment(
    aes(x = 7.5, xend = 13, y = 0.86 - 0.097-0.065+yshift_, yend = 0.86 - 0.097-0.065+yshift_),
    inherit.aes = FALSE,
    linetype = "dashed",  # or "solid"
    colour = "grey40"
  ) +
  geom_text(
    aes(
      x = 4.5+xshift, y = 0.86-0.13-0.065-0.01+yshift_,
      label = paste0("E*group('[',Y*phantom('oii')*'|'*Y>0, ']')  == ", sprintf("%.2f", m3HQ+m3P+m3HG))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )



p4 <- p4 + theme(legend.position = c(0.85-0.4,0.85))+
  geom_text(
    aes(
      x = 4.5+xshift, y = 0.86+yshift_,
      label = paste0("E*group('[', Y[H[Q]]*'|'*Y>0, ']') == ", sprintf("%.2f", m4HQ))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )+
  geom_text(
    aes(
      x = 4.5+xshift, y = 0.86-0.065+yshift_,
      label = paste0("E*group('[', Y[P]*phantom('o')*'|'*Y>0, ']')== ", sprintf("%.2f", m4P))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )+
  geom_text(
    aes(
      x = 4.5+xshift, y = 0.86-0.13+yshift_,
      label = paste0("E*group('[', Y[HG]*'|'*Y>0, ']')== ", sprintf("%.2f", m4HG))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )+
  geom_segment(
    aes(x = 7.5, xend = 13, y = 0.86 - 0.097-0.065+yshift_, yend = 0.86 - 0.097-0.065+yshift_),
    inherit.aes = FALSE,
    linetype = "dashed",  # or "solid"
    colour = "grey40"
  ) +
  geom_text(
    aes(
      x = 4.5+xshift, y = 0.86-0.13-0.065-0.01+yshift_,
      label = paste0("E*group('[',Y*phantom('oii')*'|'*Y>0, ']')  == ", sprintf("%.2f", m4HQ+m4P+m4HG))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )


p5 <- p5 + theme(legend.position = c(0.85-0.4,0.85))+
  geom_text(
    aes(
      x = 4.5+xshift, y = 0.86+yshift_,
      label = paste0("E*group('[', Y[H[Q]]*'|'*Y>0, ']') == ", sprintf("%.2f", m5HQ))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )+
  geom_text(
    aes(
      x = 4.5+xshift, y = 0.86-0.065+yshift_,
      label = paste0("E*group('[', Y[P]*phantom('o')*'|'*Y>0, ']')== ", sprintf("%.2f", m5P))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )+
  geom_text(
    aes(
      x = 4.5+xshift, y = 0.86-0.13+yshift_,
      label = paste0("E*group('[', Y[HG]*'|'*Y>0, ']')== ", sprintf("%.2f", m5HG))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )+
  geom_segment(
    aes(x = 7.5, xend = 13, y = 0.86 - 0.097-0.065+yshift_, yend = 0.86 - 0.097-0.065+yshift_),
    inherit.aes = FALSE,
    linetype = "dashed",  # or "solid"
    colour = "grey40"
  ) +
  geom_text(
    aes(
      x = 4.5+xshift, y = 0.86-0.13-0.065-0.01+yshift_,
      label = paste0("E*group('[',Y*phantom('oii')*'|'*Y>0, ']')  == ", sprintf("%.2f", m5HQ+m5P+m5HG))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )



p6 <- p6 + theme(legend.position = c(0.85-0.4,0.85))+
  geom_text(
    aes(
      x = 4.5+xshift, y = 0.86+yshift_,
      label = paste0("E*group('[', Y[H[Q]]*'|'*Y>0, ']') == ", sprintf("%.2f", m6HQ))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )+
  geom_text(
    aes(
      x = 4.5+xshift, y = 0.86-0.065+yshift_,
      label = paste0("E*group('[', Y[P]*phantom('o')*'|'*Y>0, ']')== ", sprintf("%.2f", m6P))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )+
  geom_text(
    aes(
      x = 4.5+xshift, y = 0.86-0.13+yshift_,
      label = paste0("E*group('[', Y[HG]*'|'*Y>0, ']')== ", sprintf("%.2f", m6HG))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )+
  geom_segment(
    aes(x = 7.5, xend = 13, y = 0.86 - 0.097-0.065+yshift_, yend = 0.86 - 0.097-0.065+yshift_),
    inherit.aes = FALSE,
    linetype = "dashed",  # or "solid"
    colour = "grey40"
  ) +
  geom_text(
    aes(
      x = 4.5+xshift, y = 0.86-0.13-0.065-0.01+yshift_,
      label = paste0("E*group('[',Y*phantom('oii')*'|'*Y>0, ']')  == ", sprintf("%.2f", m6HQ+m6P+m6HG))),
    parse = TRUE,
    hjust=0,
    inherit.aes = FALSE
  )


combined_plot <- (p1 | p2) / (p3 | p4) / (p5 | p6)

R <- combined_plot + plot_annotation(
  theme = theme(
    plot.margin = margin(b = 20, l = 50,t=25,r=15) # 'b' for bottom space, 'l' for left space
  )
)

R_lab <- ggdraw(R) +
  draw_label("Fast detection", x = 0.275, y = .99, angle = 0, size = 18, fontface="bold") +
  draw_label(
    expression(paste("(", 1/gamma == 1, " day)")),
    x = 0.4625, y = .99, angle = 0, size = 18, fontface="bold"
  ) +
  draw_label("Slow detection", x = 0.275 + 0.4125, y = .99, angle = 0, size = 18, fontface="bold") +
  draw_label(
    expression(paste("(", 1/gamma == 5, " days)")),
    x = 0.4625 + 0.4125, y = .99, angle = 0, size = 18, fontface="bold"
  ) +
  draw_label("No cohorting", x = 0.03, y = .825, angle = 90, size = 20, fontface="bold") +
  draw_label("Strong cohorting", x = 0.03, y = .4625, angle = 90, size = 20, fontface="bold") +
  draw_label(
    expression(paste("(", epsilon == 1 - r[C], ")")),
    x = 0.03, y = .6125, angle = 90, size = 20, fontface="bold"
  ) +
  draw_label("Strong cohorting", x = 0.03, y = .15 + 0.02, angle = 90, size = 20, fontface="bold") +
  draw_label(
    expression(paste("(", epsilon == 1, ")")),
    x = 0.03, y = .285 + 0.02, angle = 90, size = 20, fontface="bold"
  ) +
  draw_label(expression(n),
             x = 0.565, y = 0.01, angle = 0, size = 25, fontface="plain")


print(R_lab)







