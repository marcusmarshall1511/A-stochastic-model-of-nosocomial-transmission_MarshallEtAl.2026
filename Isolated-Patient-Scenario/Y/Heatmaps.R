library(parallel)
library(ggplot2)
library(dplyr)
library(Rcpp)
library(rstudioapi)


library(gridExtra)
library(metR)
library(cowplot)
library(latex2exp)



current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))

# 1. Compile the C++ engine
message("Compiling C++ simulation engine...")
sourceCpp("Y_engine_truncated.cpp")
message("Compilation successful. Starting simulations...")


progress_file <- "logs/progress.txt"

if (!file.exists(progress_file)) {
    writeLines("", progress_file) 
}
# 2. Wrapper function to run all three cohorts
run_cohorts <- function(n, gamma, rC, beta, zeta, eps) {
  S <- fast_trial_cpp(n, gamma, rC, "S", beta, zeta, eps)
  M <- fast_trial_cpp(n, gamma, rC, "M", beta, zeta, eps)
  N <- fast_trial_cpp(n, gamma, rC, "N", beta, zeta, eps)
  
  df <- rbind(
    data.frame(x = 0:(length(S)-1), y = S, category = "S"),
    data.frame(x = 0:(length(M)-1), y = M, category = "M"),
    data.frame(x = 0:(length(N)-1), y = N, category = "N")
  )
  return(df)
}

### PARAMETERS
n    <- 100
zeta <- 1/3
widthheight <- 2 #resolution of heatmaps

Xpl      <- list()
range_pl <- numeric(0)

# Determine the number of cores to use based on Slurm allocation
# If Slurm doesn't set it, default to 1 to be safe
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
message(sprintf("Running on %d cores...", n_cores))

# 3. The Main Grid Loop
for (rC in c(0.9, 0.5, 0)) {
    eps <- 1-rC
    message(sprintf("Starting rC=%.2f eps=%.2f at %s", rC, eps, Sys.time()))
    
    # Capture loop variables explicitly for workers
    rC_  <- rC
    eps_ <- eps
    
    Rvals     <- seq(1, 10, length.out = widthheight)
    gammavals <- 3 / seq(1, 21, length.out = widthheight)
    grid      <- expand.grid(R = Rvals, gamma = gammavals)
    
    # The parallel loop
    results <- mclapply(1:nrow(grid), function(i) {
      cat(sprintf("[%d/%d] rC=%.2f eps=%.2f R=%.2f gamma=%.4f\n", 
                  i, nrow(grid), rC_, eps_, grid$R[i], grid$gamma[i]),
          file = progress_file, append = TRUE)
      
      R     <- grid$R[i]
      gamma <- grid$gamma[i]
      beta  <- R * 1/(12*7)
      
      # Call the C++ wrapper
      df <- run_cohorts(n = n, gamma = gamma, rC = rC_, beta = beta, zeta = zeta, eps = eps_)
      df$category <- factor(df$category, levels = c("S", "M", "N"))
      
      mean_ <- function(X) {
        X <- X[2:length(X)] / (1 - X[1])
        sum(X * seq_along(X))
      }
      
      list(
        S = mean_(df$y[df$category == "S"]),
        M = mean_(df$y[df$category == "M"]),
        N = mean_(df$y[df$category == "N"])
      )
    }, mc.cores = n_cores)
    
    # Reassemble results
    grid$S <- sapply(results, `[[`, "S")
    grid$M <- sapply(results, `[[`, "M")
    grid$N <- sapply(results, `[[`, "N")
    
    # Plotting
    pltnameS <- paste0("p", rC_, eps_, "S")
    pltnameM <- paste0("p", rC_, eps_, "M")
    pltnameN <- paste0("p", rC_, eps_, "N")
    
    Spl <- ggplot(grid, aes(x = R, y = 1/gamma, fill = S)) +
      geom_tile() + labs(title = pltnameS)
    
    Mpl <- ggplot(grid, aes(x = R, y = 1/gamma, fill = M)) +
      geom_tile() + labs(title = pltnameM)
    
    Npl <- ggplot(grid, aes(x = R, y = 1/gamma, fill = N)) +
      geom_tile() + labs(title = pltnameN)
    
    range_pl <- range(c(range_pl, grid$S, grid$M, grid$N), na.rm = TRUE)
    
    Xpl[[pltnameS]] <- Spl
    Xpl[[pltnameM]] <- Mpl
    Xpl[[pltnameN]] <- Npl
    
    saveRDS(Spl, paste0(pltnameS, ".rds"))
    saveRDS(Mpl, paste0(pltnameM, ".rds"))
    saveRDS(Npl, paste0(pltnameN, ".rds"))
    
    message(sprintf("Done rC=%.2f eps=%.2f at %s", rC_, eps_, Sys.time()))
}

# Final plot formatting
for (i in seq_along(Xpl)) {
    Xpl[[i]] <- Xpl[[i]] + scale_fill_viridis_c(limits = range_pl, option = "turbo")
}

saveRDS(Xpl, "Xpl_list_careful.rds")
message("All (careful) simulations complete!")


###################################################

for (rC in c(0.9, 0.5, 0)) {
  eps <- 1
  message(sprintf("Starting rC=%.2f eps=%.2f at %s", rC, eps, Sys.time()))
  
  # Capture loop variables explicitly for workers
  rC_  <- rC
  eps_ <- eps
  
  
  # The parallel loop
  results <- mclapply(1:nrow(grid), function(i) {
    cat(sprintf("[%d/%d] rC=%.2f eps=%.2f R=%.2f gamma=%.4f\n", 
                i, nrow(grid), rC_, eps_, grid$R[i], grid$gamma[i]),
        file = progress_file, append = TRUE)
    
    R     <- grid$R[i]
    gamma <- grid$gamma[i]
    beta  <- R * 1/(12*7)
    
    # Call the C++ wrapper
    df <- run_cohorts(n = n, gamma = gamma, rC = rC_, beta = beta, zeta = zeta, eps = eps_)
    df$category <- factor(df$category, levels = c("S", "M", "N"))
    
    mean_ <- function(X) {
      X <- X[2:length(X)] / (1 - X[1])
      sum(X * seq_along(X))
    }
    
    list(
      S = mean_(df$y[df$category == "S"]),
      M = mean_(df$y[df$category == "M"]),
      N = mean_(df$y[df$category == "N"])
    )
  }, mc.cores = n_cores)
  
  # Reassemble results
  grid$S <- sapply(results, `[[`, "S")
  grid$M <- sapply(results, `[[`, "M")
  grid$N <- sapply(results, `[[`, "N")
  
  # Plotting
  pltnameS <- paste0("p", rC_, eps_, "S")
  pltnameM <- paste0("p", rC_, eps_, "M")
  pltnameN <- paste0("p", rC_, eps_, "N")
  
  Spl <- ggplot(grid, aes(x = R, y = 1/gamma, fill = S)) +
    geom_tile() + labs(title = pltnameS)
  
  Mpl <- ggplot(grid, aes(x = R, y = 1/gamma, fill = M)) +
    geom_tile() + labs(title = pltnameM)
  
  Npl <- ggplot(grid, aes(x = R, y = 1/gamma, fill = N)) +
    geom_tile() + labs(title = pltnameN)
  
  range_pl <- range(c(range_pl, grid$S, grid$M, grid$N), na.rm = TRUE)
  
  Xpl[[pltnameS]] <- Spl
  Xpl[[pltnameM]] <- Mpl
  Xpl[[pltnameN]] <- Npl
  
  saveRDS(Spl, paste0(pltnameS, ".rds"))
  saveRDS(Mpl, paste0(pltnameM, ".rds"))
  saveRDS(Npl, paste0(pltnameN, ".rds"))
  
  message(sprintf("Done rC=%.2f eps=%.2f at %s", rC_, eps_, Sys.time()))
}

# Final plot formatting
for (i in seq_along(Xpl)) {
  Xpl[[i]] <- Xpl[[i]] + scale_fill_viridis_c(limits = range_pl, option = "turbo")
}

saveRDS(Xpl, "Xpl_list_uncareful.rds")
message("All (uncareful) simulations complete!")



################################################################################





pl0d90_rds <- readRDS(file = "p0.90.1S.rds")
pl0d90uc_rds <- readRDS(file = "p0.91S.rds")
# pl0d75_rds <- readRDS(file = "p0.750.25S.rds")
# pl0d75uc_rds <- readRDS(file = "p0.751S.rds")
pl0d50_rds <- readRDS(file = "p0.50.5S.rds")
pl0d50uc_rds <- readRDS(file = "p0.51S.rds")
# pl0d25_rds <- readRDS(file = "p0.250.75S.rds")
# pl0d25uc_rds <- readRDS(file = "p0.251S.rds")
pl0d0_rds <- readRDS(file = "p01S.rds")

pl0d90 <- pl0d90_rds$data

pl0d90uc <- pl0d90uc_rds$data
pl0d90uc$N <- pl0d90$N

# pl0d75<- pl0d75_rds$data
# 
# pl0d75uc <- pl0d75uc_rds$data
# pl0d75uc$N <- pl0d75$N

pl0d50 <- pl0d50_rds$data

pl0d50uc <- pl0d50uc_rds$data
pl0d50uc$N <- pl0d50$N

# pl0d25 <- pl0d25_rds$data
# 
# pl0d25uc <- pl0d25uc_rds$data
# pl0d25uc$N <- pl0d25$N

pl0d0<- pl0d0_rds$data

pl0d90$Category <- "0d90e0d10"
pl0d90uc$Category <- "0d90e1d00"
# pl0d75$Category <- "0d75e0d25"
# pl0d75uc$Category <- "0d75e1d00"
pl0d50$Category <- "0d50e0d50"
pl0d50uc$Category <- "0d50e1d00"
# pl0d25$Category <- "0d25e0d75"
# pl0d25uc$Category <- "0d25e1d00"
pl0d0$Category <- "0d0e1d00"

pl <- rbind(pl0d90,pl0d90uc,pl0d50,pl0d50uc,pl0d0)

# range_ <- range(pl$S,pl$M,pl$N)
range_ <- range(pl$N)
range_diff <- range((pl$S-pl$N)/pl$N,(pl$M-pl$N)/pl$N)
# range_diff <- c(-1,1)

pl$NS <- (pl$S - pl$N) / pl$N
pl$NM <- (pl$M - pl$N) / pl$N
strs <- c("0d90e0d10","0d90e1d00","0d50e0d50","0d50e1d00","0d0e1d00")

Xpl_local<- list()

# Extract unique rC values for the top row (one N plot per rC)
rc_values <- c(0.90,0.50,0)

# Create the 4 unique N plots for top row
N_plots <- list()

for (rc in rc_values) {
  rc_str <- sprintf("0d%d", rc * 100)
  matching_str <- strs[grep(rc_str, strs)][1]  # Get first match (eps: 1-rC variant)
  
  breaksn <- seq(1, 10, 1)
  
  N_plot <- ggplot(pl[pl$Category == matching_str,], aes(x=R, y=1/gamma, fill=N)) +
    theme_bw()+
    geom_raster(interpolate=TRUE) +
    geom_contour2(aes(z=N), breaks=breaksn) +
    geom_text_contour(aes(z = N),
                      color = "black",
                      size=4.5,
                      fontface='bold',
                      stroke=0.1,
                      rotate=TRUE,
                      label.placer = label_placer_fraction(frac=0.5),
                      breaks=breaksn,
                      skip=0) +
    scale_fill_viridis_c(option="turbo", limits = range(pl$N)) +
    scale_y_continuous(breaks=c(1/3, 1:7),
                       labels=expression(1/3, 1,2,3,4,5,6,7)) +
    labs(title = paste0("rC: ", rc)) +
    theme(text=element_text(size=16),
          legend.position="right",
          # axis.title.x = element_blank(),
          # axis.title.y = element_blank()
    )
  
  N_plots[[as.character(rc)]] <- N_plot
}

# Create NS and NM plots for each of the 8 category strings
NS_NM_plots <- list()

for (str in strs) {
  rC_t <- as.numeric(substr(str,3,4))/100
  eps_t <- 1-rC_t
  
  if (substr(str,6,9) == "1d00"){
    str_title <- paste0("rC: ",rC_t," eps: 1",sep='')
  } else {
    str_title <- paste0("rC: ",rC_t," eps: 1-rC",sep='')
  }
  
  breaks <- c(-0.390,-0.375, round(c(seq(-0.350,-0.05,0.05),seq(0.05,0.350,0.05)),2))
  breakss <- breaks
  breaksm <- breaks
  
  if (str=="0d90e0d10"){
    breaksm <- c(breaks, -0.26)
  } else if (str=="0d50e0d50"){
    breakss <- seq(-0.025,-10,-0.025)
    breaksm <- seq(-0.025,-10,-0.025)
  } else if (str=="0d50e1d00"){
    breakss <- seq(0.025,10,0.025)
    breaksm <- seq(0.025,10,0.025)
  }
  
  # NS plot
  NS_plot <- ggplot(pl[pl$Category==str,], aes(x=R, y=1/gamma, fill=NS)) +
    theme_bw()+
    geom_raster(interpolate=TRUE) +
    geom_contour2(aes(z=NS), breaks=breakss) +
    geom_text_contour(aes(z = NS),
                      color = "black",
                      size=4.5,
                      fontface='bold',
                      stroke=0.1,
                      rotate=TRUE,
                      label.placer = label_placer_fraction(frac=0.5),
                      breaks=breakss,
                      skip=1) +
    scale_fill_gradientn(colours=c("#053061", "#4575b4", "#abd9e9", "white", "#fee090", "#f46d43", "#67001F"),
                         limits=range_diff,
                         values=scales::rescale(c(min(range_diff), -max(abs(range_diff))/2, -max(abs(range_diff))/10, 0, max(abs(range_diff))/10, max(abs(range_diff))/2, max(range_diff))),
                         oob=scales::squish) +
    scale_y_continuous(breaks=c(1/3, 1:7),
                       labels=expression(1/3, 1,2,3,4,5,6,7)) +
    # labs(title=paste0(str_title, " (NS)")) +
    theme(text=element_text(size=16),
          legend.position="none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  # NM plot
  NM_plot <- ggplot(pl[pl$Category==str,], aes(x=R, y=1/gamma, fill=NM)) +
    theme_bw()+
    geom_raster(interpolate=TRUE) +
    geom_contour2(aes(z=NM), breaks=breaksm) +
    geom_text_contour(aes(z = NM),
                      color = "black",
                      size=4.5,
                      fontface='bold',
                      stroke=0.1,
                      rotate=TRUE,
                      label.placer = label_placer_fraction(frac=0.5),
                      breaks=breaksm,
                      skip=1) +
    scale_fill_gradientn(colours=c("#053061", "#4575b4", "#abd9e9", "white", "#fee090", "#f46d43", "#67001F"),
                         limits=range_diff,
                         values=scales::rescale(c(min(range_diff), -max(abs(range_diff))/2, -max(abs(range_diff))/10, 0, max(abs(range_diff))/10, max(abs(range_diff))/2, max(range_diff))),
                         oob=scales::squish) +
    scale_y_continuous(
      breaks=c(1/3, 1:7),
      labels=expression(1/3, 1,2,3,4,5,6,7)
      # labels = NULL
    ) +
    # labs(title=paste0(str_title, " (NM)")) +
    theme(text=element_text(size=16),
          legend.position="none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          # axis.ticks.y = element_blank()
    )
  # LAYER 2: Force the -0.26 label (only if it exists for this string)
  if (-0.26 %in% breaksm) {
    NM_plot <- NM_plot +
      geom_text_contour(aes(z = NM),
                        color = "black", size = 4.5, fontface = 'bold',
                        stroke = 0.1, rotate = TRUE,
                        label.placer = label_placer_fraction(frac=0.5),
                        breaks = -0.26, 
                        skip = 0,      # Don't skip this one
                        min.size = 0)  # Force it even on short lines
  }
  
  
  if (-0.23 %in% breakss) {
    
    NS_plot <- ggplot(pl[pl$Category==str,], aes(x=R, y=1/gamma, fill=NS)) +
      theme_bw()+
      geom_raster(interpolate=TRUE) +
      geom_contour2(aes(z=NS), breaks=breakss) +
      geom_text_contour(aes(z = NS),
                        color = "black",
                        size=4.5,
                        fontface='bold',
                        stroke=0.1,
                        rotate=TRUE,
                        label.placer = label_placer_fraction(frac=0.5),
                        breaks= c(breaks,-0.24),
                        skip=1) +
      scale_fill_gradientn(colours=c("#053061", "#4575b4", "#abd9e9", "white", "#fee090", "#f46d43", "#67001F"),
                           limits=range_diff,
                           values=scales::rescale(c(min(range_diff), -max(abs(range_diff))/2, -max(abs(range_diff))/10, 0, max(abs(range_diff))/10, max(abs(range_diff))/2, max(range_diff))),
                           oob=scales::squish) +
      ylim(1/3,8)+
      scale_y_continuous(breaks=c(1/3, 1:7),
                         labels=expression(1/3, 1,2,3,4,5,6,7)) +
      
      # labs(title=paste0(str_title, " (NS)")) +
      theme(text=element_text(size=16),
            legend.position="none",
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
      geom_text_contour(aes(z = NS),
                        color = "black", size = 4.5, fontface = 'bold',
                        stroke = 0.1, rotate = TRUE,
                        label.placer = label_placer_fraction(frac=0.5),
                        breaks = -0.24, 
                        skip = 0,      # Don't skip this one
                        min.size = 0)  # Force it even on short lines
  }
  
  
  
  # if (substr(str,6,9) != "1d00"){
  #  NS_plot <- NS_plot +scale_x_continuous(labels=NULL) 
  #  NM_plot <- NM_plot +scale_x_continuous(labels=NULL) 
  # }
  
  NS_NM_plots[[str]] <- list(NS=NS_plot, NM=NM_plot)
}


library(patchwork)


#Nplots

design <- "
AB
CD
"

N_plots_no_legend <- list()

N_plots_no_legend[["0.9"]] <-
  N_plots[["0.9"]] +
  labs(
    x = expression("Reproduction number, "~R),
    y = expression("Average time until detection, "~1/gamma),
    title = expression(r[C] == 0.9)
  ) +
  theme(plot.title = element_text(hjust = 0))



N_plots_no_legend[["0.9"]] <- N_plots[["0.9"]]+
  theme(legend.position="none",plot.title=element_text(hjust=0),
        axis.title.y = element_text(size=12,angle=90))+
  labs(x=expression("Reproduction number, "~R),
       y=expression("Average time until detection, "~1/gamma),
       title=expression(r[C] == 0.9))+
  scale_y_continuous(breaks = c(1/3, 1:7),
                     labels = expression(scriptstyle(1/3), 1, 2, 3, 4, 5, 6, 7))
N_plots_no_legend[["0.9_"]] <- N_plots[["0.9"]]+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        legend.position="none",
        plot.title=element_blank())+
  labs(x=expression("Reproduction number, "~R),
       y=expression("Average time until detection, "~1/gamma),
       title=expression(r[C] == 0.9))+
  scale_y_continuous(labels=NULL)

N_plots_no_legend[["0.5"]] <- N_plots[["0.5"]]+
  theme(legend.position="none",plot.title=element_text(hjust=0),
        axis.title.y = element_text(size=12,angle=90))+
  labs(x=expression("Reproduction number, "~R),
       y=expression("Average time until detection, "~1/gamma), 
       title=expression(r[C] == 0.5))+
  scale_y_continuous(breaks = c(1/3, 1:7),
                     labels = expression(scriptstyle(1/3), 1, 2, 3, 4, 5, 6, 7))
N_plots_no_legend[["0"]] <- N_plots[["0"]]+
  theme(legend.position="none",
        plot.title=element_text(hjust=0),
        axis.title.y = element_text(size=12,angle=90))+
  labs(x=expression("Reproduction number, "~R),
       y=expression("Average time until detection, "~1/gamma),
       title=expression(r[C] == 0))+
  scale_y_continuous(breaks = c(1/3, 1:7))


N_plots_no_legend[["0.9"]] <- N_plots_no_legend[["0.9"]]+
  theme(axis.title.x = element_blank())

# N_plots_no_legend[["0.75"]] <- N_plots_no_legend[["0.75"]]

N_plots_no_legend[["0"]] <- N_plots_no_legend[["0"]] + 
  theme(axis.title.y = element_blank()) +
  scale_y_continuous(breaks = c(1/3, 1:7), labels=NULL)






Ns <- N_plots_no_legend[["0.9"]] +
  plot_spacer() +
  N_plots_no_legend[["0.5"]] +
  N_plots_no_legend[["0"]] +
  plot_layout(
    design = design,
    guides = "collect"
  ) +
  labs(fill = expression(bold(E)*"[Y | Y > 0, no cohorting]")) +
  theme(
    legend.position = "right",
    legend.key.height = unit(3, "cm"),
    legend.key.width  = unit(0.8, "cm"),
    legend.title = element_text(
      angle = -90,
      vjust = 0.5,
      hjust = 0.5,
      size = 16
    ),
    legend.title.position = "right"
  )

show(Ns)




ggsave("Ns_plot_new.pdf", width=11*0.75,height=10*0.75,device=cairo_pdf, dpi=300)






design <- "
AB
CD
EF
GH
"


A <- NS_NM_plots$"0d90e0d10"$NS + scale_x_continuous(labels=NULL)+
  scale_y_continuous(breaks = c(1/3, 1:7),
                     labels = expression(scriptstyle(1/3), 1, 2, 3, 4, 5, 6, 7))
# labs(y=expression("Average time until detection, "~1/gamma))+
# theme(
#       axis.title.y = element_text(angle=90),
#       # axis.ticks.y = element_blank()
# )
B <- NS_NM_plots$"0d90e0d10"$NM + scale_x_continuous(labels=NULL)+
  scale_y_continuous(labels=NULL, breaks=c(1/3, 1:7))
C <- NS_NM_plots$"0d90e1d00"$NS +labs(y=expression("Average time until detection, "~1/gamma))+
  scale_y_continuous(breaks = c(1/3, 1:7),
                     labels = expression(scriptstyle(1/3), 1, 2, 3, 4, 5, 6, 7))
# theme(
#   axis.title.y = element_text(angle=90),
#   # axis.ticks.y = element_blank()
# )
D <- NS_NM_plots$"0d90e1d00"$NM + scale_y_continuous(labels=NULL,breaks=c(1/3, 1:7))

E <- NS_NM_plots$"0d50e0d50"$NS + scale_x_continuous(labels=NULL)+
  scale_y_continuous(breaks = c(1/3, 1:7),
                     labels = expression(scriptstyle(1/3), 1, 2, 3, 4, 5, 6, 7))
# labs(y=expression("Average time until detection, "~1/gamma))+
# theme(
#   axis.title.y = element_text(angle=90),
#   # axis.ticks.y = element_blank()
# )
F <- NS_NM_plots$"0d50e0d50"$NM + scale_x_continuous(labels=NULL)+
  scale_y_continuous(labels=NULL,breaks=c(1/3, 1:7))
G <- NS_NM_plots$"0d50e1d00"$NS +labs(x=expression("Reproduction number, "~R),
                                      y=expression("Average time until detection, "~1/gamma))+
  scale_y_continuous(breaks = c(1/3, 1:7),
                     labels = expression(scriptstyle(1/3), 1, 2, 3, 4, 5, 6, 7))
# theme(
#   axis.title.y = element_text(angle=90),
#   axis.title.x = element_text()
#   # axis.ticks.y = element_blank()
# )
H <- NS_NM_plots$"0d50e1d00"$NM + scale_y_continuous(labels=NULL,breaks=c(1/3, 1:7))
# labs(x=expression("Reproduction number, "~R))+
# theme(
#   axis.title.y = element_text(angle=90),
#   axis.title.x = element_text()
#   # axis.ticks.y = element_blank()
# )


R <- A +
  B+
  C+
  D+
  E+
  F+
  G+
  H+
  plot_layout(
    design = design,
    guides = "collect"
  ) +
  labs(fill = expression(
    phantom("Relative differenceoo , ") ~ frac(bold(E)*"[Y | Y > 0, cohorting]", bold(E)*"[Y | Y > 0, no cohorting]") - 1
  )) +
  theme(
    legend.position = "right",
    # legend.box.just = "center",
    # legend.justification = "center",
    legend.key.height = unit(4, "cm"),
    legend.key.width  = unit(0.8, "cm"),
    legend.title = element_text(
      angle = -90,
      vjust = 0.5,
      hjust = 0.5,
      size = 12),
    legend.title.position = "right",
    
    
  )

R <- R + plot_annotation(
  theme = theme(
    plot.margin = margin(b = 20, l = 65,t=22) # 'b' for bottom space, 'l' for left space
  )
)


R_lab <- ggdraw(R) +
  draw_label(expression("Average time until detection, "~1/gamma), x = 0.1, y = 0.28, angle = 90, size = 14) +
  draw_label(expression("Average time until detection, "~1/gamma), x = 0.1, y = 0.75, angle = 90, size = 14) +
  draw_label(expression("Reproduction number, "~R),  x = 0.45,  y = 0.02, size = 14)+
  draw_label("Strong Cohorting",  x = 0.315,  y = 0.98, size = 15, fontface="bold")+
  draw_label("Medium Cohorting",  x = 0.64,  y = 0.98, size = 15, fontface="bold")+
  draw_label(expression(r[C]==0.9),  x = 0.02,  y = 0.74, size = 14, fontface="bold",angle=90)+
  draw_label(expression(epsilon== 1-r[C]),  x = 0.05,  y = 0.84, size = 14, fontface="bold",angle=90)+
  draw_label(expression(epsilon== 1),  x = 0.05,  y = 0.64, size = 14, fontface="bold",angle=90)+
  draw_label(expression(r[C]==0.5),  x = 0.02,  y = 0.27, size = 14, fontface="bold",angle=90)+
  draw_label(expression(epsilon== 1-r[C]),  x = 0.05,  y = 0.37, size = 14, fontface="bold",angle=90)+
  draw_label(expression(epsilon== 1),  x = 0.05,  y = 0.17, size = 14, fontface="bold",angle=90)+
  draw_label("Relative difference,",  x = 0.95,  y = 0.65, size = 17, fontface="plain",angle=-90)
R_lab
# ggsave("Reduction_new.pdf", plot=R_lab,width=8*0.9,height=11*0.9,device="pdf", dpi=300)
