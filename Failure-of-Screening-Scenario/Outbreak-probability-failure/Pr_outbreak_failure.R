library(parallel)
library(ggplot2)
library(dplyr)
library(Rcpp)
library(metR)


library(patchwork)
library(gridExtra)
library(cowplot)
library(latex2exp)


library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))

# 1. Compile the C++ engine
# message("Compiling C++ simulation engine...")
sourceCpp("Y_engine_failure_truncated.cpp")
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
n    <- 0
zeta <- 1/3

Xpl      <- list()
range_pl <- numeric(0)

# Determine the number of cores to use based on Slurm allocation
# If Slurm doesn't set it, default to 1 to be safe
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
message(sprintf("Running on %d cores...", n_cores))

# 3. The Main Grid Loop
for (rC in c(0.9, 0.75, 0.5, 0.25,0)) {
    for (eps in c(1, 1-rC)) {
    message(sprintf("Starting rC=%.2f eps=%.2f at %s", rC, eps, Sys.time()))
    
    # Capture loop variables explicitly for workers
    rC_  <- rC
    eps_ <- eps
    
    Rvals     <- seq(1, 10, length.out = 50)
    gammavals <- 3 / seq(1, 21, length.out = 50)
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
        return(1-X[1])
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
    
    saveRDS(list(
  data  = grid,          # has R, gamma, S, M, N
  plots = list(S = Spl, M = Mpl, N = Npl)
), paste0("p", rC_, eps_, ".rds"))
    
    message(sprintf("Done rC=%.2f eps=%.2f at %s", rC_, eps_, Sys.time()))
}
}

# Final plot formatting
for (i in seq_along(Xpl)) {
    Xpl[[i]] <- Xpl[[i]] + scale_fill_viridis_c(limits = range_pl, option = "turbo")
}

saveRDS(Xpl, "Xpl_list_careful.rds")
message("All simulations complete!")



################################################################################

# split into datasets
pl0d90_rds  <- readRDS("p0.90.1.rds")
pl0d90uc_rds <- readRDS("p0.91.rds")
pl0d75_rds  <- readRDS("p0.750.25.rds")
pl0d75uc_rds <- readRDS("p0.751.rds")
pl0d50_rds  <- readRDS("p0.50.5.rds")
pl0d50uc_rds <- readRDS("p0.51.rds")
pl0d25_rds  <- readRDS("p0.250.75.rds")
pl0d25uc_rds <- readRDS("p0.251.rds")
pl0d0_rds   <- readRDS("p01.rds")

pl0d90  <- pl0d90_rds$data
pl0d90uc <- pl0d90uc_rds$data;  pl0d90uc$N <- pl0d90$N
pl0d75  <- pl0d75_rds$data
pl0d75uc <- pl0d75uc_rds$data;  pl0d75uc$N <- pl0d75$N
pl0d50  <- pl0d50_rds$data
pl0d50uc <- pl0d50uc_rds$data;  pl0d50uc$N <- pl0d50$N
pl0d25  <- pl0d25_rds$data
pl0d25uc <- pl0d25uc_rds$data;  pl0d25uc$N <- pl0d25$N
pl0d0   <- pl0d0_rds$data

pl0d90$Category  <- "0d90e0d10";  pl0d90uc$Category <- "0d90e1d00"
pl0d75$Category  <- "0d75e0d25";  pl0d75uc$Category <- "0d75e1d00"
pl0d50$Category  <- "0d50e0d50";  pl0d50uc$Category <- "0d50e1d00"
pl0d25$Category  <- "0d25e0d75";  pl0d25uc$Category <- "0d25e1d00"
pl0d0$Category   <- "0d0e1d00"

pl <- rbind(pl0d90, pl0d90uc, pl0d75, pl0d75uc,
            pl0d50, pl0d50uc, pl0d25, pl0d25uc, pl0d0)

range_ <- range(pl$S, pl$M, pl$N, na.rm = TRUE)

# shared scale (applied to all plots)
shared_scale <- scale_fill_viridis_c(
  option = "turbo",
  limits = range_,
  name   = "Probability of outbreak originating from general patient",
  guide  = guide_colorbar(
    direction      = "vertical",
    title.position = "right",
    title.hjust    = 0.5,
    barwidth       = unit(0.5, "cm"),
    barheight      = unit(15, "cm"),
    title.theme    = element_text(angle=-90)
  )
)

# plot builder function
breaksn <- seq(0, 1, 0.1)

make_plot <- function(data, var) {
  ggplot(data, aes(x = R, y = 1/gamma, fill = .data[[var]])) +
    theme_bw() +
    geom_raster(interpolate = TRUE) +
    geom_contour2(aes(z = .data[[var]]), breaks = breaksn) +
    geom_text_contour(aes(z = .data[[var]]),
                      color = "black", size = 4.5, fontface = "bold",
                      stroke = 0.1, rotate = TRUE,
                      label.placer = label_placer_fraction(frac = 0.5),
                      breaks = breaksn, skip = 0) +
    shared_scale +
    scale_y_continuous(breaks = c(1/3, 1:7),
                       labels = expression(1/3, 1, 2, 3, 4, 5, 6, 7)) +
    theme(text = element_text(size = 16),
          legend.position = "right",
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
}

# bBuild plots
d90 <- pl[pl$Category == "0d90e1d00", ]
d50 <- pl[pl$Category == "0d50e1d00", ]

# Top row: rC=0.9
p_90N <- make_plot(d90, "N") + theme(axis.text.x = element_blank())
p_90M <- make_plot(d90, "M") + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
p_90S <- make_plot(d90, "S") + theme(axis.text.x = element_blank(), axis.text.y = element_blank())

# Bottom row: rC=0.5
p_50N <- make_plot(d50, "N")
p_50M <- make_plot(d50, "M") + theme(axis.text.y = element_blank())
p_50S <- make_plot(d50, "S") + theme(axis.text.y = element_blank())

# assemble
design <- "
ABC
DEF
"

p <- (p_90N + p_90M + p_90S +
        p_50N + p_50M + p_50S +
        plot_layout(design = design, guides = "collect")) +
  plot_annotation(theme = theme(
    legend.position = "right",
    plot.margin     = margin(b = 20, l = 50,t=30)
  ))

# add axis labels via cowplot
R_lab <- ggdraw(p) +
  draw_label(expression("Average time until detection, " ~ 1/gamma),
             x = 0.05, y = 0.51, angle = 90, size = 20) +
  draw_label(expression(r[C]==0.9),
             x = 0.01, y = 0.75, angle = 90, size = 18, fontface = "bold") +
  draw_label(expression(r[C]==0.5),
             x = 0.01, y = 0.3, angle = 90, size = 18, fontface = "bold") +
  draw_label(expression("Reproduction number, " ~ R),
             x = 0.48, y = 0.02, size = 20)+
  draw_label("No cohorting",
             x = 0.225, y = 0.98, size = 18, fontface='bold')+
  draw_label("Medium cohorting",
             x = 0.5, y = 0.98, size = 18, fontface='bold')+
  draw_label("Strong cohorting",
             x = 0.775, y = 0.98, size = 18, fontface='bold')

show(R_lab)


ggsave("p_out_failure.pdf",plot=R_lab, width=15*0.75,height=10*0.75,device="pdf", dpi=300)