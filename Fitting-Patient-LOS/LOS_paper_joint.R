library(ggplot2)
library(dplyr)
library(expm)
library(DEoptim)
library(parallel) 
set.seed(100)


library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))



output_file <- "minimum_of_two.txt"
writeLines("", output_file)

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


rhypoexp <- function(n, par){
  samples <- numeric(0)
  for (i in 1:n){
    sample <- 0
    for (rate in par){
      sample <- sample+rexp(n=1,rate=rate)
    }
    samples <- c(samples, sample)}
  return(samples)
}

hypo_par <- c(0.0703227823411154, 8.3198836017103, 8.30400755726981, 8.64881497333471, 8.72619937561971, 8.7971877006993, 10.3122997921632, 9.51076894142005 )


minimum_of_two <- function(n,par,meanlog,sdlog){
  # samples <- pmin(rhypoexp(n,par),rlnorm(n=n,meanlog=meanlog,sdlog=sdlog))
  samples <- pmin(rhypoexp(n,par),rgamma(n=1000, shape=5, scale=1.4))
}

raw_data <- minimum_of_two(n=1000,hypo_par,meanlog=1.63, sdlog=0.5)
# raw_data

################################################################################

dhypoexp <- function(x, par) {
  alpha <- c(1, numeric(length(par) - 1))
  Theta <- matrix(0, nrow = length(par), ncol = length(par))
  for (i in 1:length(par)) {
    Theta[i, i] <- -par[i]
    if (i < length(par)) {
      Theta[i, i + 1] <- par[i]
    }
  }
  
  ones <- rep(1, length(par))
  sapply(x, function(xx) {
    as.numeric(-alpha %*% expm(xx * Theta) %*% Theta %*% ones)
  })
}


domain <- seq(0,80,0.1)

dens <- dhypoexp(domain, c(6.69470955813761, 0.555473683554152, 0.55547824700609, 0.555478430962926 
  ))



df_dens <- data.frame(x=domain, y=dens, Category = 'Hypoexponential')

# domain <- seq(0,150,0.1)

# dens <- dhypoexp(domain, c(0.070445, 4.926052, 4.926058, 4.926068, 4.92609
# ))

# df_dens <- rbind(df_dens, data.frame(x=domain, y=dens, Category = 'LOS of COVID'))

dens <- dgamma(domain, shape=5, scale = 1.4)

# df_dens <- rbind(df_dens, data.frame(x=domain, y=dens, Category = 'LOI of COVID'))

dens <- dexp(domain, rate=0.178985825471926)

df_dens <- rbind(df_dens, data.frame(x=domain, y=dens, Category = 'Exponential'))

levels = c("Data", "Exponential", "Hypoexponential")
           # , "LOI of COVID", "LOS of COVID")

df_dens$Category <- factor(
  df_dens$Category,
  levels = levels
)

breaks <- seq(0.5, 20.5, by = 1)
labels <- 1:20  # centers of bins
binned <- cut(
  raw_data,
  breaks = breaks,
  labels = labels,
  right = FALSE,  # so each bin is [a, b)
  include.lowest = TRUE
)
df <- as.data.frame(table(factor(binned, levels = as.character(labels))))
colnames(df) <- c("value", "count")
df$prob <- df$count / sum(df$count)
df$Category <- "Data"

df$Category <- factor(
  df$Category,
  levels = levels
)
df$value <- as.numeric(as.character(df$value))


# p1 <- ggplot(df_dens, aes(x=x, y=y,fill= Category, color=Category))+
#   geom_area(df_dens[df_dens$Category == "Hypoexponential" | df_dens$Category == "Exponential", ], alpha = 0.4, position = 'identity')+
#   scale_fill_manual(values= c('LOI of COVID'="#00BA38", 'LOS of COVID' = '#F8766D' , 'Hypoexponential' = "#C49A00", 'Exponential' = "#00B6EB"))+
#   scale_color_manual(values=c('LOI of COVID'="#00BA38", 'LOS of COVID' = '#F8766D' , 'Hypoexponential' = "#C49A00", 'Exponential' = "#00B6EB"  ))+
#   theme_bw()+
#   theme(text = element_text(size=20), legend.title = element_blank(), legend.position = c(0.7, 0.85),
#         # axis.text.x = element_blank()
#   )+
#   labs(x = 'Days', y = 'Probability density')+
#   ylim(0,0.16)+
#   xlim(0,80)
# 
# show(p1)




p1 <- ggplot(df, aes(x=value))+
  geom_col(aes(y = prob, fill = Category, color = Category), alpha = 0.4, position = "identity")+
  labs(y = 'Probability density', x = 'Isolated patient removal process (Days)')+
  geom_area(data=df_dens,aes(x=x,y=y, fill=Category,color=Category),alpha = 0.4, position = "identity")+
  scale_fill_manual(values= c('Data'="grey", 'Exponential' = '#F8766D' , 'Hypoexponential' = '#00BFC4' ))+
  scale_color_manual(values=c('Data'="black", 'Exponential' = '#F8766D' , 'Hypoexponential' = '#00BFC4'  ))+
  # guides(
  #   color = guide_legend(
  #     override.aes = list(
  #       fill = c("grey", '#F8766D', '#00BFC4'),
  #       alpha = c(0.4, 0.4, 0.4),
  #       linewidth = c(.5, 1, 1),
  #       linetype = c(0, 1, 1),
  #       # color= c(NA, '#F8766D', '#00BFC4'),
  #       size=c(1,1,0.5)
  #     )
    # )
  # ) +
  ylim(0,0.16)+
  # xlim(0,80)+
  theme_bw()+
  theme(text = element_text(size=20), legend.title = element_blank(), legend.position = c(0.75, 0.85),
        # axis.text.x = element_blank()
  )+
  xlim(0,20)


show(p1)
ggsave(paste0("LOS_paper_joint.pdf"), p1, width = 6, height = 6, units = "in", dpi = 300,
       device = "pdf")
