library(ggplot2)
library(readxl)
library(dplyr)
library(expm)
#LOI

library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))


data <- read_excel("COVID-19_ICU_cases.xlsx")
raw_data <- data$`Before april 2021` %>%
  na.omit %>%
  as.vector()

raw_data

df_raw <- data.frame(Raw = raw_data, Category = 'Data')

# p1 <- ggplot(df_raw, aes(x=Raw))+
#   stat_count(aes(y = after_stat(count) / sum(after_stat(count))), fill='grey', color='black')+
#   labs(y = 'Probability density')
# 
# show(p1)

# dhypoexp <- function(x, par) {
#   alpha <- c(1, numeric(length(par) - 1))
#   Theta <- matrix(0, nrow = length(par), ncol = length(par))
#   for (i in 1:length(par)) {
#     Theta[i, i] = -par[i]
#     if (i < length(par)) {
#       Theta[i, i + 1] = par[i]
#     }
#   }
#   ones <- rep(1, length(par))
#   as.numeric(-alpha %*% expm(x * Theta) %*% Theta %*% ones)
# }


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


domain <- seq(0,150,0.1)

dens <- dhypoexp(domain, c(0.0703227823411154, 8.3198836017103, 8.30400755726981, 8.64881497333471, 8.72619937561971, 8.7971877006993, 10.3122997921632, 9.51076894142005 
))

df_dens <- data.frame(x=domain, y=dens, Category = 'Hypoexponential')

dens_exp <- dexp(domain, rate=0.0667245565690008)
  
  # dexp(domain, c(0.066715))

df_dens <- rbind(df_dens, data.frame(x=domain, y=dens_exp, Category = 'Exponential'))

# df_dens <- rbind(df_dens, data.frame(x=0, y=1/0.066715, Category = 'Exponential'))

p1 <- ggplot(df_raw, aes(x=Raw))+
  stat_count(aes(y = after_stat(count) / sum(after_stat(count)), fill = Category, color= Category), alpha=0.4)+
  labs(y = 'Probability density', x = 'LOS of COVID patients in ICU (Days)')+
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
  #       color= c("black", '#F8766D', '#00BFC4'),
  #       size=c(1,1,0.5)
  #     )
  #   )
  # ) +
  ylim(0,0.16)+
  xlim(0,80)+
  theme_bw()+
  theme(text = element_text(size=20), legend.title = element_blank(), legend.position = c(0.75, 0.85),
        # axis.text.x = element_blank()
  )
  

show(p1)

# 
ggsave(paste0("LOS_fit_COVID.pdf"), p1, width = 6, height = 6, units = "in", dpi = 300,
       device = "pdf")


