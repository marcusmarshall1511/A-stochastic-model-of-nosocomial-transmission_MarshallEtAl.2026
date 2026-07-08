library(future)
library(future.apply)
library(progressr)
library(ggplot2)
library(dplyr)
library(scales)
library(metR)



library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))




options(progressr.enable = TRUE)
sink(stdout(), type = "message")  # redirect messages to stdout

P_out <- function(r_I,R,distr="Hyp"){
  
  if (distr == "Hyp"){
  delta_iso <- I(list(c(0,6.69470955813761, 0.555473683554152, 0.55547824700609, 0.555478430962926)))
  }
  if (distr == "Exp"){
  delta_iso <- I(list(c(0,0.178985825471926)))
  }
  
  LOR <- 7
  beta <- R/(LOR*12)
  beta_0 <- 2 * (1-r_I) * beta
  prod <- 1
  for (l in 2:length(delta_iso[[1]])){
    prod <- prod*delta_iso[[1]][l]/(delta_iso[[1]][l]+beta_0)
  }
  # for (l in 2:length(delta_iso[[1]])){
  #   prod <- prod*delta_iso[[1]][l]/(delta_iso[[1]][l]+beta*(1-r_I)*NHQ)
  # }
  return(1-prod)
}


rIvals <- seq(0,1,length.out=100)
Rvals <- seq(1,10,length.out=100)
grid <- expand.grid(rI = rIvals, R = Rvals)


grid$Pout_hyp <- mapply(P_out, grid$rI, grid$R, "Hyp")

grid$Pout_exp <- mapply(P_out, grid$rI, grid$R, "Exp")



grid$ExRE <- (grid$Pout_exp-grid$Pout_hyp)/grid$Pout_hyp
grid$ExRE[grid$rI==1] <- 0


range_ <- range(grid$Pout_hyp, grid$Pout_exp)
  



pHyp <- ggplot(grid, aes(x=R, y=rI, fill=Pout_hyp))+
  geom_raster()+
  geom_contour2(aes(z=Pout_hyp)) +
  geom_text_contour(aes(z = Pout_hyp),
                    color = "black",
                    size=4.5,
                    fontface='bold',
                    stroke=0.1,
                    rotate=TRUE,
                    label.placer = label_placer_fraction(frac=0.5),
                    skip=0) +
  theme_bw()+
  scale_fill_viridis_c(
    option = "turbo",
    guide = guide_colorbar(
      barwidth = 20,        # long horizontal bar
      barheight = 1.5,      # thin height
      title.position = "top",
      title.hjust = 0.5,
      direction = "horizontal"
    )
  ) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )+
  labs(fill = "Outbreak probability", y = bquote("Isolation efficacy, " ~ r[italic(I)]), x = bquote("Reproduction number, " ~ R))+
  theme(
    legend.title = element_text(angle = 0),  # rotate text vertically
    text=element_text(size=20)
  )

show(pHyp)

ggsave("Pout_Hyp.pdf", pHyp, width = 8*0.6, height = 9*0.6, units = "in", dpi = 300,
       device = "pdf")




  
pExRE <- ggplot(grid, aes(x=R, y=rI, fill=100*ExRE))+
  geom_raster()+
  geom_contour2(aes(z=100*ExRE)) +
  geom_text_contour(aes(z = 100*ExRE),
                    color = "black",
                    size=4.5,
                    fontface='bold',
                    stroke=0.1,
                    rotate=TRUE,
                    label.placer = label_placer_fraction(frac=0.5),
                    skip=0) +
  theme_bw()+
  scale_fill_viridis_c(
    option = "inferno",
    direction=-1,
    guide = guide_colorbar(
      barwidth = 20,        # long horizontal bar
      barheight = 1.5,      # thin height
      title.position = "top",
      title.hjust = 0.5,
      direction = "horizontal"
    )
  ) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )+
  labs(
    fill = "Relative difference (%)",
    y = bquote("Isolation efficacy, " ~ r[italic(I)]),
    x = bquote("Reproduction number, " ~ R)
  )+
  theme(
    legend.title = element_text(angle = 0),  # rotate text vertically
    text=element_text(size=20)
  )

show(pExRE)

ggsave("Pout_ExRE.pdf", pExRE, width = 8*0.6, height = 9*0.6, units = "in", dpi = 300,
       device = "pdf")


saveRDS(grid, "P_out_list.rds")