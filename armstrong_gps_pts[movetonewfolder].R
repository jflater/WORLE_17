install.packages("sf")
library(sf)
library(ggplot2)

armstrong_shape <- st_read("../../Armstrong_RainSim_2018_plots_grid.shx")

ggplot() + 
  geom_sf(data = armstrong_shape, size = 3, color = "black", fill = "cyan1") + 
  ggtitle("Armstrong Plots") + 
  coord_sf()
