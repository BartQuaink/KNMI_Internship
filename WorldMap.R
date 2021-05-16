# Required functions
require(pacman)

# install required libs
libs <- c("raster", "ggplot2", "rgdal", "lubridate", "readr", "sf", "dplyr", "tidyverse" , 
          "leaflet", "rnaturalearth", "rnaturalearthdata")
pacman::p_load(char = libs)

# World Mapping -----------------------------------------------------------

#' This script has all the plotting for geographical data
#' 
#' 
#' 
#' 

lon_all <- c(ncvar_get(nc, "lon"))
lat_all <- c(ncvar_get(nc, "lat"))
g_matrix <- matrix(c(lon_all, lat_all, 1:length(lat_all)), ncol = 3)

df <- data.frame(lon_all, lat_all, idx = 1:length(lat_all))

m <- leaflet(data = g_matrix) %>%
  addMarkers(g_matrix[, 1], g_matrix[, 2])

# Filter locations
lat_all <- ncvar_get(nc, "lat")
lat_all[lat_all < -55 & lat_all > -60] <- NA

lon_all <- ncvar_get(nc, "lon")
lon_all[lon_all < -80 & lon_all > -90] <- NA

all_clouds <- ncvar_get(nc, "geum_avhrr_1b_cloud_fraction")
all_clouds[all_clouds != 0] <- NA


## ggplot
world_map <- map_data("world")

world_2 <- ne_countries(scale = "medium", returnclass = "sf")

# ( , X , ) -> lower is more west
ggplot() +
  geom_sf(data = world_2, fill = "grey90") +
  geom_point(aes(x = as.numeric(lon_all), y = as.numeric(lat_all), colour = clouds), 
             size = 1, 
             alpha = 1) +
  scale_colour_gradient(low = "blue", high = "yellow", 
                        name = "",
                        guide=guide_colourbar(reverse = TRUE),
                        limits = c(0, 15000), 
                        breaks = c(0, 2500, 5000, 7500, 10000, 12500, 15000)) +
  labs(title = "IASI observations with cloud cover = 0",
       subtitle = "Date : 2011-09-25",
       x = "longitude",
       y = "latitude") +
  coord_sf(xlim = c(-90,-60), ylim = c(-65,-30), expand = FALSE) +
  theme_bw()

# Single point
ggplot() +
  geom_sf(data = world_2, fill = "grey90") +
  geom_point(aes(x = as.numeric(lon), y = as.numeric(lat)), 
             size = 2) +
  labs(title = "All IASI observations on 2011-08-01",
       x = "longitude",
       y = "latitude") +
  coord_sf(xlim = c(-95,-60), ylim = c(-65,-30), expand = FALSE) +
  theme_bw()

