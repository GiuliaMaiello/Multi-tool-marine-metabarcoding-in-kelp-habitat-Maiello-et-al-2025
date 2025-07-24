##### Libraries #####
library("dplyr")
library("ggmap")
library("ggplot2")
library("ggrepel")
library("ggspatial")
library("grid")
library("readr")
library("stadiamaps")

mp_Staz =  read.csv("CHNMS_coordinates_all_Sept2022.csv", sep = ";", header = T)
is.num <- sapply(mp_Staz, is.numeric)
mp_Staz[is.num] <- lapply(mp_Staz[is.num], round, 2)

register_stadiamaps(key= "7271dcb5-8253-4f4d-a3c7-6e1ce32244ce")

map6 <- get_stadiamap(bbox = c(left = -121.5,
                               bottom = 34.2,
                               right =  -119.5,
                               top = 35.65),
                      zoom = 12,
                      maptype = "stamen_terrain_background")  

colnames(mp_Staz)[c(2,3)] = c("lat","lon")

map = ggmap(map6) + 
  geom_point(data = mp_Staz, aes(x = lon, y = lat, col = SampleType, shape = X), size = 4) + 
  scale_color_manual(values = c("Green4","darkorange4","red4", "orangered2")) +
  geom_text_repel(data = mp_Staz, aes(x=lon, y=lat, label=Site, color = SampleType), size = 4, force = 15,
                  fontface = "bold", inherit.aes = F, max.overlaps = 20) + 
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.2, "in"), pad_y = unit(6.5, "in"), 
                         style = north_arrow_fancy_orienteering) +
  xlab("Longitude (W)") +
  ylab("Latitude (N)") + 
  theme(legend.position = "right", legend.key.size = unit(1,"cm"))
