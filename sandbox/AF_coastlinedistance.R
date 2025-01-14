library(sf)
library(tidyverse)
library(rnaturalearth)
library(smoothr)
library(raster)
library(here)


# doing this with a high-resolution map of the coastline gives too much resolution around the Chesapeake, Long Island, etc. I want a really coarse shape (similar to the 10m isobath approach)

xmin=-125
xmax=-106
ymin=20
ymax=42


usamap <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")
valid_usamap <- usamap[st_is_valid(usamap), ]
usamap_un <- st_union(valid_usamap)[1] %>%
  st_cast("MULTILINESTRING")


bbox1 <- st_set_crs(st_as_sf(as(raster::extent(xmin, xmax, ymin, ymax), "SpatialPolygons")), st_crs(usamap))

wusmap <-  usamap_un %>% 
  st_intersection(bbox1) # can replace with st_crop when CRAN version of sf updates 

# neusmap <- usamap %>% 
#   st_intersection(bbox1) %>% # can replace with st_crop when CRAN version of sf updates 
#   st_difference(bbox2) # get rid of extra non coastal line 

smoothmap <- wusmap %>% 
  smoothr::smooth(method="ksmooth", smoothness=8)
# smoother was applied incrementally more until the Chesapeake went away 
# https://cran.r-project.org/web/packages/smoothr/vignettes/smoothr.html

st_length(smoothmap)

smoothgeom <- smoothmap %>% 
  as("Spatial") %>% 
  geom()

geomdists <- pointDistance(smoothgeom[-nrow(smoothgeom), c("x", "y")], smoothgeom[-1, c("x", "y")], lonlat=TRUE)
coastdistdat <- data.frame(smoothgeom[, c('x','y')], seglength=c(0, geomdists))
coastdistdat$lengthfromhere <- rev(cumsum(rev(coastdistdat[,"seglength"])))
# first row should match st_length(smoothmap)

st_write(smoothmap, here("processed-data","westcoast_coastline.shp"))
write_rds(coastdistdat, here("processed-data","westcoast_coastdistdat.rds"))
rm(list=ls())

shape_data <- st_read("./processed-data/coastline.shp")
# visualization
ggplot() +
  geom_sf(data = smoothmap, fill = "lightblue", color = "darkblue") +  # Customize fill and border colors
  theme_minimal() +  # Use a minimal theme
  labs(title = "Shapefile Visualization")  # Add a title

ggplot(data=coastdistdat, aes(x=x, y=y)) +
  geom_point()

