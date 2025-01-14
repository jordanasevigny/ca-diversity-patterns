# map out the sites surveyed and color points by year - should add in jitter to see all years surveyed
# date: 10/28/2024
# author: Jordana Sevigny
# contact: jordana.sevigny@gmail.com


# libraries
library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library(dplyr)
library(tidyverse)

# load a world map
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

# load in data
survey_taxa_dmy <- data.frame(read.csv("./data/MARINe_survey_taxa_dmy.csv", header=TRUE))

# add a column of the unique survey ID
survey_taxa_dmy$survey_ID <- paste(survey_taxa_dmy$site_code, survey_taxa_dmy$survey_rep, sep="_")

# make new dataframe with unique survey_ID, lat, lon, Year, and Date
unique_survey <- survey_taxa_dmy %>%
  select(c("survey_ID", "latitude", "longitude", "year", "sample_date")) %>%
  distinct(survey_ID, latitude, longitude, year, sample_date) %>%
  separate(survey_ID, into = c("site_code", "survey_rep"), sep = "_")

# need to add some longitude jitter into the points so they all show up
# will not add jitter for survey 1, will add rep*jitter to each latitude after (survey 2 = rep 1)
jitter <- 1.25

unique_survey$lon_jitter <- NA
for (i in 1:nrow(unique_survey)) {
  if (unique_survey$survey_rep[i] == 1) {
    unique_survey$lon_jitter[i] = unique_survey$longitude[i]
  } else {
    unique_survey$lon_jitter[i] = unique_survey$longitude[i] + (as.numeric(unique_survey$survey_rep[i])-1)*jitter
  }
}

# map with survey points colored by year and jittered
ggplot(data = world) +
  geom_sf() +
  geom_point(data = unique_survey, aes(x = lon_jitter, y = latitude, color = year, group = year), size = 2.5) +
  scale_color_gradient(name = "year", low="blue", high="red") +
  xlab("Longitude") + ylab("Latitude") +
  labs(color="Sample date") +
  ggtitle("MARINe surveys by year") +
  coord_sf(xlim = c(-110, -130), ylim = c(25, 45), expand = FALSE)
ggsave("MARINe_site_map_and_count.png")

# map with survey point size determined by reps
maxrep_unique_survey <- unique_survey %>%
  group_by(site_code) %>%
  slice_max(survey_rep)

ggplot(data = world) +
  geom_sf() +
  geom_point(data = maxrep_unique_survey, aes(x = longitude, y = latitude, size=survey_rep), color="red", alpha=0.3) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("MARINe surveys sized by reps") +
  coord_sf(xlim = c(-110, -130), ylim = c(25, 45), expand = FALSE)

# weird that there are decimals on reps... must have been retroactively named maybe..

# histogram of sites & sample freq
# make new dataframe with unique location name, survey_ID, lat, lon, Year, and Date
unique_survey_hist <- survey_taxa_dmy %>%
  select(c("marine_site_name", "survey_ID", "latitude", "longitude", "year", "sample_date")) %>%
  distinct(marine_site_name, survey_ID, latitude, longitude, year, sample_date)

ggplot(data=unique_survey_hist, aes(marine_site_name, fill=latitude)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=3)) +
  ggtitle("MARINe sites' survey counts")
ggsave("MARINe_site_survey_count.png") 

# map only summer surveys