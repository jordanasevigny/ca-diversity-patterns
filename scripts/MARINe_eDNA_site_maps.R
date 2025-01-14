# eDNA sites vs 2001-2003 MARINe sites

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

# 2001 to 2003

# load in data
MARINe_df <- data.frame(read.csv("./data/filtered_MARINe_survey_taxa_dmy.csv", header=TRUE)) %>%
  filter(year %in% c(2001, 2002, 2003)) %>%
  select(c("marine_site_name", "latitude", "longitude")) %>%
  distinct() %>%
  rename(Sites = marine_site_name, Latitude = latitude, Longitude = longitude) %>%
  mutate(Source = "MARINe")
#write.csv(MARINe_df, 'filtered_MARINe_sites_2001-2003.csv')

eDNA_site_df <- data.frame(read.csv("./data/eDNA_sites.csv", header=TRUE)) %>%
  select(c("Sites", "Latitude", "Longitude")) %>%
  mutate(Source = "eDNA")

combo_df <- rbind(MARINe_df, eDNA_site_df)

# plot map
ggplot(data = world) +
  geom_sf() +
  geom_point(data = combo_df, aes(x = Longitude, y = Latitude, color = Source), size = 2.5) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("eDNA Summer 2024 and MARINe 2001-2003 survey sites") +
  coord_sf(xlim = c(-110, -130), ylim = c(25, 45), expand = FALSE) +
  facet_wrap(~ Source) 
ggsave("filtered_MARINe_2001-2003_eDNA_site_map.png")

# 2001 only

# load in data
MARINe_df <- data.frame(read.csv("./data/MARINe_survey_taxa_dmy.csv", header=TRUE)) %>%
  filter(year %in% 2001) %>%
  select(c("marine_site_name", "latitude", "longitude")) %>%
  distinct() %>%
  rename(Sites = marine_site_name, Latitude = latitude, Longitude = longitude) %>%
  mutate(Source = "MARINe")
eDNA_site_df <- data.frame(read.csv("./data/eDNA_sites.csv", header=TRUE)) %>%
  select(c("Sites", "Latitude", "Longitude")) %>%
  mutate(Source = "eDNA")

combo_df <- rbind(MARINe_df, eDNA_site_df)

# plot map
ggplot(data = world) +
  geom_sf() +
  geom_point(data = combo_df, aes(x = Longitude, y = Latitude, color = Source), size = 2.5) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("MARINe 2001 and eDNA survey sites") +
  coord_sf(xlim = c(-110, -130), ylim = c(25, 45), expand = FALSE) +
  facet_wrap(~ Source) 
ggsave("MARINe_2001_eDNA_site_map.png")

# 2002 only

# load in data
MARINe_df <- data.frame(read.csv("./data/MARINe_survey_taxa_dmy.csv", header=TRUE)) %>%
  filter(year %in% 2002) %>%
  select(c("marine_site_name", "latitude", "longitude")) %>%
  distinct() %>%
  rename(Sites = marine_site_name, Latitude = latitude, Longitude = longitude) %>%
  mutate(Source = "MARINe")
eDNA_site_df <- data.frame(read.csv("./data/eDNA_sites.csv", header=TRUE)) %>%
  select(c("Sites", "Latitude", "Longitude")) %>%
  mutate(Source = "eDNA")

combo_df <- rbind(MARINe_df, eDNA_site_df)

# plot map
ggplot(data = world) +
  geom_sf() +
  geom_point(data = combo_df, aes(x = Longitude, y = Latitude, color = Source), size = 2.5) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("MARINe 2002 and eDNA survey sites") +
  coord_sf(xlim = c(-110, -130), ylim = c(25, 45), expand = FALSE) +
  facet_wrap(~ Source) 
ggsave("MARINe_2002_eDNA_site_map.png")

# 2003 only

# load in data
MARINe_df <- data.frame(read.csv("./data/MARINe_survey_taxa_dmy.csv", header=TRUE)) %>%
  filter(year %in% 2003) %>%
  select(c("marine_site_name", "latitude", "longitude")) %>%
  distinct() %>%
  rename(Sites = marine_site_name, Latitude = latitude, Longitude = longitude) %>%
  mutate(Source = "MARINe")
eDNA_site_df <- data.frame(read.csv("./data/eDNA_sites.csv", header=TRUE)) %>%
  select(c("Sites", "Latitude", "Longitude")) %>%
  mutate(Source = "eDNA")

combo_df <- rbind(MARINe_df, eDNA_site_df)

# plot map
ggplot(data = world) +
  geom_sf() +
  geom_point(data = combo_df, aes(x = Longitude, y = Latitude, color = Source), size = 2.5) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("MARINe 2003 and eDNA survey sites") +
  coord_sf(xlim = c(-110, -130), ylim = c(25, 45), expand = FALSE) +
  facet_wrap(~ Source) 
ggsave("MARINe_2003_eDNA_site_map.png")


# are any of the 2001 to 2003 sites duplicates?
MARINe_df <- data.frame(read.csv("./data/MARINe_survey_taxa_dmy.csv", header=TRUE)) %>%
  filter(year %in% c(2001, 2002, 2003)) %>%
  select(c("marine_site_name", "latitude", "longitude", "year", "sample_date")) %>%
  distinct()
MARINe_2001_2003_survey_count <- MARINe_df %>%
  count(marine_site_name) %>%
  left_join(MARINe_df, by="marine_site_name") %>%
  select(-sample_date) %>%
  distinct()

# plot map
ggplot(data = world) +
  geom_sf() +
  geom_point(data = MARINe_2001_2003_survey_count, aes(x = longitude, y = latitude, color = n), size = 2.5) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("MARINe 2001-2003 survey sites colored by survey #") +
  coord_sf(xlim = c(-110, -130), ylim = c(25, 45), expand = FALSE)
ggsave("MARINe_count_2001-2003_eDNA_site_map.png")


# marine map of each year's surveys

# load in data
MARINe_df <- data.frame(read.csv("./data/MARINe_survey_taxa_dmy.csv", header=TRUE)) %>%
  select(c("marine_site_name", "latitude", "longitude", "year")) %>%
  distinct() %>%
  arrange(year) %>%
  drop_na(longitude, latitude, year)
survey_years <- unique(MARINe_df$year)
# plot maps
for(yr in survey_years) {
  MARINe_by_year <- MARINe_df %>% filter(year==yr)
  p <- ggplot(data = world) +
    geom_sf() +
    geom_point(data = MARINe_by_year, aes(x = longitude, y = latitude), size = 2.5) +
    xlab("Longitude") + ylab("Latitude") +
    ggtitle(paste("MARINe survey sites", yr)) +
    coord_sf(xlim = c(-110, -130), ylim = c(25, 45), expand = FALSE)
  
  print(p)
  ggsave(filename = paste0("MARINe_Sites_map_", yr, ".png"), 
         plot = p)
}

# edna only map

eDNA_site_df <- data.frame(read.csv("./data/eDNA_sites.csv", header=TRUE)) %>%
  select(c("Sites", "Latitude", "Longitude")) %>%
  mutate(Source = "eDNA")

ggplot(data = world) +
  geom_sf() +
  geom_point(data = eDNA_site_df, aes(x = Longitude, y = Latitude), size = 2.5) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("eDNA survey sites") +
  coord_sf(xlim = c(-110, -130), ylim = c(25, 45), expand = FALSE)
ggsave("eDNA_site_map.png")