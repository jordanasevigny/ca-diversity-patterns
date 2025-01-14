# plot sites by years to get a feel for sampling 
# date: 11/7/2024
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

# load in MARINe data
survey_taxa_dmy <- data.frame(read.csv("./data/MARINe_survey_taxa_dmy.csv", header=TRUE))

# load in eDNA site data
eDNA <- read.csv("./data/eDNA_sites.csv", header=TRUE)

# years - site data
years_sites <- survey_taxa_dmy %>%
  select(c(year, marine_site_name, latitude)) %>%
  distinct() %>%
  arrange(latitude) %>%
  mutate(lat_idx = as.integer(factor(latitude)))

# get all combinations of sites and years
all_combinations <- expand.grid(
  marine_site_name = unique(years_sites$marine_site_name),
  year = seq(min(years_sites$year), max(years_sites$year)))

# Mark surveyed years with 1, others with 0
all_combinations$surveyed <- 0
for (i in 1:nrow(all_combinations)) {
  for (j in 1:nrow(years_sites)) {
    if (all_combinations$marine_site_name[i] == years_sites$marine_site_name[j] & 
      all_combinations$year[i] == years_sites$year[j]) {
      all_combinations$surveyed[i] <- 1
      break
    }
    else {
    }
  }
}

# add latitude ID (smallest to largest lat) to all_combinations of site survey MARINe data
all_combinations$lat_ID <- 0
for (i in 1:nrow(all_combinations)){
  for (j in 1:nrow(years_sites)) {
    if(all_combinations$marine_site_name[i] == years_sites$marine_site_name[j]) {
      all_combinations$lat_ID[i] <- years_sites$lat_idx[j]
    }
  }
}

# add eDNA site column to MARINe site data
all_combinations$eDNA_site_samp <- 0
for (i in 1:nrow(all_combinations)) {
  for (j in 1:nrow(eDNA)) {
    if (all_combinations$marine_site_name[i] == eDNA$eDNA_MARINe_Sites[j]) {
      all_combinations$eDNA_site_samp[i] <- 2
    }
  }
}

all_combinations$e_m_surveyed <- all_combinations$surveyed + all_combinations$eDNA_site_samp
all_combinations$e_m_surveyed <- as.factor(all_combinations$e_m_surveyed)
all_combinations$surveyed <- as.factor(all_combinations$surveyed)

frames <- all_combinations %>%
  select(c(year, lat_ID, eDNA_site_samp)) %>%
  filter(eDNA_site_samp==2) %>%
  select(-eDNA_site_samp)

ggplot(all_combinations, aes(x=year, y=lat_ID)) +
  geom_tile(aes(fill=surveyed)) +
  scale_fill_manual(values=c("0"="white", "1"="pink3")) +
  geom_rect(data=frames, fill=NA, color="black",
            aes(xmin=as.numeric(year)-0.5, xmax=as.numeric(year)+0.5, ymin=lat_ID-0.5, ymax=lat_ID+0.5), size=0.25)




ggplot(all_combinations, aes(x=year, y=lat_ID, color=e_m_surveyed)) +
  geom_point(size=0.8) +
  scale_color_manual(values = c("0" = "white", "1" = "black", "2" = "white", "3" = "red")) +
  theme_minimal() +
  ggtitle("site x year surveyed (red = eDNA & MARINe site, black = MARINe only)")

year_most_eDNA_sites <- all_combinations %>%
  filter(e_m_surveyed == 3) %>%
  group_by(year) %>%
  summarize(count_3 = n()) %>%
  arrange(desc(count_3)) %>%
  slice_head(n=5)
year_most_eDNA_sites
