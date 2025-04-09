# distance decay beta divrsity 2001-2003 Species-Genus Jaccard MARINe data
# author: Jordana Sevigny, jordana.sevigny@gmail.com
# 11/15/2024


# load libraries
library(tidyverse)
library(betapart)
library(geosphere)
library(here)
library(dplyr)

# load data
survey_taxa_dates <- data.frame(read.csv(here("data", "processed_data", "biodiversity", "marine_species_20241025_dates_20241028_merged.csv")))
# add survey ID column
survey_taxa_dates$survey_ID <- paste(survey_taxa_dates$site_code, survey_taxa_dates$survey_rep, sep="_")

# Taxanomic levels to use in order from narrowest to broadest!!
taxa <- c("Species", "Genus")

# clean the data to only use species-resolved data
survey_tax_min_dmy <- survey_taxa_dates %>%
  filter(lowest_taxonomic_resolution %in% taxa) %>%
  filter(year %in% c(2001, 2002, 2003)) %>%
  mutate(sample_date = as.Date(sample_date)) %>%
  group_by(marine_site_name) %>%
  arrange(sample_date) %>%
  filter(sample_date == min(sample_date)) %>%
  ungroup()

# make new column of combo species or genus depending on what is lowest tax. Species_lump seemed messier
survey_tax_min_dmy$low_tax_name <- ifelse(survey_tax_min_dmy$Species == "NULL", survey_tax_min_dmy$Genus, survey_tax_min_dmy$Species)

# drop all columns but survey_ID and lowest taxonomy and deduplicate (a genus issue)
survey_low_tax <- survey_tax_min_dmy %>%
  dplyr::select(c("survey_ID", "low_tax_name")) %>%
  distinct()



# add presence column then pivot and fill remaining cells with 0
presence_absence_matrix <- survey_low_tax %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = low_tax_name, values_from = presence, values_fill = 0)
presence_absence_matrix <- as.data.frame(presence_absence_matrix)

# should add a warning it that ensures the col and row len = unique species and sites/surveys
# set row names to site and remove site as a column
rownames(presence_absence_matrix) <- presence_absence_matrix$survey_ID
presence_absence_matrix <- presence_absence_matrix %>% 
  dplyr::select(-survey_ID)

# next step is to load betapart and do calculation then average each site!
# will probably need to delete the within site betas?

# make betapart object for expedited computation
betapart_object <- betapart.core(presence_absence_matrix)
# calculate pairwise jaccard and sorensen beta diversity
beta_pair_jaccard <- beta.pair(betapart_object, "jaccard")


# separate dist matrices of jaccard (site x site dimensions)
jacc_turnover <- as.matrix(beta_pair_jaccard[[1]])
jacc_nestedness <- as.matrix(beta_pair_jaccard[[2]])
jacc_totaldiss <- as.matrix(beta_pair_jaccard[[3]])

# switch matrix into long data
jacc_total_ave_df <- as.data.frame(as.table(jacc_totaldiss)) %>%
  separate(Var1, into = c("site1", "survey_rep_1"), sep = "_") %>%
  separate(Var2, into = c("site2", "survey_rep_2"), sep = "_") %>%
  dplyr::select(-c("survey_rep_1", "survey_rep_2")) %>%
  mutate(Freq = 1-Freq)

## Distances


# find the latitude for each survey site to put sites into order with latitude
# not sure i need those 2 mutate lines....
site_lon_lat <- survey_tax_min_dmy %>%
  dplyr::select(c(site_code, longitude, latitude)) %>%
  distinct() %>%
  mutate(site1 = as.character(site_code)) %>%
  mutate(site2 = as.character(site_code))


# haversine geom_dist in km
geo_dist <- distm(site_lon_lat[, c("longitude", "latitude")], fun=distHaversine) / 1000
rownames(geo_dist) <- colnames(geo_dist) <- site_lon_lat$site_code
geo_dist_df <- as.data.frame(as.table(as.matrix(geo_dist))) %>%
  rename(site1=Var1, site2=Var2, geographic_distance=Freq)

# merge jaccard diss and geo distance dataframes
distance_decay <- jacc_total_ave_df %>%
  inner_join(geo_dist_df, by = c("site1", "site2"))

# plot distance decay plot
ggplot(distance_decay, aes(x=geographic_distance, y=Freq)) +
  geom_point() +
  labs(
    title="Distance-decay of beta diversity (Jaccard Index)",
    x = "Geographic distance (km)",
    y = "Jaccard Similarity"
  )

# add colors to the points based on within bioregion pairs and between bioregion pairs
distance_decay_lat <- distance_decay %>%
  left_join(site_lon_lat, by="site1") %>%
  rename(site2 = site2.x) %>%
  dplyr::select(-c("site2.y", "site_code")) %>%
  left_join(site_lon_lat, by="site2") %>%
  rename(site1=site1.x) %>%
  dplyr::select(-c("site_code", "site1.y"))

distance_decay_lat$interregion <- 0
for(i in 1:nrow(distance_decay_lat)) {
  if (distance_decay_lat$latitude.x[i] > 34.4486 & distance_decay_lat$latitude.y[i] > 34.4486 | 
      distance_decay_lat$latitude.x[i] < 34.4486 & distance_decay_lat$latitude.y[i] < 34.4486) {
    distance_decay_lat$interregion[i] <- 0
  } else {
    distance_decay_lat$interregion[i] <- 1
  }
}

# plot distance decay plot
ggplot(distance_decay_lat, aes(x=geographic_distance, y=Freq, color=as.factor(interregion))) +
  geom_point() +
  labs(
    title="Distance-decay of beta diversity (Jaccard Index); 2001-2003 Species & Genus",
    x = "Geographic distance (km)",
    y = "Jaccard Similarity",
    color = "Inter-region pair (1); 
    Intra-region pair (0)"
  )
#ggsave("distance_decay_V1.png")

rm(list=ls())
