# need to use alexa's coastal distance code and get it matched up with the site points then use it to calculate distance between sites
# then recreate the distance decay plot with updated distances
# Example dataframes


# load libraries
library(here)
library(tidyverse)
library(sf)
library(betapart)
library(nlme)


# load data
survey_taxa_dates <- data.frame(read.csv(here("data", "processed_data", "biodiversity", "marine_species_20241025_dates_20241028_merged.csv")))
# add survey ID column
survey_taxa_dates$survey_ID <- paste(survey_taxa_dates$site_code, survey_taxa_dates$survey_rep, sep="_")

# Taxonomic levels to use in order from narrowest to broadest!!
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


# load data
coastdistdat <- readRDS(here("processed-data-alexa","westcoast_coastdistdat.rds"))

# get length to point from reference point (furthest south - I think)
get_length <- function(lon, lat, distdf) {
  tmp <- distdf %>% 
    mutate(abs.diff.x2 = abs(x-lon)^2,
           abs.diff.y2 = abs(y-lat)^2,
           abs.diff.xy = sqrt(abs.diff.x2 + abs.diff.y2
           )) %>% 
    filter(abs.diff.xy == min(abs.diff.xy)) %>% 
    dplyr::select(lengthfromhere) %>% 
    pull()
  return(tmp)
}

site_dist_from_ref <- site_lon_lat %>% 
  rowwise() %>% 
  mutate(coastdist_km = (get_length(lon=longitude, lat=latitude, distdf = coastdistdat))/1000) %>% 
  ungroup()

dist_mat <- matrix(
  nrow = nrow(site_dist_from_ref),
  ncol = nrow(site_dist_from_ref),
  dimnames = list(site_dist_from_ref$site_code, site_dist_from_ref$site_code)
)

# fill matrix with distances
for (i in 1:nrow(site_dist_from_ref)) {
  for (j in 1:nrow(site_dist_from_ref)) {
    dist_mat[i,j] <- abs(site_dist_from_ref$coastdist_km[i] - site_dist_from_ref$coastdist_km[j])
  }
}

# distance matrix into table format
geo_dist_df <- as.data.frame(as.table(as.matrix(dist_mat))) %>%
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
    distance_decay_lat$interregion[i] <- "intra-region"
  } else {
    distance_decay_lat$interregion[i] <- "inter-region"
  }
}

distance_decay_lat$interregion <- as.factor(distance_decay_lat$interregion)

# plot distance decay plot
ggplot(distance_decay_lat, aes(x=geographic_distance, y=Freq, color=interregion)) +
  geom_point() +
  labs(
    title="Distance-decay of beta diversity (Jaccard Index); 2001-2003 Species & Genus",
    x = "Geographic distance (km)",
    y = "Jaccard Similarity"
  )
#ggsave("distance_decay_coastdist_V1.png")

# model to compare inter and intra region pair-wise comparisons
model <- lme(Freq ~ poly(geographic_distance, 1), random = ~ geographic_distance | interregion, data=distance_decay_lat)

# plot distance decay plot with model
distance_decay_lat$predicted <- predict(model)
ggplot(distance_decay_lat, aes(x=geographic_distance, y=Freq, color=interregion)) +
  geom_point() +
  geom_line(aes(y=predicted),size=1) +
  labs(
    title="Observed and predicted similarity vs. distance; 2001-2003 Species & Genus",
    x = "Geographic distance (km)",
    y = "Jaccard Similarity"
  )




rm(list=ls())