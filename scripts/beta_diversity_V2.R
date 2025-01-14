# using betapart to calculate beta diversity for MARINe sites 2001-2003 (using the earliest survey for each location in these years) 
# for Genus & Species minimum taxanomic level (genus-genus and species-species comparisons only)
# author: Jordana Sevigny, jordana.sevigny@gmail.com
# 11/14/2024


# load libraries
library(tidyverse)
library(betapart)

# load data
# MARINe filtered excludes islands & inland Baja site
survey_taxa_dmy <- data.frame(read.csv("./data/filtered_MARINe_survey_taxa_dmy.csv", header=TRUE))

## the data needs to be cleaned (some are species, some genus, some appear to have multiple names)
## this version will use species & genus and just 2001-2003

########### beta diversity ##########
## each unique sample/time species/genus is a column
## each site is a row
## each call is a 1 or 0 depicting presence or absence

## I guess each survey will be considered a 'site' then the pairwise beta div attached to a single site could be averaged
# add a column of the unique survey ID
survey_taxa_dmy$survey_ID <- paste(survey_taxa_dmy$site_code, survey_taxa_dmy$survey_rep, sep="_")

# Taxanomic levels to use in order from narrowest to broadest!!
taxa <- c("Species", "Genus")

# clean the data to only use 2001-2003 data
survey_tax_min_dmy <- survey_taxa_dmy %>%
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
  select(c("survey_ID", "low_tax_name")) %>%
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
  select(-survey_ID)

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


##### Find average of surveys for each site ######

# find the latitude for each survey site to put sites into order with latitude
site_lat <- survey_tax_min_dmy %>%
  select(c(site_code, latitude, marine_site_name)) %>%
  distinct() %>%
  mutate(site_code_1 = as.character(site_code))

#(1)
# jaccard turnover
jacc_turnover_ave_df <- as.data.frame(as.table(jacc_turnover)) %>%
  separate(Var1, into = c("site_code_1", "survey_rep_1"), sep = "_") %>%
  separate(Var2, into = c("site_code_2", "survey_rep_2"), sep = "_") %>%
  select(-c("survey_rep_1", "survey_rep_2")) %>%
  group_by(site_code_1, site_code_2) %>%
  summarize(avg_diss = mean(Freq)) %>%
  ungroup() %>%
  full_join(site_lat, by="site_code_1") %>%
  select(-site_code) %>%
  mutate(
    site_code_1 = fct_reorder(site_code_1, latitude),
    site_code_2 = factor(site_code_2, levels = levels(site_code_1)))

# PC is 6016
# plot Jaccard turnover
ggplot(jacc_turnover_ave_df, aes(site_code_1, site_code_2, fill = avg_diss)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", name = "Replacement") +
  labs(x = "Sites", y = "Sites", title = "Jaccard Species Turnover Heatmap 2001-2003 Species & Genus") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),
        axis.text.y = element_text(size=4))
ggsave("J_turnover_V2.png") 

#(2)
# jaccard nestedness
jacc_nest_ave_df <- as.data.frame(as.table(jacc_nestedness)) %>%
  separate(Var1, into = c("site_code_1", "survey_rep_1"), sep = "_") %>%
  separate(Var2, into = c("site_code_2", "survey_rep_2"), sep = "_") %>%
  select(-c("survey_rep_1", "survey_rep_2")) %>%
  group_by(site_code_1, site_code_2) %>%
  summarize(avg_diss = mean(Freq)) %>%
  ungroup() %>%
  full_join(site_lat, by="site_code_1") %>%
  select(-site_code) %>%
  mutate(
    site_code_1 = fct_reorder(site_code_1, latitude),
    site_code_2 = factor(site_code_2, levels = levels(site_code_1)))

# plot Jaccard nestedness
ggplot(jacc_nest_ave_df, aes(site_code_1, site_code_2, fill = avg_diss)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", name = "Nestedness") +
  labs(x = "Sites", y = "Sites", title = "Jaccard Species Nestedness Heatmap 2001-2003 Species & Genus") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),
        axis.text.y = element_text(size=6))
ggsave("J_nest_V2.png") 

#(3)
# jaccard composite
jacc_total_ave_df <- as.data.frame(as.table(jacc_totaldiss)) %>%
  separate(Var1, into = c("site_code_1", "survey_rep_1"), sep = "_") %>%
  separate(Var2, into = c("site_code_2", "survey_rep_2"), sep = "_") %>%
  select(-c("survey_rep_1", "survey_rep_2")) %>%
  group_by(site_code_1, site_code_2) %>%
  summarize(avg_diss = mean(Freq)) %>%
  ungroup() %>%
  full_join(site_lat, by="site_code_1") %>%
  select(-site_code) %>%
  mutate(
    site_code_1 = fct_reorder(site_code_1, latitude),
    site_code_2 = factor(site_code_2, levels = levels(site_code_1)))

# plot Jaccard composite
ggplot(jacc_total_ave_df, aes(site_code_1, site_code_2, fill = avg_diss)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", name = "Dissimilarity") +
  labs(x = "Sites", y = "Sites", title = "Jaccard Composite Heatmap 2001-2003 Species & Genus") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),
        axis.text.y = element_text(size=6))
ggsave("J_composite_V2.png") 











#(3)
# jaccard composite but force the site ID's into numerical integer order
jacc_total_ave_df <- as.data.frame(as.table(jacc_totaldiss)) %>%
  separate(Var1, into = c("site_code_1", "survey_rep_1"), sep = "_") %>%
  separate(Var2, into = c("site_code_2", "survey_rep_2"), sep = "_") %>%
  select(-c("survey_rep_1", "survey_rep_2")) %>%
  group_by(site_code_1, site_code_2) %>%
  summarize(avg_diss = mean(Freq)) %>%
  ungroup() %>%
  mutate(avg_diss = 1-avg_diss) %>%
  full_join(site_lat, by="site_code_1") %>%
  select(-site_code) %>%
  arrange(latitude) %>%
  mutate(site_code_1 = factor(site_code_1, levels=unique(site_code_1))) %>%
  group_by(site_code_1) %>%
  mutate(site_code_1ID = as.integer(site_code_1)) %>%
  ungroup()

key <- jacc_total_ave_df %>%
  select(c('site_code_1', 'site_code_1ID')) %>%
  distinct()
key$site_code_1 <- as.character(key$site_code_1)
jacc_total_ave_df$site_code_2ID <- NA
for (i in 1:nrow(jacc_total_ave_df)) {
  for (k in 1:nrow(key)) {
    if (jacc_total_ave_df$site_code_2[i] == key$site_code_1[k]){
      jacc_total_ave_df$site_code_2ID[i] <- key$site_code_1ID[k]
    }
  }
}
jacc_total_ave_df

# sites_ids <- jacc_total_ave_df %>%
#   select(c("marine_site_name", "site_code_1ID")) %>%
#   distinct()
# write.csv(sites_ids, "marine_site_to_ID_key.csv", row.names=FALSE)
#   
# plot Jaccard composite
ggplot(jacc_total_ave_df, aes(site_code_1ID, site_code_2ID, fill = avg_diss)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", name = "Similarity") +
  labs(x = "Sites", y = "Sites", title = "Jaccard Composite Heatmap 2001-2003 Species & Genus") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),
        axis.text.y = element_text(size=6)) +
  geom_vline(xintercept = 6.5) +
  geom_vline(xintercept = 20, linetype="dotdash") +
  geom_hline(yintercept = 6.5) +
  geom_hline(yintercept = 20, linetype="dotdash")

ggsave("J_composite_V2_filter_sites.png") 
