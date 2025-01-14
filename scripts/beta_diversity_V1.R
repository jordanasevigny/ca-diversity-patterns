# using betapart to calculate beta diversity for MARINe sites
# author: Jordana Sevigny, jordana.sevigny@gmail.com
# 10/31/2024

# load libraries
library(tidyverse)
library(betapart)

# load data
survey_taxa_dmy <- data.frame(read.csv("./data/MARINe_survey_taxa_dmy.csv", header=TRUE))

## the data needs to be cleaned (some are species, some genus, some appear to have multiple names)
## will try to write the code to format the data for betapart without being cleaned first (will use species only data first)

########### beta diversity ##########
## each species is a column
## each site is a row
## each call is a 1 or 0 depicting presence or absence

## I guess each survey will be considered a 'site' then the pairwise beta div attached to a single site could be averaged
# add a column of the unique survey ID
survey_taxa_dmy$survey_ID <- paste(survey_taxa_dmy$site_code, survey_taxa_dmy$survey_rep, sep="_")

# clean the data to only use species-resolved data
survey_species_dmy <- survey_taxa_dmy %>%
  filter(lowest_taxonomic_resolution=="Species")
# drop all columns but survey_ID and Species
survey_species <- survey_species_dmy %>%
  select(c("survey_ID", "Species"))

# add presence column then pivot and fill remaining cells with 0
presence_absence_matrix <- survey_species %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = Species, values_from = presence, values_fill = list(presence = 0))
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
beta_pair_sorensen <- beta.pair(betapart_object, "sorensen")

# separate dist matrices of jaccard (site x site dimensions)
jacc_turnover <- as.matrix(beta_pair_jaccard[[1]])
jacc_nestedness <- as.matrix(beta_pair_jaccard[[2]])
jacc_totaldiss <- as.matrix(beta_pair_jaccard[[3]])

# separate dist matrices of sorensen
sor_turnover <- as.matrix(beta_pair_sorensen[[1]])
sor_nestedness <- as.matrix(beta_pair_sorensen[[2]])
sor_totaldiss <- as.matrix(beta_pair_sorensen[[3]])


##### Find average of surveys for each site ######

# find the latitude for each survey site to put sites into order with latitude
site_lat <- survey_taxa_dmy %>%
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
  labs(x = "Sites", y = "Sites", title = "Jaccard Species Turnover Heatmap") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),
        axis.text.y = element_text(size=4))
ggsave("J_turnover.png") 

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
  scale_fill_gradient(low = "blue", high = "red", name = "Nestedness (1=full nest, 0=no nest)") +
  labs(x = "Sites", y = "Sites", title = "Jaccard Species Nestedness Heatmap") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),
        axis.text.y = element_text(size=6))
ggsave("J_nest.png") 

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
  labs(x = "Sites", y = "Sites", title = "Jaccard Composite Heatmap") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),
        axis.text.y = element_text(size=6))
ggsave("J_composite.png") 

#(4)
# sorensen turnover
sor_turnover_ave_df <- as.data.frame(as.table(sor_turnover)) %>%
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

# plot sorensen turnover
ggplot(sor_turnover_ave_df, aes(site_code_1, site_code_2, fill = avg_diss)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", name = "Replacement") +
  labs(x = "Sites", y = "Sites", title = "Sorensen Species Turnover Heatmap") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),
        axis.text.y = element_text(size=6))
ggsave("S_turnover.png") 

#(5)
# sorensen nest
sor_nest_ave_df <- as.data.frame(as.table(sor_nestedness)) %>%
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

# plot sorensen nest
ggplot(sor_nest_ave_df, aes(site_code_1, site_code_2, fill = avg_diss)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", name = "Nestedness (1=full nest, 0=no nest)") +
  labs(x = "Sites", y = "Sites", title = "Sorensen Species Nestedness Heatmap") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),
        axis.text.y = element_text(size=6))
ggsave("S_nest.png") 

# (6)
# sorensen composite
sor_total_ave_df <- as.data.frame(as.table(sor_totaldiss)) %>%
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

# plot sorensen composite
ggplot(sor_total_ave_df, aes(site_code_1, site_code_2, fill = avg_diss)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", name = "Dissimilarity") +
  labs(x = "Sites", y = "Sites", title = "Sorensen Composite Heatmap") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),
        axis.text.y = element_text(size=6))
ggsave("S_composite.png") 
