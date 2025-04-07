# using betapart to calculate beta diversity for MARINe sites 2001-2003 (using the earliest survey for each location in these years) 
# for Genus & Species minimum taxanomic level (genus-genus and species-species comparisons only)
# author: Jordana Sevigny, jordana.sevigny@gmail.com
# 11/14/2024


# load libraries
library(tidyverse)
library(betapart)

# load data
# MARINe filtered excludes islands & inland Baja site
survey_taxa_dmy <- data.frame(read.csv("./data/preliminary_processed_data/filtered_MARINe_survey_taxa_dmy.csv", header=TRUE))


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
df2 <- survey_tax_min_dmy %>%
  select(c('marine_site_name', 'site_id')) %>%
  distinct()
# make new column of combo species or genus depending on what is lowest tax. Species_lump seemed messier
survey_tax_min_dmy$low_tax_name <- ifelse(survey_tax_min_dmy$Species == "NULL", survey_tax_min_dmy$Genus, survey_tax_min_dmy$Species)

# drop all columns but site_id and lowest taxonomy and deduplicate (a genus issue)
n_survey_low_tax <- survey_tax_min_dmy %>%
  select(c("site_id", "low_tax_name")) %>%
  distinct()

# add presence column then pivot and fill remaining cells with 0
n_presence_absence_matrix <- n_survey_low_tax %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = low_tax_name, values_from = presence, values_fill = 0)
n_presence_absence_matrix <- as.data.frame(n_presence_absence_matrix)

# should add a warning it that ensures the col and row len = unique species and sites/surveys
# set row names to site and remove site as a column
rownames(n_presence_absence_matrix) <- n_presence_absence_matrix$site_id
n_presence_absence_matrix <- n_presence_absence_matrix %>%
  select(-site_id)

# next step is to load betapart and do calculation then average each site!
# will probably need to delete the within site betas?

# make betapart object for expedited computation
n_betapart_object <- betapart.core(n_presence_absence_matrix)
# calculate pairwise jaccard and sorensen beta diversity
n_beta_pair_jaccard <- beta.pair(n_betapart_object, "jaccard")


# separate dist matrices of jaccard (site x site dimensions)
#jacc_turnover <- as.matrix(beta_pair_jaccard[[1]])
#jacc_nestedness <- as.matrix(beta_pair_jaccard[[2]])
n_jacc_totaldiss <- as.matrix(n_beta_pair_jaccard[[3]])


# jaccard composite
n_jacc_total_ave_df <- as.data.frame(as.table(n_jacc_totaldiss)) %>%
  mutate(Freq = 1-Freq) %>%
  mutate(Var1 = factor(Var1, levels = sort(unique(as.numeric(as.character(Var1)))))) %>%
  mutate(Var2 = factor(Var2, levels = sort(unique(as.numeric(as.character(Var2))))))
  

# plot Jaccard composite
ggplot(n_jacc_total_ave_df, aes(x = Var1, y = Var2, fill = Freq)) +
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
ggsave("J_composite_V2.png") 
