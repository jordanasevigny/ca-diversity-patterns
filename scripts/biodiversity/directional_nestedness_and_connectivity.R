# Combines directional nestedness and connectivity
# Jordana Sevigny
# Created: 12-2024 (around there - forgot to put this in)


# load libraries
library(tidyverse)
library(betapart)

# load data
survey_taxa_dates <- data.frame(read.csv(here("data", "processed_data", "biodiversity", "marine_species_20241025_dates_20241028_merged.csv")))

## the data needs to be cleaned (some are species, some genus, some appear to have multiple names)
## this version will use species & genus and just 2001-2003

########### beta diversity ##########
## each unique sample/time species/genus is a column
## each site is a row
## each call is a 1 or 0 depicting presence or absence

## I guess each survey will be considered a 'site' then the pairwise beta div attached to a single site could be averaged
# add a column of the unique survey ID
survey_taxa_dates$survey_ID <- paste(survey_taxa_dates$site_code, survey_taxa_dates$survey_rep, sep="_")

# Taxanomic levels to use in order from narrowest to broadest!!
taxa <- c("Species", "Genus")

# # clean the data to only use 2001-2003 data
# survey_tax_min_dmy <- survey_taxa_dates %>%
#   filter(lowest_taxonomic_resolution %in% taxa) %>%
#   filter(year %in% c(2001, 2002, 2003)) %>%
#   mutate(sample_date = as.Date(sample_date)) %>%
#   group_by(marine_site_name) %>%
#   arrange(sample_date) %>%
#   filter(sample_date == min(sample_date)) %>%
#   ungroup() %>%
#   filter(latitude > 32.5)

# clean the data to use all but post 2020 data
survey_tax_min_dmy <- survey_taxa_dates %>%
  filter(lowest_taxonomic_resolution %in% taxa) %>%
  filter(year <= 2020) %>%
  mutate(sample_date = as.Date(sample_date)) %>%
  group_by(marine_site_name) %>%
  arrange(sample_date) %>%
  filter(sample_date == min(sample_date)) %>%
  ungroup() %>%
  filter(latitude > 32.5)

# make new column of combo species or genus depending on what is lowest tax. Species_lump seemed messier
survey_tax_min_dmy$low_tax_name <- ifelse(survey_tax_min_dmy$Species == "NULL", survey_tax_min_dmy$Genus, survey_tax_min_dmy$Species)

# drop all columns but site_ID and lowest taxonomy and deduplicate (a genus issue)
survey_low_tax <- survey_tax_min_dmy %>%
  select(c("site_id", "low_tax_name")) %>%
  distinct()

# add presence column then pivot and fill remaining cells with 0
presence_absence_matrix <- survey_low_tax %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = low_tax_name, values_from = presence, values_fill = 0)
presence_absence_matrix <- as.data.frame(presence_absence_matrix)

# should add a warning it that ensures the col and row len = unique species and sites/surveys
# set row names to site and remove site as a column
rownames(presence_absence_matrix) <- presence_absence_matrix$site_id
presence_absence_matrix <- presence_absence_matrix %>% 
  select(-site_id)

comm_data <- presence_absence_matrix
# Initialize a matrix to store directional nestedness values
nestedness_matrix <- matrix(0, nrow = nrow(comm_data), ncol = nrow(comm_data),
                            dimnames = list(rownames(comm_data), rownames(comm_data)))



# Loop through pairs of sites to calculate directional nestedness
for (site_a in 1:nrow(comm_data)) {
  for (site_b in 1:nrow(comm_data)) {
    # Intersection: Number of shared species between the two sites
    intersection <- sum(comm_data[site_a,] & comm_data[site_b,])
    
    # Number of species in each site (denominator for the normalization)
    species_site_a <- sum(comm_data[site_a,])
    species_site_b <- sum(comm_data[site_b,])
    
    # Directional nestedness: Intersection divided by the number of species in each site
    directional_nestedness <- intersection / species_site_a  # For Site A's perspective
    # Alternatively, you can calculate it as intersection / species_site_b for Site B's perspective
    
    # Store the value in the matrix
    nestedness_matrix[site_a, site_b] <- directional_nestedness
    
  }
}


nestedness_df <- as.data.frame(as.table(nestedness_matrix))


# find the latitude for each survey site to put sites into order with latitude
site_lat <- survey_tax_min_dmy %>%
  select(c(site_code, latitude, marine_site_name)) %>%
  distinct() %>%
  mutate(site_code_1 = as.character(site_code))


# force the site ID's into numerical integer order
nestedness_df2 <- nestedness_df %>%
  mutate(Var1 = factor(Var1, levels = sort(unique(as.numeric(as.character(Var1)))))) %>%
  mutate(Var2 = factor(Var2, levels = sort(unique(as.numeric(as.character(Var2))))))


ggplot(nestedness_df2, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", name = "Nestedness") +
  labs(x = "Sites (intersection / # spp. in Site X)", y = "Sites (intersection / # spp. in Site Y)", title = "Unidirectional Nestedness 2001-2003 Species & Genus above 32.5 lat") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),
        axis.text.y = element_text(size=6))


# convert back to matrix
nest3 <- nestedness_df2 %>%
  select(c("Var1", "Var2", "Freq")) %>%
  pivot_wider(names_from = Var1, values_from = Freq)
rownames(nest3) <- nest3$Var2

# load in connectivity matrix (make sure orientation is the same)
con_df <- read.csv('./data/preliminary_processed_data/filtered_MARINe_2001_through_2003_count_to_from_poly.csv')

con_df2 <- con_df[,4:ncol(con_df)]

######### make a split connectivity plot
con_df2_symmetric <- con_df2 %>%
  rename(X = log10_count) %>%  # Rename for clarity
  inner_join(con_df2, 
             by = c("origin_polygon_up" = "destination_polygon_up", 
                    "destination_polygon_up" = "origin_polygon_up")) %>%
  rename(Y = log10_count)  # Rename the second count column

ggplot(con_df2_symmetric, aes(x = X, y = Y)) +
  geom_point() +
  labs(x = "log10(Number of Floats (Origin → Destination))",
       y = "log10(Number of Floats (Destination → Origin))",
       title = "Symmetric Float Transitions") +
  theme_minimal()

################ combine connectivity with nestedness

colnames(con_df2)[colnames(con_df2) == "origin_polygon_up"] <- "Var1"
colnames(con_df2)[colnames(con_df2) == "destination_polygon_up"] <- "Var2"
# Merge the two tables based on Row and Column
combined_df <- merge(nestedness_df2, con_df2, by = c("Var1", "Var2"))


lm_model <- lm(Freq ~ log10_count, data = combined_df)

# Get p-value from the summary of the linear model
summary_model <- summary(lm_model)
p_value <- summary_model$coefficients[2, 4] 

ggplot(combined_df, aes(x = log10_count, y = Freq)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x="connectivity (log10 of count)", y="nestedness (1 being high)", title = "connectivity vs nestedness")



combined_df$count <- 10^(combined_df$log10_count)

ggplot(combined_df, aes(x = count, y = Freq)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x="connectivity (count)", y="nestedness (1 being high)", title = "connectivity vs nestedness")

