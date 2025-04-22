# Combines directional nestedness and connectivity
# Jordana Sevigny
# Created: 12-2024 (around there - forgot to put this in)


# load libraries
library(tidyverse)
library(betapart)
library(here)
library(patchwork)

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

# clean the data to only use 2001-2003 data
survey_tax_min_dmy <- survey_taxa_dates %>%
  filter(lowest_taxonomic_resolution %in% taxa) %>%  # filter for only the rows where the desired taxa are available
  filter(year %in% c(2001, 2002, 2003)) %>% # filter for the years we care about
  mutate(sample_date = as.Date(sample_date)) %>% # reformat date
  group_by(marine_site_name) %>%# group data but location
  arrange(sample_date) %>%# reorganize each location's data by sample date
  filter(sample_date == min(sample_date)) %>% # select the earliest sample
  ungroup() # ungroup
  # %>% filter(latitude > 32.5) # input a latitude cutoff

# # clean the data to use all but post 2020 data
# survey_tax_min_dmy <- survey_taxa_dates %>%
#   filter(lowest_taxonomic_resolution %in% taxa) %>%
#   filter(year <= 2020) %>%
#   mutate(sample_date = as.Date(sample_date)) %>%
#   group_by(marine_site_name) %>%
#   arrange(sample_date) %>%
#   filter(sample_date == min(sample_date)) %>%
#   ungroup() %>%
#   filter(latitude > 32.5)

# Make new column of combo species or genus depending on what is lowest tax. Species_lump seemed messier
survey_tax_min_dmy$low_tax_name <- ifelse(survey_tax_min_dmy$Species == "NULL", survey_tax_min_dmy$Genus, survey_tax_min_dmy$Species)

# Load in the site ID key to pair the nestedness with connectivity
marine_site_id <- data.frame(read.csv(here("data", "processed_data", "pre_QE_archive", "MARINe_site_names_ids.csv"))) %>%
  select(-X)

survey_tax_min_dmy_siteid <- left_join(survey_tax_min_dmy, marine_site_id, by="marine_site_name") %>%
  drop_na(site_id) # drop data from anywhere without an id (e.g. channel islands)

# Drop all columns but site_ID and lowest taxonomy and deduplicate (a genus issue)
survey_low_tax <- survey_tax_min_dmy_siteid %>%
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
# 
# 
# ggplot(nestedness_df2, aes(Var1, Var2, fill = Freq)) +
#   geom_tile() +
#   scale_fill_gradient(low = "blue", high = "red", name = "Nestedness") +
#   labs(x = "Sites (intersection / # spp. in Site X)", y = "Sites (intersection / # spp. in Site Y)", title = "Unidirectional Nestedness 2001-2003 Species & Genus above 32.5 lat") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),
#         axis.text.y = element_text(size=6))
# 
# 
# # convert back to matrix
# nest3 <- nestedness_df2 %>%
#   select(c("Var1", "Var2", "Freq")) %>%
#   pivot_wider(names_from = Var1, values_from = Freq)
# rownames(nest3) <- nest3$Var2

# load in connectivity matrix (make sure orientation is the same)
# filtered_MARINe_2001_through_2003_count_to_from_poly.csv
con_df2 <- read.csv('./data/connectivity_data/qe/04-07-2025_final_count_cleaned_connectivity_30-60_1_13_2001_2003.csv')

#con_df2 <- con_df[,4:ncol(con_df)]

######### make a split connectivity plot
# con_df2_symmetric <- con_df2 %>%
#   rename(X = log10_count) %>%  # Rename for clarity
#   inner_join(con_df2, 
#              by = c("origin_polygon_up" = "destination_polygon_up", 
#                     "destination_polygon_up" = "origin_polygon_up")) %>%
#   rename(Y = log10_count)  # Rename the second count column
# 
# ggplot(con_df2_symmetric, aes(x = X, y = Y)) +
#   geom_point() +
#   labs(x = "log10(Number of Floats (Origin → Destination))",
#        y = "log10(Number of Floats (Destination → Origin))",
#        title = "Symmetric Float Transitions") +
#   theme_minimal()

################ combine connectivity with nestedness

colnames(con_df2)[colnames(con_df2) == "origin_polygon_up"] <- "Var1"
colnames(con_df2)[colnames(con_df2) == "destination_polygon_up"] <- "Var2"
# Merge the two tables based on Row and Column
combined_df <- merge(nestedness_df2, con_df2, by = c("Var1", "Var2"))

lm_model <- lm(Freq ~ log10(count), data = combined_df)

# Get p-value from the summary of the linear model
summary_model <- summary(lm_model)
p_value <- summary_model$coefficients[2, 4] 
r_squared <- summary_model$r.squared 

p1 <- ggplot(combined_df, aes(x = log10(count), y = Freq)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x="Connectivity (log10[count])", 
    y="Nestedness (1 being more nested)", 
    title = "Nestedness (MARINe) vs Connectivity (raw count) 2001-2003\nIncludes same site comparisons") +
  annotate("text", x = min(log10(combined_df$count), na.rm = TRUE), 
           y = max(combined_df$Freq, na.rm = TRUE), 
           label = paste0("R² = ", round(r_squared, 2), 
                          ", p = ", format.pval(p_value, digits = 3, eps = 0.001)),
           hjust = 0, vjust = 1)


# Drop rows where site A is compared to site A (100% sim)
combined_df2 <- combined_df %>%
  filter(Var1!=Var2)

lm_model2 <- lm(Freq ~ log10(count), data = combined_df2)

# Get p-value from the summary of the linear model
summary_model2 <- summary(lm_model2)
p_value2 <- summary_model2$coefficients[2, 4] 
r_squared2 <- summary_model2$r.squared 

p2 <- ggplot(combined_df2, aes(x = log10(count), y = Freq)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x="Connectivity (log10[count])", 
    y="Nestedness (1 being more nested)", 
    title = "Nestedness (MARINe) vs Connectivity (raw count) 2001-2003\nExcludes same site comparisons") +
  annotate("text", x = min(log10(combined_df_norm2$count), na.rm = TRUE), 
           y = max(combined_df2$Freq, na.rm = TRUE), 
           label = paste0("R² = ", round(r_squared2, 2), 
                          ", p = ", format.pval(p_value2, digits = 3, eps = 0.001)),
           hjust = 0, vjust = 1)
# 
# combined_df$count <- 10^(combined_df$log10_count)
# 
# ggplot(combined_df, aes(x = count, y = Freq)) +
#   geom_point() +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(x="Connectivity (count)", y="Nestedness (1 being high)", title = "connectivity vs nestedness")
# 

#####################
# Normalized connectivity by release X Nestedness

con_release_counts <- read.csv('./data/connectivity_data/qe/04-07-2025_release_counts_connectivity_30-60_1_13_2001_2003.csv')

# X0 column is release counts by release site
con_release_counts$Var1 <- sort(unlist(unique(con_df2$Var1))) # Var 1 here is release site ID - why I named it this... beats me
con_df2_norm <- left_join(con_df2, con_release_counts, by="Var1")
con_df2_norm$count_norm <- con_df2_norm$count / con_df2_norm$X0

# Merge the two tables based on Row and Column
combined_df_norm <- merge(nestedness_df2, con_df2_norm, by = c("Var1", "Var2"))

lm_model_norm <- lm(Freq ~ log10(count_norm), data = combined_df_norm)

# Get p-value from the summary of the linear model
summary_model_norm <- summary(lm_model_norm)
p_value_norm <- summary_model_norm$coefficients[2, 4] 
r_squared_norm <- summary_model_norm$r.squared 

#plot(combined_df_norm$count_norm, combined_df_norm$Freq)

p3 <- ggplot(combined_df_norm, aes(x = log10(count_norm), y = Freq)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x="Connectivity (log10[count/total floats released by release cell])", 
    y="Nestedness (1 being more nested)", 
    title = "Nestedness (MARINe) vs Connectivity (normalized count) 2001-2003\nIncludes same site comparisons") +
  annotate("text", x = min(log10(combined_df_norm$count_norm), na.rm = TRUE), 
           y = max(combined_df_norm$Freq, na.rm = TRUE), 
           label = paste0("R² = ", round(r_squared_norm, 2), 
                          ", p = ", format.pval(p_value_norm, digits = 3, eps = 0.001)),
           hjust = 0, vjust = 1)

# Drop rows where site A is compared to site A (100% sim)
combined_df_norm2 <- combined_df_norm %>%
  filter(Var1!=Var2)

lm_model_norm2 <- lm(Freq ~ log10(count_norm), data = combined_df_norm2)

# Get p-value from the summary of the linear model
summary_model_norm2 <- summary(lm_model_norm2)
p_value_norm2 <- summary_model_norm2$coefficients[2, 4] 
r_squared_norm2 <- summary_model_norm2$r.squared 

p4 <- ggplot(combined_df_norm2, aes(x = log10(count_norm), y = Freq)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x="Connectivity (log10[count/total floats released by release cell])", 
    y="Nestedness (1 being more nested)", 
    title = "Nestedness (MARINe) vs Connectivity (normalized count) 2001-2003\nExcludes same site comparisons") +
  annotate("text", x = min(log10(combined_df_norm2$count_norm), na.rm = TRUE), 
           y = max(combined_df_norm2$Freq, na.rm = TRUE), 
           label = paste0("R² = ", round(r_squared_norm2, 2), 
                          ", p = ", format.pval(p_value_norm2, digits = 3, eps = 0.001)),
           hjust = 0, vjust = 1)

(p1 | p2) / (p3 | p4)
