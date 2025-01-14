

# load libraries
library(tidyverse)
library(betapart)

# load data
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
  ungroup() %>%
  filter(latitude > 32.5)

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
  separate(Var1, into = c("site_code_1", "survey_rep_1"), sep = "_") %>%
  separate(Var2, into = c("site_code_2", "survey_rep_2"), sep = "_") %>%
  select(-c("survey_rep_1", "survey_rep_2")) %>%
#  group_by(site_code_1, site_code_2) %>%
#  summarize(avg_diss = mean(Freq)) %>%
  ungroup() %>%
  full_join(site_lat, by="site_code_1") %>%
  select(-site_code) %>%
  arrange(latitude) %>%
  mutate(site_code_1 = factor(site_code_1, levels=unique(site_code_1))) %>%
  group_by(site_code_1) %>%
  mutate(site_code_1ID = as.integer(site_code_1)) %>%
  ungroup()

key <- nestedness_df2 %>%
  select(c('site_code_1', 'site_code_1ID')) %>%
  distinct()
key$site_code_1 <- as.character(key$site_code_1)
nestedness_df2$site_code_2ID <- NA
for (i in 1:nrow(nestedness_df2)) {
  for (k in 1:nrow(key)) {
    if (nestedness_df2$site_code_2[i] == key$site_code_1[k]){
      nestedness_df2$site_code_2ID[i] <- key$site_code_1ID[k]
    }
  }
}
nestedness_df2


ggplot(nestedness_df2, aes(site_code_1ID, site_code_2ID, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", name = "Nestedness") +
  labs(x = "Sites (intersection / # spp. in Site X)", y = "Sites (intersection / # spp. in Site Y)", title = "Unidirectional Nestedness 2001-2003 Species & Genus above 32.5 lat") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),
        axis.text.y = element_text(size=6))


# convert back to matrix
nest3 <- nestedness_df2 %>%
  select(c("site_code_1ID", "site_code_2ID", "Freq")) %>%
  pivot_wider(names_from = site_code_1ID, values_from = Freq)

colnames(nest3) <- paste0("X", colnames(nest3))
nest3 <-  nest3 %>%
  rename(rownames_df = Xsite_code_2ID)
nest3 <- as.data.frame(nest3)
rownames(nest3) <- nest3[, 1]
nest3 <- nest3[, -1]



# load in connectivity matrix (make sure orientation is the same)
con_df <- read.csv('./data/filtered_MARINe_samp5000_2001_through_2003_count_to_from_poly.csv')

# Define the full set of row names (1 to 64)
all_rows <- 1:nrow(con_df)
con_df <- as.data.frame(con_df)
rownames(con_df) <- con_df[, 1]
rownames_df <- con_df$polygon_id_origin

con_df <- con_df %>%
  rename(rownames_df = polygon_id_origin)
# Add missing rows and fill all columns with 0
con_df_w <- con_df %>%
  complete(rownames_df = all_rows, fill = as.list(rep(0, ncol(con_df) - 1)))

# Add missing columns

all_columns <- paste0("X", 1:64)  # Define full column names
missing_cols <- setdiff(all_columns, colnames(con_df_w))  # Find missing columns
con_df_w[missing_cols] <- 0  # Add missing columns with 0


# Reorder the columns for consistency
con_df_w <- con_df_w %>%
  select(rownames_df, sort(names(con_df_w)[-1]))
con_df_w <- con_df_w[, -1]
numeric_order <- as.numeric(gsub("X", "", colnames(con_df_w)))
con_df_w <- con_df_w[, order(numeric_order)]
con_df_w[is.na(con_df_w)] <- 0

con_df_w <- as.data.frame(con_df_w)


# combine dataframes
# Convert the tables to data frames in long format

nest3 <- as.matrix(nest3)
con_df_w <- as.matrix(con_df_w)

table1_long <- as.data.frame(as.table(nest3))
colnames(table1_long) <- c("Row", "Column", "Value1")
table2_long <- as.data.frame(as.table(con_df_w))
table2_long$Row <- as.numeric(table2_long$Var1)

table2_long <- table2_long[, c("Row", "Var2", "Freq")]
colnames(table2_long) <- c("Row", "Column", "Value2")

# Merge the two tables based on Row and Column
combined_df <- merge(table1_long, table2_long, by = c("Row", "Column"))


ggplot(combined_df, aes(x = Value2, y = Value1)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x="connectivity", y="nestedness (1 being high)", title = "connectivity vs nestedness")


