
# Load required packages
library(tidyverse)
library(betapart)
library(vegan)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(here)

# THE eDNA FILE IS TOO BIG FOR GIT
# Read in pre-processed eDNA data
edna <- data.frame(read.csv(here("data","processed_data", "biodiversity", "eDNA_summer2024_allprimers_allmeta.csv")))

# Find average hit number of blanks
threshold <- edna %>%
  filter(BagID=="blank") %>%
  group_by(ESVId, Location) %>%
  summarise(Total_Hits = sum(Hits, na.rm = TRUE), .groups = "drop") %>%
  summarize(Ave_Hits = mean(Total_Hits, na.rm = TRUE)) %>%
  pull(Ave_Hits)

# Filter by primer
primer1 <- '18Sv9_89'
primer2 <- 'MiFishU'
primer3 <- 'UniCOI'

edna_single_p <- edna %>%
  filter(TestId==primer1)

# Drop any samples with meta data, Group by ESVId and Location then sum Hits, drop inland location
edna_grouped <- edna_single_p %>%
  filter(!is.na(Location)) %>%
  group_by(ESVId, Location) %>%
  summarise(Total_Hits = sum(Hits, na.rm = TRUE), .groups = "drop") %>%
  filter(Total_Hits > threshold) %>% # set minimum threshold (currently it is average blank sample copy number)
  filter(Location != "coastal estuary") %>%
  filter(Location != "GuerroNegro") %>%
  filter(Location != "Catavina") %>%
  filter(Location != "PlayaCoyote")

# Create presence/absence matrix
presence_matrix <- edna_grouped %>%
  filter(Total_Hits > 1) %>%
  mutate(present = 1) %>%
  select(ESVId, Location, present) %>%
  pivot_wider(names_from = Location, values_from = present, values_fill = 0) %>%
  column_to_rownames("ESVId") %>%
  as.matrix()

# Transpose to site x ESV
site_by_esv <-  as.data.frame(t(presence_matrix))
# # After transposing, make sure rownames are locations
# rownames(site_by_esv) <- colnames(presence_matrix)

comm_data <- site_by_esv

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


#nestedness_df <- as.data.frame(as.table(nestedness_matrix))



# Define a Mode function
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

site_lat <- edna %>%
  mutate(Latitude = round(Latitude, digits = 2)) %>%
  filter(!is.na(Latitude)) %>%
  group_by(Location) %>%
  summarize(Latitude = Mode(Latitude)) %>%  # most common latitude per location
  ungroup() %>%
  filter(Location %in% rownames(nestedness_matrix)) %>%
  arrange(desc(Latitude))  # use arrange(Latitude) for south to north

# Reorder matrix
dn_sorted <- nestedness_matrix[site_lat$Location, site_lat$Location]
my_breaks <- c(seq(0, 0.24, length.out = 20),
               seq(0.25, 0.74, length.out = 10),
               seq(0.75, 1, length.out = 20))
# Plot heatmap
pheatmap(as.matrix(dn_sorted),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = FALSE,
         main = "Directional Nestedness PCE 18Sv9_89 (Y in X)",
         color = gray.colors(length(my_breaks) - 1, start = 1, end = 0),
         breaks = my_breaks)

