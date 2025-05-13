# Betadiversity heatmap for eDNA ESVs per primer
# By Jordana Sevigny
# 04/24/25


# Load required packages
library(tidyverse)
library(betapart)
library(vegan)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(here)

# Read in pre-processed eDNA data
edna <- data.frame(read.csv(here("data","processed_data", "biodiversity", "eDNA_summer204_allprimers_allmeta.csv")))

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

# Drop any samples with meta data, Group by ESVId and Location then sum Hits
edna_grouped <- edna_single_p %>%
  filter(!is.na(Location)) %>%
  group_by(ESVId, Location) %>%
  summarise(Total_Hits = sum(Hits, na.rm = TRUE), .groups = "drop") %>%
  filter(Total_Hits > threshold) # set minimum threshold (currently it is average blank sample copy number)

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
# After transposing, make sure rownames are locations
rownames(site_by_esv) <- colnames(presence_matrix)

# Compute beta diversity
beta_obj <- betapart.core(site_by_esv)
#beta_sor <- beta.pair(beta_obj)$beta.sor
beta_jac <- beta.pair(beta_obj, index.family = "jaccard")$beta.jac

# Convert dist to full matrix
beta_matrix <- as.matrix(beta_jac)

# Now assign row/col names from your site_by_esv
site_names <- rownames(site_by_esv)
rownames(beta_matrix) <- site_names
colnames(beta_matrix) <- site_names


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
  filter(Location %in% rownames(beta_matrix)) %>%
  arrange(desc(Latitude))  # use arrange(Latitude) for south to north

# Reorder matrix
beta_sor_sorted <- beta_matrix[site_lat$Location, site_lat$Location]

# Step 6: Plot heatmap
pheatmap(as.matrix(beta_sor_sorted),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = FALSE,
         main = "Jaccard Beta Diversity (Sorted by Latitude) for Primer 18Sv9_89 - Dissimilarity",
         color = viridis(100))

# Plot to more closely match MARINe plotting
# Calculate similarity
similarity_matrix <- 1 - as.matrix(beta_sor_sorted)

# Convert matrix to long data frame
similarity_df <- as.data.frame(similarity_matrix) %>%
  rownames_to_column(var = "Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "Freq")

# Now you can use your plotting code directly!
ggplot(similarity_df, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", name = "Similarity") +
  labs(x = "Sites", y = "Sites", title = "Jaccard Beta Diversity (Sorted by Latitude) for Primer UniCOI") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6))
