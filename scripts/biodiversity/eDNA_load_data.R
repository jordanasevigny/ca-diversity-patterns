# Script to load in eDNA data
# By Jordana Sevigny: jordana.sevigny@gmail.com
# April 22, 2025


# load libraries
library(tidyverse)
library(betapart)
library(here)


# load metadata
meta_df <- data.frame(read.csv(here("data", "edna_jv_data", "summer_2024", "Summer_Exped_And_Central_Coast_All_Metadata_2024_JKS.csv"))) %>%
  select(-X.)
#meta_fall_df <- data.frame(read.csv(here("data", "edna_jv_data", "fall_2023", "JVB2844-samples.csv")))


# Define cleaning function
clean_primer_df <- function(file_path, meta_data, start_col = 13) {
  df <- data.frame(read.csv(file_path))
  
  # Get original colnames starting from start_col
  raw_colnames <- colnames(df)[start_col:ncol(df)]
  
  # Remove ".0" or ".1"
  final_names <- gsub("\\.(0|1)$", "", raw_colnames)
  colnames(df)[start_col:ncol(df)] <- final_names
  
  return(df)
}

# Apply to each primer file
p18Sv9_df <- clean_primer_df(here("data", "edna_jv_data", "summer_2024", "JVB3735-18Sv9_89-read-data.csv"), meta_df)
pMiFishU_df <- clean_primer_df(here("data", "edna_jv_data", "summer_2024", "JVB3735-MiFishU-read-data.csv"), meta_df)
pUniCOI_df <- clean_primer_df(here("data", "edna_jv_data", "summer_2024", "JVB3735-UniCOI-read-data.csv"), meta_df)

# Combine the three dataframes by column name
combined_df <- bind_rows(p18Sv9_df, pMiFishU_df, pUniCOI_df)

# Pivot into long format
df_long <- combined_df %>%
  pivot_longer(
    cols = 13:ncol(.),          # Columns with sample IDs and hits
    names_to = "SampleID",     # New column name for sample identifiers
    values_to = "Hits")         # New column name for counts/hits


df_joined <- df_long %>%
  mutate(SampleID_clean = ifelse(
    sub("^X", "", SampleID) %in% meta_df$SampleID,
    sub("^X", "", SampleID),
    SampleID
  )) %>%
  filter(SampleID_clean != "X") %>%  # remove rows where ID is just "X" (there was an extraneous column)
  left_join(meta_df, by = c("SampleID_clean" = "SampleID"))

# THERE APPEARS TO BE ADDITIONAL SAMPLES IN THE EDNA THAT HAVE NO META DATA
df_joined <- df_joined %>%
  select(-SampleID) %>%
  rename(SampleID = SampleID_clean)

# setdiff(meta_df$SampleID, df_joined$SampleID_clean)
# meta_df$SampleID
# df_joined$SampleID_clean
# unique(df_joined$SampleID[is.na(df_joined$Location) | df_joined$Location == ""])
# any(df_joined$SampleID=="ZOTLLVPY")
# unique(df_joined$Location)
# Save as a new dataframe for further analysis
# THE eDNA FILE IS TOO BIG FOR GIT
write.csv(df_joined, here("data","processed_data", "biodiversity", "eDNA_summer2024_allprimers_allmeta.csv"), row.names = FALSE)
