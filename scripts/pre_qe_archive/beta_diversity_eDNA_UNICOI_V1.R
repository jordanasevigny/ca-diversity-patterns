
library(tidyverse)
library(betapart)

metadat <- read.csv("./data/eDNA_JV_processed/summer_2024/Summer_Exped_And_Central_Coast_All_Metadata_2024.csv") %>%
  select(-1)
unicoi_df <- data.frame(read.csv("./data/eDNA_JV_processed/summer_2024/JVB3735-UniCOI-read-data.csv", header=TRUE)) %>%
  select(-last_col()) # empty column

# # number species is equal to number unique ESVID
# sum(unicoi_df$Species != "")
# length(unique(unicoi_df$ESVId[unicoi_df$Species != ""]))

long_data <- unicoi_df %>%
  pivot_longer(cols = 13:ncol(unicoi_df), names_to = "SampleID") %>%
  separate(SampleID, c("SampleID", "Subsample"))

# need to add meta into esv dataframe
ESV_meta_combo <- long_data %>%
  left_join(metadat, by='SampleID')
