library(tidyverse)
library(here)
metadat <- read_csv(here("./data/eDNA JV processed/Summer_Exped_And_Central_Coast_All_Metadata_2024.csv"))
dat <- read_csv(here("./data/eDNA JV processed/JVB2844-16S-read-data.csv"))

metadat_samples <- unique(metadat$SampleID)
dat_samples <- dat |> 
  pivot_longer(cols = 13:ncol(dat), names_to = "SampleID") |> 
  separate(SampleID, c("SampleID", "Subsample")) |> 
  pull(SampleID) |> 
  unique()

intersect(dat_samples, metadat_samples)