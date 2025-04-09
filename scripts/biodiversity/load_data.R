# load in MARINe datasets and join them
# date: 10.28.2024
# author: Jordana Sevigny
# contact: jordana.sevigny@gmail.com

# libraries
library(tidyverse)
library(here)

# read in species presence data from MARINe surveys
survey_taxa <- data.frame(read.csv(here("data", "marine_data", "jordana_sevigny_cbs_species_list_20241025.csv"), header=TRUE))

# read in MARINe survey dates & remove replicate columns
survey_dates <- data.frame(read.csv(here("data", "marine_data", "jordana_sevigny_cbs_survey_dates_20241028.csv"), header=TRUE)) %>%
  select(-c("marine_sort_order", "year", "pc_point_type"))

# merge datasets - add date column to survey_taxa (dmy = day, month, year)
survey_taxa_dates <- merge(survey_taxa, survey_dates, by=c("marine_site_name", "survey_rep"), all=TRUE) %>%
  rename("sample_date" = "Min.sample_date.")

# write the dataset to data folder
write.csv(x=survey_taxa_dates, file=here("data", "processed_data", "biodiversity", "marine_species_20241025_dates_20241028_merged.csv"))


# Trial loading in the merged dataframe
survey_taxa_dates <- data.frame(read.csv(here("data", "processed_data", "biodiversity", "marine_species_20241025_dates_20241028_merged.csv")))
# check for empty dates
any(is.na(survey_taxa_dates$sample_date))
# no missing dates

# barplot of taxanomic resolution across all identifications for all surveys
ggplot(data=survey_taxa_dates, aes(x=lowest_taxonomic_resolution)) +
  geom_bar()
#ggsave("taxa_resolution_count.png")
# could do everything at genus level... 
