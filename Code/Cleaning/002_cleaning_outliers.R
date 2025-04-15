# Script to validate fluorescence data by:
#  * Checking for outliers
#  * Remove outliers


# Looking for outliers (based on four technical replicates)
# Inspiration on how to do the procedure comes from a blog-post from Jorik Bot, 
# https://jorikbot.com/post/an-r-function-to-removing-qpcr-outliers-within-technical-replicates/ [last update Sep 15, 2020, accessed May 3 2024]

# Load packages ####
library(tidyverse)

# Criteria for removing technical replicates
# * A  replicate will be removed if the distance between its value and the median value of all technical replicates is greater than a set threshold. Here, I use the standard deviation (sd) * 1.5 as the threshold. 
# * All technical replicates will be removed if all values are more than a set threshold apart (too much variation) 
#THIS STILL NEEDS TO BE SET UP!! #####

# Read dataset ####
data_raw <- read_csv("Data/Clean/merged_raw.csv")


# Check for outliers ####
descriptive_stats <- data_raw |> 
  #group data by sample, enzyme and type
  group_by(sample_ID, enzyme, type) |> 
  summarise(median_fluo = median(fluorescense, na.rm = TRUE),
            sd = sd(fluorescense, na.rm = TRUE))

dev_test <- data_raw |> 
  full_join(descriptive_stats, by = c("sample_ID", "enzyme", "type")) |> 
  # calculate absolute distance to mean
  mutate(diff_median = abs(fluorescense - median_fluo)) |>
  # calculate threshold value (relationship to sd can be changed)
  mutate(threshold = sd * 1.5) |> 
  # check whether replicate falls within the threshold value
  mutate(keep = if_else(diff_median < threshold, TRUE, FALSE))

# Count number of TRUEs in keep column (i.e. samples that pass QC)
count_true <- dev_test |> 
  group_by(sample_ID, enzyme, type) |> 
  summarise(count_keep = sum(keep, na.rm = TRUE)) |> 
  ungroup()

# Filter out all outliers
clean_data <- dev_test |> 
  full_join(count_true, by = c("sample_ID", "enzyme", "type")) |> 
  # keeps technical replicates that pass QC and samples where at least 2 replicates pass QC
  filter(keep == TRUE, count_keep > 1) |> 
  # remove unnecessary columns
  select(-(median_fluo:count_keep))

# Export cleaned dataset
write.csv(clean_data, "Data/Clean/clean_QC.csv", row.names = FALSE)





