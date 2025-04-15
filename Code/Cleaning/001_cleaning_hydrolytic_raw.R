# Script to clean raw data for hydrolytic enzyme activities

# Load packages
library(tidyverse) # ver 2.0.0
library(writexl) # ver 4.3.1
library(readxl)

#########################################################################################
#### Read raw data ####
# Loop through all raw data files from TecanSpark. They will in the end be merged into one large dataset. 

# Step 1: Set up file paths
# Set the path to the folder containing the .xlsx raw-data files
folder_path <- "./Data/Raw/Raw_platereader" # this is where the raw data from the platereader is stored
folder_path_clean <- "./Data/Clean/Intermediate" # this is where the files end up after the initial data cleaning step

# Step 2: Generate a list of all raw data files in "folder path"
xlsx_files <- list.files(folder_path, pattern = "\\.xlsx$", full.names = TRUE)

# Step 3: make lists to produce plate-layout with 12 columns and 8 rows (96 well plate)
cols <- c(1:12)
rows <- LETTERS[1:8]

# Step 4: Loop through each raw-data file 
for (file_path in xlsx_files) {
  # Load the .xlsx file as a data frame
  df <- read_excel(file_path, range = "B45:M52", col_names = FALSE)
  
  # Modifications of data
  # Change column and row names and transpose data to long format.
  colnames(df)[1:12] = cols 
  df$row <- rows
  df <- df |> 
    pivot_longer(cols = 1:12, names_to = "num", values_to = "fluo") |> 
    unite("well_ID", c(num, row))
  
  # Add column with unique dataset name to enable downstream merging with metadata
  df <- df %>% mutate(file_name = basename(file_path))
  
  # Generate a new file name for the cleaned df
  new_file_name <- sub(".xlsx$", "_clean.csv", basename(file_path))
  
  # New file path for cleaned data
#  file_path_clean <- file.path("./Data/Clean/Intermediate", new_file_name)
  
  # Save the df as a CSV file in the folder Clean/Intermediate 
  write.csv(df, file.path(folder_path_clean, new_file_name), row.names = FALSE)
}

# Merge all clean data in one long-format file
# Extract file names from first cleaning step
clean_file_list <- list.files(path = "./Data/Clean/Intermediate", full.names = TRUE)
# Make tibble with individual datasets
data_list <- map(clean_file_list, read_csv)
# Bind individual datasets
all_cleaned_step1 <- bind_rows(data_list)
# Clean file column to only keep relevant ID-information
all_cleaned_step1$file_name <- substr(all_cleaned_step1$file_name, nchar(all_cleaned_step1$file_name) - 19, nchar(all_cleaned_step1$file_name) - 5)

fluo_data <- all_cleaned_step1 |> 
  separate(file_name, into = c("date", "runID"), sep = "_", remove = FALSE)
  
#### Read sample placement ####
# This code extracts informaion about the plate layouts, and which samples are located in which plate position. 

# Set the path to the folder containing the .xlsx files
folder_path2 <- "./Data/Raw/Plate_setups"
folder_path_clean2 <- "./Data/Clean/Intermediate_setup"

# Generate a list of all plate setups
xlsx_files <- list.files(folder_path2, pattern = "\\.xlsx$", full.names = TRUE)

for (file_path in xlsx_files) {
  # Load the .xlsx file as a data frame
  df <- read_excel(file_path, range = "B2:M9", col_names = FALSE)
  
  # Modifications of data
  # Change column and row names and transpose data to long format.
  colnames(df)[1:12] = cols 
  df$row <- rows
  df <- df |> 
    pivot_longer(cols = 1:12, names_to = "num", values_to = "sample_ID") |> 
    unite("well_ID", c(num, row))
  
  # Add column with origin of dataset to make merge with sample ID easy
  df <- df %>% mutate(file_name = basename(file_path))
  
  # Generate a new file name for the cleaned df
  new_file_name <- sub(".xlsx$", "_clean.csv", basename(file_path))
  
  # New file path for cleaned data
  #file_path_clean <- file.path("./Data/Clean/Intermediate_setup", new_file_name)
  
  # Save the df as a CSV file in a new folder "Cleaned_data"
  write.csv(df, file.path(folder_path_clean2, new_file_name), row.names = FALSE)
}

# Merge all clean sample-location data in one long-format file

clean_file_list2 <- list.files(path = "./Data/Clean/Intermediate_setup", full.names = TRUE)
data_list2 <- map(clean_file_list2, read_csv)
all_cleaned2 <- bind_rows(data_list2)
# Clean file column to only keep relevant ID-information
all_cleaned2$file_name <- substr(all_cleaned2$file_name, nchar(all_cleaned2$file_name) - 11, nchar(all_cleaned2$file_name) - 5)
samples <- all_cleaned2 |> # rename df "samples"
  rename(analysis_round = file_name)
samples <- samples |>
  separate(analysis_round, into = c("round", "analysis_round"), sep = "_") |> 
  select(- round)
  


#### Analysis type ####
# Information about which enzyme and incubation length added. 
plate_content <- read_xlsx("./Data/Raw/Hydrolytic_filenames.xlsx")

plate_content <- plate_content |> 
  separate(file_name, into = c("date", "runID"), remove = FALSE) 

#### Merge all data ####
# Merges sample data, analysis types and fluorescence data into one df

merged_data <- full_join(fluo_data, plate_content, by = "runID")
merged_data$analysis_round <- as.character(merged_data$analysis_round)
merged_data <- full_join(merged_data, samples, by = c("well_ID", "analysis_round")) |> 
  rename(date = date.x) |> 
  select(-date.y) |>
  rename(file_name = file_name.x) |> 
  select(-file_name.y) |> 
  relocate(fluorescense = fluo, .after = last_col())
         

# Export merged data as .csv file
write.csv(merged_data, "Data/Clean/merged_raw.csv", row.names = FALSE)
