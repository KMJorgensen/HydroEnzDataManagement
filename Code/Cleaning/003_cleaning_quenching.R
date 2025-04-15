# Script to check degree of quenching

# Check quenching of aP, CBH and BG against MUB.
# Check quenching of LAP against MUC. 

# Load packages ####
library(tidyverse)

# Read data ####
data <- read_csv("Data/Clean/clean_QC.csv")

# Split data based on substrate
# MUB is used to assay all hydrolytic enzymes except LAP
MUB_q <- data |> 
  filter(enzyme == "MUB_quench")

MUB_q <- MUB_q |> 
  group_by(sample_ID) |> 
  summarise(mean = mean(fluorescense),
            n = n()) 

MUB_standard <- MUB_q[MUB_q$sample_ID == "MUB_ctrl", "mean"] 
MUB_standard_mean <- MUB_standard[[1,1]]

MUB_q <- MUB_q |> 
  mutate("MUB_standard" = MUB_standard_mean) |> 
  mutate("Q_ratio" = mean/MUB_standard_mean)

view(MUB_q) # View the dataset and have a look at the quenching ration (Q_ratio). If it is below 1 it means there is some quenching. Previously, we have allowed samples with a Q_ratio of >0.6 but this threshold should be discussed. If all samples pass are above the threshold, you don't have to do anything more. Samples below the threshold may have to be diluted and re-analysed. 



