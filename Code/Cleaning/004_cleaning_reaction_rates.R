# Code to calculate reaction rates based on the difference between incubated and t0 samples
# Before running this script, you have to make sure that the quenching ratios are not too low

# Calculations are based on the cleaned data generated from script 002_cleaning_outliers and is called:
# clean_QC.csv (exported in the end of the outlier cleaning script)

enzyme_data <- read_csv("Data/Clean/clean_QC.csv")
soil_data <- read_xlsx("Data/Raw/soil_weight_moisture.xlsx")

# Filter data to only contain your samples (no positive controls, and no quanching tests)
data_reactions <- enzyme_data |> 
  group_by(sample_ID, runID, enzyme, type) |> 
  summarise(mean = mean(fluorescense)) |> 
  filter(sample_ID != "MUB_ctrl") |> #removes the positive control from the data
  filter(enzyme != "MUB_quench") # removes the quenching test from the data

data_inc <- data_reactions |> 
  filter(type == "inc")

data_t0 <- data_reactions |> 
  filter(type == "t0")


# Define incubation times (in hours)
# In this example aP was incubated for 15 min (time1) and CBH and BG for 45 min (time2)
time1 <- 15/60
time2 <- 45/60

enzyme_data_rates <- left_join(data_inc, data_t0, by = c("sample_ID", "enzyme")) |> #joins incubated and t0 data together
  rename(c(mean_inc = mean.x, mean_t0 = mean.y)) |> 
  select(sample_ID, enzyme, mean_inc, mean_t0) |> 
  mutate(inc_time = case_when(
    enzyme == "aP" ~ time1, #all rows where aP is the enzyme, incubation time is time1
    enzyme != "aP" ~ time2)) |> # all rows where aP is NOT the enzyme, incubation time is time2
  mutate(net_fluo = mean_inc - mean_t0) |> # calculate net fluorescence "created" during incubation
  mutate(net_fluo_h = net_fluo / inc_time) |> # express it per hour
  mutate(net_fluo_h = if_else(net_fluo_h < 0, 0, net_fluo_h))  # change negative values to 0 (i.e. no activity)

# Calculate soil concentrations in the reaction
soil_data <- soil_data |> 
  mutate(dw_sample = dw_cont * fw_g) # calculates dry weight content of the sample that you weighed in

# Calculate soil concentrations in the reaction (CHANGE BASED ON YOUR OWN DILUTIONS)
# First dilution was done with approx 5 g soil (FW) in 10 ml acetate buffer (conc1)
# Second dilution was 100µl of conc1 in 10 ml acetate buffer (total volume 10.1 ml)
# Third dilution (reaction concentration) was 200 µl sample (conc2)with 50 µl enzyme substrate and 10 µl NaOH

soil_data <- soil_data |> 
  mutate(conc1 = dw_sample / 10) |> 
  mutate(conc2 = (conc1 * 0.1)/10.1) |> 
  mutate(conc3 = (conc2 * 0.2)/0.26)

# Join reaction rates with soil data to do final conversions of fluorescence to concentrations of reaction product

joined_data <- full_join(enzyme_data_rates, soil_data, by = "sample_ID")

joined_data2 <- joined_data |> 
  select(sample_ID, enzyme, net_fluo_h, conc3) # select variables needed to finalise dataset

finalising_step1 <- joined_data2 |> 
  mutate(standard_slope = 2.2879*10^-6) |> # adds the slope from standard curve 
  mutate(fluo_gDW_h = net_fluo_h/conc3) |> # epress activities per gDW of soil
  mutate(nmol_MUB_gDW_h = fluo_gDW_h * standard_slope) # sample reaction rate: converts fluorescence to nmol MUB per gDW of soil and hour.

finalising_step2 <- finalising_step1 |> 
  select(sample_ID, enzyme, nmol_MUB_gDW_h)

# export finished dataset
write_csv(finalising_step2, "Data/Clean/clean_enzyme_data.csv")
    