
##############################################################################
########################## Gain and Loss #####################################
######################     Species level      ###############################
##############################################################################
setwd("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species")

# Load necessary library
library(dplyr)

# Read the data
data <- read.table("merged_abundance_table_filtered.tsv", header = TRUE, sep = "\t")
dim(data )
#[1]  327 1815
data[1:5, 1:6]
#  SampleID Treatment Time SubjectID Diet t__SGB10068
#1     A004      aFMT    A         4    3     0.00000
#2     A013   Placebo    A        13    2     0.00000
#3     A015   Placebo    A        15    3     0.00000
#4     A017   Placebo    A        17    3     0.00000
#5     A019   Placebo    A        19    2     0.00031
data[is.na(data)] <- 0
# Filter data for Time points A and B
data_A <- data %>% filter(Time == "A")
data_B <- data %>% filter(Time == "B")
data_D <- data %>% filter(Time == "D")
data_C <- data %>% filter(Time == "C")

# Merge data for A and B by SubjectID
merged_data_A_B <- merge(data_A, data_B, by = "SubjectID", suffixes = c("_species_A", "_species_B"))
merged_data_D_C <- merge(data_D, data_C, by = "SubjectID", suffixes = c("_species_D", "_species_C"))
species_merged_data <- full_join(merged_data_A_B, merged_data_D_C, by="SubjectID" )
species_merged_data[1:5, 1:6]
#  SubjectID SampleID_species_A Treatment_species_A Time_species_A Diet_species_A t__SGB10068_species_A
#1         4               A004                aFMT              A              3               0.00000
#2        13               A013             Placebo              A              2               0.00000
#3        15               A015             Placebo              A              3               0.00000
#4        17               A017             Placebo              A              3               0.00000
#5        19               A019             Placebo              A              2               0.00031

dim(species_merged_data)
#[1]   87 7257


##############################################################################
######################### Swap and Persistence ###############################
#########################     Strain level     ###############################
##############################################################################
setwd("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/5.matrix/")

# Load necessary library
library(dplyr)
library(tidyr) # Ensure tidyr is loaded for pivot_longer

# Define function to reformat a single file
reformat_file <- function(file_path) {
  # Read the file
  data <- read.table(file_path, header = TRUE, sep = "\t")

  # Gather the data into long format
  reformatted <- data %>% 
    pivot_longer(cols = -1, names_to = "TimePoint", values_to = "Value") %>% 
    rename(SubjectID = 1) %>% 
    mutate(FileName = basename(file_path))

  # Return the reformatted data
  reformatted
}

# Path to the folder containing the files
folder_path <- "~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/5.matrix/"  # Update with your folder path

# List all files in the folder
file_list <- list.files(folder_path, full.names = TRUE)

# Initialize an empty data frame to store combined results
combined_results <- data.frame()

# Loop through each file and reformat
for (file in file_list) {
  # Reformat the current file
  reformatted <- reformat_file(file)

  # Combine results
  combined_results <- bind_rows(combined_results, reformatted)
}

# Pivot wider to match required output structure
final_results <- combined_results %>% 
  pivot_wider(names_from = FileName, values_from = Value) %>% 
  arrange(SubjectID, TimePoint)

# Replace NA with 0 in final_results
final_results[is.na(final_results)] <- 0

names(final_results) <- gsub("_nGD_dynamic_matrix.tsv", "", names(final_results))

# Identify t__ names in species_columns that are not in strain_columns
species_columns <- grep("^t__", colnames(data), value = TRUE)
strain_columns <- grep("^t__", colnames(final_results), value = TRUE)

#new_cols <- species_columns[grepl("^t__", species_columns)]
new_cols <- setdiff(species_columns, strain_columns)

# Add these new columns with 0 to final_results
for(col in new_cols) {
  final_results[[col]] <- 0
}

# Check the result: print first few rows and columns
print(final_results[1:5, 1:7])

# Filter rows for TimePoint A, B, C, and D
final_results_filtered <- final_results %>% 
  filter(TimePoint %in% c("A", "B", "D", "C"))

# Filter data for Time points A and B
final_results_A <- final_results_filtered %>% filter(TimePoint == "A")
final_results_B <- final_results_filtered %>% filter(TimePoint == "B")
final_results_D <- final_results_filtered %>% filter(TimePoint == "D")
final_results_C <- final_results_filtered %>% filter(TimePoint == "C")

# Merge data for A and B by SubjectID
final_results_merged_data_A_B <- merge(final_results_A, final_results_B, by = "SubjectID", suffixes = c("_strain_A", "_strain_B"))
final_results_merged_data_D_C <- merge(final_results_D, final_results_C, by = "SubjectID", suffixes = c("_strain_D", "_strain_C"))
strain_merged_data <- full_join(final_results_merged_data_A_B, final_results_merged_data_D_C, by="SubjectID" )
dim(strain_merged_data)
#[1]   90 7245
strain_merged_data[1:5, 1:6]
#  SubjectID TimePoint_strain_A t__SGB10068_strain_A t__SGB10115_group_strain_A t__SGB10130_strain_A
#1         4                  A                    0                          0                    0
#2        13                  A                    0                          0                    0
#3        15                  A                    0                          0                    0
#4        17                  A                    0                          0                    0
#5        19                  A                    0                          0                    0


##############################################################################
##################### Merge Gain Loss Swap Persistence #######################
##############################################################################
full_merged_data <- merge(strain_merged_data, species_merged_data, by="SubjectID", all = TRUE )
dim(full_merged_data)
#[1]    90 9477

# Identify species SGB columns
species_columns <- grep("^t__", colnames(data), value = TRUE)
strain_columns <- grep("^t__", colnames(final_results), value = TRUE)


# Initialize result data frame
results <- full_merged_data %>% 
  rowwise() %>% 
  mutate(
    # Species level
    Total_species_A_B = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_A")) > 0 | get(paste0(species, "_species_B")) > 0
    })),
    Total_species_B_D = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_B")) > 0 | get(paste0(species, "_species_D")) > 0
    })),
    Total_species_D_C = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_D")) > 0 | get(paste0(species, "_species_C")) > 0
    })),
    
    Total_species_A = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_A")) > 0
    })),
    Total_species_B = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_B")) > 0
    })),
    Total_species_D = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_D")) > 0
    })),
    Total_species_C = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_C")) > 0
    })),
    
    # Species level LOSS
    Loss_count_A_B = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_A")) > 0 & get(paste0(species, "_species_B")) == 0
    })),
    Loss_count_B_D = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_B")) > 0 & get(paste0(species, "_species_D")) == 0
    })),
    Loss_count_D_C = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_D")) > 0 & get(paste0(species, "_species_C")) == 0
    })),
    
     # Species level GAIN
    Gain_count_A_B = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_A")) == 0 & get(paste0(species, "_species_B")) > 0
    })),
    Gain_count_B_D = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_B")) == 0 & get(paste0(species, "_species_D")) > 0
    })),
    Gain_count_D_C = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_D")) == 0 & get(paste0(species, "_species_C")) > 0
    })),
    
    # Strain level SWAP
    Swap_count_A_B = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_strain_A")) > 0 & get(paste0(species, "_strain_B")) > 0 &
      get(paste0(species, "_strain_A")) != get(paste0(species, "_strain_B"))
    })),
    Swap_count_B_D = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_strain_B")) > 0 & get(paste0(species, "_strain_D")) > 0 &
       get(paste0(species, "_strain_B")) != get(paste0(species, "_strain_D"))
    })),
    Swap_count_D_C = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_strain_D")) > 0 & get(paste0(species, "_strain_C")) > 0 &
      get(paste0(species, "_strain_D")) != get(paste0(species, "_strain_C"))
    })),
    
    # Strain level PERSISTENCE
    Persistence_count_A_B = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_strain_A")) > 0 & 
      get(paste0(species, "_strain_A")) == get(paste0(species, "_strain_B"))
    })),
    Persistence_count_B_D = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_strain_B")) > 0 &
      get(paste0(species, "_strain_B")) == get(paste0(species, "_strain_D"))
    })),
    Persistence_count_D_C = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_strain_D")) > 0 &
      get(paste0(species, "_strain_D")) == get(paste0(species, "_strain_C"))
    })),
    
    # Strain level Questioned Gain
    QGain_count_A_B = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_A")) > 0 & get(paste0(species, "_species_B")) > 0 &
      get(paste0(species, "_strain_A")) == 0 & get(paste0(species, "_strain_B")) > 0
    })),
    QGain_count_B_D = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_B")) > 0 & get(paste0(species, "_species_D")) > 0 &
      get(paste0(species, "_strain_B")) == 0 & get(paste0(species, "_strain_D")) > 0
    })),
    QGain_count_D_C = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_D")) > 0 & get(paste0(species, "_species_C")) > 0 &
          get(paste0(species, "_strain_D")) == 0 & get(paste0(species, "_strain_C")) > 0
    })),
    
    # Strain level Questioned Loss
    QLoss_count_A_B = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_A")) > 0 & get(paste0(species, "_species_B")) > 0 &
          get(paste0(species, "_strain_A")) > 0 & get(paste0(species, "_strain_B")) == 0
    })),
    QLoss_count_B_D = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_B")) > 0 & get(paste0(species, "_species_D")) > 0 &
          get(paste0(species, "_strain_B")) > 0 & get(paste0(species, "_strain_D")) == 0
    })),
    QLoss_count_D_C = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_D")) > 0 & get(paste0(species, "_species_C")) > 0 &
          get(paste0(species, "_strain_D")) > 0 & get(paste0(species, "_strain_C")) == 0
    })),
    
    # Strain level Questioned Absence
    QAbsence_count_A_B = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_A")) > 0 & get(paste0(species, "_species_B")) > 0 &
          get(paste0(species, "_strain_A")) == 0 & get(paste0(species, "_strain_B")) == 0
    })),
    QAbsence_count_B_D = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_B")) > 0 & get(paste0(species, "_species_D")) > 0 &
          get(paste0(species, "_strain_B")) == 0 & get(paste0(species, "_strain_D")) == 0
    })),
    QAbsence_count_D_C = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_D")) > 0 & get(paste0(species, "_species_C")) > 0 &
          get(paste0(species, "_strain_D")) == 0 & get(paste0(species, "_strain_C")) == 0
    })),
    
    # Strain level Gain-Persistence 
    Gain_Persistence_count = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_strain_A")) == 0 & get(paste0(species, "_strain_B")) == 1 & get(paste0(species, "_strain_D")) == 1
    })),
    # Strain level Swap-Persistence
    Swap_Persistence_count = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_strain_A")) == 1 & get(paste0(species, "_strain_B")) == 2 & get(paste0(species, "_strain_D")) == 2
    })),
    # Species level Loss-Persistence
    Loss_Persistence_count = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_species_A")) > 0 & get(paste0(species, "_species_B")) == 0 & get(paste0(species, "_species_D")) == 0
    })),
    
     # Strain level Succession rate = Succeeded / Successional
     # Succeeded = (0-1-1, 1-2-2), Successional = (0-1, 1-2)
     Succeeded_v1_count = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_strain_A")) == 0 & get(paste0(species, "_strain_B")) == 1 & get(paste0(species, "_strain_D")) == 1
    })),
    Succeeded_v2_count = sum(sapply(species_columns, function(species) {
      get(paste0(species, "_strain_A")) == 1 & get(paste0(species, "_strain_B")) == 2 & get(paste0(species, "_strain_D")) == 2
    })),
    Succeeded_count = Succeeded_v1_count + Succeeded_v2_count,
    Successional_count = Gain_count_A_B + Swap_count_A_B,
    

    # % Species level Gain (divided by 2nd time point)
    Gain_A_B_percentage = Gain_count_A_B / Total_species_B * 100,
    Gain_B_D_percentage = Gain_count_B_D / Total_species_D * 100,
    Gain_D_C_percentage = Gain_count_D_C / Total_species_C * 100,
    # % Species level Loss (divided by 1st time point)
    Loss_A_B_percentage = Loss_count_A_B / Total_species_A * 100,
    Loss_B_D_percentage = Loss_count_B_D / Total_species_B * 100,
    Loss_D_C_percentage = Loss_count_D_C / Total_species_D * 100,
    # % Strain level Swap (divided by 2nd time point)
    Swap_A_B_percentage = Swap_count_A_B / Total_species_B * 100,
    Swap_B_D_percentage = Swap_count_B_D / Total_species_D * 100,
    Swap_D_C_percentage = Swap_count_D_C / Total_species_C * 100,
    # % Strain level Persistence (divided by 2nd time point)
    Persistence_A_B_percentage = Persistence_count_A_B / Total_species_B * 100,
    Persistence_B_D_percentage = Persistence_count_B_D / Total_species_D * 100,
    Persistence_D_C_percentage = Persistence_count_D_C / Total_species_C * 100,
    # % Strain level Questioned Gain (divided by 2nd time point)    
    QGain_A_B_percentage = QGain_count_A_B / Total_species_B * 100,
    QGain_B_D_percentage = QGain_count_B_D / Total_species_D * 100,
    QGain_D_C_percentage = QGain_count_D_C / Total_species_C * 100,
    # % Strain level Questioned Loss (divided by 2nd time point)    
    QLoss_A_B_percentage = QLoss_count_A_B / Total_species_B * 100,
    QLoss_B_D_percentage = QLoss_count_B_D / Total_species_D * 100,
    QLoss_D_C_percentage = QLoss_count_D_C / Total_species_C * 100,
    # % Strain level Questioned Absence (divided by 2nd time point)            
    QAbsence_A_B_percentage = QAbsence_count_A_B / Total_species_B * 100,
    QAbsence_B_D_percentage = QAbsence_count_B_D / Total_species_D * 100,
    QAbsence_D_C_percentage = QAbsence_count_D_C / Total_species_C * 100,
    
    # % Strain level Gain_Persistence (divided by 3rd time point D)                     
    Gain_Persistence_percentage = Gain_Persistence_count / Total_species_D * 100,
    # Strain level Swap_Persistence
    Swap_Persistence_percentage = Swap_Persistence_count / Total_species_D * 100,
    # % Strain level Loss_Persistence (divided by 1st time point A)
    Loss_Persistence_percentage = Loss_Persistence_count / Total_species_A * 100,
    
    # % Strain level Succession rate
    Succession_rate = Succeeded_count / Successional_count * 100
  
  ) %>% 
  ungroup() %>% 
  select(SubjectID, #Diet = Diet_A, Treatment = Treatment_A,
         Total_species_A, Total_species_B, Total_species_D, Total_species_C,
         Gain_A_B_percentage, Gain_B_D_percentage, Gain_D_C_percentage,
         Loss_A_B_percentage,Loss_B_D_percentage, Loss_D_C_percentage,
         #Swap_count_A_B, Swap_count_B_D, Swap_count_D_C, 
         Swap_A_B_percentage, Swap_B_D_percentage, Swap_D_C_percentage,
         Persistence_A_B_percentage, Persistence_B_D_percentage, Persistence_D_C_percentage,
         QGain_A_B_percentage, QGain_B_D_percentage, QGain_D_C_percentage,
         QLoss_A_B_percentage, QLoss_B_D_percentage, QLoss_D_C_percentage,
         QAbsence_A_B_percentage, QAbsence_B_D_percentage, QAbsence_D_C_percentage,
         Gain_Persistence_percentage, Loss_Persistence_percentage, Swap_Persistence_percentage,
         Succession_rate
)

# Print or write the results
# write.table(results, "results.txt", sep = "\t", row.names = FALSE, quote = FALSE)

dim(results)
#[1] 90 30

results
# A tibble: 90 × 30
#   SubjectID Total_species_A Total_species_B Total_species_D Total_species_C Gain_A_B_percentage Gain_B_D_percentage Gain_D_C_percentage
#       <int>           <int>           <int>           <int>           <int>               <dbl>               <dbl>               <dbl>
# 1         4             226             220             288             256               22.7                 31.2               12.1 
# 2        13             135             147             137             146               32.0                 21.2               26.0 
# 3        15             168             172             145             181               25                   22.1               30.4 
# 4        17             136             129              NA              NA               21.7                 NA                 NA   
# 5        19             447             368             440             374                7.88                20.9                6.95
# 6        22             373             283             375             318               12.7                 33.6               11.9 
# 7        27             144             147              NA              NA               29.9                 NA                 NA   
# 8        31             202              74             143             155               18.9                 62.9               36.8 
# 9        37              NA              NA             186             227               NA                   NA                 26.0 
#10        38             237             213             253             179               18.8                 24.1                3.35
# ℹ 80 more rows
# ℹ 22 more variables: Loss_A_B_percentage <dbl>, Loss_B_D_percentage <dbl>, Loss_D_C_percentage <dbl>, Swap_A_B_percentage <dbl>,
#   Swap_B_D_percentage <dbl>, Swap_D_C_percentage <dbl>, Persistence_A_B_percentage <dbl>, Persistence_B_D_percentage <dbl>,
#   Persistence_D_C_percentage <dbl>, QGain_A_B_percentage <dbl>, QGain_B_D_percentage <dbl>, QGain_D_C_percentage <dbl>,
#   QLoss_A_B_percentage <dbl>, QLoss_B_D_percentage <dbl>, QLoss_D_C_percentage <dbl>, QAbsence_A_B_percentage <dbl>,
#   QAbsence_B_D_percentage <dbl>, QAbsence_D_C_percentage <dbl>, Gain_Persistence_percentage <dbl>, Gain_Swap_percentage <dbl>,
#   Loss_Persistence_percentage <dbl>, Succession_rate <dbl>
# ℹ Use `print(n = ...)` to see more rows




###############################################################################
setwd("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species")
# Select data
results3 <- results[, c("SubjectID", "Gain_Persistence_percentage", "Loss_Persistence_percentage", "Swap_Persistence_percentage", "Succession_rate")]
# Read metadata and merge them
md <- read_tsv(file = "metadata2.tsv", col_names = TRUE, show_col_types = FALSE)
results3_merged <- merge(md, results3, by="SubjectID", all = TRUE )
head(results3_merged)

###############################################################################
########################## Covariates into quartiles ##########################
results3_merged2 <- results3_merged %>% filter(Time!= "18")

# Weight change into quartiles
library(dplyr)
# Filter for time 0 and time 6
time_0 <- results3_merged2 %>% filter(Time == 0)
time_6 <- results3_merged2 %>% filter(Time == 6)
# Merge the data frames on SubjectID
merged_data <- merge(time_0, time_6, by = "SubjectID", suffixes = c("_0", "_6"))
# Calculate weight change
weight_change_df <- merged_data %>%
  mutate(Weight_change = Weight_6 - Weight_0) %>%
  select(SubjectID, Weight_change) %>% distinct()
results4 <- merge(results3_merged, weight_change_df, by = "SubjectID")
results4$Weight_change_Quartile <- cut(
  results4$Weight_change,
  breaks = quantile(results4$Weight_change, probs = seq(0, 1, 0.25), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("1", "2", "3", "4"))
results4$Weight_change_Quartile <- as.factor(results4$Weight_change_Quartile)

# Categorize Waist Circumference change into quartiles
WC_change_df <- merged_data %>%
  mutate(WC_change = WC_6 - WC_0) %>%
  select(SubjectID, WC_change) %>% distinct()
results4 <- merge(results4, WC_change_df, by = "SubjectID")  
results4$WC_change_Quartile <- cut(
  results4$WC_change,
  breaks = quantile(results4$WC_change, probs = seq(0, 1, 0.25), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("1", "2", "3", "4"))
results4$WC_change_Quartile <- as.factor(results4$WC_change_Quartile)

# Categorize Age into quartiles
results4$Age_Quartile <- cut(
  results4$Age,
  breaks = quantile(results4$Age, probs = seq(0, 1, 0.25), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("1", "2", "3", "4"))
results4$Age_Quartile <- as.factor(results4$Age_Quartile)

# Categorize Richness into quartiles
results4$Richness_Quartile <- cut(
  results4$Richness,
  breaks = quantile(results4$Richness, probs = seq(0, 1, 0.25), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("1", "2", "3", "4"))
results4$Richness_Quartile <- as.factor(results4$Richness_Quartile)

results4$Gender <- as.factor(results4$Gender)
results4$Diet <- as.factor(results4$Diet)
results4$Diet <- relevel(results4$Diet, ref = "2")
results4$Treatment <- as.factor(results4$Treatment)
results4$Treatment <- relevel(results4$Treatment, ref = "Placebo")



###############################################################################
##  Adjusted covariates for the overall level
###############################################################################

# Re-center adjusted overall values
# Fit a linear model
# Using emmeans (Successor of lsmeans)
# install.packages("emmeans")  # If not installed
library(emmeans)

# colonization_rate
succession_rate <- results4[ , c("SubjectID", "Diet", "Treatment", 
                "Age_Quartile", "Gender", "WC_change_Quartile", "Weight_change_Quartile", "Richness_Quartile", 
                "Succession_rate")] %>% distinct()
succession_rate_model <- lm(Succession_rate ~ Treatment + Age_Quartile + Gender + WC_change_Quartile + 
                     Weight_change_Quartile + Richness_Quartile, succession_rate, na.action = na.exclude) 
succession_rate$predicted_value <- predict(succession_rate_model)
lsmeans_succession_rate <- emmeans(succession_rate_model, ~ Treatment)
lsmeans_succession_rate <- as.data.frame(lsmeans_succession_rate)[, c("Treatment", "emmean")]
# Merge the LS‐Means Back into Your Data:
succession_rate <- succession_rate %>% left_join(lsmeans_succession_rate, by = "Treatment")
# Calculate the Re-centered (Adjusted) Raw Values:
succession_rate <- succession_rate %>% mutate(adjusted_raw = Succession_rate - (predicted_value - emmean))


# Gain-Persistence
persistent_gain <- results4[ , c("SubjectID", "Diet", "Treatment", 
                "Age_Quartile", "Gender", "WC_change_Quartile", "Weight_change_Quartile", "Richness_Quartile", 
                "Gain_Persistence_percentage")] %>% distinct()
persistent_gain_model <- lm(Gain_Persistence_percentage ~ Treatment + Age_Quartile + Gender + WC_change_Quartile + 
                     Weight_change_Quartile + Richness_Quartile, persistent_gain, na.action = na.exclude) 
persistent_gain$predicted_value <- predict(persistent_gain_model)
# Get the least-squares means for 'aFMT'
lsmeans_persistent_gain <- emmeans(persistent_gain_model, ~ Treatment )
lsmeans_persistent_gain <- as.data.frame(lsmeans_persistent_gain)[, c("Treatment", "emmean")]
# Merge the LS‐Means Back into Your Data:
persistent_gain <- persistent_gain %>% left_join(lsmeans_persistent_gain, by = "Treatment")
# Calculate the Re-centered (Adjusted) Raw Values:
persistent_gain <- persistent_gain %>% mutate(adjusted_raw = Gain_Persistence_percentage - (predicted_value - emmean))


# persistent_loss
persistent_loss <- results4[ , c("SubjectID", "Diet", "Treatment", 
                "Age_Quartile", "Gender", "WC_change_Quartile", "Weight_change_Quartile", "Richness_Quartile", 
                "Loss_Persistence_percentage")] %>% distinct()
persistent_loss_model <- lm(Loss_Persistence_percentage ~ Treatment + Age_Quartile + Gender + WC_change_Quartile + 
                     Weight_change_Quartile + Richness_Quartile, persistent_loss, na.action = na.exclude) 
persistent_loss$predicted_value <- predict(persistent_loss_model)
lsmeans_persistent_loss <- emmeans(persistent_loss_model, ~ Treatment )
lsmeans_persistent_loss <- as.data.frame(lsmeans_persistent_loss)[, c("Treatment", "emmean")]
# Merge the LS‐Means Back into Your Data:
persistent_loss <- persistent_loss %>% left_join(lsmeans_persistent_loss, by = "Treatment")
# Calculate the Re-centered (Adjusted) Raw Values:
persistent_loss <- persistent_loss %>% mutate(adjusted_raw = Loss_Persistence_percentage - (predicted_value - emmean))


# Swap_Persistence
swap_persistence <- results4[ , c("SubjectID", "Diet", "Treatment", 
                "Age_Quartile", "Gender", "WC_change_Quartile", "Weight_change_Quartile", "Richness_Quartile", 
                "Swap_Persistence_percentage")] %>% distinct()
swap_persistence_model <- lm(Swap_Persistence_percentage ~ Treatment + Age_Quartile + Gender + WC_change_Quartile + 
                            Weight_change_Quartile + Richness_Quartile, swap_persistence, na.action = na.exclude) 
swap_persistence$predicted_value <- predict(swap_persistence_model)
lsmeans_swap_persistence <- emmeans(swap_persistence_model, ~ Treatment )
lsmeans_swap_persistence <- as.data.frame(lsmeans_swap_persistence)[, c("Treatment", "emmean")]
# Merge the LS‐Means Back into Your Data:
swap_persistence <- swap_persistence %>% left_join(lsmeans_swap_persistence, by = "Treatment")
# Calculate the Re-centered (Adjusted) Raw Values:
swap_persistence <- swap_persistence %>% mutate(adjusted_raw = Swap_Persistence_percentage - (predicted_value - emmean))

# Add the metric name
succession_rate$Metric <- "Succession Rate"
persistent_gain$Metric <- "Gain-Persistence"
persistent_loss$Metric <- "Loss-Persistence"
swap_persistence$Metric <- "Swap-persistence"
# Rename the percentage column to Value to merge them
succession_rate <- succession_rate %>% rename(Value = Succession_rate)
persistent_gain <- persistent_gain %>% rename(Value = Gain_Persistence_percentage)
persistent_loss <- persistent_loss %>% rename(Value = Loss_Persistence_percentage)
swap_persistence <- swap_persistence %>% rename(Value = Swap_Persistence_percentage)


persistent_gain_loss_colonization <- rbind(#as.data.frame(persistent_gain), 
                                           as.data.frame(persistent_loss), 
                                           as.data.frame(succession_rate)
                                           #as.data.frame(swap_persistence)
                                           )
head(persistent_gain_loss_colonization)
# Reverse the x and y axis
persistent_gain_loss_colonization$Metric <- factor(persistent_gain_loss_colonization$Metric, 
                levels = c("Succession Rate", "Loss-Persistence"))

ggplot(persistent_gain_loss_colonization, aes(x = Metric, y = predicted_value , fill = factor(Treatment))) +
    geom_boxplot(outlier.color = "grey", outlier.shape = 16, width = 0.7) +
    scale_fill_manual(
        values = c("Placebo" = "#7796A7", "aFMT" = "#B9772B")
    ) +
    labs(
        title = "Strain Change induced by diets and aFMT",
        x = "",
        y = "% of species"
    ) +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 16, size = 0, face = "bold"),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black", face = "bold", angle = -20),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 20, color = "black"),
        legend.position = "none",
        strip.text = element_text(size = 20, face = "bold")
    )
# Output 10*3 for horizontal 
# Output portrait 8*4 for vertical

library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "Fig 2d Strain Persist Overall")
writeData(wb, "Fig 2d Strain Persist Overall", persistent_gain_loss_colonization)
#saveWorkbook(wb, "Source Data Fig 2c.xlsx", overwrite = TRUE)



###############################################################################
##  Adjusted covariates for the diet-stratified level
###############################################################################

# colonization_rate
succession_rate <- results4[ , c("SubjectID", "Diet", "Treatment", 
                "Age_Quartile", "Gender", "WC_change_Quartile", "Weight_change_Quartile", "Richness_Quartile", 
                "Succession_rate")] %>% distinct()


# Using emmeans (Successor of lsmeans)
# install.packages("emmeans")  # If not installed
library(emmeans)

succession_rate_diet1 <- succession_rate[succession_rate$Diet=="1",]
succession_rate_diet1_model <- lm(Succession_rate ~ Treatment + Age_Quartile + Gender + WC_change_Quartile + Weight_change_Quartile + Richness_Quartile, succession_rate_diet1,  na.action = na.exclude) 
succession_rate_diet1$predicted_value <- predict(succession_rate_diet1_model)
lsmeans_succession_rate_diet1 <- emmeans(succession_rate_diet1_model, ~ Treatment)
lsmeans_succession_rate_diet1_df <- as.data.frame(lsmeans_succession_rate_diet1)[, c("Treatment", "emmean")]
# Merge the LS‐Means Back into Your Data:
succession_rate_diet1 <- succession_rate_diet1 %>% left_join(lsmeans_succession_rate_diet1_df, by = "Treatment")
# Calculate the Re-centered (Adjusted) Raw Values:
succession_rate_diet1 <- succession_rate_diet1 %>% mutate(adjusted_raw = Succession_rate - (predicted_value - emmean))

succession_rate_diet2 <- succession_rate[succession_rate$Diet=="2",]
succession_rate_diet2_model <- lm(Succession_rate ~ Treatment + Age_Quartile + Gender + WC_change_Quartile + Weight_change_Quartile + Richness_Quartile, succession_rate_diet2,  na.action = na.exclude) 
succession_rate_diet2$predicted_value <- predict(succession_rate_diet2_model)
lsmeans_succession_rate_diet2 <- emmeans(succession_rate_diet2_model, ~ Treatment)
lsmeans_succession_rate_diet2_df <- as.data.frame(lsmeans_succession_rate_diet2)[, c("Treatment", "emmean")]
# Merge the LS‐Means Back into Your Data:
succession_rate_diet2 <- succession_rate_diet2 %>% left_join(lsmeans_succession_rate_diet2_df, by = "Treatment")
# Calculate the Re-centered (Adjusted) Raw Values:
succession_rate_diet2 <- succession_rate_diet2 %>% mutate(adjusted_raw = Succession_rate - (predicted_value - emmean))

succession_rate_diet3 <- succession_rate[succession_rate$Diet=="3",]
succession_rate_diet3_model <- lm(Succession_rate ~ Treatment + Age_Quartile + Gender + WC_change_Quartile + Weight_change_Quartile + Richness_Quartile, succession_rate_diet3,  na.action = na.exclude) 
succession_rate_diet3$predicted_value <- predict(succession_rate_diet3_model)
lsmeans_succession_rate_diet3 <- emmeans(succession_rate_diet3_model, ~ Treatment)
lsmeans_succession_rate_diet3_df <- as.data.frame(lsmeans_succession_rate_diet3)[, c("Treatment", "emmean")]
# Merge the LS‐Means Back into Your Data:
succession_rate_diet3 <- succession_rate_diet3 %>% left_join(lsmeans_succession_rate_diet3_df, by = "Treatment")
# Calculate the Re-centered (Adjusted) Raw Values:
succession_rate_diet3 <- succession_rate_diet3 %>% mutate(adjusted_raw = Succession_rate - (predicted_value - emmean))

succession_rate_3diets <- rbind(succession_rate_diet1,
                                          succession_rate_diet2,
                                          succession_rate_diet3)
succession_rate_3diets$Metric <- "Succession Rate"
succession_rate_3diets <- succession_rate_3diets %>% rename(Value = Succession_rate)


# Gain-Persistence
persistent_gain <- results4[ , c("SubjectID", "Diet", "Treatment", 
                "Age_Quartile", "Gender", "WC_change_Quartile", "Weight_change_Quartile", "Richness_Quartile", 
                "Gain_Persistence_percentage")] %>% distinct()


persistent_gain_diet1 <- persistent_gain[persistent_gain$Diet=="1",]
persistent_gain_diet1_model <- lm(Gain_Persistence_percentage ~ Treatment + Age_Quartile + Gender + WC_change_Quartile + Weight_change_Quartile + Richness_Quartile, persistent_gain_diet1, na.action = na.exclude) 
persistent_gain_diet1$predicted_value <- predict(persistent_gain_diet1_model)
lsmeans_persistent_gain_diet1 <- emmeans(persistent_gain_diet1_model, ~ Treatment )
lsmeans_persistent_gain_diet1 <- as.data.frame(lsmeans_persistent_gain_diet1)[, c("Treatment", "emmean")]
# Merge the LS‐Means Back into Your Data:
persistent_gain_diet1 <- persistent_gain_diet1 %>% left_join(lsmeans_persistent_gain_diet1, by = "Treatment")
# Calculate the Re-centered (Adjusted) Raw Values:
persistent_gain_diet1 <- persistent_gain_diet1 %>% mutate(adjusted_raw = Gain_Persistence_percentage - (predicted_value - emmean))

persistent_gain_diet2 <- persistent_gain[persistent_gain$Diet=="2",]
persistent_gain_diet2_model <- lm(Gain_Persistence_percentage ~ Treatment + Age_Quartile + Gender + WC_change_Quartile + Weight_change_Quartile + Richness_Quartile, persistent_gain_diet2, na.action = na.exclude) 
lsmeans_persistent_gain_diet2 <- emmeans(persistent_gain_diet2_model, ~ Treatment)
persistent_gain_diet2$predicted_value <- predict(persistent_gain_diet2_model)
lsmeans_persistent_gain_diet2 <- as.data.frame(lsmeans_persistent_gain_diet2)[, c("Treatment", "emmean")]
# Merge the LS‐Means Back into Your Data:
persistent_gain_diet2 <- persistent_gain_diet2 %>% left_join(lsmeans_persistent_gain_diet2, by = "Treatment")
# Calculate the Re-centered (Adjusted) Raw Values:
persistent_gain_diet2 <- persistent_gain_diet2 %>% mutate(adjusted_raw = Gain_Persistence_percentage - (predicted_value - emmean))


persistent_gain_diet3 <- persistent_gain[persistent_gain$Diet=="3",]
persistent_gain_diet3_model <- lm(Gain_Persistence_percentage ~ Treatment + Age_Quartile + Gender + WC_change_Quartile + Weight_change_Quartile + Richness_Quartile, persistent_gain_diet3, na.action = na.exclude) 
lsmeans_persistent_gain_diet3 <- emmeans(persistent_gain_diet3_model, ~ Treatment)
persistent_gain_diet3$predicted_value <- predict(persistent_gain_diet3_model)
lsmeans_persistent_gain_diet3 <- as.data.frame(lsmeans_persistent_gain_diet3)[, c("Treatment", "emmean")]
# Merge the LS‐Means Back into Your Data:
persistent_gain_diet3 <- persistent_gain_diet3 %>% left_join(lsmeans_persistent_gain_diet3, by = "Treatment")
# Calculate the Re-centered (Adjusted) Raw Values:
persistent_gain_diet3 <- persistent_gain_diet3 %>% mutate(adjusted_raw = Gain_Persistence_percentage - (predicted_value - emmean))

persistent_gain_3diets <- rbind(persistent_gain_diet1,
                                        persistent_gain_diet2,
                                        persistent_gain_diet3)
persistent_gain_3diets$Metric <- "Gain-Persistence"
persistent_gain_3diets <- persistent_gain_3diets %>% rename(Value = Gain_Persistence_percentage)


# persistent_loss
persistent_loss <- results4[ , c("SubjectID", "Diet", "Treatment", 
                "Age_Quartile", "Gender", "WC_change_Quartile", "Weight_change_Quartile", "Richness_Quartile", 
                "Loss_Persistence_percentage")] %>% distinct()

persistent_loss_diet1 <- persistent_loss[persistent_loss$Diet=="1",]
persistent_loss_diet1_model <- lm(Loss_Persistence_percentage ~ Treatment + Age_Quartile + Gender + WC_change_Quartile + Weight_change_Quartile + Richness_Quartile, persistent_loss_diet1, na.action = na.exclude) 
persistent_loss_diet1$predicted_value <- predict(persistent_loss_diet1_model)
lsmeans_persistent_loss_diet1 <- emmeans(persistent_loss_diet1_model, ~ Treatment )
lsmeans_persistent_loss_diet1 <- as.data.frame(lsmeans_persistent_loss_diet1)[, c("Treatment", "emmean")]
# Merge the LS‐Means Back into Your Data:
persistent_loss_diet1 <- persistent_loss_diet1 %>% left_join(lsmeans_persistent_loss_diet1, by = "Treatment")
# Calculate the Re-centered (Adjusted) Raw Values:
persistent_loss_diet1 <- persistent_loss_diet1 %>% mutate(adjusted_raw = Loss_Persistence_percentage - (predicted_value - emmean))


persistent_loss_diet2 <- persistent_loss[persistent_loss$Diet=="2",]
persistent_loss_diet2_model <- lm(Loss_Persistence_percentage ~ Treatment + Age_Quartile + Gender + WC_change_Quartile + Weight_change_Quartile + Richness_Quartile, persistent_loss_diet2, na.action = na.exclude) 
persistent_loss_diet2$predicted_value <- predict(persistent_loss_diet2_model)
lsmeans_persistent_loss_diet2 <- emmeans(persistent_loss_diet2_model, ~ Treatment )
lsmeans_persistent_loss_diet2 <- as.data.frame(lsmeans_persistent_loss_diet2)[, c("Treatment", "emmean")]
# Merge the LS‐Means Back into Your Data:
persistent_loss_diet2 <- persistent_loss_diet2 %>% left_join(lsmeans_persistent_loss_diet2, by = "Treatment")
# Calculate the Re-centered (Adjusted) Raw Values:
persistent_loss_diet2 <- persistent_loss_diet2 %>% mutate(adjusted_raw = Loss_Persistence_percentage - (predicted_value - emmean))


persistent_loss_diet3 <- persistent_loss[persistent_loss$Diet=="3",]
persistent_loss_diet3_model <- lm(Loss_Persistence_percentage ~ Treatment + Age_Quartile + Gender + WC_change_Quartile + Weight_change_Quartile + Richness_Quartile, persistent_loss_diet3, na.action = na.exclude) 
persistent_loss_diet3$predicted_value <- predict(persistent_loss_diet3_model)
lsmeans_persistent_loss_diet3 <- emmeans(persistent_loss_diet3_model, ~ Treatment )
lsmeans_persistent_loss_diet3 <- as.data.frame(lsmeans_persistent_loss_diet3)[, c("Treatment", "emmean")]
# Merge the LS‐Means Back into Your Data:
persistent_loss_diet3 <- persistent_loss_diet3 %>% left_join(lsmeans_persistent_loss_diet3, by = "Treatment")
# Calculate the Re-centered (Adjusted) Raw Values:
persistent_loss_diet3 <- persistent_loss_diet3 %>% mutate(adjusted_raw = Loss_Persistence_percentage - (predicted_value - emmean))

persistent_loss_3diets <- rbind(persistent_loss_diet1,
                                        persistent_loss_diet2,
                                        persistent_loss_diet3)
persistent_loss_3diets$Metric <- "Loss-Persistence"
persistent_loss_3diets <- persistent_loss_3diets %>% rename(Value = Loss_Persistence_percentage)


# Swap_Persistence
swap_persistence <- results4[ , c("SubjectID", "Diet", "Treatment", 
                "Age_Quartile", "Gender", "WC_change_Quartile", "Weight_change_Quartile", "Richness_Quartile", 
                "Swap_Persistence_percentage")] %>% distinct()

swap_persistence_diet1 <- swap_persistence[swap_persistence$Diet=="1",]
swap_persistence_diet1_model <- lm(Swap_Persistence_percentage ~ Treatment + Age_Quartile + Gender + WC_change_Quartile + Weight_change_Quartile + Richness_Quartile, swap_persistence_diet1, na.action = na.exclude) 
swap_persistence_diet1$predicted_value <- predict(swap_persistence_diet1_model)
lsmeans_swap_persistence_diet1 <- emmeans(swap_persistence_diet1_model, ~ Treatment )
lsmeans_swap_persistence_diet1 <- as.data.frame(lsmeans_swap_persistence_diet1)[, c("Treatment", "emmean")]
# Merge the LS‐Means Back into Your Data:
swap_persistence_diet1 <- swap_persistence_diet1 %>% left_join(lsmeans_swap_persistence_diet1, by = "Treatment")
# Calculate the Re-centered (Adjusted) Raw Values:
swap_persistence_diet1 <- swap_persistence_diet1 %>% mutate(adjusted_raw = Swap_Persistence_percentage - (predicted_value - emmean))


swap_persistence_diet2 <- swap_persistence[swap_persistence$Diet=="2",]
swap_persistence_diet2_model <- lm(Swap_Persistence_percentage ~ Treatment + Age_Quartile + Gender + WC_change_Quartile + Weight_change_Quartile + Richness_Quartile, swap_persistence_diet2, na.action = na.exclude) 
swap_persistence_diet2$predicted_value <- predict(swap_persistence_diet2_model)
lsmeans_swap_persistence_diet2 <- emmeans(swap_persistence_diet2_model, ~ Treatment )
lsmeans_swap_persistence_diet2 <- as.data.frame(lsmeans_swap_persistence_diet2)[, c("Treatment", "emmean")]
# Merge the LS‐Means Back into Your Data:
swap_persistence_diet2 <- swap_persistence_diet2 %>% left_join(lsmeans_swap_persistence_diet2, by = "Treatment")
# Calculate the Re-centered (Adjusted) Raw Values:
swap_persistence_diet2 <- swap_persistence_diet2 %>% mutate(adjusted_raw = Swap_Persistence_percentage - (predicted_value - emmean))


swap_persistence_diet3 <- swap_persistence[swap_persistence$Diet=="3",]
swap_persistence_diet3_model <- lm(Swap_Persistence_percentage ~ Treatment + Age_Quartile + Gender + WC_change_Quartile + Weight_change_Quartile + Richness_Quartile, swap_persistence_diet3, na.action = na.exclude) 
swap_persistence_diet3$predicted_value <- predict(swap_persistence_diet3_model)
lsmeans_swap_persistence_diet3 <- emmeans(swap_persistence_diet3_model, ~ Treatment )
lsmeans_swap_persistence_diet3 <- as.data.frame(lsmeans_swap_persistence_diet3)[, c("Treatment", "emmean")]
# Merge the LS‐Means Back into Your Data:
swap_persistence_diet3 <- swap_persistence_diet3 %>% left_join(lsmeans_swap_persistence_diet3, by = "Treatment")
# Calculate the Re-centered (Adjusted) Raw Values:
swap_persistence_diet3 <- swap_persistence_diet3 %>% mutate(adjusted_raw = Swap_Persistence_percentage - (predicted_value - emmean))

swap_persistence_3diets <- rbind(swap_persistence_diet1,
                                        swap_persistence_diet2,
                                        swap_persistence_diet3)
swap_persistence_3diets$Metric <- "Swap-Persistence"
swap_persistence_3diets <- swap_persistence_3diets %>% rename(Value = Swap_Persistence_percentage)



###############################################################################
# Merge the persistence scenarios
###############################################################################

persistent_gain_loss_colonization <- rbind(persistent_gain_3diets,
                                          persistent_loss_3diets,
                                          #swap_persistence_3diets,
                                          succession_rate_3diets)

# Remove rows with NA or NaN values in the column "adjusted_raw" 
#persistent_gain_loss_colonization <- persistent_gain_loss_colonization[!is.na(persistent_gain_loss_colonization$adjusted_raw) & !is.nan(persistent_gain_loss_colonization$adjusted_raw), ]
#persistent_gain_loss_colonization$adjusted_raw[persistent_gain_loss_colonization$adjusted_raw<0] <- 0
persistent_gain_loss_colonization <- persistent_gain_loss_colonization %>% filter(predicted_value > 0 | predicted_value == 0)

#persistent_gain_loss_colonization2 <- persistent_gain_loss_colonization %>% filter(Metric == "Loss-Persistence"| Metric == "Succession Rate")

# Create a separate annotation data frame with p values of diet-afmt interaction effects
interaction_df <- data.frame(
    Diet = factor("2", levels = c("1", "2", "3")),  # Must match the original factor levels
    #predicted_value = c(45, 15, 45),  # Adjust this y-value if needed
    predicted_value = c(45, 45),
    #Metric = c("Succession Rate", "Gain-Persistence", "Loss-Persistence"),
    Metric = c("Succession Rate", "Loss-Persistence"),
    #label = c("p.int = 9.46e-04", "p.int = 0.55", "p.int = 1.12e-03")
    label = c("p.int = 9.46e-04", "p.int = 1.12e-03")
)

# Create a separate annotation data frame with p values of afmt effects
afmt_succesionrate_df <- data.frame(
    Diet = factor(c("2", "3"), levels = c("1", "2", "3")),  # Must match the original factor levels
    predicted_value = c(38, 40),  # Adjust this y-value if needed
    Metric = c("Succession Rate"),
    label = c("p = 1.68e-03", "p = 6.76e-05")
)
# Create a separate annotation data frame with p values of afmt effects
afmt_losspersistence_df <- data.frame(
    Diet = factor(c("3"), levels = c("1", "2", "3")),  # Must match the original factor levels
    predicted_value = c(27),  # Adjust this y-value if needed
    Metric = c("Loss-Persistence"),
    label = c("p = 2.09e-03")
)
interaction_df$Metric <- factor(interaction_df$Metric, 
                                levels = c("Succession Rate","Gain-Persistence", "Loss-Persistence"))
afmt_succesionrate_df$Metric <- factor(afmt_succesionrate_df$Metric, 
                                levels = c("Succession Rate","Gain-Persistence", "Loss-Persistence"))
afmt_losspersistence_df$Metric <- factor(afmt_losspersistence_df$Metric, 
                                levels = c("Succession Rate","Gain-Persistence", "Loss-Persistence"))                                
persistent_gain_loss_colonization$Metric <- factor(persistent_gain_loss_colonization$Metric, 
                                levels = c("Succession Rate","Gain-Persistence", "Loss-Persistence"))
persistent_gain_loss_colonization$Diet <- factor(persistent_gain_loss_colonization$Diet, levels = c("1","2","3"))
persistent_gain_loss_colonization$Treatment <- factor(persistent_gain_loss_colonization$Treatment, levels = c("Placebo","aFMT"))


persistent_gain_loss_colonization2 <- persistent_gain_loss_colonization %>% filter(Metric == "Succession Rate" | Metric == "Loss-Persistence")
# Stratified by Diet and Treatment
ggplot(persistent_gain_loss_colonization2, 
       aes(x = factor(Diet), y = predicted_value, fill = factor(Treatment))) +
    #geom_point(position = position_jitterdodge(jitter.width = 0.2, 
    #                                           dodge.width = 0.7), 
    #           size = 3, alpha=0.1) +
    geom_boxplot(outlier.color = "grey", outlier.shape = 16, width = 0.7, 
                 position = position_dodge(width = 0.7)) +
    scale_fill_manual(
        values = c("Placebo" = "#7796A7", "aFMT" = "#B9772B")
        #labels = c("0" = "Placebo", "1" = "aFMT")
    ) +
    scale_x_discrete(labels = c("1" = "HDG", "2"= "MED", "3" = "GreenMED") ) +
    labs(
        title = "Species/Strain Change by Diet and aFMT",
        x = "",
        y = "Percentage (%)",
        color = "Diet"  # Optional if you want a title for the color legend
    ) +
    facet_wrap(~ Metric, scales = "free", nrow = 1,
               labeller = as_labeller(c(
                   "Succession Rate" = "Succession Rate (strain)\n(0-6-14)",
                   #"Gain-Persistence" = "Gain-Persistence (strain)\n(0-6-14)", 
                   "Loss-Persistence" = "Loss-Persistence (species)\n(0-6-14)"
                   #"Swap-Persistence" = "Swap-Persistence (strain)\n(0-6-14)",
               ))
    ) +
    geom_text(
        data = interaction_df,
        aes(x = Diet, y = predicted_value, label = label),
        inherit.aes = FALSE,
        size = 6
        #fontface = "bold"
    ) +
    geom_text(
        data = afmt_succesionrate_df,
        aes(x = Diet, y = predicted_value, label = label),
        inherit.aes = FALSE,
        size = 4
        #fontface = "bold"
    ) +
    geom_text(
        data = afmt_losspersistence_df,
        aes(x = Diet, y = predicted_value, label = label),
        inherit.aes = FALSE,
        size = 4
        #fontface = "bold"
    ) +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 0, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black", face = "bold"),
        axis.text.y = element_text(size = 20, color = "black", face = "bold"),
        legend.title = element_blank(),  # Remove if you don't want a legend title
        legend.text = element_text(size = 16, color = "black"),
        legend.position = "right",
        strip.text = element_text(size = 16, face = "bold")
    )


# Output 13*8
# Output 9*8 for two scenarios


addWorksheet(wb, "Fig 2d Strain Persist Diets")
writeData(wb, "Fig 2d Strain Persist Diets", persistent_gain_loss_colonization2)
saveWorkbook(wb, "Source Data Fig 2d.xlsx", overwrite = TRUE)


