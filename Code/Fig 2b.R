
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


# Select data
#results2 <- results[, c("SubjectID", "Loss_A_B_percentage","Gain_A_B_percentage", "Swap_A_B_percentage", "Persistence_A_B_percentage", 
#                                "QGain_A_B_percentage", "QLoss_A_B_percentage", "QAbsence_A_B_percentage")]
#results2 <- results[, c(1, 6:26)]

results_A_B <- results[, c("SubjectID", "Loss_A_B_percentage","Gain_A_B_percentage", "Swap_A_B_percentage", "Persistence_A_B_percentage", 
                        "QGain_A_B_percentage", "QLoss_A_B_percentage", "QAbsence_A_B_percentage")] %>% filter(Gain_A_B_percentage != "NA")
results_B_D <- results[, c("SubjectID", "Loss_B_D_percentage","Gain_B_D_percentage", "Swap_B_D_percentage", "Persistence_B_D_percentage", 
                        "QGain_B_D_percentage", "QLoss_B_D_percentage", "QAbsence_B_D_percentage")] %>% filter(Gain_B_D_percentage != "NA")
results_D_C <- results[, c("SubjectID", "Loss_D_C_percentage","Gain_D_C_percentage", "Swap_D_C_percentage", "Persistence_D_C_percentage", 
                        "QGain_D_C_percentage", "QLoss_D_C_percentage", "QAbsence_D_C_percentage")] %>% filter(Gain_D_C_percentage != "NA")   

library(reshape2)
library(dplyr)
library(stringr)
#########################################################################################################
results_melted <- melt(results_A_B, id="SubjectID")
head(results_melted)
#  SubjectID            variable    value
#1         4 Loss_A_B_percentage 24.77876
#2        13 Loss_A_B_percentage 25.92593
#3        15 Loss_A_B_percentage 23.21429
#4        17 Loss_A_B_percentage 25.73529
#5        19 Loss_A_B_percentage 24.16107
#6        22 Loss_A_B_percentage 33.78016

# extract A_B 
results_melted2 <- results_melted %>%
  mutate(new_column = sapply(str_split(variable, "_"), function(x) {
    if(length(x) >= 3) {
      # Paste together elements from the second element until the second-to-last element
      paste(x[3:(length(x))-1], collapse = "_")
    } else {
      NA_character_
    }
  }))

# extract the text before the first underscore (for example, "gain" from "gain_A_B_percentage")    
results_melted2 <- results_melted2 %>% mutate(prefix = str_extract(variable, "^[^_]+")) 
names(results_melted2) <- c("SubjectID", "Scenario", "Percentage", "Time", "Change")
head(results_melted2)

mean_by_scenario_A_B <- results_melted2 %>%
  group_by(Scenario, Time, Change) %>%
  summarise(mean_percentage = mean(Percentage, na.rm = TRUE),
            .groups = "drop")

mean_by_scenario_A_B
# A tibble: 7 × 4
#  Scenario                   Time  Change      mean_percentage
#  <fct>                      <chr> <chr>                 <dbl>
#1 Loss_A_B_percentage        A_B   Loss                  27.5 
#2 Gain_A_B_percentage        A_B   Gain                  21.8 
#3 Swap_A_B_percentage        A_B   Swap                   4.26
#4 Persistence_A_B_percentage A_B   Persistence           18.4 
#5 QGain_A_B_percentage       A_B   QGain                  6.52
#6 QLoss_A_B_percentage       A_B   QLoss                 11.0 
#7 QAbsence_A_B_percentage    A_B   QAbsence              38.1 

#########################################################################################################
# repeat the calculation for B_D
results_melted <- melt(results_B_D, id="SubjectID")
head(results_melted)

# extract B_D 
results_melted2 <- results_melted %>%
  mutate(new_column = sapply(str_split(variable, "_"), function(x) {
    if(length(x) >= 3) {
      # Paste together elements from the second element until the second-to-last element
      paste(x[3:(length(x))-1], collapse = "_")
    } else {
      NA_character_
    }
  }))

# extract the text before the first underscore (for example, "gain" from "gain_A_B_percentage")    
results_melted2 <- results_melted2 %>% mutate(prefix = str_extract(variable, "^[^_]+")) 
names(results_melted2) <- c("SubjectID", "Scenario", "Percentage", "Time", "Change")
head(results_melted2)

mean_by_scenario_B_D <- results_melted2 %>%
  group_by(Scenario, Time, Change) %>%
  summarise(mean_percentage = mean(Percentage, na.rm = TRUE),
            .groups = "drop")

mean_by_scenario_B_D
# A tibble: 7 × 4
#  Scenario                   Time  Change      mean_percentage
#  <fct>                      <chr> <chr>                 <dbl>
#1 Loss_B_D_percentage        B_D   Loss                  18.7 
#2 Gain_B_D_percentage        B_D   Gain                  29.9 
#3 Swap_B_D_percentage        B_D   Swap                   3.87
#4 Persistence_B_D_percentage B_D   Persistence           17.0 
#5 QGain_B_D_percentage       B_D   QGain                  9.91
#6 QLoss_B_D_percentage       B_D   QLoss                  4.86
#7 QAbsence_B_D_percentage    B_D   QAbsence              34.5 

##############################################################################################################################
# repeat the calculation for D_C
results_melted <- melt(results_D_C, id="SubjectID")
head(results_melted)

# extract D_C 
results_melted2 <- results_melted %>%
  mutate(new_column = sapply(str_split(variable, "_"), function(x) {
    if(length(x) >= 3) {
      # Paste together elements from the second element until the second-to-last element
      paste(x[3:(length(x))-1], collapse = "_")
    } else {
      NA_character_
    }
  }))

# extract the text before the first underscore (for example, "gain" from "gain_A_B_percentage")    
results_melted2 <- results_melted2 %>% mutate(prefix = str_extract(variable, "^[^_]+")) 
names(results_melted2) <- c("SubjectID", "Scenario", "Percentage", "Time", "Change")
head(results_melted2)

mean_by_scenario_D_C <- results_melted2 %>%
  group_by(Scenario, Time, Change) %>%
  summarise(mean_percentage = mean(Percentage, na.rm = TRUE),
            .groups = "drop")

mean_by_scenario_D_C
# A tibble: 7 × 4
#  Scenario                   Time  Change      mean_percentage
#  <fct>                      <chr> <chr>                 <dbl>
#1 Loss_D_C_percentage        D_C   Loss                  21.9 
#2 Gain_D_C_percentage        D_C   Gain                  18.5 
#3 Swap_D_C_percentage        D_C   Swap                   8.79
#4 Persistence_D_C_percentage D_C   Persistence           14.8 
#5 QGain_D_C_percentage       D_C   QGain                  5.62
#6 QLoss_D_C_percentage       D_C   QLoss                 10.5 
#7 QAbsence_D_C_percentage    D_C   QAbsence              41.7 

##############################################################################################################################
# Stack plot for gain loss swap and persistence
library(dplyr)
library(ggplot2)
mean_by_scenario <- rbind(mean_by_scenario_A_B, mean_by_scenario_B_D, mean_by_scenario_D_C)

# Transform data so that loss values are negative (if not already done)
mean_by_scenario <- mean_by_scenario %>%
  mutate(mean_percentage = if_else(Change == "Loss", -abs(mean_percentage), mean_percentage))

# Ensure the Change variable is a factor with the desired order 
mean_by_scenario <- mean_by_scenario %>%mutate(Change = factor(Change, 
                    levels = c("QAbsence", "QLoss", "QGain","Persistence", "Swap", "Gain", "Loss")))


# Create the side-by-side bar plot with revised x-axis labels.
ggplot(mean_by_scenario, aes(x = Time, y = mean_percentage, fill = Change)) +
    geom_hline(yintercept = 0, color = "grey", linetype = "solid", size = 0.5) +
    geom_bar(stat = "identity", position = "stack", width = 0.5, alpha=0.7) +
    labs(x = "Month", y = "Average % of strain/species",
         title = "Average % of Strain/Species Dynamics") +
    #scale_fill_manual(values = c("Persistence" = "#D95F02", "Swap" = "#404080", "Gain" = "#69b3a2", "Loss" = "grey")) +
    # Change x-axis labels to the desired intervals
    scale_x_discrete(labels = c("0-6", "6-14", "14-18")) +
    scale_y_continuous(limits = c(-30, 101), breaks = seq(-25, 100, by = 25)) +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0, size = 12, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black", face = "bold"),
        axis.text.y = element_text(size = 16, color = "black", face = "bold"),
        legend.text = element_text(size = 16, color = "black"),
        legend.title = element_blank(),
        legend.position = "right"
    ) +
    scale_fill_manual(
        values = c(
            "QAbsence"    = "grey80",  # very light grey
            "QLoss"       = "grey60",  # medium grey
            "QGain"       = "grey40",  # grey
            "Persistence" = "#D95F02",  # orange
            "Swap"        = "#404080",  # dark blue-gray
            "Gain"        = "#69b3a2",  # teal
            "Loss"        = "grey20"      # dark grey
        ), 
        labels = c("QAbsence" = "?-?",
                   "QLoss" = "1-?",
                   "QGain" = "?-1",
                   "Persistence"  = "Persistence",
                   "Swap" = "Swap",
                   "Gain" = "Gain (species)", 
                   "Loss" = "Loss (species)")
    )

# Output 6*5




library(openxlsx)
wb <- createWorkbook()

addWorksheet(wb, "Fig 2b Microbial Turnover")
writeData(wb, "Fig 2b Microbial Turnover", final_data)

saveWorkbook(wb, "Source Data Fig 2b Microbial Turnover.xlsx", overwrite = TRUE)
