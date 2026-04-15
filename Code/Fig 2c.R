##############################################################################
# Gain and Loss: Species level 
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



##############################################################################
# Gain Loss Swap during M0–M6
##############################################################################

setwd("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species")
# Select data
results3 <- results[, c("SubjectID", "Gain_Persistence_percentage", "Loss_Persistence_percentage", "Swap_Persistence_percentage", "Succession_rate")]
# Read metadata and merge them
library(readr)
md <- read_tsv(file = "metadata2.tsv", col_names = TRUE, show_col_types = FALSE)
results3_merged <- merge(md, results3, by="SubjectID", all = TRUE )

# Covariates into quartiles
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

head(results4)


###############################################################################

# Select data
gain_loss_swap_persistnce <- results[, c("SubjectID", "Loss_A_B_percentage","Gain_A_B_percentage", "Swap_A_B_percentage", "Persistence_A_B_percentage")]

gain_loss_swap_persistnce2 <- full_join(results4, gain_loss_swap_persistnce, by = "SubjectID")
head(gain_loss_swap_persistnce2)


library(emmeans)
# Loss
loss <- gain_loss_swap_persistnce2[ , c("SubjectID", "Diet", 
                "Age_Quartile", "Gender", "WC_change_Quartile", "Weight_change_Quartile", "Richness_Quartile", 
                "Loss_A_B_percentage")] %>% distinct()
loss_model <- lm(Loss_A_B_percentage ~ Diet + Age_Quartile + Gender + WC_change_Quartile + 
                     Weight_change_Quartile + Richness_Quartile, loss, na.action = na.exclude)
loss$predicted_value <- predict(loss_model)                     
lsmeans_loss <- emmeans(loss_model, ~ Diet)
lsmeans_loss <- as.data.frame(lsmeans_loss)[, c("Diet", "emmean")]
loss_lsmeans_loss <- loss %>% left_join(lsmeans_loss, by = "Diet") %>% distinct()
loss_lsmeans_loss <- loss_lsmeans_loss %>% mutate(adjusted_raw = Loss_A_B_percentage - (predicted_value - emmean))
loss_lsmeans_loss <- as.data.frame(loss_lsmeans_loss)
loss_lsmeans_loss <- loss_lsmeans_loss %>% rename(Value = Loss_A_B_percentage)

# Gain
gain <- gain_loss_swap_persistnce2[ , c("SubjectID", "Diet", 
                "Age_Quartile", "Gender", "WC_change_Quartile", "Weight_change_Quartile", "Richness_Quartile", 
                "Gain_A_B_percentage")] %>% distinct()
gain_model <- lm(Gain_A_B_percentage ~ Diet + Age_Quartile + Gender + WC_change_Quartile + 
    Weight_change_Quartile + Richness_Quartile, gain, na.action = na.exclude) 
gain$predicted_value <- predict(gain_model)
lsmeans_gain <- emmeans(gain_model, ~ Diet)
lsmeans_gain <- as.data.frame(lsmeans_gain)[, c("Diet", "emmean")]
gain_lsmeans_gain <- gain %>% left_join(lsmeans_gain, by = "Diet") %>% distinct()
gain_lsmeans_gain <- gain_lsmeans_gain %>% mutate(adjusted_raw = Gain_A_B_percentage - (predicted_value - emmean))
gain_lsmeans_gain <- as.data.frame(gain_lsmeans_gain)
gain_lsmeans_gain <- gain_lsmeans_gain %>% rename(Value = Gain_A_B_percentage)

# Swap
swap <- gain_loss_swap_persistnce2[ , c("SubjectID", "Diet", 
                "Age_Quartile", "Gender", "WC_change_Quartile", "Weight_change_Quartile", "Richness_Quartile", 
                "Swap_A_B_percentage")] %>% distinct()
swap_model <- lm(Swap_A_B_percentage ~ Diet + Age_Quartile + Gender + WC_change_Quartile + 
                     Weight_change_Quartile + Richness_Quartile, swap, na.action = na.exclude)
swap$predicted_value <- predict(swap_model)                   
lsmeans_swap <- emmeans(swap_model, ~ Diet)
lsmeans_swap <- as.data.frame(lsmeans_swap)[, c("Diet", "emmean")]
swap_lsmeans_swap <- swap %>% left_join(lsmeans_swap, by = "Diet") %>% distinct()
swap_lsmeans_swap <- swap_lsmeans_swap %>% mutate(adjusted_raw = Swap_A_B_percentage - (predicted_value - emmean))
swap_lsmeans_swap <- as.data.frame(swap_lsmeans_swap)
swap_lsmeans_swap <- swap_lsmeans_swap %>% rename(Value = Swap_A_B_percentage)

# Persistence
persistence <- gain_loss_swap_persistnce2[ , c("SubjectID", "Diet", 
                "Age_Quartile", "Gender", "WC_change_Quartile", "Weight_change_Quartile", "Richness_Quartile", 
                "Persistence_A_B_percentage")] %>% distinct()
persistence_model <- lm(Persistence_A_B_percentage ~ Diet + Age_Quartile + Gender + WC_change_Quartile + 
                     Weight_change_Quartile + Richness_Quartile, persistence, na.action = na.exclude)
persistence$predicted_value <- predict(persistence_model)                   
lsmeans_persistence <- emmeans(persistence_model, ~ Diet)
lsmeans_persistence <- as.data.frame(lsmeans_persistence)[, c("Diet", "emmean")]
persistence_lsmeans_persistence <- persistence %>% left_join(lsmeans_persistence, by = "Diet") %>% distinct()
persistence_lsmeans_persistence <- persistence_lsmeans_persistence %>% mutate(adjusted_raw = Persistence_A_B_percentage - (predicted_value - emmean))
persistence_lsmeans_persistence <- as.data.frame(persistence_lsmeans_persistence)
persistence_lsmeans_persistence <- persistence_lsmeans_persistence %>% rename(Value = Persistence_A_B_percentage)

loss_lsmeans_loss$Metric <- "Loss"
gain_lsmeans_gain$Metric <- "Gain"
swap_lsmeans_swap$Metric <- "Swap"
persistence_lsmeans_persistence$Metric <- "Persistence"

lsmeans_gain_loss_swap_persistence <- rbind(as.data.frame(gain_lsmeans_gain), 
                                as.data.frame(loss_lsmeans_loss), 
                                as.data.frame(swap_lsmeans_swap)
                                #as.data.frame(persistence_lsmeans_persistence)
                                )


lsmeans_gain_loss_swap_persistence$Diet <- factor(lsmeans_gain_loss_swap_persistence$Diet, levels = c("1","2","3"))
lsmeans_gain_loss_swap_persistence$Metric <- factor(lsmeans_gain_loss_swap_persistence$Metric, 
                                   levels = c("Gain", "Loss","Swap"))

plot_results <- lsmeans_gain_loss_swap_persistence[,c("SubjectID", "Diet", "Value", "predicted_value", "Metric")]


###############################################################################
# ggplot
library(ggplot2)

# Create a separate annotation data frame with p values of diet effects
annotation_df <- data.frame(
    Diet = factor("2", levels = c("1", "2", "3")),  # Must match the original factor levels
    predicted_value = c(50, 50, 7),  # Adjust this y-value if needed
    Metric = c("Gain", "Loss", "Swap"),
    label = c("p = 0.06", "p = 1.6e-04", "p = 0.19")
)

ggplot(plot_results, aes(x = factor(Diet), y = predicted_value, fill = factor(Diet))) +
    geom_boxplot(outlier.color = "grey", outlier.shape = 16, width = 0.7) +
    scale_fill_manual(
        values = c("1" = "#EED6B7", "2" = "#00788B", "3" = "#049B86")
    ) +
    scale_x_discrete(labels = c("1" = "HDG", "2" ="MED", "3" = "GreenMED")) +
    labs(
        title = "Strain Change induced by diets, pre-aFMT",
        x = "",
        y = "% of species"
    ) +
    facet_wrap(
        ~ Metric, scales = "free", nrow = 1,
        labeller = as_labeller(c(
            "Gain" = "Gain (species)\n(0-6)",
            "Loss" = "Loss (species)\n(0-6)",
            "Swap" = "Swap (strain)\n(0-6)"
        ))
    ) +
    geom_text(
        data = annotation_df,
        aes(x = Diet, y = predicted_value, label = label),
        inherit.aes = FALSE,
        size = 7
        #fontface = "bold"
    ) +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black", face = "bold"),
        axis.text.y = element_text(size = 20, color = "black", face = "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 18, face = "bold")
    )

 
# Output 10*6
library(openxlsx)
wb <- createWorkbook()

addWorksheet(wb, "Fig 2c Strain Gain Loss Swap")
writeData(wb, "Fig 2c Strain Gain Loss Swap", plot_results)

saveWorkbook(wb, "Source Data Fig 2c.xlsx", overwrite = TRUE)


