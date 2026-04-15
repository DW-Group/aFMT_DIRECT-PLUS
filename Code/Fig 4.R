
#################
# Load packages #
#################

options(warn=-1)
library(caret)
library(cowplot)
library(patchwork)
library(iml)
library(tidyverse)
#library(ggembl)
library(vegan)
#library(ggpubr)
#library(ggsignif)
library(mlr3)
library(mlr3learners)
library("mlr3tuning")
library(ranger)
library(mlr3viz)
library(pROC)
library("paradox")
library('missForest')
library(scales)


# Read the data
################################################################################################
# Dataset 1: Load biomarker data, Load the base characteristics of the participants
################################################################################################
setwd("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/8.phenotype_prediction")
base_biomarker_data <- read.table("metadata_biomarkers2.tsv", header=T, sep="\t", check.names = FALSE)
base_biomarker_data <- base_biomarker_data %>% select("sno_time","SubjectID","Time","Age","Weight","Gender","Diet","Treatment") 

# Filter data for Time points A and B
base_biomarker_data_A <- base_biomarker_data %>% filter(Time == "0")
base_biomarker_data_B <- base_biomarker_data %>% filter(Time == "6")
base_biomarker_data_D <- base_biomarker_data %>% filter(Time == "14")
base_biomarker_data_C <- base_biomarker_data %>% filter(Time == "18")

# Merge data for A and B by SubjectID
base_merged_biomarker_data_A_B <- merge(base_biomarker_data_A, base_biomarker_data_B, by = "SubjectID", suffixes = c("_A", "_B"), all = TRUE)
base_merged_biomarker_data_D_C <- merge(base_biomarker_data_D, base_biomarker_data_C, by = "SubjectID", suffixes = c("_D", "_C"), all = TRUE)
base_biomarker_data_reformat <- full_join(base_merged_biomarker_data_A_B, base_merged_biomarker_data_D_C, by="SubjectID" )

# Calculate weight change (weight at time 6 - weight at time 0)
base_biomarker_data_reformat$Weight_Baseline <- 0
base_biomarker_data_reformat <- base_biomarker_data_reformat %>%
    group_by(SubjectID) %>%
    mutate(Weight_Change0_6 = Weight_B - Weight_A, 
           Weight_Change6_14 = Weight_D - Weight_B,
           Weight_Change6_18 = Weight_C - Weight_B, 
           Weight_regain6_14 = Weight_Change6_14*100/abs(Weight_Change0_6), 
           Weight_regain6_18 = Weight_Change6_18*100/abs(Weight_Change0_6) ) %>%
           #filter(SubjectID != "124" & SubjectID != "125") %>% #Remove the 2 outliers
    ungroup()
dim(base_biomarker_data_reformat)
#[1]  90 35


################################################################################################
# Dataset 2: Load the microbiome species-level data
################################################################################################
species_data <- read.table("merged_abundance_table_filtered.tsv", header = TRUE, sep = "\t")
species_data <- species_data %>%
  mutate(
    Month = recode(Time, "A" = 0, "B" = 6, "C" = 14, "D" = 18),
    sno_time = paste(SubjectID, Month, sep = "_")
  ) %>%
  select(Month, sno_time, everything())

# Remove species (columns) with mean abundance < 0.01 from your data frame where rows are samples and columns are species
species_to_keep <- names(species_data)[8:ncol(species_data)][colMeans(species_data[, 8:ncol(species_data)]) >= 0.01]
microbiome_data <- cbind(species_data[, 1:7], species_data[, species_to_keep])
  
# Filter data for Time points A and B
microbiome_data_A <- microbiome_data %>% filter(Time == "A")
microbiome_data_B <- microbiome_data %>% filter(Time == "B")
microbiome_data_D <- microbiome_data %>% filter(Time == "D")
microbiome_data_C <- microbiome_data %>% filter(Time == "C")
dim(microbiome_data_D)
#[1]   77 529
dim(microbiome_data_B)
#[1]   79 529

# Find common row names
common_ids <- intersect(rownames(microbiome_data_D), rownames(microbiome_data_B))
# Subset both data frames to the common rows
D_common <- microbiome_data_D[common_ids, ]
B_common <- microbiome_data_B[common_ids, ]
# Now subtract elementwise
microbiome_change_BD <- D_common[, 8:length(microbiome_data_D)] - B_common[, 8:length(microbiome_data_B)]
microbiome_change_BD$"SubjectID" <- B_common$SubjectID
#merged_data_changeBD_C <- merge(change_B_D, microbiome_data_C, by = "SubjectID", suffixes = c("_changeBD", "_C"), all = TRUE)

# Rename all columns in microbiome_change_BD except "SubjectID"
microbiome_change_BD_renamed <- microbiome_change_BD %>%
  rename_with(~ paste0(.x, "_change_BD_B"), -SubjectID)
dim(microbiome_change_BD_renamed)
#[1]   77 523
# Now merge; since the columns from microbiome_change_BD already have the suffix, they won't be merged without the suffix.
predict_data2 <- merge(base_biomarker_data_reformat, microbiome_change_BD_renamed, by = "SubjectID")
dim(predict_data2)
#[1]   77 557


################################################################################################
# Dataset 3: Load the fecal metabolites data
################################################################################################
fecal_metabolites_data <- read.table("Batch_normalizedData_90p.tsv", header=T, sep="\t", check.names = FALSE)
#fecal_metabolites_data <- read.table("Fecal_metabolites_Batch-norm_Imputed_Data_90p.tsv", header=T, sep="\t", check.names = FALSE)
names(fecal_metabolites_data)[1] <- "SampleID"
fecal_metabolites_data[is.na(fecal_metabolites_data)] <- 0

# Remove metabolites (columns) with mean abundance < 1 from your data frame where rows are samples and columns are metabolites
#metabolites_to_keep <- names(fecal_metabolites_data)[9:ncol(fecal_metabolites_data)][colMeans(fecal_metabolites_data[, 9:ncol(fecal_metabolites_data)]) >= 1]
#fecal_metabolites_data <- cbind(fecal_metabolites_data[, 1:8], fecal_metabolites_data[, metabolites_to_keep])
dim(fecal_metabolites_data)
fecal_metabolites_data[1:5, 1:18]

###################################
# assign fecal chemical real names
library(readr)
fecal_chem_annotation <- read_tsv("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/Fecal_chem_annotation.tsv",
  na = c("", "NA", "N/A"), locale = locale(encoding = "UTF-8"), guess_max = 10000)
fecal_chem_annotation <- as.data.frame(fecal_chem_annotation)
# Make a copy of the original data just in case
fecal_metabolites_renamed <- fecal_metabolites_data
# Get the current column names from column 9 onward (which are CHEM_IDs)
chem_ids <- colnames(fecal_metabolites_renamed)[9:ncol(fecal_metabolites_renamed)]
# Create a named vector to map CHEM_ID to CHEMICAL_NAME
id_to_name <- setNames(fecal_chem_annotation$CHEMICAL_NAME, 
                       as.character(fecal_chem_annotation$CHEM_ID))
# Replace column names using the mapping
new_names <- sapply(chem_ids, function(id) {
  # Return CHEMICAL_NAME if it exists, otherwise keep the original
  if (!is.null(id_to_name[[id]])) id_to_name[[id]] else id
})
# Use grepl to identify column names that look like "X-12345" or are purely numeric
# These typically represent unnamed metabolites
cols_to_remove <- grepl("^X-\\d+$", unname(new_names))
# Keep columns that are not matches OR are part of the metadata (first 8 columns)
fecal_metabolites_cleaned <- fecal_metabolites_renamed[, !cols_to_remove | seq_along(cols_to_remove) <= 8]
dim(fecal_metabolites_cleaned)
#[1] 331 1188

# Filter data for Time points A, B, D, and C.
fecal_metabolites_data_A <- fecal_metabolites_cleaned %>% filter(GROUP_NUMBER == "A")
fecal_metabolites_data_B <- fecal_metabolites_cleaned %>% filter(GROUP_NUMBER == "B")
fecal_metabolites_data_D <- fecal_metabolites_cleaned %>% filter(GROUP_NUMBER == "D")
fecal_metabolites_data_C <- fecal_metabolites_cleaned %>% filter(GROUP_NUMBER == "C")
dim(fecal_metabolites_data_B)
#[1]   80 1188
dim(fecal_metabolites_data_D)
#[1]   77 1188

# Find common row names
common_ids <- intersect(rownames(fecal_metabolites_data_D), rownames(fecal_metabolites_data_B))
# Subset both data frames to the common rows
fecal_metabolites_D_common <- fecal_metabolites_data_D[common_ids, ]
fecal_metabolites_B_common <- fecal_metabolites_data_B[common_ids, ]
# Now subtract elementwise
fecal_metabolites_change_BD <- fecal_metabolites_D_common[, 9:length(fecal_metabolites_data_D)] - fecal_metabolites_B_common[, 9:length(fecal_metabolites_data_B)]
fecal_metabolites_change_BD$"SubjectID" <- fecal_metabolites_B_common$SubjectID
#merged_data_changeBD_C <- merge(fecal_metabolites_change_B_D, fecal_metabolites_data_C, by = "SubjectID", suffixes = c("_changeBD", "_C"), all = TRUE)

# Rename all columns in microbiome_change_BD except "SubjectID"
fecal_metabolites_change_BD_renamed <- fecal_metabolites_change_BD %>%
  rename_with(~ paste0(.x, "_change_BD_B"), -SubjectID)
dim(fecal_metabolites_change_BD_renamed)

# Merge predict_data3 and fecal_metabolites_data
predict_data3 <- merge(predict_data2, fecal_metabolites_change_BD_renamed, by = "SubjectID")
dim(predict_data3)
#[1]    74 1737



################################################################################################
# Dataset 4: Load the strain-level change data
################################################################################################
##############################  Gain and Loss (Species level)     ###############################

#setwd("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species")
# Load necessary library
library(dplyr)

# Read the data
data <- read.table("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/merged_abundance_table_filtered.tsv", header = TRUE, sep = "\t")
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

# Remove species (columns) with mean abundance < 0.01 from your data frame where rows are samples and columns are species
species_to_keep <- names(data)[8:ncol(data)][colMeans(data[, 8:ncol(data)]) >= 0.01]
data <- cbind(data[, 1:7], data[, species_to_keep])
dim(data)
#[1] 327 528

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
#[1]   87 2109


##########################################################################################
######################### Strain level: Swap and Persistence ###################
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
# Use the subset of strain profiling data which contains at least 20% (18) subjects
#folder_path <- "~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/5.matrix/"  # Update with your folder path
folder_path <- "~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/5.matrix2_18subjects/"  # Update with your folder path

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
# Filter data for Time points A and B
final_results_A <- final_results %>% filter(TimePoint == "A")
final_results_B <- final_results %>% filter(TimePoint == "B")
final_results_D <- final_results %>% filter(TimePoint == "D")
final_results_C <- final_results %>% filter(TimePoint == "C")

library(dplyr)
library(tidyr)
library(stringr)

# Filter to only the desired timepoints: B, C, and D
bdc_data <- final_results %>% 
  filter(TimePoint %in% c("B", "C", "D")) %>% 
  arrange(SubjectID, TimePoint)

# Pivot wider: one row per SubjectID, and columns for each t__SGB variable by TimePoint
wide_df <- bdc_data %>% 
  pivot_wider(id_cols = SubjectID, 
              names_from = TimePoint, 
              values_from = starts_with("t__SGB"))

# Get the base names of the SGB columns by removing the trailing _B, _C, or _D suffix.
sgb_base_names <- wide_df %>% 
  select(-SubjectID) %>% 
  names() %>% 
  sub("_[BCD]$", "", .) %>% 
  unique()

# Loop over each base name and calculate the differences for B-D and B-C
change_df <- wide_df
for(base in sgb_base_names) {
  colB <- paste0(base, "_B")
  colC <- paste0(base, "_C")
  colD <- paste0(base, "_D")
  
  # Compute change_BD_B: assign 1 if B != D, otherwise 0.
  if(all(c(colB, colD) %in% names(wide_df))) {
    #change_df[[paste0(base, "_strain_change_BD_B")]] <- ifelse(wide_df[[colB]] == wide_df[[colD]], 0, 1)
    change_df[[paste0(base, "_strain_change_BD_B")]] <- ifelse(
    wide_df[[colB]] > 0 & wide_df[[colD]] > 0,
    ifelse(wide_df[[colB]] == wide_df[[colD]], 0, 1),
    ifelse(
      species_merged_data[[paste0(base, "_species_B")]] > 0 & 
      species_merged_data[[paste0(base, "_species_D")]] > 0, 2, 1)
      )
   }
  # Compute change_BC_B: assign 1 if B != C, otherwise 0.
  #if(all(c(colB, colC) %in% names(wide_df))) {
  #  change_df[[paste0(base, "_strain_change_BC_B")]] <- ifelse(wide_df[[colB]] != wide_df[[colC]], 1, 0)
  #}
}

# Replace all NA values in change_df with 1
change_df[is.na(change_df)] <- 1

# Keep only the SubjectID and the newly computed change columns
strain_change_df <- change_df %>% 
  select(SubjectID, ends_with("_strain_change_BD_B"))

dim(strain_change_df)
#[1]  90 203

########## Strain_change_df will be used for prediction ###################



##############################################################################
##################### Merge Gain Loss Swap Persistence #######################
##############################################################################
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
final_results_merged_data_A_B <- merge(final_results_A[,-2], final_results_B[,-2], by = "SubjectID", suffixes = c("_strain_A", "_strain_B"))
final_results_merged_data_D_C <- merge(final_results_D[,-2], final_results_C[,-2], by = "SubjectID", suffixes = c("_strain_D", "_strain_C"))

strain_merged_data <- full_join(final_results_merged_data_A_B, final_results_merged_data_D_C, by="SubjectID" )
dim(strain_merged_data)
#[1]   90 7241
strain_merged_data[1:5, 1:6]
#  SubjectID t__SGB10068_strain_A t__SGB10115_group_strain_A t__SGB10130_strain_A
#1         4                    0                          0                    0
#2        13                    0                          0                    0
#3        15                    0                          0                    0
#4        17                    0                          0                    0
#5        19                    0                          0                    0


##################### Merge Gain Loss Swap Persistence #######################
full_merged_data <- merge(strain_merged_data, species_merged_data, by="SubjectID", all = TRUE )
dim(full_merged_data)
#[1]    90 14497

# Identify species SGB columns
species_columns <- grep("^t__", colnames(data), value = TRUE)
strain_columns <- grep("^t__", colnames(final_results), value = TRUE)


# Initialize result data frame
strain_results <- full_merged_data %>% 
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
         #Total_species_A, Total_species_B, Total_species_D, Total_species_C,
         #Gain_A_B_percentage, Gain_B_D_percentage, Gain_D_C_percentage,
         Gain_B_D_percentage,
         #Loss_A_B_percentage,Loss_B_D_percentage, Loss_D_C_percentage,
         Loss_B_D_percentage, 
         #Swap_count_A_B, Swap_count_B_D, Swap_count_D_C, 
         #Swap_A_B_percentage, Swap_B_D_percentage, Swap_D_C_percentage,
         Swap_B_D_percentage, 
         #Persistence_A_B_percentage, Persistence_B_D_percentage, Persistence_D_C_percentage,
         #Persistence_B_D_percentage,
         #QGain_A_B_percentage, QGain_B_D_percentage, QGain_D_C_percentage,
         #QGain_B_D_percentage,
         #QLoss_A_B_percentage, QLoss_B_D_percentage, QLoss_D_C_percentage,
         #QAbsence_A_B_percentage, QAbsence_B_D_percentage, QAbsence_D_C_percentage,
         #QAbsence_B_D_percentage, 
         #Gain_Persistence_percentage, Loss_Persistence_percentage, Swap_Persistence_percentage,
         Succession_rate
)

strain_results <- as.data.frame(strain_results)
dim(strain_results)
#[1] 90 5




################################################################################################
# Model 1: RF based on Biomarker data
################################################################################################
# MRMR is a feature selection method that chooses features with maximum relevance to the target while minimizing redundancy among them.
# Set a seed for reproducibility 
#[1] 115 120 195 225 255 285
set.seed(285)
#set.seed(666)
#set.seed(888888)

library(mRMRe)

# Example data: replace 'df' with your data frame.
# Ensure all columns are numeric. If not, you might need to convert or remove non-numeric columns.
# Also, mRMRe expects a matrix with numeric values.
predict_data0 <- base_biomarker_data_reformat[, !names(base_biomarker_data_reformat) %in% c("SubjectID", "sno_time_A", "Time_A", "sno_time_B", "Time_B", #"Gender_B", "Diet_B", "Treatment_B",
                                                               "sno_time_D", "Time_D", #"Gender_D", "Diet_D", "Treatment_D",
                                                               "sno_time_C", "Time_C", #"Gender_C", "Diet_C", "Treatment_C", 
                                                               "Weight_Baseline", "Weight_Change0_6",
                                                               "Weight_Change6_14", "Weight_regain6_14",
                                                               "Weight_Change6_18")] #"Weight_regain6_18")]

dim(predict_data0) 
#[1] 90 61

predictors <- predict_data0 %>% select(Weight_regain6_18, ends_with("_B")) 
predictors$Treatment_B <- ifelse(predictors$Treatment_B == "aFMT", 1, 0)
predictors <- predictors %>% filter(Weight_regain6_18!= "NA")
dim(predictors)
#[1] 88 13


# Initialize an empty data frame to store the performance metrics
results_df <- data.frame(feature_count = integer(),
                         Group = character(),
                         Model = character(),
                         Accuracy = numeric(),
                         AUC = numeric(),
                         Sensitivity = numeric(),
                         Specificity = numeric(),
                         stringsAsFactors = FALSE)

groups <- c("1", "0")                   
for (grp in groups) {
    # 1) group data
    predictors_grp <- predictors %>% filter(Treatment_B == grp)
    
    # Create an mRMRe Data object. 
    # IMPORTANT: mRMRe requires column names and numeric data.
    data_matrix <- as.data.frame(lapply(predictors_grp, as.numeric))
    data_matrix$Weight_regain6_18 <- ifelse(data_matrix$Weight_regain6_18 > median(data_matrix$Weight_regain6_18), 1, 0)
    mrmr_data <- mRMR.data(data = data_matrix)
    
    # Identify the index of the target variable column (e.g., column named "target")
    target_index <- which(colnames(data_matrix) == "Weight_regain6_18")
    
    # Perform mRMR feature selection: select the top # features.
    # Note: adjust feature_count as needed.
                             
    ########################################################
    # --- Random Forest Regression with mlr3 using the Selected Features ---
    
    # Uses mlr3 and the mlr3learners package to perform Random Forest Regression (using the ranger learner) to predict weight regain percentage
    # Load required packages
    library(mlr3)
    library(mlr3learners)
    library(mlr3measures)
    library(ggplot2)
    
    # Define the range of feature counts and number of repeats
    feature_counts <- seq(2, 5, by = 1)
    
    # Loop over each feature count
    for(fc in feature_counts) {
        # Perform mRMR feature selection: select the top # features.
        # Note: adjust feature_count as needed.
        feature_selection <- mRMR.classic(data = mrmr_data, target_indices = target_index, feature_count = fc)
        
        # Extract selected feature names:
        selected_features <- solutions(feature_selection)[[1]]
        #print(selected_features)
        
        data_matrix_selected <- data_matrix %>% select(Weight_regain6_18, selected_features[,1])
        
        # Classification task to predict weight regain as binary (0/1)
        taskClas <- TaskClassif$new(
          id = "rf_classification",
          backend = data_matrix_selected %>% 
            filter(!is.na(Weight_regain6_18)) %>% 
            #mutate(Weight_regain6_18 = ifelse(Weight_regain6_18 > median(Weight_regain6_18), "1", "0")) %>% 
            mutate(Weight_regain6_18 = as.factor(Weight_regain6_18)),
          target = "Weight_regain6_18"
        )
        
        # --- Initialize and Train Models ---
        # Initialize a Random Forest Classification learner (using ranger)
        n_pred <- ncol(data_matrix_selected) - 1  # excluding target
        rfModel_class <- lrn("classif.ranger", 
                         importance = "permutation",
                         mtry = round(sqrt(n_pred)),
                         min.node.size = 5,
                         sample.fraction = 0.8,
                         num.threads = parallel::detectCores())
        rfModel_class$predict_type <- "prob"
        
        # --- Set up 5-fold cross-validation ---
        cv5 <- rsmp("cv", folds = 5)
        cv5$instantiate(taskClas)
        
        # --- Run 10-fold CV resampling ---
        resample_result <- resample(taskClas, rfModel_class, cv5)
        
        # --- Aggregate performance metrics ---
        accuracy    <- resample_result$aggregate(msr("classif.acc"))
        auc         <- resample_result$aggregate(msr("classif.auc"))
        sensitivity <- resample_result$aggregate(msr("classif.sensitivity"))
        specificity <- resample_result$aggregate(msr("classif.specificity"))
        
        #cat("10-fold CV results:\n")
        cat("Accuracy:", round(accuracy, 3), "\n")
        cat("AUC:", round(auc, 3), "\n")
        cat("Sensitivity:", round(sensitivity, 3), "\n")
        cat("Specificity:", round(specificity, 3), "\n")
        
        # --- Save the metrics for this repetition and feature count ---
        results_df <- rbind(results_df,
                            data.frame(feature_count = fc,
                                       Group = grp,
                                       Model = "Model 1: Metadata",
                                       Accuracy = accuracy,
                                       AUC = auc,
                                       Sensitivity = sensitivity,
                                       Specificity = specificity))
        
        }
    }
print(results_df)



################################################################################################
# Model 2: RF based on Metadata + Microbiome change
################################################################################################

library(mRMRe)

# Example data: replace 'df' with your data frame.
# Ensure all columns are numeric. If not, you might need to convert or remove non-numeric columns.
# Also, mRMRe expects a matrix with numeric values.
#predict_data <- merge(base_biomarker_data_reformat, microbiome_change_BD_renamed, by = "SubjectID")
predict_data <- predict_data2
dim(predict_data)
# 77 982
predict_data <- predict_data[, !names(base_biomarker_data_reformat) %in% c("SubjectID", "sno_time_A", "Time_A", 
                                                               "sno_time_B", "Time_B", #"Gender_B", "Diet_B", "Treatment_B",
                                                               "sno_time_D", "Time_D", #"Gender_D", "Diet_D", "Treatment_D",
                                                               "sno_time_C", "Time_C", #"Gender_C", "Diet_C", "Treatment_C", 
                                                               "Weight_Baseline", "Weight_Change0_6",
                                                               "Weight_Change6_14", "Weight_regain6_14",
                                                               "Weight_Change6_18")]#, "Weight_regain6_18")]

dim(predict_data) 
#[1] 77 2171

predictors <- predict_data %>% select(Weight_regain6_18, ends_with("_B"))
predictors$Treatment_B <- ifelse(predictors$Treatment_B == "aFMT", 1, 0)
predictors <- predictors %>% filter(Weight_regain6_18!= "NA")
dim(predictors)
#[1] 75 1816

groups <- c("1", "0")                   
for (grp in groups) {
    # 1) group data
    predictors_grp <- predictors %>% filter(Treatment_B == grp)
        
    # Create an mRMRe Data object. 
    # IMPORTANT: mRMRe requires column names and numeric data.
    data_matrix <- as.data.frame(lapply(predictors_grp, as.numeric))
    data_matrix$Weight_regain6_18 <- ifelse(data_matrix$Weight_regain6_18 > median(data_matrix$Weight_regain6_18), 1, 0)
    mrmr_data <- mRMR.data(data = data_matrix)
    
    # Identify the index of the target variable column (e.g., column named "target")
    target_index <- which(colnames(data_matrix) == "Weight_regain6_18")
    
                      
    ########################################################
    # --- Random Forest Regression with mlr3 using the Selected Features ---
    
    # Uses mlr3 and the mlr3learners package to perform Random Forest Regression (using the ranger learner) to predict weight regain percentage
    # Load required packages
    library(mlr3)
    library(mlr3learners)
    library(mlr3measures)
    library(ggplot2)
    
    # Set a seed for reproducibility
    #set.seed(111)
    # Define the range of feature counts and number of repeats
    feature_counts <- seq(10, 100, by = 10)
    
    # Loop over each feature count
    for(fc in feature_counts) {
        # Perform mRMR feature selection: select the top # features.
        # Note: adjust feature_count as needed.
        feature_selection <- mRMR.classic(data = mrmr_data, target_indices = target_index, feature_count = fc)
        
        # Extract selected feature names:
        selected_features <- solutions(feature_selection)[[1]]
        #print(selected_features)
        
        data_matrix_selected <- data_matrix %>% select(Weight_regain6_18, selected_features[,1])
        
        # Classification task to predict weight regain as binary (0/1)
        taskClas <- TaskClassif$new(
          id = "rf_classification",
          backend = data_matrix_selected %>% 
            filter(!is.na(Weight_regain6_18)) %>% 
            #mutate(Weight_regain6_18 = ifelse(Weight_regain6_18 > median(Weight_regain6_18), "1", "0")) %>% 
            mutate(Weight_regain6_18 = as.factor(Weight_regain6_18)),
          target = "Weight_regain6_18"
        )
        
        # --- Initialize and Train Models ---
        # Initialize a Random Forest Classification learner (using ranger)
        n_pred <- ncol(data_matrix_selected) - 1  # excluding target
        rfModel_class <- lrn("classif.ranger", 
                         importance = "permutation",
                         mtry = round(sqrt(n_pred)),
                         min.node.size = 5,
                         sample.fraction = 0.8,
                         num.threads = parallel::detectCores())
        rfModel_class$predict_type <- "prob"
        
        # --- Set up 5-fold cross-validation ---
        cv5 <- rsmp("cv", folds = 5)
        cv5$instantiate(taskClas)
        
        # --- Run 10-fold CV resampling ---
        resample_result <- resample(taskClas, rfModel_class, cv5)
        
        # --- Aggregate performance metrics ---
        accuracy    <- resample_result$aggregate(msr("classif.acc"))
        auc         <- resample_result$aggregate(msr("classif.auc"))
        sensitivity <- resample_result$aggregate(msr("classif.sensitivity"))
        specificity <- resample_result$aggregate(msr("classif.specificity"))
        
        #cat("10-fold CV results:\n")
        cat("Accuracy:", round(accuracy, 3), "\n")
        cat("AUC:", round(auc, 3), "\n")
        cat("Sensitivity:", round(sensitivity, 3), "\n")
        cat("Specificity:", round(specificity, 3), "\n")
        
        # --- Save the metrics for this repetition and feature count ---
        results_df <- rbind(results_df,
                            data.frame(feature_count = fc,
                                       Group = grp,
                                       Model = "Model 2: 1 + Species change (6 vs 14)",
                                       Accuracy = accuracy,
                                       AUC = auc,
                                       Sensitivity = sensitivity,
                                       Specificity = specificity))
        
        }
    }

print(results_df)


################################################################################################
# Model 3: RF based on Metadata + Microbiome change + Fecal Metabolite change
################################################################################################

library(mRMRe)

# Example data: replace 'df' with your data frame.
# Ensure all columns are numeric. If not, you might need to convert or remove non-numeric columns.
# Also, mRMRe expects a matrix with numeric values.
predict_data <- predict_data3 #, fecal_metabolites_change_BD_renamed, by = "SubjectID")
dim(predict_data)
# 77 2162
predict_data <- predict_data[, !names(base_biomarker_data_reformat) %in% c("SubjectID", "sno_time_A", "Time_A", 
                                                               "sno_time_B", "Time_B", #"Gender_B", "Diet_B", "Treatment_B",
                                                               "sno_time_D", "Time_D", #"Gender_D", "Diet_D", "Treatment_D",
                                                               "sno_time_C", "Time_C", #"Gender_C", "Diet_C", "Treatment_C", 
                                                               "Weight_Baseline", "Weight_Change0_6",
                                                               "Weight_Change6_14", "Weight_regain6_14",
                                                               "Weight_Change6_18")]#, "Weight_regain6_18")]

dim(predict_data) 
#[1] 74 1709

predictors <- predict_data %>% select(Weight_regain6_18, ends_with("_B"))
predictors$Treatment_B <- ifelse(predictors$Treatment_B == "aFMT", 1, 0)
predictors <- predictors %>% filter(Weight_regain6_18!= "NA")
dim(predictors)
#[1] 72 1462

groups <- c("1", "0")                   
for (grp in groups) {
    # 1) group data
    predictors_grp <- predictors %>% filter(Treatment_B == grp)
        
    # Create an mRMRe Data object. 
    # IMPORTANT: mRMRe requires column names and numeric data.
    data_matrix <- as.data.frame(lapply(predictors_grp, as.numeric))
    data_matrix$Weight_regain6_18 <- ifelse(data_matrix$Weight_regain6_18 > median(data_matrix$Weight_regain6_18), 1, 0)
    
    mrmr_data <- mRMR.data(data = data_matrix)
    
    # Identify the index of the target variable column (e.g., column named "target")
    target_index <- which(colnames(data_matrix) == "Weight_regain6_18")
    
                      
    ########################################################
    # --- Random Forest Regression with mlr3 using the Selected Features ---
    
    # Uses mlr3 and the mlr3learners package to perform Random Forest Regression (using the ranger learner) to predict weight regain percentage
    # Load required packages
    library(mlr3)
    library(mlr3learners)
    library(mlr3measures)
    library(ggplot2)
    
    # Set a seed for reproducibility
    #set.seed(111)
    # Define the range of feature counts and number of repeats
    feature_counts <- seq(10, 100, by = 10)
    
    # Loop over each feature count
    for(fc in feature_counts) {
        # Perform mRMR feature selection: select the top # features.
        # Note: adjust feature_count as needed.
        feature_selection <- mRMR.classic(data = mrmr_data, target_indices = target_index, feature_count = fc)
        
        # Extract selected feature names:
        selected_features <- solutions(feature_selection)[[1]]
        #print(selected_features)
        
        data_matrix_selected <- data_matrix %>% select(Weight_regain6_18, selected_features[,1])
        
        # Classification task to predict weight regain as binary (0/1)
        taskClas <- TaskClassif$new(
          id = "rf_classification",
          backend = data_matrix_selected %>% 
            filter(!is.na(Weight_regain6_18)) %>% 
            #mutate(Weight_regain6_18 = ifelse(Weight_regain6_18 > median(Weight_regain6_18), "1", "0")) %>% 
            mutate(Weight_regain6_18 = as.factor(Weight_regain6_18)),
          target = "Weight_regain6_18"
        )
        
        # --- Initialize and Train Models ---
        # Initialize a Random Forest Classification learner (using ranger)
        n_pred <- ncol(data_matrix_selected) - 1  # excluding target
        rfModel_class <- lrn("classif.ranger", 
                         importance = "permutation",
                         mtry = round(sqrt(n_pred)),
                         min.node.size = 5,
                         sample.fraction = 0.8,
                         num.threads = parallel::detectCores())
        rfModel_class$predict_type <- "prob"
        
        # --- Set up 5-fold cross-validation ---
        cv5 <- rsmp("cv", folds = 5)
        cv5$instantiate(taskClas)
        
        # --- Run 10-fold CV resampling ---
        resample_result <- resample(taskClas, rfModel_class, cv5)
        
        # --- Aggregate performance metrics ---
        accuracy    <- resample_result$aggregate(msr("classif.acc"))
        auc         <- resample_result$aggregate(msr("classif.auc"))
        sensitivity <- resample_result$aggregate(msr("classif.sensitivity"))
        specificity <- resample_result$aggregate(msr("classif.specificity"))
        
        #cat("10-fold CV results:\n")
        cat("Accuracy:", round(accuracy, 3), "\n")
        cat("AUC:", round(auc, 3), "\n")
        cat("Sensitivity:", round(sensitivity, 3), "\n")
        cat("Specificity:", round(specificity, 3), "\n")
        
        # --- Save the metrics for this repetition and feature count ---
        results_df <- rbind(results_df,
                            data.frame(feature_count = fc,
                                       Group = grp,
                                       Model = "Model 3: 2 + Metabolites change (6 vs 14)",
                                       Accuracy = accuracy,
                                       AUC = auc,
                                       Sensitivity = sensitivity,
                                       Specificity = specificity))
        
        }
    }
print(results_df)


################################################################################################
# Model 4: RF based on Metadata + Microbiome change + Strain change data
################################################################################################

library(mRMRe)
# Ensure all columns are numeric. If not, you might need to convert or remove non-numeric columns.
# Also, mRMRe expects a matrix with numeric values.
 
predict_data4 <- merge(predict_data2, strain_change_df, by = "SubjectID")
dim(predict_data4)
#[1]    77 1184

# Rename all columns in strain_results except "SubjectID" by appending the suffix "_strain_B"
strain_results_renamed <- strain_results %>%
  rename_with(~ paste0(.x, "_strain_B"), -SubjectID)
# Merge the data frames; now the strain_results columns already have the suffix.
predict_data6 <- merge(predict_data4, strain_results_renamed, by = "SubjectID")

dim(predict_data6)
#[1]    90 674
rownames(predict_data6) <- predict_data6$SubjectID


predict_data7 <- predict_data6[, !names(predict_data6) %in% c("SubjectID", "sno_time_A", "Time_A", 
                                                               "sno_time_B", "Time_B", #"Gender_B", "Diet_B", "Treatment_B",
                                                               "sno_time_D", "Time_D", #"Gender_D", "Diet_D", "Treatment_D",
                                                               "sno_time_C", "Time_C", #"Gender_C", "Diet_C", "Treatment_C", 
                                                               "Weight_Baseline", "Weight_Change0_6",
                                                               "Weight_Change6_14", "Weight_regain6_14",
                                                               "Weight_Change6_18")]#, "Weight_regain6_18")]

predict_data7$Treatment_A <- ifelse(predict_data7$Treatment_A == "aFMT", 1,
                              ifelse(predict_data7$Treatment_A == "Placebo", 0, predict_data7$Treatment_A))
predict_data7$Treatment_B <- ifelse(predict_data7$Treatment_B == "aFMT", 1,
                              ifelse(predict_data7$Treatment_B == "Placebo", 0, predict_data7$Treatment_B))
predict_data7$Treatment_D <- ifelse(predict_data7$Treatment_D == "aFMT", 1,
                              ifelse(predict_data7$Treatment_D == "Placebo", 0, predict_data7$Treatment_D))                              
predict_data7$Treatment_C <- ifelse(predict_data7$Treatment_C == "aFMT", 1,
                              ifelse(predict_data7$Treatment_C == "Placebo", 0, predict_data7$Treatment_C))
predict_data7$Treatment_A <- as.numeric(predict_data7$Treatment_A)
predict_data7$Treatment_B <- as.numeric(predict_data7$Treatment_B)
predict_data7$Treatment_D <- as.numeric(predict_data7$Treatment_D)
predict_data7$Treatment_C <- as.numeric(predict_data7$Treatment_C)
predict_data7[is.na(predict_data7)] <- 0
dim(predict_data7)
#[1]    90 660

predict_data <- predict_data7 %>% select(Weight_regain6_18, ends_with("_B"))
predictors$Treatment_B <- ifelse(predictors$Treatment_B == "aFMT", 1, 0)
dim(predict_data) 
#[1] 90 353

predictors <- predict_data %>% filter(Weight_regain6_18!= "NA")
dim(predictors)
#[1] 90 353


groups <- c("1", "0")                   
for (grp in groups) {
    # 1) group data
    predictors_grp <- predictors %>% filter(Treatment_B == grp)
        
    # Create an mRMRe Data object. 
    # IMPORTANT: mRMRe requires column names and numeric data.
    data_matrix <- as.data.frame(lapply(predictors_grp, as.numeric))
    data_matrix$Weight_regain6_18 <- ifelse(data_matrix$Weight_regain6_18 > median(data_matrix$Weight_regain6_18), 1, 0)
    
    mrmr_data <- mRMR.data(data = data_matrix)
    
    # Identify the index of the target variable column (e.g., column named "target")
    target_index <- which(colnames(data_matrix) == "Weight_regain6_18")
    
    
    ########################################################
    # --- Random Forest Regression with mlr3 using the Selected Features ---
    
    # Uses mlr3 and the mlr3learners package to perform Random Forest Regression (using the ranger learner) to predict weight regain percentage
    # Load required packages
    library(mlr3)
    library(mlr3learners)
    library(mlr3measures)
    library(ggplot2)
    
    # Set a seed for reproducibility
    #set.seed(111)
    # Define the range of feature counts and number of repeats
    feature_counts <- seq(10, 100, by = 10)
    
    #results_df <- results_df %>% filter(Model != "Model6: Model2 + Strain change (6 vs 14)")
    # Loop over each feature count
    for(fc in feature_counts) {
        # Perform mRMR feature selection: select the top # features.
        # Note: adjust feature_count as needed.
        feature_selection <- mRMR.classic(data = mrmr_data, target_indices = target_index, feature_count = fc)
        
        # Extract selected feature names:
        selected_features <- solutions(feature_selection)[[1]]
        #print(selected_features)
        
        data_matrix_selected <- data_matrix %>% select(Weight_regain6_18, selected_features[,1])
        
        # Classification task to predict weight regain as binary (0/1)
        taskClas <- TaskClassif$new(
          id = "rf_classification",
          backend = data_matrix_selected %>% 
            filter(!is.na(Weight_regain6_18)) %>% 
            #mutate(Weight_regain6_18 = ifelse(Weight_regain6_18 > median(Weight_regain6_18), "1", "0")) %>% 
            mutate(Weight_regain6_18 = as.factor(Weight_regain6_18)),
          target = "Weight_regain6_18"
        )
        
        # Initialize a Random Forest Classification learner with some tuned parameters:
        # For example, mtry is set to the square root of the number of predictors,
        # min.node.size to 5, sample.fraction to 0.8, and using available cores.
        n_pred <- ncol(data_matrix_selected) - 1  # excluding target
        rfModel_class <- lrn("classif.ranger", 
                         importance = "permutation",
                         mtry = round(sqrt(n_pred)),
                         min.node.size = 5,
                         sample.fraction = 0.8,
                         num.threads = parallel::detectCores())
        rfModel_class$predict_type <- "prob"
        
        # --- Set up 5-fold cross-validation ---
        cv5 <- rsmp("cv", folds = 5)
        cv5$instantiate(taskClas)
        
        # --- Run 10-fold CV resampling ---
        resample_result <- resample(taskClas, rfModel_class, cv5)
        
        # --- Aggregate performance metrics ---
        accuracy    <- resample_result$aggregate(msr("classif.acc"))
        auc         <- resample_result$aggregate(msr("classif.auc"))
        sensitivity <- resample_result$aggregate(msr("classif.sensitivity"))
        specificity <- resample_result$aggregate(msr("classif.specificity"))
        
        #cat("10-fold CV results:\n")
        cat("Accuracy:", round(accuracy, 3), "\n")
        cat("AUC:", round(auc, 3), "\n")
        cat("Sensitivity:", round(sensitivity, 3), "\n")
        cat("Specificity:", round(specificity, 3), "\n")
        
        # --- Save the metrics for this repetition and feature count ---
        results_df <- rbind(results_df,
                            data.frame(feature_count = fc,
                                       Group = grp,
                                       Model = "Model 4: 2 + Strain type change (6 vs 14)",
                                       Accuracy = accuracy,
                                       AUC = auc,
                                       Sensitivity = sensitivity,
                                       Specificity = specificity))
        
        }
    }
print(results_df)



################################################################################################
# Model 5- ALL: RF based on Metadata + Microbiome + Fecal Metabolites + Strain data
################################################################################################

library(mRMRe)
# Ensure all columns are numeric. If not, you might need to convert or remove non-numeric columns.
# Also, mRMRe expects a matrix with numeric values.
 
predict_data5 <- merge(predict_data3, strain_change_df, by = "SubjectID")
dim(predict_data5)
#[1]    76 2294

# Rename all columns in strain_results except "SubjectID" by appending the suffix "_strain_B"
strain_results_renamed <- strain_results %>%
  rename_with(~ paste0(.x, "_strain_B"), -SubjectID)
# Merge the data frames; now the strain_results columns already have the suffix.
predict_data6 <- merge(predict_data5, strain_results_renamed, by = "SubjectID")

dim(predict_data6)
#[1]    76 2298
rownames(predict_data6) <- predict_data6$SubjectID


predict_data7 <- predict_data6[, !names(predict_data6) %in% c("SubjectID", "sno_time_A", "Time_A", 
                                                               "sno_time_B", "Time_B", #"Gender_B", "Diet_B", "Treatment_B",
                                                               "sno_time_D", "Time_D", #"Gender_D", "Diet_D", "Treatment_D",
                                                               "sno_time_C", "Time_C", #"Gender_C", "Diet_C", "Treatment_C", 
                                                               "Weight_Baseline", "Weight_Change0_6",
                                                               "Weight_Change6_14", "Weight_regain6_14",
                                                               "Weight_Change6_18")]#, "Weight_regain6_18")]

predict_data7$Treatment_A <- ifelse(predict_data7$Treatment_A == "aFMT", 1,
                              ifelse(predict_data7$Treatment_A == "Placebo", 0, predict_data7$Treatment_A))
predict_data7$Treatment_B <- ifelse(predict_data7$Treatment_B == "aFMT", 1,
                              ifelse(predict_data7$Treatment_B == "Placebo", 0, predict_data7$Treatment_B))
predict_data7$Treatment_D <- ifelse(predict_data7$Treatment_D == "aFMT", 1,
                              ifelse(predict_data7$Treatment_D == "Placebo", 0, predict_data7$Treatment_D))                              
predict_data7$Treatment_C <- ifelse(predict_data7$Treatment_C == "aFMT", 1,
                              ifelse(predict_data7$Treatment_C == "Placebo", 0, predict_data7$Treatment_C))
predict_data7$Treatment_A <- as.numeric(predict_data7$Treatment_A)
predict_data7$Treatment_B <- as.numeric(predict_data7$Treatment_B)
predict_data7$Treatment_D <- as.numeric(predict_data7$Treatment_D)
predict_data7$Treatment_C <- as.numeric(predict_data7$Treatment_C)
predict_data7[is.na(predict_data7)] <- 0
dim(predict_data7)
#[1]    74 2284

predict_data <- predict_data7 %>% select(Weight_regain6_18, ends_with("_B"))
dim(predict_data) 
#[1] 74 1977

predictors <- predict_data %>% filter(Weight_regain6_18!= "NA")
dim(predictors)
#[1] 74 1977


roc_df <- data.frame(
              Group = character(),
              FPR = numeric(),
              TPR = numeric(),
              AUC = numeric(),
              FeatureCount = character(),
            stringsAsFactors = FALSE)
            
groups <- c("1", "0")                   
for (grp in groups) {
    # 1) group data
    predictors_grp <- predictors %>% filter(Treatment_B == grp)
    
    # Create an mRMRe Data object. 
    # IMPORTANT: mRMRe requires column names and numeric data.
    data_matrix <- as.data.frame(lapply(predictors_grp, as.numeric))
    data_matrix$Weight_regain6_18 <- ifelse(data_matrix$Weight_regain6_18 > median(data_matrix$Weight_regain6_18), 1, 0)
    #data_matrix$Weight_regain6_18 <- as.factor(data_matrix$Weight_regain6_18)
    
    mrmr_data <- mRMR.data(data = data_matrix)
    
    # Identify the index of the target variable column (e.g., column named "target")
    target_index <- which(colnames(data_matrix) == "Weight_regain6_18")
    
    
    ########################################################
    # --- Random Forest Regression with mlr3 using the Selected Features ---
    
    # Uses mlr3 and the mlr3learners package to perform Random Forest Regression (using the ranger learner) to predict weight regain percentage
    # Load required packages
    library(mlr3)
    library(mlr3learners)
    library(mlr3measures)
    library(ggplot2)
    
    # Set a seed for reproducibility
    #set.seed(123456)
    
    # Define the range of feature counts and number of repeats
    feature_counts <- seq(10, 100, by = 10)
    
    #results_df <- results_df %>% filter(Model != "Model 6: All")
    predictions_list <- list()
    roc_curve_data <- list()
    
    for(fc in feature_counts) {
        cat("Running for feature count:", fc, "\n")
        
        # Feature selection
        feature_selection <- mRMR.classic(data = mrmr_data, target_indices = target_index, feature_count = fc)
        selected_features <- solutions(feature_selection)[[1]]
        
        data_matrix_selected <- data_matrix %>%  select(Weight_regain6_18, selected_features[,1])
    
        taskClas <- TaskClassif$new(
          id = "rf_classification",
          backend = data_matrix_selected %>% 
            filter(!is.na(Weight_regain6_18)) %>% 
            mutate(Weight_regain6_18 = as.factor(Weight_regain6_18)),
          target = "Weight_regain6_18"
        )
        
        # Define RF model
        n_pred <- ncol(data_matrix_selected) - 1
        rfModel_class <- lrn("classif.ranger", 
                         importance = "permutation",
                         mtry = round(sqrt(n_pred)),
                         min.node.size = 5,
                         sample.fraction = 0.8,
                         num.threads = parallel::detectCores())
        rfModel_class$predict_type <- "prob"
        
        # 5-fold CV
        cv5 <- rsmp("cv", folds = 5)
        cv5$instantiate(taskClas)
        resample_result <- resample(taskClas, rfModel_class, cv5, store_models = TRUE)
    
        # Store predictions for plotting
        pred_dt <- as.data.table(resample_result$prediction())
        pred_dt$feature_count <- fc
        predictions_list[[as.character(fc)]] <- pred_dt
    
        # Extract ROC data using pROC
        try({
            roc_obj <- roc(pred_dt$truth, pred_dt$prob.1)
            roc_df <- rbind(roc_df,
                data.frame(
                  Group = grp,
                  FPR = 1 - roc_obj$specificities,
                  TPR = roc_obj$sensitivities,
                  AUC = round(roc_obj$auc, 3),
                  FeatureCount = as.factor(fc)
                ))
            roc_curve_data[[as.character(fc)]] <- roc_df
        })
    
        # Metrics
        accuracy    <- resample_result$aggregate(msr("classif.acc"))
        auc         <- resample_result$aggregate(msr("classif.auc"))
        sensitivity <- resample_result$aggregate(msr("classif.sensitivity"))
        specificity <- resample_result$aggregate(msr("classif.specificity"))
    
        cat("Accuracy:", round(accuracy, 3), "\n")
        cat("AUC:", round(auc, 3), "\n")
        cat("Sensitivity:", round(sensitivity, 3), "\n")
        cat("Specificity:", round(specificity, 3), "\n")
    
        results_df <- rbind(results_df,
                            data.frame(feature_count = fc,
                                       Group = grp,
                                       Model = "Model 5: All",
                                       Accuracy = accuracy,
                                       AUC = auc,
                                       Sensitivity = sensitivity,
                                       Specificity = specificity))
        }
    }
print(results_df)


################################################################################################
# Plot
################################################################################################
results_df <- results_df %>% mutate(Model_main = trimws(sub(":.*$", "", as.character(Model))))

model_label_text <- results_df %>%
    distinct(Model) %>%
    mutate(model_num = as.integer(sub(".*?(\\d+).*", "\\1", Model))) %>%
    arrange(model_num) %>%
    pull(Model) %>%
    paste(collapse = "\n")
# Add left-aligned text block inside the panel (top-left); tweak y if needed
#x_left <- if (is.factor(results_df$Model)) levels(results_df$Model)[1] else unique(results_df$Model)[1]

  
p <- ggplot(results_df, aes(x = Model_main, y = AUC, fill = Group)) +
    geom_boxplot(outlier.colour = "grey", alpha=0.9) +
    ylab("AUC") +
    ggtitle("Weight Regain (Month 6 vs 18)") +
    theme_minimal() +
    theme(
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 20, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 17, angle=-0, vjust = 0.5,hjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 20, face = "bold", color = "black"),
        axis.ticks.x = element_blank(),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),   # <- control legend label size
        legend.position = c(0.62, 0.28),          # <- position in normalized (x,y)
        legend.justification = c(1, 0)           # <- align the corner
    ) +
    scale_fill_manual(values = c("0" = "#7796A7", "1" = "#B9772B"),labels = c("0"="Placebo", "1"="aFMT")) +
    #scale_y_continuous(breaks = seq(0.5, 1, by = 0.1), limits = c(0.5, 1)) +
    #guides(fill = guide_legend(ncol = 2)) 
    annotate("text",
             x = 2.8, y = 0.55,              # right、above bottom
             label = model_label_text,
             hjust = 0, vjust = 1,
             size = 5, lineheight = 1.05)

p
#Output 8*7



# Output source data
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "AUC of RF models")
writeData(wb, "AUC of RF models", results_df)
saveWorkbook(wb, "Source Data Fig 4a.xlsx", overwrite = TRUE)



####################################################################################
# AUROC curves of the full integrated model 
####################################################################################

############# Do the calculation for Placebo group #############
roc_curve_data_fixed_group0 <- lapply(roc_curve_data, function(df) {
    df[df$Group == "0", , drop = FALSE] })

all_roc_df <- bind_rows(roc_curve_data_fixed_group0$`100`)
head(all_roc_df)
# Show row count for each FeatureCount
#all_roc_df %>% count(FeatureCount)
table(all_roc_df$FeatureCount)

# Check roc_curve_data to see which variable(s) don't have the same length as others.
#target_n <- max(all_roc_df %>% count(FeatureCount) %>% select(n))
#target_n
#roc_curve_data_fixed <- lapply(roc_curve_data_fixed_group0, function(df) {
#    if (nrow(df) != target_n) {
#         diff <- target_n - nrow(df)
#        df <- rbind(df, df[rep(nrow(df), diff), , drop = FALSE])
#    }
#    df
#})

# Extract FPR, TPR, and AUC columns for each element and bind row-wise
#mean_df <- do.call(cbind, lapply(all_roc_df, function(df) df[, c("Group","FPR", "TPR", "AUC")]))
library(dplyr)
library(tidyr)

mean_df <- all_roc_df %>%
  group_by(FeatureCount) %>%
  mutate(rowid = row_number()) %>%   # row index within each FeatureCount
  ungroup() %>%
  select(rowid, FeatureCount, Group, FPR, TPR, AUC) %>%
  mutate(FeatureCount = as.character(FeatureCount)) %>%  # for clean column names
  pivot_wider(
    id_cols     = rowid,
    names_from  = FeatureCount,
    values_from = c(Group, FPR, TPR, AUC),
    names_glue  = "{FeatureCount}.{.value}",
    names_vary  = "slowest",          # yields 10.Group 10.FPR 10.TPR 10.AUC, then 20.*, etc.
    values_fill = NA
  ) %>%
  select(-rowid)

mean_df

# Compute row-wise means for FPR, TPR, AUC across different feature counts
mean_df2 <- mean_df %>%
  mutate(
    Group = "0",
    FPR = rowMeans(select(., matches("\\.FPR$")), na.rm = TRUE),
    TPR = rowMeans(select(., matches("\\.TPR$")), na.rm = TRUE),
    AUC = rowMeans(select(., matches("\\.AUC$")), na.rm = TRUE)
  )  %>% select(Group, FPR, TPR, AUC) %>% 
  mutate( FeatureCount = "Mean")
  
all_roc_df_group0 <- rbind(all_roc_df, mean_df2)


############# Do the calculation for aFMT group #############
roc_curve_data_fixed_group1 <- lapply(roc_curve_data, function(df) {
    df[df$Group == "1", , drop = FALSE] })

all_roc_df <- bind_rows(roc_curve_data_fixed_group1$`100`)
head(all_roc_df)
# Show row count for each FeatureCount
#all_roc_df %>% count(FeatureCount)

# Check roc_curve_data to see which variable(s) don't have the same length as others.
#target_n <- max(all_roc_df %>% count(FeatureCount) %>% select(n))
#target_n
#roc_curve_data_fixed <- lapply(roc_curve_data_fixed_group0, function(df) {
#    if (nrow(df) != target_n) {
#         diff <- target_n - nrow(df)
#        df <- rbind(df, df[rep(nrow(df), diff), , drop = FALSE])
#    }
#    df
#})

# Extract FPR, TPR, and AUC columns for each element and bind row-wise
#mean_df <- do.call(cbind, lapply(all_roc_df, function(df) df[, c("Group","FPR", "TPR", "AUC")]))
library(dplyr)
library(tidyr)
library(matrixStats)


mean_df <- all_roc_df %>%
  group_by(FeatureCount) %>%
  mutate(rowid = row_number()) %>%   # row index within each FeatureCount
  ungroup() %>%
  select(rowid, FeatureCount, Group, FPR, TPR, AUC) %>%
  mutate(FeatureCount = as.character(FeatureCount)) %>%  # for clean column names
  pivot_wider(
    id_cols     = rowid,
    names_from  = FeatureCount,
    values_from = c(Group, FPR, TPR, AUC),
    names_glue  = "{FeatureCount}.{.value}",
    names_vary  = "slowest",          # yields 10.Group 10.FPR 10.TPR 10.AUC, then 20.*, etc.
    values_fill = NA
  ) %>%
  select(-rowid)

mean_df

# Compute row-wise means for FPR, TPR, AUC across different feature counts
mean_df2 <- mean_df %>%
  mutate(
    Group = "1",
    FPR = rowMeans(select(., matches("\\.FPR$")), na.rm = TRUE),
    TPR = rowMeans(select(., matches("\\.TPR$")), na.rm = TRUE),
    AUC = rowMeans(select(., matches("\\.AUC$")), na.rm = TRUE)
  )  %>% select(Group, FPR, TPR, AUC) %>% 
  mutate( FeatureCount = "Mean")
  
all_roc_df_group1 <- rbind(all_roc_df, mean_df2)

####### merge the 2 data ############

all_roc_df2 <- rbind(all_roc_df_group0, all_roc_df_group1)

####### Alternative: smoothed AUROCs ############
library(zoo)
smooth_roc <- all_roc_df2 %>%
  group_by(Group, FeatureCount) %>%
  arrange(FPR) %>%
  mutate(
    TPR_smooth = na.approx(TPR, x = FPR, na.rm = FALSE),
    FPR_smooth = na.approx(FPR, x = FPR, na.rm = FALSE)
  )
custom_colors <- c(
    "#aec7e8", "#1f77b4",  "#2ca02c", "#17becf",
    "#bcbd22", "#9467bd", "#8c564b", "#e377c2",  "#ff7f0e",  # Continued
    "#d62728", "Black")

smooth_roc_non_mean_group0 <- smooth_roc %>% filter(FeatureCount != "Mean" & Group == "0")
smooth_roc_non_mean_group1 <- smooth_roc %>% filter(FeatureCount != "Mean" & Group == "1")
smooth_roc_mean <- smooth_roc %>% filter(FeatureCount == "Mean")

results_df_model5 <- results_df %>% filter(Model == "Model 5: All") %>% filter(AUC != "NaN")
tapply(results_df_model5$AUC, results_df_model5$Group, mean, na.rm = TRUE)
#       0         1 
#0.9133690 0.9438272 



p <- ggplot() +
    geom_line(data = smooth_roc_non_mean_group0, aes(x = FPR_smooth, y = TPR_smooth, group = FeatureCount), color = "grey") +
    geom_line(data = smooth_roc_non_mean_group1, aes(x = FPR_smooth, y = TPR_smooth, group = FeatureCount), color = "grey") +
    geom_line(data = smooth_roc_mean, aes(x = FPR_smooth, y = TPR_smooth, color = Group), size = 2, alpha = 1) +
    annotate("text", x = 0.58, y = 0.53, label = "Mean AUC",
             hjust = 0, vjust = 0, size = 5, fontface = "bold", color = "grey20") +
    annotate("text", x = 0.58, y = 0.38, label = "aFMT = 0.94\nPlacebo =0.91",
             hjust = 0, vjust = 0, size = 5, color = "grey20") +
    labs(title = "AUROC Curves of Model 5",
         x = "FPR", 
         y = "TPR",
         size = "Mean AUC") +
    scale_color_manual(values = c("0" = "#7796A7", "1" = "#B9772B"),labels = c("0"="Placebo", "1"="aFMT")) +
    theme(
        axis.line = element_line(color = "black"),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 16, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        legend.position = c(0.87, 0.1),
        legend.justification = c(1, 0),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.1, "cm")
    ) +
    guides(size = guide_legend(ncol = 1))
p
# Output 5*5
#ggsave("../Figure4b. Weight_regain AUROCs smoothed mean aFMT vs Placebo.pdf", plot = p, width = 5, height = 5)


# Output source data
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "AUROC of Model 5")
writeData(wb, "AUROC of Model 5", smooth_roc)
saveWorkbook(wb, "Source Data Fig 4b.xlsx", overwrite = TRUE)



####################################################################################
# Extract importance of the top 30 features
####################################################################################

################################################################################################
# Model 5- ALL: RF based on Metadata + Microbiome + Fecal Metabolites + Strain data
################################################################################################
set.seed(285)
library(matrixStats)
library(mRMRe)

# Ensure all columns are numeric. If not, you might need to convert or remove non-numeric columns.
# Also, mRMRe expects a matrix with numeric values.
 
predict_data5 <- merge(predict_data3, strain_change_df, by = "SubjectID")
dim(predict_data5)
#[1]    76 2294

# Rename all columns in strain_results except "SubjectID" by appending the suffix "_strain_B"
strain_results_renamed <- strain_results %>%
  rename_with(~ paste0(.x, "_strain_B"), -SubjectID)
# Merge the data frames; now the strain_results columns already have the suffix.
predict_data6 <- merge(predict_data5, strain_results_renamed, by = "SubjectID")

dim(predict_data6)
#[1]    76 2298
rownames(predict_data6) <- predict_data6$SubjectID


predict_data7 <- predict_data6[, !names(predict_data6) %in% c("SubjectID", "sno_time_A", "Time_A", 
                                                               "sno_time_B", "Time_B", #"Gender_B", "Diet_B", "Treatment_B",
                                                               "sno_time_D", "Time_D", #"Gender_D", "Diet_D", "Treatment_D",
                                                               "sno_time_C", "Time_C", #"Gender_C", "Diet_C", "Treatment_C", 
                                                               "Weight_Baseline", "Weight_Change0_6",
                                                               "Weight_Change6_14", "Weight_regain6_14",
                                                               "Weight_Change6_18")]#, "Weight_regain6_18")]

predict_data7$Treatment_A <- ifelse(predict_data7$Treatment_A == "aFMT", 1,
                              ifelse(predict_data7$Treatment_A == "Placebo", 0, predict_data7$Treatment_A))
predict_data7$Treatment_B <- ifelse(predict_data7$Treatment_B == "aFMT", 1,
                              ifelse(predict_data7$Treatment_B == "Placebo", 0, predict_data7$Treatment_B))
predict_data7$Treatment_D <- ifelse(predict_data7$Treatment_D == "aFMT", 1,
                              ifelse(predict_data7$Treatment_D == "Placebo", 0, predict_data7$Treatment_D))                              
predict_data7$Treatment_C <- ifelse(predict_data7$Treatment_C == "aFMT", 1,
                              ifelse(predict_data7$Treatment_C == "Placebo", 0, predict_data7$Treatment_C))
predict_data7$Treatment_A <- as.numeric(predict_data7$Treatment_A)
predict_data7$Treatment_B <- as.numeric(predict_data7$Treatment_B)
predict_data7$Treatment_D <- as.numeric(predict_data7$Treatment_D)
predict_data7$Treatment_C <- as.numeric(predict_data7$Treatment_C)
predict_data7[is.na(predict_data7)] <- 0
dim(predict_data7)
#[1]    74 2284

predict_data <- predict_data7 %>% select(Weight_regain6_18, ends_with("_B"))
dim(predict_data) 
#[1] 74 1977

predictors <- predict_data %>% filter(Weight_regain6_18!= "NA")
dim(predictors)
#[1] 74 1977


roc_df <- data.frame(
              #Group = character(),
              FPR = numeric(),
              TPR = numeric(),
              AUC = numeric(),
              FeatureCount = character(),
            stringsAsFactors = FALSE)
            
# Initialize an empty data frame to store the performance metrics
results_df <- data.frame(feature_count = integer(),
                         #Group = character(),
                         Model = character(),
                         Accuracy = numeric(),
                         AUC = numeric(),
                         Sensitivity = numeric(),
                         Specificity = numeric(),
                         stringsAsFactors = FALSE)

# Create an mRMRe Data object. 
# IMPORTANT: mRMRe requires column names and numeric data.
data_matrix <- as.data.frame(lapply(predictors, as.numeric))
data_matrix$Weight_regain6_18 <- ifelse(data_matrix$Weight_regain6_18 > median(data_matrix$Weight_regain6_18), 1, 0)
#data_matrix$Weight_regain6_18 <- as.factor(data_matrix$Weight_regain6_18)

mrmr_data <- mRMR.data(data = data_matrix)

# Identify the index of the target variable column (e.g., column named "target")
target_index <- which(colnames(data_matrix) == "Weight_regain6_18")


########################################################
# --- Random Forest Regression with mlr3 using the Selected Features ---

# Uses mlr3 and the mlr3learners package to perform Random Forest Regression (using the ranger learner) to predict weight regain percentage
# Load required packages
library(mlr3)
library(mlr3learners)
library(mlr3measures)
library(ggplot2)

# Set a seed for reproducibility
#set.seed(123456)

# Define the range of feature counts and number of repeats
feature_counts <- seq(10, 100, by = 10)

results_df <- results_df %>% filter(Model != "Model 5: All")
predictions_list <- list()
roc_curve_data <- list()

for(fc in feature_counts) {
    cat("Running for feature count:", fc, "\n")
    
    # Feature selection
    feature_selection <- mRMR.classic(data = mrmr_data, target_indices = target_index, feature_count = fc)
    selected_features <- solutions(feature_selection)[[1]]
    
    data_matrix_selected <- data_matrix %>%  select(Weight_regain6_18, selected_features[,1])

    taskClas <- TaskClassif$new(
      id = "rf_classification",
      backend = data_matrix_selected %>% 
        filter(!is.na(Weight_regain6_18)) %>% 
        mutate(Weight_regain6_18 = as.factor(Weight_regain6_18)),
      target = "Weight_regain6_18"
    )
    
    # Define RF model
    n_pred <- ncol(data_matrix_selected) - 1
    rfModel_class <- lrn("classif.ranger", 
                     importance = "permutation",
                     mtry = round(sqrt(n_pred)),
                     min.node.size = 5,
                     sample.fraction = 0.8,
                     num.threads = parallel::detectCores())
    rfModel_class$predict_type <- "prob"
    
    # 5-fold CV
    cv5 <- rsmp("cv", folds = 5)
    cv5$instantiate(taskClas)
    resample_result <- resample(taskClas, rfModel_class, cv5, store_models = TRUE)

    # Store predictions for plotting
    pred_dt <- as.data.table(resample_result$prediction())
    pred_dt$feature_count <- fc
    predictions_list[[as.character(fc)]] <- pred_dt

    # Extract ROC data using pROC
    try({
        roc_obj <- roc(pred_dt$truth, pred_dt$prob.1)
        roc_df <- rbind(roc_df,
            data.frame(
              #Group = grp,
              FPR = 1 - roc_obj$specificities,
              TPR = roc_obj$sensitivities,
              AUC = round(roc_obj$auc, 3),
              FeatureCount = as.factor(fc)
            ))
        roc_curve_data[[as.character(fc)]] <- roc_df
    })

    # Metrics
    accuracy    <- resample_result$aggregate(msr("classif.acc"))
    auc         <- resample_result$aggregate(msr("classif.auc"))
    sensitivity <- resample_result$aggregate(msr("classif.sensitivity"))
    specificity <- resample_result$aggregate(msr("classif.specificity"))

    cat("Accuracy:", round(accuracy, 3), "\n")
    cat("AUC:", round(auc, 3), "\n")
    cat("Sensitivity:", round(sensitivity, 3), "\n")
    cat("Specificity:", round(specificity, 3), "\n")

    results_df <- rbind(results_df,
                        data.frame(feature_count = fc,
                                   #Group = grp,
                                   Model = "Model 5: All",
                                   Accuracy = accuracy,
                                   AUC = auc,
                                   Sensitivity = sensitivity,
                                   Specificity = specificity))
        }
print(results_df)

########### Continued as above ############
# Ensure that your resampling stored models:
#resample_result <- resample(taskClas, rfModel_class, cv5, store_models = TRUE)

# Extract importance from each fold:
importance_list <- map(resample_result$learners, function(lrn_obj) {
    # Extract the importance vector from each learner
    lrn_obj$importance()
})

# Combine all importance vectors into a single dataframe
importance_df <- bind_rows(lapply(seq_along(importance_list), function(i) {
  tibble(
    feature = names(importance_list[[i]]),
    importance = unname(importance_list[[i]]),
    run = paste0("Run_", i)
  )
}))

# Pivot to wide format with features as rows and runs as columns
importance_wide <- importance_df %>%
  pivot_wider(names_from = run, values_from = importance) %>%
  column_to_rownames("feature")

# Compute summary statistics
num_cols <- sapply(importance_wide, is.numeric)
importance_wide$mean_importance <- rowMedians(
  as.matrix(importance_wide[, num_cols, drop = FALSE]), na.rm = TRUE)
#importance_wide$mean_importance <- rowMeans(importance_wide, na.rm = TRUE)
importance_wide$sd_importance <- apply(importance_wide, 1, sd, na.rm = TRUE)

# Sort by mean importance
importance_wide <- importance_wide[order(-importance_wide$mean_importance), ]

# View top features
head(importance_wide, 10)

importance_df  <- importance_wide %>%
    select(mean_importance) %>%
    mutate(Importance = mean_importance,
           Feature = rownames(importance_wide))%>%
    select(Feature, Importance)

# View the top features ranked by importance
head(importance_df, 20)


# Set all column names of all datasets as a dataframe
biomarker_names_df <- data.frame(
  Feature = colnames(base_biomarker_data_reformat),
  Dataset = "Biomarkers",
  row.names = NULL)
# foods_names_df <- data.frame(
#   Feature = colnames(nuts_merged_data_A_B_C)[!colnames(nuts_merged_data_A_B_C) %in% "SubjectID"],  
#   Dataset = "Foods",
#   row.names = NULL)
species_names_df <- data.frame(
  Feature = colnames(microbiome_change_BD_renamed)[!colnames(microbiome_change_BD_renamed) %in% "SubjectID"],
  Dataset = "Species",
  row.names = NULL)
fecalmetabolites_names_df <- data.frame(
  Feature = colnames(fecal_metabolites_change_BD_renamed)[!colnames(fecal_metabolites_change_BD_renamed) %in% "SubjectID"],
  Dataset = "Fecal Metabolites",
  row.names = NULL)
fecalmetabolites_names_df2 <- fecalmetabolites_names_df %>% mutate(Feature = paste0("X", Feature))

strain_names_df <- data.frame(
  Feature = c(colnames(strain_change_df)[-1], paste0(colnames(strain_results)[-1], "_strain_B")), 
  Dataset = "Strain",
  row.names = NULL)

datasets_names_df <- rbind(biomarker_names_df, species_names_df, fecalmetabolites_names_df2, strain_names_df)

library(dplyr)
importance_df2 <- importance_df %>%
  left_join(datasets_names_df %>% select(Feature, Dataset), by = "Feature")
table(importance_df2[1:50,]$Dataset)

head(importance_df2, 20)


#biomarker, foods, host genetics, species/strain, fecal metabolites 
########################################################################################
# Rename Species/Strain name to SGB
#setwd("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/8.phenotype_prediction")
library(dplyr) 
library(stringr)
taxonomy <- read.table("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/taxonomy_all_SGB.txt", header=T, sep="\t", check.names = FALSE)
# Subset the importance_df2 for Species
species_imp <- importance_df2 %>% 
  filter(Dataset == "Species") %>%
  # Remove the suffix "_change_BD_B" if present, and create a new column "Feature_clean"
  mutate(Feature_clean = str_remove(Feature, "_change_BD_B$"))
# Now join with taxonomy using the cleaned feature names matching taxonomy$t__SGB
new_species_dataset <- species_imp %>% left_join(taxonomy[-1], by = c("Feature_clean" = "t__SGB"))

new_species_dataset <- new_species_dataset %>%
  mutate(Species_SGB = paste0(Species, " (", SGB, ")"))
# Remove the leading "s__" from the values in the Species_SGB column using the gsub()
new_species_dataset$Species_SGB <- gsub("^s__", "", new_species_dataset$Species_SGB)
new_species_dataset2 <- new_species_dataset[, c(1,2,3,12)]
names(new_species_dataset2)[4] <- "New_name"


# New dataframe for strains
strain_imp <- importance_df2 %>% 
  filter(Dataset == "Strain") %>%
  # Remove the suffix "_change_BD_B" if present, and create a new column "Feature_clean"
  mutate(Feature_clean = str_remove(Feature, "_strain_change_BD_B$"))
# Now join with taxonomy using the cleaned feature names matching taxonomy$t__SGB
new_strain_dataset <- strain_imp %>% left_join(taxonomy[-1], by = c("Feature_clean" = "t__SGB"))

new_strain_dataset <- new_strain_dataset %>%
  mutate(Species_SGB = paste0(Species, " (", SGB, ")"))
# Remove the leading "s__" from the values in the Species_SGB column using the gsub()
new_strain_dataset$Species_SGB <- gsub("^s__", "", new_strain_dataset$Species_SGB)
new_strain_dataset2 <- new_strain_dataset[, c(1,2,3,12)]
names(new_strain_dataset2)[4] <- "New_name"


metrics <- c("Gain_B_D_percentage", "Loss_B_D_percentage", "Swap_B_D_percentage", "Succession_rate")
new_strain_dataset2 <- new_strain_dataset2 %>%
  mutate(New_name = case_when(
    grepl("Gain_B_D_percentage", Feature) ~ "Gain_rate_6_14 (%)",
    grepl("Loss_B_D_percentage", Feature) ~ "Loss_rate_6_14 (%)",
    grepl("Swap_B_D_percentage", Feature) ~ "Swap_rate_6_14 (%)",
    grepl("Succession_rate", Feature) ~ "Succession rate",
    TRUE ~ New_name
  ))


########################################################################################
# assign fecal chemical real names
library(readr)
fecal_chem_annotation <- read_tsv("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/Fecal_chem_annotation.tsv",
  na = c("", "NA", "N/A"), locale = locale(encoding = "UTF-8"), guess_max = 10000)
fecal_chem_annotation <- as.data.frame(fecal_chem_annotation)
# Subset the importance_df2 for Metabolites
Fecal_Metabolites_imp <- importance_df2 %>% 
  filter(Dataset == "Fecal Metabolites") %>%
  # Remove a leading "X" if present, and remove the suffix "_change_BD_B"
  mutate(Feature_clean = Feature %>% 
           str_remove("^X") %>% 
           str_remove("_change_BD_B$"))

# Now join with taxonomy using the cleaned feature names matching taxonomy$t__SGB
new_fecal_metabolite_dataset <- Fecal_Metabolites_imp %>% 
  left_join(fecal_chem_annotation %>% mutate(CHEM_ID = as.character(CHEM_ID)), 
            by = c("Feature_clean" = "CHEM_ID"))

new_fecal_metabolite_dataset2 <- new_fecal_metabolite_dataset[, c(1,2,3,7)]
names(new_fecal_metabolite_dataset2)[4] <- "New_name"


########################################################################################
# Assign Foods names
foods <- read.table("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/food_names.txt", header=T, sep="\t", check.names = FALSE)
foods_imp <- importance_df2 %>% 
  filter(Dataset == "Foods") %>%
  # Remove the suffix "_B" if present, and create a new column "Feature_clean"
  mutate(Feature_clean = str_remove(Feature, "_B$"))
new_foods_dataset <- foods_imp %>% left_join(foods, by = c("Feature_clean" = "nutNum"))
new_foods_dataset2 <- new_foods_dataset[, c(1,2,3,5)]
names(new_foods_dataset2)[4] <- "New_name"

########################################################################################
#biomarker, 
biomarker_imp <- importance_df2 %>% filter(Dataset == "Biomarkers")%>%
  # Remove the suffix "_B" if present, and create a new column "Feature_clean"
  mutate(Feature_clean = str_remove(Feature, "_B$"))
names(biomarker_imp)[4] <- "New_name"


########################################################################################
top_importance_df <- rbind(biomarker_imp, new_foods_dataset2, 
        new_species_dataset2, new_strain_dataset2, new_fecal_metabolite_dataset2)
# Reorder the data frame by the Importance column.
top_importance_df_ordered <- top_importance_df %>% arrange(desc(Importance))

head(top_importance_df_ordered, 20)
#                  Feature   Importance           Dataset                        New_name
#1                  PC22_B 0.0018890290     Host genetics              Host_genetics_PC22
#2 t__SGB15089_change_BD_B 0.0010803907           Species Vescimonas_coprocola (SGB15089)
#3                Leptin_B 0.0009073100        Biomarkers                          Leptin
#4  X100001313_change_BD_B 0.0007829267 Fecal Metabolites        gamma-glutamylmethionine
#5        X535_change_BD_B 0.0007725938 Fecal Metabolites                         uridine
#6  X100000774_change_BD_B 0.0007148969 Fecal Metabolites             phenyllactate (PLA)

########################################################################################
top_importance_df_ordered <- top_importance_df_ordered %>%
  mutate(Dataset = if_else(Dataset == "Fecal Metabolites", "Metabolites", Dataset))


library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "Features importance of Model 5")
writeData(wb, "Features importance of Model 5", top_importance_df_ordered)

saveWorkbook(wb, "Supplementary Table 7.xlsx", overwrite = TRUE)



#top_importance_df_ordered[39,4] <- "Anaerobutyricum_soehngenii strain (SGB4537)"
library(dplyr)
plot_data <- top_importance_df_ordered[1:30, ]

# if New_name is a factor, convert it to character
if (is.factor(plot_data$New_name)) {
  plot_data$New_name <- as.character(plot_data$New_name)
}
plot_data <- plot_data %>%
  mutate(
    .dup_group  = New_name %in% New_name[duplicated(New_name)],  # 该名称是否有重复
    .is_first   = !duplicated(New_name),                         # 是否为该名称的首个出现
    New_name    = ifelse(.dup_group & .is_first, paste0(New_name, "*"), New_name)
  ) %>%
  select(-.dup_group, -.is_first)

library(dplyr)
library(stringr)
plot_data <- plot_data %>%
  mutate(New_name = str_replace_all(New_name, c("alpha" = "α", "beta" = "β", "gamma" = "γ")) )

library(showtext)
library(sysfonts)

# 2) Add a Unicode font (Google Noto Sans works great)
#font_add_google("Noto Sans", "Noto Sans")
#showtext_auto()  # use text shaping for plots

p <- ggplot(plot_data, aes(x = Importance, y = reorder(New_name, Importance), fill = Dataset)) +
    geom_bar(stat = "identity", alpha=0.9) +
    #coord_flip() +
    labs(title = "Top 30 Features of Model 5",
       x = "Importance",
       y = "Feature") +
    #scale_x_continuous(limits = c(0, 0.002)) +
    #theme_minimal() +
    scale_fill_manual(values = c(
        "Biomarkers" = "#e7298a", # Magenta pink
        "Foods" = "#1b9e77",  # Teal green
        "Host genetics" = "#7570b3",  # Dusty purple
        "Species" = "#66a61e",  # Olive green
        "Metabolites" = "#d95f02", # Burnt orange
        "Strain" = "#e6ab02"  # Golden yellow
    )) +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 14, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 12, hjust = 0.9, color = "black", family="sans"),
        axis.text.y = element_text(size = 12, color = "black", family="sans"),
        #axis.ticks.x = element_blank()
        legend.title = element_text( face = "bold", size = 14),
        legend.text = element_text(size = 14),   # <- control legend label size
        legend.position = c(0.97, 0.02),          # <- position in normalized (x,y)
        legend.justification = c(1, 0),           # <- align the corner
        #legend.position = "bottom"
    ) + theme(text = element_text(family = "DejaVu Sans"))
p
# Output portrait: 10*9

#ggsave("../Figure4d. Weight_regain Top 30 Features aFMT vs Placebo.pdf", plot = p, width = 7, height = 7)


# Output source data
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "Importance of top 30 features")
writeData(wb, "Importance of top 30 features", plot_data)
saveWorkbook(wb, "Source Data Fig 4d.xlsx", overwrite = TRUE)


####################################################
#  Reformat the names of the top 30 features
####################################################
library(dplyr)
library(stringr)
library(ggplot2)

#--------------------------------------------------
# Function to generate plotmath expressions for y-axis labels
#--------------------------------------------------
make_species_expr <- function(x, dataset) {
  
  if (is.na(x)) return('""')
  
  # Global replacement for the swap-rate label
  if (x == "Swap_rate_6_14 (%)") {
    return("'Strain swap rate (%) M6–M14'")
  }
  
  # For non-Species entries: keep as plain text, with underscores replaced by spaces
  if (is.na(dataset) || dataset != "Species") {
    x_clean <- gsub("_", " ", x)
    x_clean <- gsub("'", "\\\\'", x_clean)
    x_clean = sub("^([a-z])", "\\U\\1", x_clean, perl = TRUE)
    return(paste0("'", x_clean, "'"))
  }
  
  # Extract the main name and the trailing bracket part
  main_part <- str_trim(str_extract(x, "^[^(]+"))
  bracket_part <- str_extract(x, "\\([^()]+\\)$")
  
  # If there is no underscore, it is unlikely to be a formatted species label
  if (!str_detect(main_part, "_")) {
    x_clean <- gsub("_", " ", x)
    x_clean <- gsub("'", "\\\\'", x_clean)
    return(paste0("'", x_clean, "'"))
  }
  
  tokens <- unlist(str_split(main_part, "_"))
  
  # Remove internal SGB-like tokens if there is already a final bracketed SGB
  # Example: GGB9627_SGB15081 (SGB15081) -> GGB9627 (SGB15081)
  if (!is.na(bracket_part)) {
    bracket_id <- str_extract(bracket_part, "SGB\\d+")
    if (!is.na(bracket_id)) {
      tokens <- tokens[tokens != bracket_id]
    }
  }
  
  # Rebuild main part after cleanup
  main_part_cleaned <- paste(tokens, collapse = "_")
  
  # If fewer than 2 tokens remain, treat as plain text
  if (length(tokens) < 2) {
    plain_label <- gsub("_", " ", main_part_cleaned)
    if (!is.na(bracket_part)) {
      plain_label <- paste0(plain_label, " ", bracket_part)
    }
    plain_label <- gsub("'", "\\\\'", plain_label)
    return(paste0("'", plain_label, "'"))
  }
  
  genus <- tokens[1]
  second <- tokens[2]
  
  # Case 1: "Genus_sp_xxx" style -> italicize genus only, render "sp." in plain text
  if (str_detect(genus, "^[A-Z][A-Za-z-]+$") && second == "sp") {
    
    remainder <- tokens[-c(1, 2)]
    remainder_text <- if (length(remainder) > 0) paste(remainder, collapse = " ") else ""
    remainder_text <- gsub("'", "\\\\'", remainder_text)
    
    if (!is.na(bracket_part) && remainder_text != "") {
      bracket_clean <- gsub("^\\(|\\)$", "", bracket_part)
      bracket_clean <- gsub("'", "\\\\'", bracket_clean)
      return(
        paste0(
          "italic('", genus, "')~'sp.'~'",
          remainder_text, "'~'(", bracket_clean, ")'"
        )
      )
    } else if (!is.na(bracket_part) && remainder_text == "") {
      bracket_clean <- gsub("^\\(|\\)$", "", bracket_part)
      bracket_clean <- gsub("'", "\\\\'", bracket_clean)
      return(
        paste0(
          "italic('", genus, "')~'sp.'~'(", bracket_clean, ")'"
        )
      )
    } else if (is.na(bracket_part) && remainder_text != "") {
      return(
        paste0(
          "italic('", genus, "')~'sp.'~'",
          remainder_text, "'"
        )
      )
    } else {
      return(
        paste0(
          "italic('", genus, "')~'sp.'"
        )
      )
    }
  }
  
  # Case 2: full Latin-style name -> italicize the full main part
  genus_ok <- str_detect(genus, "^[A-Z][A-Za-z-]+$")
  second_ok <- !str_detect(second, "^[0-9]+$") &&
    !str_detect(second, "^(SGB|GGB|GB|MGYG)\\d+$") &&
    str_detect(second, "^[A-Za-z][A-Za-z0-9.-]*$")
  
  if (genus_ok && second_ok) {
    main_clean <- gsub("_", " ", main_part_cleaned)
    main_clean <- gsub("'", "\\\\'", main_clean)
    
    if (!is.na(bracket_part)) {
      bracket_clean <- gsub("^\\(|\\)$", "", bracket_part)
      bracket_clean <- gsub("'", "\\\\'", bracket_clean)
      return(paste0("italic('", main_clean, "')~'(", bracket_clean, ")'"))
    } else {
      return(paste0("italic('", main_clean, "')"))
    }
  }
  
  # Case 3: not a real Latin-style species name -> plain text
  plain_label <- gsub("_", " ", main_part_cleaned)
  if (!is.na(bracket_part)) {
    plain_label <- paste0(plain_label, " ", bracket_part)
  }
  plain_label <- gsub("'", "\\\\'", plain_label)
  return(paste0("'", plain_label, "'"))
}

#--------------------------------------------------
# Prepare plotting data
#--------------------------------------------------
plot_data2 <- plot_data %>%
  mutate(
    New_name_expr = mapply(make_species_expr, New_name, Dataset)
  ) %>%
  arrange(Importance) %>%
  mutate(
    New_name_expr = factor(New_name_expr, levels = New_name_expr)
  )

#--------------------------------------------------
# Plot
#--------------------------------------------------
library(showtext)
library(sysfonts)
p <- ggplot(plot_data2, aes(x = Importance, y = New_name_expr, fill = Dataset)) +
  geom_col(alpha = 1) +
  labs(
    title = "Top 30 Features of Model 5",
    x = "Importance",
    y = NULL
  ) +
  scale_y_discrete(labels = function(x) parse(text = x)) +
  scale_x_continuous(
    labels = function(x) ifelse(x == 0, "0", x)
  ) +
  scale_fill_manual(values = c(
    "Biomarkers" = "#e7298a",
    "Foods" = "#1b9e77",
    "Host genetics" = "#7570b3",
    "Species" = "#54A24B",
    "Metabolites" = "#d95f02",
    "Strain" = "#e6ab02"
  )) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12, color = "black", family = "sans",hjust = 0.9),
    axis.text.y = element_text(size = 12, color = "black", family = "sans", hjust = 1),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14),
    legend.position = c(0.97, 0.02),
    legend.justification = c(1, 0),
    text = element_text(family = "DejaVu Sans"),
    plot.margin = margin(5.5, 10, 5.5, 5.5)
  ) +
  coord_cartesian(clip = "off")

p
#ggsave("../Figure4d. Weight_regain Top 30 Features aFMT vs Placebo.pdf", plot = p, width = 6, height = 7)




#######################################################################################
# Sum importance of ALL features
#######################################################################################

importance_sum <- top_importance_df_ordered %>%
  group_by(Dataset) %>%
  summarise(total_importance = sum(Importance, na.rm = TRUE))

print(importance_sum)
# A tibble: 5 × 2
#  Dataset           total_importance
#  <chr>                        <dbl>
#1 Biomarkers                 0.00474
#2 Fecal Metabolites          0.0419 
#3 Foods                      0.00376
#4 Species                    0.0214 
#5 Strain                     0.00931

#######################################################################################
# Plot
# Rank the datasets by importance
importance_sum_ordered <- importance_sum %>%
  arrange(total_importance)
importance_sum_ordered <- importance_sum_ordered %>%
  mutate(Dataset = if_else(Dataset == "Fecal Metabolites", "Metabolites", Dataset))

p <- ggplot(importance_sum_ordered, aes(x = total_importance, y = reorder(Dataset, total_importance), fill = Dataset)) +
    geom_bar(stat = "identity", alpha=1) +
    #coord_flip() +
    labs(title = "Total Importance of Model 5",
       x = "Importance",
       y = "Dataset") +
    #scale_x_continuous(limits = c(0, 0.005)) +
    theme_minimal() +
    scale_fill_manual(values = c(
        "Biomarkers" = "#e7298a", # Magenta pink
        "Foods" = "#1b9e77",  # Teal green
        "Host genetics" = "#7570b3",  # Dusty purple
        "Species" = "#54A24B",  # Olive green
        "Metabolites" = "#d95f02", # Burnt orange
        "Strain" = "#e6ab02"  # Golden yellow
    )) +
    theme(
        #panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 14, face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, face = "bold", color = "black"),
        #axis.ticks.x = element_blank()
        #legend.title = element_text( face = "bold"),
        legend.position = "none"
    ) 
p
# Output: 5*4
#ggsave("../Figure4c. Weight_regain Total Importance aFMT vs Placebo.pdf", plot = p, width = 4, height = 5)


# Output source data
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "Importance sum of each dataset")
writeData(wb, "Importance sum of each dataset", importance_sum_ordered)
saveWorkbook(wb, "Source Data Fig 4c.xlsx", overwrite = TRUE)






#########################################################################
# Correlation between top features and weight regain M6–M18
#########################################################################

head(importance_df)
#                                        Feature  Importance
#t__SGB4953_change_BD_B   t__SGB4953_change_BD_B 0.006220616
#t__SGB15322_change_BD_B t__SGB15322_change_BD_B 0.004044112
#LDLc_B                                   LDLc_B 0.003551747
#X100004275_change_BD_B   X100004275_change_BD_B 0.003482872
#t__SGB4644_change_BD_B   t__SGB4644_change_BD_B 0.003382953
#Insulin_B                             Insulin_B 0.002674002



############################ t__SGB4953_change_BD_B ############################
top_feature <- "t__SGB4953_change_BD_B"
# Extract abundance data at Month 6 and 14
t__SGB4953_B <- microbiome_data_B[,c("SubjectID","Month","SampleID","Diet","Treatment","t__SGB4953")]
t__SGB4953_D <- microbiome_data_D[,c("SubjectID","Month","SampleID","Diet","Treatment","t__SGB4953")]
head(t__SGB4953_B)
#  SubjectID Month SampleID Diet Treatment t__SGB4953
#1         4     6     B004    3      aFMT    0.00000
#2        13     6     B013    2   Placebo    0.00000
#3        15     6     B015    3   Placebo    0.18154
#4        17     6     B017    3   Placebo    0.00000
#5        19     6     B019    2   Placebo    0.00490
#6        22     6     B022    3   Placebo    0.15649
names(t__SGB4953_B)[6] <- "Abundance_B"
names(t__SGB4953_D)[6] <- "Abundance_D"

predictors_top_df <- predict_data2[, c("SubjectID","Diet_A","Treatment_A","Weight_regain6_18", top_feature)]
head(predictors_top_df)
names(predictors_top_df) <- c("SubjectID","Diet","Treatment","Weight_regain6_18", "AbundanceChange")

# Use an inner_join() on SubjectID so that only subjects present in all three tables are kept, 
# then append Abundance_B and Abundance_D to predictors_top_df.
library(dplyr)
predictors_top_df_merged <- predictors_top_df %>%
  inner_join(
    t__SGB4953_B %>%
      dplyr::select(SubjectID, Abundance_B) %>%
      distinct(SubjectID, .keep_all = TRUE),
    by = "SubjectID"
  ) %>%
  inner_join(
    t__SGB4953_D %>%
      dplyr::select(SubjectID, Abundance_D) %>%
      distinct(SubjectID, .keep_all = TRUE),
    by = "SubjectID"
  )

head(predictors_top_df_merged)
#  SubjectID Diet Treatment Weight_regain6_18 AbundanceChange Abundance_B Abundance_D
#1         4    3      aFMT          41.77215         0.00000     0.00000     0.00000
#2        13    2   Placebo          39.58333         0.00000     0.00000     0.00000
#3        15    3   Placebo          61.90476         0.14980     0.18154     0.33134
#4        19    2   Placebo         -31.95876         0.51592     0.00490     0.02554
#5        22    3   Placebo          50.63291        -0.15649     0.15649     0.52082
#6        27    2      aFMT          25.97403         0.00000     0.00000     0.00000

# remove rows where Abundance_B and Abundance_D equal to 0, which will also result in the zero abundance change
predictors_top_df_merged_nonzero <- predictors_top_df_merged %>%
  dplyr::filter(!(Abundance_B == 0 & Abundance_D == 0))
  
# Compute the median
#median_wr <- median(predictors_top_df$Weight_regain6_18, na.rm = TRUE)
# Add label column
#predictors_top_df$Regain_group <- ifelse(predictors_top_df$Weight_regain6_18 > median_wr, "High regainer", "Low regainer")
#predictors_top_df_merged_nonzero$Treatment <- ifelse(predictors_top_df_merged_nonzero$Treatment == "1", "aFMT", "Placebo")

# Remove points with very extreme high or low values
#predictors_top_df_filtered <- predictors_top_df_merged_nonzero
predictors_top_df_filtered <- predictors_top_df_merged_nonzero %>% filter(Weight_regain6_18 <= 200 & abs(AbundanceChange) <= 0.5)

# Calculate spearman correlation
cor_test_result <- cor.test(
  predictors_top_df_filtered$AbundanceChange,
  predictors_top_df_filtered$Weight_regain6_18,
  method = "spearman"
)
cor_test_result$estimate  # spearman correlation coefficient (r)
#     rho 
#0.5477941 
cor_test_result$p.value
#[1] 0.001156525



# compare the correlation between AbundanceChange and Weight_regain6_18 across two treatment groups (e.g., aFMT vs Placebo)
df_afmt <- predictors_top_df_filtered[predictors_top_df_filtered$Treatment == "aFMT", ]
df_placebo <- predictors_top_df_filtered[predictors_top_df_filtered$Treatment == "Placebo", ]
cor_afmt <- cor.test(df_afmt$AbundanceChange, df_afmt$Weight_regain6_18, method = "spearman")
cor_placebo <- cor.test(df_placebo$AbundanceChange, df_placebo$Weight_regain6_18, method = "spearman")
cor_afmt
cor_placebo

cor_afmt

#    Spearman's rank correlation rho

#data:  df_afmt$AbundanceChange and df_afmt$Weight_regain6_18
#S = 3583.3, p-value = 0.05444
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#     rho 
#0.343245 

cor_placebo

#    Spearman's rank correlation rho

#data:  df_placebo$AbundanceChange and df_placebo$Weight_regain6_18
#S = 7472.1, p-value = 0.06085
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.2990547 

library(MASS)

p <- ggplot(predictors_top_df_filtered, aes(x = AbundanceChange, y = Weight_regain6_18, colour = Treatment)) +
    geom_point(size = 2, alpha = 1) +
    #geom_smooth(method = "rlm", se = TRUE, linetype = "dashed", size = 1) +  # fit for aFMT and Placebo respectively
    geom_smooth(method = "rlm", se = TRUE, linetype = "dashed", size = 1, colour = "grey", aes(group = 1)) +  # fit for all
    #geom_hline(yintercept = median_wr, color = "grey50", linetype = "dotted") +
    labs(
        title = "Lachnospiraceae bacterium (SGB4953)",
        subtitle = paste0(
                    "Spearman's r = ", signif(cor_test_result$estimate, 2), ", p = ", signif(cor_test_result$p.value, 2)
                    #"aFMT: Spearman's r = ", signif(cor_afmt$estimate, 2), ", p = ", signif(cor_afmt$p.value, 1),
                    #"\nPlacebo: Spearman's r = ", signif(cor_placebo$estimate, 2), ", p = ", signif(cor_placebo$p.value, 1)
                    ),
        x = "Abundance Change % (6–14)", 
        y = "Weight regain % (6–18)"
    ) +
    scale_color_manual(values = c("Placebo" = "#7796A7", "aFMT" = "#B9772B")) +
    #scale_color_manual(values = c("High regainer" = "#006241", "Low regainer" = "#66CDAA")) +
    #annotate("text", x = -0.3, y = median_wr + 50, label = "High regainers", color = "#006241", fontface = "bold") +
    #annotate("text", x = -0.3, y = median_wr - 50, label = "Low regainers", color = "#66CDAA", fontface = "bold") +
    #theme_minimal() +
    theme(
        panel.background = element_blank(),
        #axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "grey", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0, size = 14),
        axis.line = element_line(color = "black"),
        axis.title.x = element_text(size = 14, face = "bold", family="sans"),
        axis.title.y = element_text(size = 14, face = "bold", family="sans"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        #axis.ticks.x = element_blank()
        #legend.title = element_text( face = "bold"),
        legend.position = "right"
    ) 
p
# Save the plot to a PDF file
# ggsave("../Figure4e. Top1 Lachnospiraceae bacterium (SGB4953) abundance_vs_weight_regain nonzero.pdf", plot = p, width = 6, height = 5)

# Output 6*5

# Output source data
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "Lachnospiraceae bacterium")
writeData(wb, "Lachnospiraceae bacterium", predictors_top_df_filtered)
#saveWorkbook(wb, "Source Data Fig 4e.xlsx", overwrite = TRUE)



################################# Top 2 ########################################
############################ t__SGB15322_change_BD_B ############################
top_feature <- "t__SGB15322_change_BD_B"
# Extract abundance data at Month 6 and 14
t__SGB15322_B <- microbiome_data_B[,c("SubjectID","Month","SampleID","Diet","Treatment","t__SGB15322")]
t__SGB15322_D <- microbiome_data_D[,c("SubjectID","Month","SampleID","Diet","Treatment","t__SGB15322")]
head(t__SGB15322_B)
#  SubjectID Month SampleID Diet Treatment t__SGB15322
#1         4     6     B004    3      aFMT     0.00000
#2        13     6     B013    2   Placebo     0.00000
#3        15     6     B015    3   Placebo     0.00000
#4        17     6     B017    3   Placebo     0.00000
#5        19     6     B019    2   Placebo     0.02556
#6        22     6     B022    3   Placebo     1.51759
names(t__SGB15322_B)[6] <- "Abundance_B"
names(t__SGB15322_D)[6] <- "Abundance_D"

predictors_top_df <- predict_data2[, c("SubjectID","Diet_A","Treatment_A","Weight_regain6_18", top_feature)]
head(predictors_top_df)
names(predictors_top_df) <- c("SubjectID","Diet","Treatment","Weight_regain6_18", "AbundanceChange")

# Use an inner_join() on SubjectID so that only subjects present in all three tables are kept, 
# then append Abundance_B and Abundance_D to predictors_top_df.
library(dplyr)
predictors_top_df_merged <- predictors_top_df %>%
  inner_join(
    t__SGB15322_B %>%
      dplyr::select(SubjectID, Abundance_B) %>%
      distinct(SubjectID, .keep_all = TRUE),
    by = "SubjectID"
  ) %>%
  inner_join(
    t__SGB15322_D %>%
      dplyr::select(SubjectID, Abundance_D) %>%
      distinct(SubjectID, .keep_all = TRUE),
    by = "SubjectID"
  )

head(predictors_top_df_merged)
#  SubjectID Diet Treatment Weight_regain6_18 AbundanceChange Abundance_B Abundance_D
#1         4    3      aFMT          41.77215         1.45378     0.00000     1.45378
#2        13    2   Placebo          39.58333         0.00000     0.00000     0.00000
#3        15    3   Placebo          61.90476         0.00000     0.00000     0.00000
#4        19    2   Placebo         -31.95876         3.31095     0.02556     0.07642
#5        22    3   Placebo          50.63291        -1.51759     1.51759     3.33651
#6        27    2      aFMT          25.97403         0.00000     0.00000     0.00000

# remove rows where Abundance_B and Abundance_D equal to 0, which will also result in the zero abundance change
predictors_top_df_merged_nonzero <- predictors_top_df_merged %>%
  dplyr::filter(!(Abundance_B == 0 & Abundance_D == 0))
  
  
# Compute the median
#median_wr <- median(predictors_top_df$Weight_regain6_18, na.rm = TRUE)
# Add label column
#predictors_top_df$Regain_group <- ifelse(predictors_top_df$Weight_regain6_18 > median_wr, "High regainer", "Low regainer")
#predictors_top_df$Treatment <- ifelse(predictors_top_df$Treatment == "1", "aFMT", "Placebo")

# Remove points with very extreme high or low values
#predictors_top_df_filtered <- predictors_top_df_merged_nonzero
#predictors_top_df_filtered <- predictors_top_df_merged_nonzero  %>% filter(Weight_regain6_18 <= 200 & AbundanceChange <= 5)
predictors_top_df_filtered <- predictors_top_df_merged_nonzero  %>% filter(Weight_regain6_18 <= 200)
# Calculate Pearson correlation
cor_test_result <- cor.test(
  predictors_top_df_filtered$AbundanceChange,
  predictors_top_df_filtered$Weight_regain6_18,
  method = "spearman"
)
cor_test_result$estimate  # spearman correlation coefficient (r)
#       rho 
# -0.4220042
cor_test_result$p.value
#[1] 0.0356155

# compare the correlation between AbundanceChange and Weight_regain6_18 across two treatment groups (e.g., aFMT vs Placebo)
df_afmt <- predictors_top_df_filtered[predictors_top_df_filtered$Treatment == "aFMT", ]
df_placebo <- predictors_top_df_filtered[predictors_top_df_filtered$Treatment == "Placebo", ]
cor_afmt <- cor.test(df_afmt$AbundanceChange, df_afmt$Weight_regain6_18, method = "spearman")
cor_placebo <- cor.test(df_placebo$AbundanceChange, df_placebo$Weight_regain6_18, method = "spearman")
cor_afmt
cor_placebo

library(MASS)
p <- ggplot(predictors_top_df_filtered, aes(x = AbundanceChange, y = Weight_regain6_18, colour = Treatment)) +
    geom_point(size = 2, alpha = 1) +
    #geom_smooth(method = "rlm", se = TRUE, linetype = "dashed", size = 1) +  # fit for aFMT and Placebo respectively
    geom_smooth(method = "rlm", se = TRUE, linetype = "dashed", size = 1, colour = "grey", aes(group = 1)) +  # fit for all
    #geom_hline(yintercept = median_wr, color = "grey50", linetype = "dotted") +
    labs(
        #title = "<i>Faecalibacterium prausnitzii</i> (SGB15322)",
        #title = expression(italic("Faecalibacterium prausnitzii") ~ "(SGB15322)"),
        title = expression(bolditalic("Faecalibacterium prausnitzii") ~ bold("(SGB15322)")),
        subtitle = paste0(
                    "Spearman's r = ", signif(cor_test_result$estimate, 2), ", p = ", signif(cor_test_result$p.value, 2)            
                    #"aFMT: Spearman's r = ", signif(cor_afmt$estimate, 2), ", p = ", signif(cor_afmt$p.value, 1),
                    #"\nPlacebo: Spearman's r = ", signif(cor_placebo$estimate, 2), ", p = ", signif(cor_placebo$p.value, 1)
                    ),
        x = "Abundance Change % (6–14)", 
        y = "Weight regain % (6–18)"
    ) +
    scale_color_manual(values = c("Placebo" = "#7796A7", "aFMT" = "#B9772B")) +
    #scale_color_manual(values = c("High regainer" = "#006241", "Low regainer" = "#66CDAA")) +
    #annotate("text", x = -0.3, y = median_wr + 50, label = "High regainers", color = "#006241", fontface = "bold") +
    #annotate("text", x = -0.3, y = median_wr - 50, label = "Low regainers", color = "#66CDAA", fontface = "bold") +
    #theme_minimal() +
    theme(
        panel.background = element_blank(),
        #axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "grey", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 14),
        #plot.title = ggtext::element_markdown(hjust = 0, size = 14, face = "plain"),
        plot.subtitle = element_text(hjust = 0, size = 14),
        axis.line = element_line(color = "black"),
        axis.title.x = element_text(size = 14, face = "bold", family="sans"),
        axis.title.y = element_text(size = 14, face = "bold", family="sans"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        #axis.ticks.x = element_blank()
        #legend.title = element_text( face = "bold"),
        legend.position = "right"
    ) 
p
# Save the plot to a PDF file
# ggsave("../Figure4e. Top2 Faecalibacterium prausnitzii (SGB15322) abundance_vs_weight_regain nonzero.pdf", plot = p, width = 6, height = 5)

# Output 6*5

addWorksheet(wb, "Faecalibacterium prausnitzii")
writeData(wb, "Faecalibacterium prausnitzii", predictors_top_df_filtered)
#saveWorkbook(wb, "Source Data Fig 4e.xlsx", overwrite = TRUE)



################################# Top 4 ########################################
############################     tryptophol         ############################
#top_feature <- importance_df$Feature[4]
top_feature <- "100004275_change_BD_B"
predictors_top_df <- predict_data7[, c("Diet_A","Treatment_A","Weight_regain6_18", top_feature)]
head(predictors_top_df)
names(predictors_top_df) <- c("Diet","Treatment","Weight_regain6_18", "AbundanceChange")

# Compute the median
#median_wr <- median(predictors_top_df$Weight_regain6_18, na.rm = TRUE)
# Add label column
#predictors_top_df$Regain_group <- ifelse(predictors_top_df$Weight_regain6_18 > median_wr, "High regainer", "Low regainer")
predictors_top_df$Treatment <- ifelse(predictors_top_df$Treatment == "1", "aFMT", "Placebo")

# Remove points with very extreme high or low values
#predictors_top_df_filtered <- predictors_top_df
predictors_top_df_filtered <- predictors_top_df %>% filter(Weight_regain6_18 <= 200)

# Calculate Pearson correlation
cor_test_result <- cor.test(
  predictors_top_df_filtered$AbundanceChange,
  predictors_top_df_filtered$Weight_regain6_18,
  method = "spearman"
)
cor_test_result$estimate  # spearman correlation coefficient (r)
cor_test_result$p.value

# compare the correlation between AbundanceChange and Weight_regain6_18 across two treatment groups (e.g., aFMT vs Placebo)
df_afmt <- predictors_top_df_filtered[predictors_top_df_filtered$Treatment == "aFMT", ]
df_placebo <- predictors_top_df_filtered[predictors_top_df_filtered$Treatment == "Placebo", ]
cor_afmt <- cor.test(df_afmt$AbundanceChange, df_afmt$Weight_regain6_18, method = "spearman")
cor_placebo <- cor.test(df_placebo$AbundanceChange, df_placebo$Weight_regain6_18, method = "spearman")
cor_afmt
cor_placebo


p <- ggplot(predictors_top_df_filtered, aes(x = AbundanceChange, y = Weight_regain6_18, colour = Treatment)) +
    geom_point(size = 2, alpha = 1) +
    #geom_smooth(method = "lm", se = TRUE, linetype = "dashed", size = 1) +   # fit for aFMT and Placebo respectively
    geom_smooth(method = "rlm", se = TRUE, linetype = "dashed", size = 1, colour = "grey", aes(group = 1)) +  # fit for all
    #geom_hline(yintercept = median_wr, color = "grey50", linetype = "dotted") +
    labs(
        title = "Tryptophol",
        subtitle = paste0(
                    "Spearman's r = ", signif(cor_test_result$estimate, 2), ", p = ", signif(cor_test_result$p.value, 2)
                    #"aFMT: Spearman's r = ", round(cor_afmt$estimate, 2), ", p = ", signif(cor_afmt$p.value, 1),
                    #"\nPlacebo: Spearman's r = ", round(cor_placebo$estimate, 2), ", p = ", signif(cor_placebo$p.value, 2)
                    ),
        x = "Abundance Change % (6–14)", 
        y = "Weight regain % (6–18)"
    ) +
    scale_color_manual(values = c("Placebo" = "#7796A7", "aFMT" = "#B9772B")) +
    #scale_color_manual(values = c("High regainer" = "#006241", "Low regainer" = "#66CDAA")) +
    #annotate("text", x = -0.3, y = median_wr + 50, label = "High regainers", color = "#006241", fontface = "bold") +
    #annotate("text", x = -0.3, y = median_wr - 50, label = "Low regainers", color = "#66CDAA", fontface = "bold") +
    #theme_minimal() +
    theme(
        panel.background = element_blank(),
        #axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "grey", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0, size = 14),
        axis.line = element_line(color = "black"),
        axis.title.x = element_text(size = 14, face = "bold", family="sans"),
        axis.title.y = element_text(size = 14, face = "bold", family="sans"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        #axis.ticks.x = element_blank()
        #legend.title = element_text( face = "bold"),
        legend.position = "right"
    ) 
p
# Save the plot to a PDF file
#ggsave("../Figure4e. Top4 Tryptophol abundance_vs_weight_regain.pdf", plot = p, width = 6, height = 5)

# Output 6*5

addWorksheet(wb, "Tryptophol")
writeData(wb, "Tryptophol", predictors_top_df_filtered)
#saveWorkbook(wb, "Source Data Fig 4e.xlsx", overwrite = TRUE)




################################# Top 4 ########################################
############################      Swap_rate_6_14 (%)       ############################
top_feature <- "Swap_B_D_percentage_strain_B"
           
predictors_top_df <- predict_data7[, c("Diet_A","Treatment_A","Weight_regain6_18", top_feature)]
head(predictors_top_df)
names(predictors_top_df) <- c("Diet","Treatment","Weight_regain6_18", "AbundanceChange")

# Compute the median
#median_wr <- median(predictors_top_df$Weight_regain6_18, na.rm = TRUE)
# Add label column
#predictors_top_df$Regain_group <- ifelse(predictors_top_df$Weight_regain6_18 > median_wr, "High regainer", "Low regainer")
predictors_top_df$Treatment <- ifelse(predictors_top_df$Treatment == "1", "aFMT", "Placebo")

# Remove points with very extreme high or low values
#predictors_top_df_filtered <- predictors_top_df
predictors_top_df_filtered <- predictors_top_df %>% filter(Weight_regain6_18 <= 200 & AbundanceChange > 0)

# Calculate Pearson correlation
cor_test_result <- cor.test(
  predictors_top_df_filtered$AbundanceChange,
  predictors_top_df_filtered$Weight_regain6_18,
  method = "spearman"
)
cor_test_result$estimate  # spearman correlation coefficient (r)
cor_test_result$p.value

# compare the correlation between AbundanceChange and Weight_regain6_18 across two treatment groups (e.g., aFMT vs Placebo)
df_afmt <- predictors_top_df_filtered[predictors_top_df_filtered$Treatment == "aFMT", ]
df_placebo <- predictors_top_df_filtered[predictors_top_df_filtered$Treatment == "Placebo", ]
cor_afmt <- cor.test(df_afmt$AbundanceChange, df_afmt$Weight_regain6_18, method = "spearman")
cor_placebo <- cor.test(df_placebo$AbundanceChange, df_placebo$Weight_regain6_18, method = "spearman")
cor_afmt
cor_placebo


p <- ggplot(predictors_top_df_filtered, aes(x = AbundanceChange, y = Weight_regain6_18, colour = Treatment)) +
    geom_point(size = 2, alpha = 1) +
    #geom_smooth(method = "lm", se = TRUE, linetype = "dashed", size = 1) +  # fit for aFMT and Placebo respectively
    geom_smooth(method = "rlm", se = TRUE, linetype = "dashed", size = 1, colour = "grey", aes(group = 1)) +  # fit for all
    #geom_hline(yintercept = median_wr, color = "grey50", linetype = "dotted") +
    labs(
        title = "Strain swap rate (6-14)",
        subtitle = paste0(
                    "Spearman's r = ", signif(cor_test_result$estimate, 2), ", p = ", signif(cor_test_result$p.value, 1)
                    #"aFMT: Spearman's r = ", round(cor_afmt$estimate, 2), ", p = ", signif(cor_afmt$p.value, 2),
                    #"\nPlacebo: Spearman's r = ", round(cor_placebo$estimate, 2), ", p = ", signif(cor_placebo$p.value, 2)
                    ),
        x = "Strain swap rate % (6–14)", 
        y = "Weight regain % (6–18)"
    ) +
    scale_color_manual(values = c("Placebo" = "#7796A7", "aFMT" = "#B9772B")) +
    #scale_color_manual(values = c("High regainer" = "#006241", "Low regainer" = "#66CDAA")) +
    #annotate("text", x = -0.3, y = median_wr + 50, label = "High regainers", color = "#006241", fontface = "bold") +
    #annotate("text", x = -0.3, y = median_wr - 50, label = "Low regainers", color = "#66CDAA", fontface = "bold") +
    #theme_minimal() +
    theme(
        panel.background = element_blank(),
        #axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "grey", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0, size = 14),
        axis.line = element_line(color = "black"),
        axis.title.x = element_text(size = 14, face = "bold", family="sans"),
        axis.title.y = element_text(size = 14, face = "bold", family="sans"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        #axis.ticks.x = element_blank()
        #legend.title = element_text( face = "bold"),
        legend.position = "right"
    ) 
p
# Save the plot to a PDF file
#ggsave("../Figure4e. Swap_rate_6_14 abundance_vs_weight_regain.pdf", plot = p, width = 6, height = 5)

# Output 6*5


addWorksheet(wb, "Strain swap rate (6-14)")
writeData(wb, "Strain swap rate (6-14)", predictors_top_df_filtered)
saveWorkbook(wb, "Source Data Fig 4e.xlsx", overwrite = TRUE)


