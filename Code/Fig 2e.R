###############################################################################
# We calcualte the engraftment rate per species
###############################################################################

setwd("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/5.matrix/")

# Load necessary library 
library(dplyr) 
library(tidyr) # Ensure tidyr is loaded for pivot_longer
library(pheatmap)
library(oddsratio)
library(stringr)
library(tibble)
library(broom)

species_data <- read.table("t__SGB14991_nGD_dynamic_matrix.tsv", row.names=1, header=T, sep="\t", check.names = FALSE)
species_data2 <- species_data %>% rownames_to_column(var = "SubjectID")

metadata <- read.table("../metadata2.tsv", header=T, sep="\t", check.names = FALSE)
metadata <- metadata %>% filter(Time != 18)



######################### Code for single SGB: Start #####################
# Use participant measurements at month 6 as covariates
# Weight change into quartiles
library(dplyr)
# Filter for time 0 and time 6
time_0 <- metadata %>% filter(Time == 0)
time_6 <- metadata %>% filter(Time == 6)
# Merge the data frames on SubjectID
merged_data <- merge(time_0, time_6, by = "SubjectID", suffixes = c("_0", "_6"))
# Calculate weight change
weight_change_df <- merged_data %>%
  mutate(Weight_change = Weight_6 - Weight_0) %>%
  select(SubjectID, Weight_change) %>% distinct()
metadata <- merge(metadata, weight_change_df, by = "SubjectID")
metadata$Weight_change_Quartile <- cut(
  metadata$Weight_change,
  breaks = quantile(metadata$Weight_change, probs = seq(0, 1, 0.25), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("1", "2", "3", "4"))
metadata$Weight_change_Quartile <- as.factor(metadata$Weight_change_Quartile)

# Categorize Age into quartiles
#time_6_age <- metadata[, c("SubjectID","Time","Age")] %>% filter(Time == 6)
metadata$Age_Quartile <- cut(
  metadata$Age,
  breaks = quantile(metadata$Age, probs = seq(0, 1, 0.25), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("1", "2", "3", "4"))
#metadata <- merge(metadata, time_6_age[, c("SubjectID","Age_Quartile")], by = "SubjectID")
metadata$Age_Quartile <- as.factor(metadata$Age_Quartile)


# Categorize WC_change into quartiles
time_0 <- metadata %>% filter(Time == 0)
time_6 <- metadata %>% filter(Time == 6)
# Merge the data frames on SubjectID
merged_data <- merge(time_0, time_6, by = "SubjectID", suffixes = c("_0", "_6"))
# Calculate weight change
WC_change_df <- merged_data %>%
  mutate(WC_change = WC_6 - WC_0) %>%
  select(SubjectID, WC_change) %>% distinct()
metadata <- merge(metadata, WC_change_df, by = "SubjectID")  
metadata$WC_change_Quartile <- cut(
  metadata$WC_change,
  breaks = quantile(metadata$WC_change, probs = seq(0, 1, 0.25), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("1", "2", "3", "4"))
metadata$WC_change_Quartile <- as.factor(metadata$WC_change_Quartile)
head(metadata)

# Categorize Richness into quartiles
#time_0 <- metadata %>% filter(Time == 0)
#time_6 <- metadata %>% filter(Time == 6)
# Merge the data frames on SubjectID
#merged_data <- merge(time_0, time_6, by = "SubjectID", suffixes = c("_0", "_6"))
# Calculate weight change
#Richness_change_df <- merged_data %>%
#  mutate(Richness_change = Richness_6 - Richness_0) %>%
#  select(SubjectID, Richness_change) %>% distinct()
#metadata <- merge(metadata, Richness_change_df, by = "SubjectID")  
#time_6_richness <- metadata[, c("SubjectID","Time","Richness")] %>% filter(Time == 6)
metadata$Richness_Quartile <- cut(
  metadata$Richness,
  breaks = quantile(metadata$Richness, probs = seq(0, 1, 0.25), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("1", "2", "3", "4"))
#metadata <- merge(metadata, time_6_richness[, c("SubjectID","Richness_Quartile")], by = "SubjectID")
metadata$Richness_Quartile <- as.factor(metadata$Richness_Quartile)
head(metadata)

metadata$Diet <- as.factor(metadata$Diet)
metadata$Gender <- as.factor(metadata$Gender)
metadata$Treatment <- as.factor(metadata$Treatment)
# Prepare metadata: select relevant columns and ensure SubjectID is available
md <- metadata %>% select(SubjectID, Gender, Diet, Treatment, Weight_change_Quartile, Age_Quartile, WC_change_Quartile, Richness_Quartile) %>% distinct()




# Merge species data with metadata by SubjectID
species_data3 <- merge(md, species_data2, by = "SubjectID")

# Create binary engraftment status: 1 if (A,B,D) meets condition, 0 otherwise.
species_data3 <- species_data3 %>%
  mutate(engrafted_status = ifelse((A == 0 & B == 1 & D == 1) | (A == 1 & B == 2 & D == 2), 1, 0))

# Create two lists: one for aFMT and one for Placebo
afmt_list <- species_data3 %>% filter(Treatment == "aFMT") %>% pull(engrafted_status)
placebo_list <- species_data3 %>% filter(Treatment == "Placebo") %>% pull(engrafted_status)
afmt_list
placebo_list


# Combine the data for logistic regression (only for aFMT and Placebo)
logistic_data <- species_data3 %>% filter(Treatment %in% c("aFMT", "Placebo"))
# Ensure Treatment is a factor and set Placebo as the reference level
logistic_data$Treatment <- factor(logistic_data$Treatment, levels = c("Placebo","aFMT"))
logistic_data$Diet <- relevel(logistic_data$Diet,ref = "2" )

model <- glm(engrafted_status ~ Treatment + Diet + Weight_change_Quartile + Age_Quartile + WC_change_Quartile + Richness_Quartile, data = logistic_data, family = "binomial")
model <- glm(engrafted_status ~ Treatment*Diet + Weight_change_Quartile + Age_Quartile + WC_change_Quartile + Richness_Quartile, data = logistic_data, family = "binomial")
summary(model)


# Extract the coefficient for TreatmentaFMT and calculate its odds ratio
coef_treatment <- coef(model)["TreatmentaFMT"]
odds_ratio <- exp(coef_treatment)

# Calculate the 95% confidence interval for the coefficient
ci <- exp(confint(model, "TreatmentaFMT"))

# Extract the p-value for the Treatment effect (comparing aFMT vs Placebo)
#p_val <- summary(model)$coefficients["TreatmentaFMT", "Pr(>|z|)"]

# Check that we have a valid row from the tidy output
LR_estimate  <- summary(model)$coefficients["TreatmentaFMT", "Estimate"]
LR_std_error <- summary(model)$coefficients["TreatmentaFMT", "Std. Error"]
LR_statistic <- summary(model)$coefficients["TreatmentaFMT", "z value"]
LR_p_value   <- summary(model)$coefficients["TreatmentaFMT", "Pr(>|z|)"]
LR_conf_low  <- round(as.numeric(ci)[1], 2)
LR_conf_high <- round(as.numeric(ci)[2], 2)

library(oddsratio)
# Calculate OR for specific increment step of continuous variable
or_glm(data = logistic_data, model = model)

######################### Code for single SGB: End #####################



# Path to the folder containing the files
# 5.matrix2 folder contains SGBs with at least 20% of all 90 participants
#folder_path <- "~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/5.matrix/"

folder_path <- "~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/5.matrix2_18subjects/"
# List all files in the folder (with full paths)
file_list <- list.files(folder_path, full.names = TRUE)

library(dplyr)
library(stringr)
library(tibble)
library(broom)

# Initialize a data frame to hold the results
ER_per_species_results <- data.frame(
  SGB = character(),
  gain_persistence = numeric(),
  swap_persistence = numeric(),
  subject_count = numeric(),
  ER_overall = numeric(),
  ER_afmt = numeric(),
  ER_placebo = numeric(),
  ER_afmt_diet1 = numeric(),
  ER_afmt_diet2 = numeric(),
  ER_afmt_diet3 = numeric(),
  ER_placebo_diet1 = numeric(),
  ER_placebo_diet2 = numeric(),
  ER_placebo_diet3 = numeric(),
  LR_odds_ratio = numeric(),
  LR_estimate = numeric(),
  LR_std_error = numeric(),
  LR_statistic = numeric(),
  LR_p_value = numeric(),
  LR_conf_low = numeric(),
  LR_conf_high = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each file
for (file in file_list) {
  # Extract the SGB code (using basename to remove path)
  SGB <- str_extract(basename(file), "SGB\\d+")
  
  # Read in the species data and convert rownames to a column
  species_data <- read.table(file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  species_data2 <- species_data %>% rownames_to_column(var = "SubjectID")
  
  # Prepare metadata: select relevant columns and ensure SubjectID is available
  md <- metadata %>% select(SubjectID, Gender, Diet, Treatment, Weight_change_Quartile, Age_Quartile, WC_change_Quartile, Richness_Quartile) %>% distinct()
  #md <- metadata %>% select(SubjectID, Diet, Treatment) %>% distinct()
  
  # Merge species data with metadata by SubjectID
  species_data3 <- merge(md, species_data2, by = "SubjectID")
  species_data3 <- species_data3 %>% drop_na()
  print(nrow(species_data3))
  # Process only if there are more than 50% participants in the merged dataset
  if (length(unique(species_data3$SubjectID)) >= 45) {
    #print(SGB)
    # Create binary engraftment status: 1 if (A, B, D) meets condition, 0 otherwise.
    species_data3 <- species_data3 %>%
      mutate(engrafted_status = ifelse((A == 0 & B == 1 & D == 1) | 
                                         (A == 1 & B == 2 & D == 2), 1, 0),
             gain_persistence = ifelse((A == 0 & B == 1 & D == 1), 1, 0),
             swap_persistence = ifelse((A == 1 & B == 2 & D == 2), 1, 0)                                                        
                                         )
    
    # Create two lists: one for aFMT and one for Placebo
    afmt_list <- species_data3 %>% select("SubjectID", "Treatment", "engrafted_status") %>% distinct() %>% filter(Treatment == "aFMT") %>% pull(engrafted_status)
    placebo_list <- species_data3 %>% select("SubjectID", "Treatment", "engrafted_status") %>% distinct() %>% filter(Treatment == "Placebo") %>% pull(engrafted_status)
    #print(afmt_list)
    #print(placebo_list)
    
    # For logistic regression, use only aFMT and Placebo data
    logistic_data <- species_data3 %>% filter(Treatment %in% c("aFMT", "Placebo"))
    #logistic_data <- full_join(species_data3[,-c(2,3)], metadata, by="SubjectID")
    # Set Placebo as the reference group
    logistic_data$Treatment <- factor(logistic_data$Treatment, levels = c("Placebo", "aFMT"))
    logistic_data$Diet <- factor(logistic_data$Diet)
    
    #print(logistic_data)
    # Initialize LR variables as NA
    LR_odds_ratio <- NA
    LR_estimate <- NA
    LR_std_error <- NA
    LR_statistic <- NA
    LR_p_value <- NA
    LR_conf_low <- NA
    LR_conf_high <- NA
    p_val <- NA
    
    # Only run the glm if each list has at least two non-zero (or non-NA) values
    if ((sum(!is.na(afmt_list) & afmt_list != 0) >= 1) && 
        (sum(!is.na(placebo_list) & placebo_list != 0) >= 1)) {
      
      # Run logistic regression: engrafted_status ~ Treatment
      model <- glm(engrafted_status ~ Treatment + Diet + Age_Quartile + Gender + WC_change_Quartile + Weight_change_Quartile + Richness_Quartile, data = logistic_data, family = "binomial")
      # Extract the coefficient for TreatmentaFMT and calculate its odds ratio
      coef_treatment <- coef(model)["TreatmentaFMT"]
      odds_ratio <- exp(coef_treatment)
      
      # Calculate the 95% confidence interval for the coefficient
      ci <- exp(confint(model, "TreatmentaFMT"))
      # log-CI
      # ci <- confint(model, "TreatmentaFMT")

      # Extract the p-value for the Treatment effect (comparing aFMT vs Placebo)
      p_val <- summary(model)$coefficients["TreatmentaFMT", "Pr(>|z|)"]
      
      # Check that we have a valid row from the tidy output
      LR_odds_ratio <-  odds_ratio
      LR_estimate  <- summary(model)$coefficients["TreatmentaFMT", "Estimate"]
      LR_std_error <- summary(model)$coefficients["TreatmentaFMT", "Std. Error"]
      LR_statistic <- summary(model)$coefficients["TreatmentaFMT", "z value"]
      LR_p_value   <- summary(model)$coefficients["TreatmentaFMT", "Pr(>|z|)"]
      LR_conf_low  <- round(as.numeric(ci)[1], 2)
      LR_conf_high <- round(as.numeric(ci)[2], 2)
    }
    
    # Calculate engraftment rates (as percentages)
    species_data3  <- species_data3 %>% select("SubjectID", "Treatment", "Diet", "engrafted_status", "gain_persistence", "swap_persistence") %>% distinct()
    gain_persistence  <- nrow(species_data3 %>% filter(gain_persistence == 1))
    swap_persistence  <- nrow(species_data3 %>% filter(swap_persistence == 1))
    engrafted_count  <- nrow(species_data3 %>% filter(engrafted_status == 1))
    engraftable_count <- nrow(species_data3)
    ER_overall <- ifelse(engraftable_count > 0, engrafted_count * 100 / engraftable_count, NA)
    
    # aFMT treatment
    engrafted_afmt   <- species_data3 %>% filter(Treatment == "aFMT", engrafted_status == 1)
    engraftable_afmt <- species_data3 %>% filter(Treatment == "aFMT")
    ER_afmt <- ifelse(nrow(engraftable_afmt) > 0, nrow(engrafted_afmt) * 100 / nrow(engraftable_afmt), NA)
    
    # Placebo treatment
    engrafted_placebo   <- species_data3 %>% filter(Treatment == "Placebo", engrafted_status == 1)
    engraftable_placebo <- species_data3 %>% filter(Treatment == "Placebo")
    ER_placebo <- ifelse(nrow(engraftable_placebo) > 0, nrow(engrafted_placebo) * 100 / nrow(engraftable_placebo), NA)
    
    # aFMT Diet groups
    engrafted_afmt_diet1   <- species_data3 %>% filter(Treatment == "aFMT", Diet == "1", engrafted_status == 1)
    engraftable_afmt_diet1 <- species_data3 %>% filter(Treatment == "aFMT", Diet == "1")
    ER_afmt_diet1 <- ifelse(nrow(engraftable_afmt_diet1) > 0, nrow(engrafted_afmt_diet1) * 100 / nrow(engraftable_afmt_diet1), NA)
    
    engrafted_afmt_diet2   <- species_data3 %>% filter(Treatment == "aFMT", Diet == "2", engrafted_status == 1)
    engraftable_afmt_diet2 <- species_data3 %>% filter(Treatment == "aFMT", Diet == "2")
    ER_afmt_diet2 <- ifelse(nrow(engraftable_afmt_diet2) > 0, nrow(engrafted_afmt_diet2) * 100 / nrow(engraftable_afmt_diet2), NA)
    
    engrafted_afmt_diet3   <- species_data3 %>% filter(Treatment == "aFMT", Diet == "3", engrafted_status == 1)
    engraftable_afmt_diet3 <- species_data3 %>% filter(Treatment == "aFMT", Diet == "3")
    ER_afmt_diet3 <- ifelse(nrow(engraftable_afmt_diet3) > 0, nrow(engrafted_afmt_diet3) * 100 / nrow(engraftable_afmt_diet3), NA)
    
    # Placebo Diet groups
    engrafted_placebo_diet1   <- species_data3 %>% filter(Treatment == "Placebo", Diet == "1", engrafted_status == 1)
    engraftable_placebo_diet1 <- species_data3 %>% filter(Treatment == "Placebo", Diet == "1")
    ER_placebo_diet1 <- ifelse(nrow(engraftable_placebo_diet1) > 0, nrow(engrafted_placebo_diet1) * 100 / nrow(engraftable_placebo_diet1), NA)
    
    engrafted_placebo_diet2   <- species_data3 %>% filter(Treatment == "Placebo", Diet == "2", engrafted_status == 1)
    engraftable_placebo_diet2 <- species_data3 %>% filter(Treatment == "Placebo", Diet == "2")
    ER_placebo_diet2 <- ifelse(nrow(engraftable_placebo_diet2) > 0, nrow(engrafted_placebo_diet2) * 100 / nrow(engraftable_placebo_diet2), NA)
    
    engrafted_placebo_diet3   <- species_data3 %>% filter(Treatment == "Placebo", Diet == "3", engrafted_status == 1)
    engraftable_placebo_diet3 <- species_data3 %>% filter(Treatment == "Placebo", Diet == "3")
    ER_placebo_diet3 <- ifelse(nrow(engraftable_placebo_diet3) > 0, nrow(engrafted_placebo_diet3) * 100 / nrow(engraftable_placebo_diet3), NA)
    
    # Append the results for this file into the results data frame
    ER_per_species_results <- rbind(ER_per_species_results, data.frame(
      SGB = SGB,
      gain_persistence = gain_persistence,
      swap_persistence = swap_persistence,
      subject_count = engraftable_count,
      ER_overall = ER_overall,
      ER_afmt = ER_afmt,
      ER_placebo = ER_placebo,
      ER_afmt_diet1 = ER_afmt_diet1,
      ER_afmt_diet2 = ER_afmt_diet2,
      ER_afmt_diet3 = ER_afmt_diet3,
      ER_placebo_diet1 = ER_placebo_diet1,
      ER_placebo_diet2 = ER_placebo_diet2,
      ER_placebo_diet3 = ER_placebo_diet3,
      LR_odds_ratio = LR_odds_ratio,
      LR_estimate = LR_estimate,
      LR_std_error = LR_std_error,
      LR_statistic = LR_statistic,
      LR_p_value = LR_p_value,
      LR_conf_low = LR_conf_low,
      LR_conf_high = LR_conf_high,
      stringsAsFactors = FALSE
    ))
  }
}

# Output the results data frame
print(ER_per_species_results)

# Extracted species with ER_overall > 0
ER_per_species_results2 <- ER_per_species_results %>% filter(ER_overall>0)
dim(ER_per_species_results)
# [1] 80 17
dim(ER_per_species_results2)
#[1] 79 17




###############################################################################
# Select species with best engraftment rate in aFMT group 
###############################################################################

# Only keep species with odds ratios > 1
ER_per_species_results4 <- ER_per_species_results3 %>% filter( LR_odds_ratio > 1)
dim(ER_per_species_results4)
#[1] 22 25

# remove species that have a higher ER in Placebo-MED than in aFMT-MED
species_remove <- c("Bacteroides_uniformis (SGB1836)", "Oscillibacter_sp_ER4 (SGB15254)", "Agathobaculum_butyriciproducens (SGB14993)", "Faecalibacterium_prausnitzii (SGB15318)", "Faecalicatena_fissicatena (SGB4871)", "Faecalibacterium_prausnitzii (SGB15316)", "Clostridiales_bacterium_KLE1615 (SGB5090)", "Gemmiger_formicilis (SGB15300)")
#Mediterraneibacter_faecis (SGB4563)
#Eubacterium_ramulus (SGB4959)
ER_per_species_results5 <- ER_per_species_results4 %>% filter(!(Species_SGB %in% species_remove))

p <- ggplot(ER_per_species_results5, aes(x = Species_SGB, 
    # Originial odds ratio
    y = LR_odds_ratio, ymin = LR_odds_ratio, ymax = LR_odds_ratio)) +
    # Log tranformation odds ratio
    #y = log(LR_odds_ratio), ymin = log(LR_conf_low), ymax = log(LR_conf_high))) +
    geom_pointrange(size = 0.5, color = "#1F78B4") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    coord_flip() +
    labs(x = "",
         y = "OR of succession\naFMT vs. Placebo",
         title = "Odds Ratio") +
    scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2) ) +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 10, face = "bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 12, color = "black", angle=0),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.title = element_blank(),  # Remove legend title
        legend.position = "right",
        strip.text = element_text(size = 16, face = "bold")  # Modify facet title style
    )
p

# Output source data
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "Fig 2e OR per species")
writeData(wb, "Fig 2e OR per species", ER_per_species_results5)
saveWorkbook(wb, "Source Data Fig 2e.xlsx", overwrite = TRUE)



###############################################################################
# Annotate species taxonomy
###############################################################################

# Create a stacked bar chart for species taxonomy
library(ggplot2)
library(cowplot)

# Only keep species with odds ratios > 1
#ER_per_species_results4 <- ER_per_species_results3 %>% filter( LR_odds_ratio > 1)
dim(ER_per_species_results5)
#[1] 15 25

unique(ER_per_species_results5$Phylum)
# [1] "p__Firmicutes"     "p__Actinobacteria"

# Remove the leading "s__" from the values in the Species_SGB column using the gsub()
ER_per_species_results5$Phylum <- gsub("^p__", "", ER_per_species_results5$Phylum)

# Create the stacked bar plot for species taxonomy
#distinct_colors <- c("#469990", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#ffd8b1", "#A65628")
distinct_colors <- c("#469990", "#377EB8")

p <- ggplot(ER_per_species_results5, aes(x = Species_SGB, fill = Phylum)) +
    geom_bar(alpha=0.8) +
    coord_flip() +
    labs(fill = "Phylum") +
    #theme_minimal() +
    scale_fill_manual(values = distinct_colors) +
    labs(title="Phylum")+
    theme(
        panel.background = element_blank(),
        #axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "white", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 11, face = "bold"),
        axis.title.x = element_text(size = 0, face = "bold"),
        axis.title.y = element_text(size = 0, face = "bold"),
        axis.text.x = element_text(size = 0, color = "white"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.ticks.x = element_blank(),
        legend.title = element_text( face = "bold"),
        legend.position = "right",
        strip.text = element_text(size = 16, face = "bold")  # Modify facet title style
    )
p

#addWorksheet(wb, "Fig 2e Strain")
#writeData(wb, "Fig 2e Strain", persistent_gain_loss_colonization)



###############################################################################
# Count the numbers of gain-persistence, swap-persistence, and participants
###############################################################################

library(ggplot2)
library(cowplot)

# Only keep species with odds ratios > 1
#ER_per_species_results4 <- ER_per_species_results3 %>% filter( LR_odds_ratio > 1)
dim(ER_per_species_results5)
#[1] 15 28

ER_per_species_results6 <- ER_per_species_results5 %>% 
            select("gain_persistence", "swap_persistence", "subject_count", "LR_odds_ratio", "Species_SGB")

library(tidyr)
library(dplyr)
# Melt the data
ER_per_species_results7 <- ER_per_species_results6 %>%
  pivot_longer(cols = c(gain_persistence, swap_persistence, subject_count), names_to = "Metric", values_to = "Value") %>%
  select(Species_SGB, LR_odds_ratio, Metric, Value)
head(ER_per_species_results7)

ER_per_species_results7$Species_SGB <- factor(ER_per_species_results7$Species_SGB)

# Create the stacked bar plot for SUPER_PATHWAY
#distinct_colors <- c("#469990", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#ffd8b1", "#A65628")
distinct_colors <- c("#4DAF4A", "#984EA3", "#FF7F00")

library(ggplot2)
library(dplyr)
library(patchwork)

# Ensure consistent species order
species_order <- levels(ER_per_species_results7$Species_SGB)

# Left panel: gain and swap persistence (plotted negative for leftward bars)
left_data <- ER_per_species_results7 %>%
  filter(Metric %in% c("gain_persistence", "swap_persistence")) %>%
  mutate(Value = -Value,
         Species_SGB = factor(Species_SGB, levels = species_order))

p1 <- ggplot(left_data, aes(x = Species_SGB, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.6, alpha=0.8) +
  coord_flip() +
  scale_y_continuous(labels = abs, limits = c(-10, 0), breaks = seq(-10, 0, by = 2) ) +
  scale_fill_manual(values = c("gain_persistence" = "#FF7F00", "swap_persistence" = "#984EA3")) +
  labs(y = "Persistence", x = NULL, title="# of events") +
  #theme_minimal(base_size = 11) +
  theme(
    panel.background = element_blank(),
    #axis.line = element_line(color = "grey"),
    panel.border = element_rect(color = "grey", fill = NA, size = 1),
    plot.title = element_text(size = 10, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 0, face = "bold"),
    axis.text.y = element_text(size = 8, hjust = 1),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(5, 2, 5, 5)  # less space on right
  )

# Right panel: subject_count
right_data <- ER_per_species_results7 %>%
  filter(Metric == "subject_count") %>%
  mutate(Species_SGB = factor(Species_SGB, levels = species_order))

p2 <- ggplot(right_data, aes(x = Species_SGB, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", width = 0.6) +
  coord_flip() +
  labs(y = "Subject Count", x = NULL, title="# of participants") +
  scale_fill_manual(values = c("subject_count" = "#4DAF4A")) +
  #theme_minimal(base_size = 11) +
  theme(
    panel.background = element_blank(),
    #axis.line = element_line(color = "grey"),
    panel.border = element_rect(color = "grey", fill = NA, size = 1),
    plot.title = element_text(size = 10, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 0, face = "bold"),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 2)  # less space on left
  )

# Combine plots with tighter spacing
(p1 | p2) + plot_layout(widths = c(1, 1), guides = "collect")




###############################################################################
# Bubble plot for succession rate per species, stratified by diets and aFMT
###############################################################################
dim(ER_per_species_results5)
#[1] 15 25

# Select SGBs in the same order of the forest plot
ER_per_species_results5$Species_SGB <- factor(ER_per_species_results5$Species_SGB,
  levels = ER_per_species_results5$Species_SGB[order(ER_per_species_results5$LR_odds_ratio, decreasing = F)])

# Reshape data into long format
df_long <- ER_per_species_results5 %>%
  select(SGB, ER_overall, ER_afmt, ER_placebo,
         ER_afmt_diet1, ER_afmt_diet2, ER_afmt_diet3,
         ER_placebo_diet1, ER_placebo_diet2, ER_placebo_diet3) %>%
  pivot_longer(cols = -SGB, names_to = "Group", values_to = "EngraftmentRate")

# Species name to SGB
taxonomy <- read.table("../taxonomy_all_SGB.txt", header=T, sep="\t", check.names = FALSE)
df_long2 <- merge(taxonomy[,-1], df_long, by="SGB")

# Define the column order
df_long2$Group <- factor(df_long2$Group, levels=c("ER_overall","ER_afmt","ER_placebo",
            "ER_afmt_diet1", "ER_placebo_diet1",
            "ER_afmt_diet2","ER_placebo_diet2",
            "ER_afmt_diet3", "ER_placebo_diet3"
            ))

# Connecting species and SGB
df_long2 <- df_long2 %>% mutate(Species_SGB = paste0(Species, " (", SGB, ")"))

# Remove the leading "s__" from the values in the Species_SGB column using the gsub()
df_long2$Species_SGB <- gsub("^s__", "", df_long2$Species_SGB)

# Filter df_long2 to include only SGBs that are in ER_per_species_results4$Species_SGB
df_long2_filtered <- df_long2 %>%
    filter(Species_SGB %in% levels(ER_per_species_results5$Species_SGB)) %>%
    # Reorder the factor levels of Species_SGB in df_long2 to match ER_per_species_results4
    mutate(Species_SGB = factor(Species_SGB, levels = levels(ER_per_species_results5$Species_SGB)))

df_long2_filtered$Group <- factor(df_long2_filtered$Group, 
        levels=c("ER_overall","ER_afmt","ER_placebo","ER_afmt_diet1","ER_placebo_diet1",
        "ER_afmt_diet2", "ER_placebo_diet2","ER_afmt_diet3","ER_placebo_diet3"))

# Define a named vector of 6 new color codes
new_colors <- c(
    "ER_overall"      = "#8B0000",  # Deep DarkRed
    "ER_afmt"         = "#D95F02",  # Strong gold-orange
    "ER_placebo"      = "#4F7CAC",  # Bold blue-gray
    "ER_afmt_diet1"   = "#E6AB02",  # Vivid bright orange/gold
    "ER_afmt_diet2"   = "#008B8B",  # Strong dark cyan
    "ER_afmt_diet3"   = "#66A61E",  # Rich teal
    "ER_placebo_diet1"= "#E6AB02",  # (Same as aFMT_diet1)
    "ER_placebo_diet2"= "#008B8B",  # (Same as aFMT_diet2)
    "ER_placebo_diet3"= "#66A61E"   # (Same as aFMT_diet3)
)

ggplot(df_long2_filtered, aes(x = Group, y = Species_SGB)) +
    geom_point(aes(size = EngraftmentRate, color = Group, alpha=EngraftmentRate/50) ) +
    #scale_size_continuous(range = c(2, 10), na.value = 0) +
    #scale_color_gradient(low = "white", high = new_colors, na.value = "grey90") +
    scale_color_manual(values = new_colors) +
    labs(title = "Within-group % of participants",
         x = "",
         y = "",
         size = "% of participants\nwith succession",
         color = "Engrafted groups") +
    guides(color = "none", alpha = "none") +
    scale_x_discrete(
        labels = c("Overall", "aFMT", "Placebo","aFMT_HDG","Placebo_HDG",
                   "aFMT_MED","Placebo_MED", "aFMT_GreenMED", "Placebo_GreenMED")) + # Set custom x-axis labels
    theme(
        panel.background = element_blank(),
        axis.line = element_line(color = "grey80"),
        panel.border = element_rect(color = "grey80", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 10, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 8, color = "black", face = "bold",vjust = 1,hjust = 0.05, angle= -30),
        axis.text.y = element_text(size = 8, color = "black", face = "bold"),
        #legend.title = element_blank(),  # Remove legend title
        legend.position = "right",
        strip.text = element_text(size = 16, face = "bold")  # Modify facet title style
    )

