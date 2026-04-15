
###############################################################################
# Calculate the decrease-persistence and increase-persistence of metabolites
###############################################################################

setwd("~/mydata/Project_aFMT/Metabolites/fecal")
# Load necessary library
library(dplyr)
library(tidyr)

# Read the dataset (assuming it is in a TSV format with one row per person)
data <- read.table("Batch_normalizedData_90p.tsv", header=T, sep="\t", check.names = FALSE)
data[is.na(data)] <- 0
# Create the new ID column
data$sno_time <- paste0(data$SubjectID, "_", data$Month)
# Move it to the first column
data <- data[, c(ncol(data), 1:(ncol(data)-1))]
head(data[, 1:10])

metadata <- read.table("../metadata.tsv", header=T, sep="\t", check.names = FALSE)
merged_data <- merge(metadata[, !(names(metadata) %in% c("SubjectID", "BMI"))], data, by = "sno_time")
merged_data[1:4, 1:18]
#  sno_time SubjectID Time   Age Weight      BMI  WC Gender Richness Diet Treatment       CLIENT_MATRIX CLIENT_SAMPLE_ID
#1    105_0       105    0 52.93  101.0 34.14008 112      0      116    2   Placebo Feces, Flash Frozen             A105
#2   105_14       105   14 54.10   88.9 30.05003  98      0      282    2   Placebo Feces, Flash Frozen             D105
#3   105_18       105   18 54.43   89.2 30.15143  99      0      255    2   Placebo Feces, Flash Frozen             C105
#4    105_6       105    6 53.43   88.6 29.94862  99      0      206    2   Placebo Feces, Flash Frozen             B105
#  CLIENT_SAMPLE_NUMBER Month PARENT_SAMPLE_NAME        30        35
#1                  785     0         BRIG-05719 6.7308097 0.8297286
#2                   35    14         BRIG-05685 0.0000000 1.0078231
#3                  261    18         BRIG-05694 0.3789397 1.2314464
#4                  334     6         BRIG-05697 0.3376543 1.5759851

merged_data <- merged_data %>% filter(Month != "18")

# Identify metabolite columns (assuming they start from column 8)
metabolite_cols <- colnames(merged_data)[18:ncol(merged_data)]

# Function to determine quartile assignment
assign_quartile <- function(value, q1, q2, q3, q4) {
  if (is.na(value)) {
    return(NA)
  } else if (value <= q1) {
    return("Q1")
  } else if (value <= q2) {
    return("Q2")
  } else if (value <= q3) {
    return("Q3")
  } else {
    return("Q4")
  }
}

# Function to categorize subjects by Gain, Loss, Persistent Gain, and Persistent Loss
perform_subject_transition_analysis <- function(data, results_df, diet_label, treatment_label) {
  
  # Initialize an empty dataframe for storing per-subject results
  subject_results <- data.frame()
  
  for (metabolite in metabolite_cols) {
    
    # Select relevant columns and remove missing values
    quartile_data <- data %>%
      select(SubjectID, GROUP_NUMBER, all_of(metabolite)) %>%
      distinct() %>%
      spread(GROUP_NUMBER, all_of(metabolite)) #%>% na.omit()
    
    # Compute quartiles for the metabolite across all subjects
    quartiles <- quantile(data[[metabolite]], probs = c(0.25, 0.5, 0.75, 1), na.rm = TRUE)
    q1 <- quartiles[1]
    q2 <- quartiles[2]
    q3 <- quartiles[3]
    q4 <- quartiles[4]  # Maximum value, ensuring full range
    
    # Assign quartiles to each subject for Time A, B, and D
    quartile_data <- quartile_data %>%
      mutate(
        Quartile_A = sapply(A, assign_quartile, q1, q2, q3, q4),
        Quartile_B = sapply(B, assign_quartile, q1, q2, q3, q4),
        Quartile_D = sapply(D, assign_quartile, q1, q2, q3, q4)
      )
    
    # Add Gain, Loss, Persistent Gain, Persistent Loss Flags
    quartile_data <- quartile_data %>%
      mutate(
        Gain = as.numeric((Quartile_A == "Q1" & (Quartile_B == "Q4" | Quartile_B == "Q3")) |
                          (Quartile_A == "Q2" & Quartile_B == "Q4")),
        
        Loss = as.numeric((Quartile_A == "Q4" & (Quartile_B == "Q1" | Quartile_B == "Q2")) |
                          (Quartile_A == "Q3" & Quartile_B == "Q1")),
        
        Persistent_Gain = as.numeric((Quartile_A == "Q1" & (Quartile_B == "Q4" | Quartile_B == "Q3") & (Quartile_D == "Q4" | Quartile_D == "Q3")) |
                                     (Quartile_A == "Q2" & Quartile_B == "Q4" & Quartile_D == "Q4")),
        
        Persistent_Loss = as.numeric((Quartile_A == "Q4" & (Quartile_B == "Q1" | Quartile_B == "Q2") & (Quartile_D == "Q1" | Quartile_D == "Q2")) |
                                     (Quartile_A == "Q3" & Quartile_B == "Q1"))
      )
    
    # Store per-metabolite subject-level results
    subject_results <- rbind(subject_results, quartile_data %>% select(SubjectID, Gain, Loss, Persistent_Gain, Persistent_Loss))
  }
  
  # Summarize data by SubjectID
  subject_summary <- subject_results %>%
    group_by(SubjectID) %>%
    summarise(
      Total_Metabolites = n(),
      Gain_Count = sum(Gain, na.rm = TRUE),
      Loss_Count = sum(Loss, na.rm = TRUE),
      Persistent_Gain_Count = sum(Persistent_Gain, na.rm = TRUE),
      Persistent_Loss_Count = sum(Persistent_Loss, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Gain_Percent = ifelse(Total_Metabolites > 0, (Gain_Count / Total_Metabolites) * 100, 0),
      Loss_Percent = ifelse(Total_Metabolites > 0, (Loss_Count / Total_Metabolites) * 100, 0),
      Persistent_Gain_Percent = ifelse(Total_Metabolites > 0, (Persistent_Gain_Count / Total_Metabolites) * 100, 0),
      Persistent_Loss_Percent = ifelse(Total_Metabolites > 0, (Persistent_Loss_Count / Total_Metabolites) * 100, 0)
    ) %>%
    mutate(Diet = diet_label, Treatment = treatment_label)  # Add grouping labels
  
  # Bind subject-level data to results
  results_df <- rbind(results_df, subject_summary)
  
  return(results_df)
}

# Initialize results data frames
results_overall <- data.frame(SubjectID = character(), Diet = character(), Treatment = character(),
                              Total_Metabolites = numeric(), Gain_Count = numeric(), Gain_Percent = numeric(),
                              Loss_Count = numeric(), Loss_Percent = numeric(),
                              Persistent_Gain_Count = numeric(), Persistent_Gain_Percent = numeric(),
                              Persistent_Loss_Count = numeric(), Persistent_Loss_Percent = numeric(),
                              stringsAsFactors = FALSE)

results_by_diet <- results_overall
results_by_diet_treatment <- results_overall

# 1. Perform transition analysis for the entire dataset (Overall)
results_overall <- perform_subject_transition_analysis(merged_data, results_overall, "Overall", "Overall")

# 2. Perform transition analysis grouped by Diet
for (diet_group in unique(merged_data$Diet)) {
  data_subset <- merged_data %>% filter(Diet == diet_group)
  results_by_diet <- perform_subject_transition_analysis(data_subset, results_by_diet, diet_group, "Overall")
}

# 3. Perform transition analysis grouped by Diet and Treatment
for (diet_group in unique(merged_data$Diet)) {
  for (treatment_group in unique(merged_data$Treatment)) {
    data_subset <- merged_data %>% filter(Diet == diet_group, Treatment == treatment_group)
    results_by_diet_treatment <- perform_subject_transition_analysis(data_subset, results_by_diet_treatment, diet_group, treatment_group)
  }
}

# Save results to separate CSV files
#write.csv(results_overall, "subject_transition_analysis_overall.csv", row.names = FALSE)
#write.csv(results_by_diet, "subject_transition_analysis_by_diet.csv", row.names = FALSE)
#write.csv(results_by_diet_treatment, "subject_transition_analysis_by_diet_treatment.csv", row.names = FALSE)

# Print first few rows of each results file
print("Overall Results:")
print(head(results_overall))

print("By Diet Results:")
print(head(results_by_diet))

print("By Diet & Treatment Results:")
print(head(results_by_diet_treatment))


library(reshape2)
#metadata2 <- metadata %>% filter(Time != "18")
results_by_diet_treatment2 <- full_join(metadata, results_by_diet_treatment[,-c(2:6, 11,12)], by="SubjectID")
head(results_by_diet_treatment2)

results_by_diet_treatment3 <- results_by_diet_treatment2 %>% filter(Time != "18")



###############################################################################
# adjusted for covariables

# Weight change into quartiles
library(dplyr)
# Filter for time 0 and time 6
time_0 <- results_by_diet_treatment2 %>% filter(Time == 0)
time_6 <- results_by_diet_treatment2 %>% filter(Time == 6)
# Merge the data frames on SubjectID
merged_data_0_6 <- merge(time_0, time_6, by = "SubjectID", suffixes = c("_0", "_6"))
# Calculate weight change
weight_change_df <- merged_data_0_6 %>%
  mutate(Weight_change = Weight_6 - Weight_0) %>%
  select(SubjectID, Weight_change) %>% distinct()
results_by_diet_treatment3 <- merge(results_by_diet_treatment3, weight_change_df, by = "SubjectID")
results_by_diet_treatment3$Weight_change_Quartile <- cut(
  results_by_diet_treatment3$Weight_change,
  breaks = quantile(results_by_diet_treatment3$Weight_change, probs = seq(0, 1, 0.25), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("1", "2", "3", "4"))
results_by_diet_treatment3$Weight_change_Quartile <- as.factor(results_by_diet_treatment3$Weight_change_Quartile)

# Categorize WC_change into quartiles
time_0 <- results_by_diet_treatment2 %>% filter(Time == 0)
time_6 <- results_by_diet_treatment2 %>% filter(Time == 6)
# Merge the data frames on SubjectID
merged_data <- merge(time_0, time_6, by = "SubjectID", suffixes = c("_0", "_6"))
# Calculate WC change
WC_change_df <- merged_data %>%
  mutate(WC_change = WC_6 - WC_0) %>%
  select(SubjectID, WC_change) %>% distinct()
results_by_diet_treatment3 <- merge(results_by_diet_treatment3, WC_change_df, by = "SubjectID")  
results_by_diet_treatment3$WC_change_Quartile <- cut(
  results_by_diet_treatment3$WC_change,
  breaks = quantile(results_by_diet_treatment3$WC_change, probs = seq(0, 1, 0.25), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("1", "2", "3", "4"))
results_by_diet_treatment3$WC_change_Quartile <- as.factor(results_by_diet_treatment3$WC_change_Quartile)


# Categorize Age into quartiles with only data at month 6
time_6_df <- results_by_diet_treatment2 %>% filter(Time == 6)
time_6_df$Age_Quartile <- cut(
  time_6_df$Age,
  breaks = quantile(time_6_df$Age, probs = seq(0, 1, 0.25), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("1", "2", "3", "4"))
results_by_diet_treatment3 <- merge(results_by_diet_treatment3, time_6_df[,c("SubjectID", "Age_Quartile")], by = "SubjectID")
results_by_diet_treatment3$Age_Quartile <- as.factor(results_by_diet_treatment3$Age_Quartile)


# Categorize Richness into quartiles with only data at month 6
time_6_df <- results_by_diet_treatment2 %>% filter(Time == 6)
time_6_df$Richness_Quartile <- cut(
  time_6_df$Richness,
  breaks = quantile(time_6_df$Richness, probs = seq(0, 1, 0.25), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("1", "2", "3", "4"))
results_by_diet_treatment3 <- merge(results_by_diet_treatment3, time_6_df[,c("SubjectID", "Richness_Quartile")], by = "SubjectID")
results_by_diet_treatment3$Richness_Quartile <- as.factor(results_by_diet_treatment3$Richness_Quartile)
head(results_by_diet_treatment3)

results_by_diet_treatment3 <- results_by_diet_treatment3 %>% select("SubjectID", "Diet", "Treatment", "Age_Quartile", "Gender", "WC_change_Quartile", "Weight_change_Quartile", "Richness_Quartile", "Gain_Percent", "Loss_Percent", "Persistent_Gain_Percent", "Persistent_Loss_Percent") %>% distinct()

# Linear model for interaction test
results_by_diet_treatment3$Diet <- as.factor(results_by_diet_treatment3$Diet)
results_by_diet_treatment3$Diet <- relevel(results_by_diet_treatment3$Diet, ref = "2")
# Changing treatment to numeric (aFMT=1, Placebo=0) in the data frame
results_by_diet_treatment3 <- results_by_diet_treatment3 %>% mutate(aFMT = ifelse(Treatment == "aFMT", 1, 0))
results_by_diet_treatment3$aFMT <- as.factor(results_by_diet_treatment3$aFMT)
results_by_diet_treatment3$Gender <- as.factor(results_by_diet_treatment3$Gender)



###############################################################################
# Decrease-Persistence and Increase-Persistence, Overall
###############################################################################

# Persistent Gain/Upregulated model
metabolite_persistent_gain <- results_by_diet_treatment3[,c("SubjectID", "Diet", "Age_Quartile", "Gender", "WC_change_Quartile","Weight_change_Quartile","Richness_Quartile","aFMT", "Persistent_Gain_Percent")]
metabolite_persistent_gain_model <- lm(Persistent_Gain_Percent ~ aFMT*Diet + Age_Quartile + Gender + WC_change_Quartile + Weight_change_Quartile + Richness_Quartile, metabolite_persistent_gain, na.action = na.exclude) 
metabolite_persistent_gain$predicted_value <- predict(metabolite_persistent_gain_model)
lsmeans_metabolite_persistent_gain <- emmeans(metabolite_persistent_gain_model, ~ aFMT)
lsmeans_metabolite_persistent_gain_df <- as.data.frame(lsmeans_metabolite_persistent_gain)[, c("aFMT", "emmean")]
# Merge the LS‐Means Back into Your Data:
metabolite_persistent_gain <- metabolite_persistent_gain %>% left_join(lsmeans_metabolite_persistent_gain_df, by = "aFMT")
# Calculate the Re-centered (Adjusted) Raw Values:
names(metabolite_persistent_gain)[9] <- "Value"
metabolite_persistent_gain <- metabolite_persistent_gain %>% mutate(adjusted_raw = Value - (predicted_value - emmean))



# Persistent Loss/Downregulated model
metabolite_persistent_loss <- results_by_diet_treatment3[,c("SubjectID", "Diet", "Age_Quartile", "Gender", "WC_change_Quartile","Weight_change_Quartile","Richness_Quartile","aFMT", "Persistent_Loss_Percent")]
metabolite_persistent_loss_model <- lm(Persistent_Loss_Percent ~ aFMT*Diet + Age_Quartile + Gender + WC_change_Quartile + Weight_change_Quartile + Richness_Quartile, metabolite_persistent_loss, na.action = na.exclude) 
metabolite_persistent_loss$predicted_value <- predict(metabolite_persistent_loss_model)
lsmeans_metabolite_persistent_loss <- emmeans(metabolite_persistent_loss_model, ~ aFMT)
lsmeans_metabolite_persistent_loss_df <- as.data.frame(lsmeans_metabolite_persistent_loss)[, c("aFMT", "emmean")]
# Merge the LS‐Means Back into Your Data:
metabolite_persistent_loss <- metabolite_persistent_loss %>% left_join(lsmeans_metabolite_persistent_loss_df, by = "aFMT")
# Calculate the Re-centered (Adjusted) Raw Values:
names(metabolite_persistent_loss)[9] <- "Value"
metabolite_persistent_loss <- metabolite_persistent_loss %>% mutate(adjusted_raw = Value - (predicted_value - emmean))


metabolite_persistent_gain$Metric <- "Increase-Persistence" 
metabolite_persistent_loss$Metric <- "Decrease-Persistence" 
metabolite_persistent_gainloss <- rbind(metabolite_persistent_gain, metabolite_persistent_loss)

metabolite_persistent_gainloss <- metabolite_persistent_gainloss %>% filter(adjusted_raw > 0 |adjusted_raw == 0)

metabolite_persistent_gainloss$Metric <- factor(metabolite_persistent_gainloss$Metric, levels = c("Increase-Persistence","Decrease-Persistence"))
# Stratified by Treatment
ggplot(metabolite_persistent_gainloss, aes(x = Metric, y = predicted_value, fill = aFMT)) +
    geom_boxplot(outliers = T, outlier.color = "grey", outlier.shape = 16, width = 0.7) +
    scale_fill_manual(
        #values = c("1" = "#EED6B7", "2" = "#00788B", "3" = "#049B86") # Custom colors
        values = c("0" = "#7796A7", "1" = "#B9772B"), # Custom colors
        labels = c("0" = "Placebo", "1" = "aFMT")
    ) +
    #scale_x_discrete(labels = c("Decrease-Persistence", "Increase-Persistence") ) +
    labs(
        title = "Metabolite change induced by diets and aFMT",
        y = "% of metabolites", x = ""
    ) +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 0, face = "bold"),
        axis.title.x = element_text(size = 24, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 20, color = "black", face = "bold", angle=30, vjust = 0.5),
        axis.text.y = element_text(size = 20, color = "black", face = "bold"),
        legend.title = element_blank(),  # Remove legend title
        legend.text = element_text(size = 22, color = "black"),
        legend.position = "right",
        strip.text = element_text(size = 16, face = "bold")  # Modify facet title style
    )
# Output 10*3


# Output source data
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "Increased-Persist Metabolites")
writeData(wb, "Increased-Persist Metabolites", metabolite_persistent_gain)

addWorksheet(wb, "Decreased-Persist Metabolites")
writeData(wb, "Decreased-Persist Metabolites", metabolite_persistent_loss)

saveWorkbook(wb, "Source Data Fig 3c.xlsx", overwrite = TRUE)




###############################################################################
# Decrease-Persistence and Increase-Persistence, stratified by Diet
###############################################################################

# Create a separate annotation data frame
p_int_df <- data.frame(
  Diet = factor("2", levels = c("1", "2", "3")),  # Must match the original factor levels
  predicted_value = c(22, 22),  # Adjust this y-value if needed
  Metric = c("Increase-Persistence", "Decrease-Persistence"),
  label = c("p.int = 0.01", "p.int = 0.04")
)

p_increase_persistence_df <- data.frame(
  Diet = factor("3", levels = c("1", "2", "3")),  # Must match the original factor levels
  predicted_value = c(10),  # Adjust this y-value if needed
  Metric = c("Increase-Persistence"),
  label = c("p = 0.015")
)

p_decrease_persistence_df <- data.frame(
  Diet = factor("3", levels = c("1", "2", "3")),  # Must match the original factor levels
  predicted_value = c(18),  # Adjust this y-value if needed
  Metric = c("Decrease-Persistence"),
  label = c("p = 0.013")
)


metabolite_persistent_gainloss$Metric <- factor(metabolite_persistent_gainloss$Metric, levels = c("Decrease-Persistence", "Increase-Persistence"))
# Stratified by Diet and Treatment
ggplot(metabolite_persistent_gainloss, aes(x = factor(Diet), y = predicted_value, fill = aFMT)) +
    geom_boxplot(outliers = T, outlier.color = "grey", outlier.shape = 16, width = 0.7) +
    #scale_y_continuous(limits = c(0, 10), breaks = seq(0, 20, by = 5) ) +
    scale_fill_manual(
        #values = c("1" = "#EED6B7", "2" = "#00788B", "3" = "#049B86") # Custom colors
        values = c("0" = "#7796A7", "1" = "#B9772B"), labels=c("0"="Placebo","1"="aFMT")
    ) +
    scale_x_discrete(
        labels = c("1" = "HDG", "2" = "MED", "3" = "GreenMED") ) +
    labs(
        title = "Metabolites change induced by diets and aFMT",
        x = "",
        y = "% of metabolites"
    ) +
    facet_wrap(~ Metric, scales = "fixed",
        labeller = as_labeller(c( "Increase-Persistence" = "Increase-Persistence\n(0-6-14)", 
                                  "Decrease-Persistence" = "Decrease-Persistence\n(0-6-14)"))
    ) +  # Separate plots by metabolites Change
    geom_text( data = p_int_df, aes(x = Diet, y = predicted_value, label = label),
        inherit.aes = FALSE, size = 7
    #fontface = "bold"
      ) +
    geom_text( data = p_increase_persistence_df, aes(x = Diet, y = predicted_value, label = label),
        inherit.aes = FALSE, size = 5
    #fontface = "bold"
      ) +
    geom_text( data = p_decrease_persistence_df, aes(x = Diet, y = predicted_value, label = label),
        inherit.aes = FALSE, size = 5
    #fontface = "bold"
      ) +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 0, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 28, vjust = 1, face = "bold"),
        axis.text.x = element_text(size = 20, color = "black", face = "bold", angle=0),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        legend.title = element_blank(),  # Remove legend title
        legend.text = element_text(size = 20, color = "black"),
        legend.position = "none",
        strip.text = element_text(size = 20, face = "bold")  # Modify facet title style
    )
# Output 8*9 


