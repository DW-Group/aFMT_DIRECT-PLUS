
####################################################################################
# Weight change
####################################################################################

setwd("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/")

# Reading the metadata
library(readr)
library(dplyr)
library(reshape2)

metadata <- read_tsv(file = "metadata2.tsv", col_names = TRUE, show_col_types = FALSE)
head(metadata)
#  SubjectID sno_time Time   Age Weight      BMI  WC Gender Richness Diet Treatment Percentage
#1         4      4_0    0 69.75   90.5 29.21617 116      0      214    3      aFMT   31.70732
#2         4     4_14   14 70.92   85.9 27.73115 105      0      276    3      aFMT   31.70732
#3         4     4_18   18 71.25   85.9 27.73115 108      0      242    3      aFMT   31.70732
#4         4      4_6    6 70.25   82.6 26.66581 106      0      209    3      aFMT   31.70732
#5        13     13_0    0 50.93   86.7 34.73001 107      1      128    2   Placebo   57.14286
#6        13    13_14   14 52.10   84.1 33.68851 107      1      132    2   Placebo   57.14286

# Calculate weight change (weight at time 6 - weight at time 0)
metadata$Weight_Baseline <- 0
metadata <- metadata %>%
  group_by(SubjectID) %>%
  mutate(Weight_Change0_6 = Weight[Time == 6] - Weight[Time == 0], 
        Weight_Change6_14 = Weight[Time == 14] - Weight[Time == 6],
        Weight_Change6_18 = Weight[Time == 18] - Weight[Time == 6], 
        Weight_regain6_14 = Weight_Change6_14*100/abs(Weight_Change0_6), 
        Weight_regain6_18 = Weight_Change6_18*100/abs(Weight_Change0_6) ) %>%
  ungroup()
  
# Ensure 'Diet' is treated as a factor for proper faceting 
metadata$Diet <- factor(metadata$Diet, levels = c(1, 2, 3), labels = c("HDG", "MED", "GreenMED")) 
# Convert 'Time' to a factor to treat time points as discrete categories 
metadata$Time <- factor(metadata$Time, levels = c(6, 14, 18), labels = c("6", "14", "18"))


weight_change <- metadata[, c("SubjectID","Time", "Diet","Treatment", "Weight_Baseline", "Weight_regain6_14", "Weight_regain6_18")]

weight_change_melted <- melt(weight_change, id=c("SubjectID","Time", "Diet","Treatment"))
# There are two outliers in HDG, subject IDs 124 and 125. try to delet them 
# Remove outliers
weight_change_melted2 <- weight_change_melted %>% filter(SubjectID != 124 & SubjectID != 125)

#data_phase1 <- weight_change_melted %>% filter(variable == "Weight_Baseline"| variable == "Weight_Change0_6")
data_phase2 <- weight_change_melted2 %>% filter(variable == "Weight_Baseline"| variable == "Weight_regain6_14")
data_phase3 <- weight_change_melted2 %>% filter(variable == "Weight_regain6_14"| variable == "Weight_regain6_18")


library(dplyr)

means_time_treatment <- weight_change_melted2 %>%
  # First, remove rows with missing Time values
  filter(!is.na(variable) & variable != "NA") %>%
  group_by(variable, Diet, Treatment) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value   = sd(value, na.rm = TRUE),
    n          = n(),
    .groups    = "drop"
  )

# Display the rows
means_time_treatment
# A tibble: 18 × 6
#   variable           Diet     Treatment mean_value sd_value     n
#   <fct>              <fct>    <chr>          <dbl>    <dbl> <int>
# 1 Weight_Baseline    HDG      Placebo       0          0       32
# 2 Weight_Baseline    HDG      aFMT          0          0       32
# 3 Weight_Baseline    MED      Placebo       0          0       72
# 4 Weight_Baseline    MED      aFMT          0          0       68
# 5 Weight_Baseline    GreenMED Placebo       0          0       80
# 6 Weight_Baseline    GreenMED aFMT          0          0       76

# Merge the individual-level and mean-level figures together 

library(ggplot2)
ggplot() +
    # Add a reference horizontal dashed line at y = 0
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 1) +
    # Phase 1: Time 0 to 6 -- all lines/points drawn in a single color (blue)
    # Connect the mean points for each Treatment group over Time
    geom_line(data = means_time_treatment,
            aes(x = variable, y = mean_value, group = Treatment, color = Treatment),
            alpha = 1, size = 1.5) +
    # Plot the mean points
    geom_point(data = means_time_treatment,
             aes(x = variable, y = mean_value, group = Treatment, color = Treatment),
             size = 3) +
    # Phase 2: Time 6 and 14 -- lines/points colored by Treatment
    geom_line(data = data_phase2,
              aes(x = variable, y = value, group = SubjectID, color = Treatment),
              alpha = 0.1, size = 1) +
    geom_point(data = data_phase2,
               aes(x = variable, y = value, group = SubjectID, color = Treatment),
               alpha = 0.1, size = 1) +
    # Phase 3: Time 14 and later -- lines/points colored by Treatment
    geom_line(data = data_phase3,
              aes(x = variable, y = value, group = SubjectID, color = Treatment),
              alpha = 0.1,  size = 1) +
    geom_point(data = data_phase3,
               aes(x = variable, y = value, group = SubjectID, color = Treatment),
               alpha = 0.1, size = 1) +
    scale_x_discrete(labels = c("6", "14", "18") ) +
    scale_color_manual(values = c("HDG" = "#D0495B", "MED" = "#4DBBD6", "GreenMED" = "#049B86", 
                                  "Placebo" = "#7796A7", "aFMT" = "#B9772B" )) +
    labs(x = "Month",y = "Weight regain % (from weight loss)") +
    # Facet by Diet so that each diet group appears in its own panel
    facet_wrap(~ Diet, scales = "fixed") +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black", face = "bold", angle = 0),
        axis.text.y = element_text(size = 16, color = "black", face = "bold"),
        legend.title = element_blank(),
        legend.position = "right",
        strip.text = element_text(size = 16, face = "bold")
    )

# Output 8*5


library(openxlsx)

wb <- createWorkbook()

addWorksheet(wb, "Fig 1d Weight individual")
writeData(wb, "Fig 1d Weight individual", weight_change_melted2)

addWorksheet(wb, "Fig 1d Weight group")
writeData(wb, "Fig 1d Weight group", means_time_treatment)

#saveWorkbook(wb, "Source Data Fig 1d Weight.xlsx", overwrite = TRUE)



####################################################################################
# Waist circumference change
####################################################################################

setwd("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/")

# Reading the metadata
library(readr)
library(dplyr)
library(reshape2)

metadata <- read_tsv(file = "metadata2.tsv", col_names = TRUE, show_col_types = FALSE)
head(metadata)
#  SubjectID sno_time Time   Age Weight      BMI  WC Gender Richness Diet Treatment Percentage
#1         4      4_0    0 69.75   90.5 29.21617 116      0      214    3      aFMT   31.70732
#2         4     4_14   14 70.92   85.9 27.73115 105      0      276    3      aFMT   31.70732
#3         4     4_18   18 71.25   85.9 27.73115 108      0      242    3      aFMT   31.70732
#4         4      4_6    6 70.25   82.6 26.66581 106      0      209    3      aFMT   31.70732
#5        13     13_0    0 50.93   86.7 34.73001 107      1      128    2   Placebo   57.14286
#6        13    13_14   14 52.10   84.1 33.68851 107      1      132    2   Placebo   57.14286

# Calculate weight change (weight at time 6 - weight at time 0)
metadata$WC_Baseline <- 0
metadata <- metadata %>%
  group_by(SubjectID) %>%
  mutate(#WC_Change0_6 = WC[Time == 6] - WC[Time == 0], 
        WC_Change6_14 = WC[Time == 14] - WC[Time == 6],
        WC_Change14_18 = WC[Time == 18] - WC[Time == 14]) %>%
  ungroup()
  
# Ensure 'Diet' is treated as a factor for proper faceting 
metadata$Diet <- factor(metadata$Diet, levels = c(1, 2, 3), labels = c("HDG", "MED", "GreenMED")) 

# Display the updated data frame
head(metadata)
# A tibble: 6 × 14
#  sno_time SubjectID  Time   Age Weight   BMI    WC Gender Richness Diet     Treatment WC_Baseline WC_Change6_14 WC_Change14_18
#  <chr>        <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl>  <dbl>    <dbl> <fct>    <chr>           <dbl>         <dbl>          <dbl>
#1 4_0              4     0  69.8   90.5  29.2   116      0      214 GreenMED aFMT                0            -1              3
#2 4_14             4    14  70.9   85.9  27.7   105      0      276 GreenMED aFMT                0            -1              3
#3 4_18             4    18  71.2   85.9  27.7   108      0      242 GreenMED aFMT                0            -1              3
#4 4_6              4     6  70.2   82.6  26.7   106      0      209 GreenMED aFMT                0            -1              3
#5 13_0            13     0  50.9   86.7  34.7   107      1      128 MED      Placebo             0             7             -8
#6 13_14           13    14  52.1   84.1  33.7   107      1      132 MED      Placebo             0             7             -8

#WC_change <- metadata[, c("SubjectID","Time", "Diet","Treatment", "WC_Baseline", "WC_Change0_6", "WC_Change6_14", "WC_Change14_18")]
WC_change <- metadata[, c("SubjectID","Time", "Diet","Treatment", "WC_Baseline", "WC_Change6_14", "WC_Change14_18")]

WC_change_melted <- melt(WC_change, id=c("SubjectID","Time", "Diet","Treatment"))
head(WC_change_melted)
#  SubjectID Time     Diet Treatment    variable value
#1         4    0 GreenMED      aFMT WC_Baseline     0
#2         4   14 GreenMED      aFMT WC_Baseline     0
#3         4   18 GreenMED      aFMT WC_Baseline     0
#4         4    6 GreenMED      aFMT WC_Baseline     0
#5        13    0      MED   Placebo WC_Baseline     0
#6        13   14      MED   Placebo WC_Baseline     0


# Remove outliers SubjectID 225
WC_change_melted2 <- WC_change_melted %>% filter(SubjectID != 225)

#data_phase1 <- WC_change_melted %>% filter(variable == "WC_Baseline"| variable == "WC_Change0_6")
data_phase2 <- WC_change_melted2 %>% filter(variable == "WC_Baseline"| variable == "WC_Change6_14")
data_phase3 <- WC_change_melted2 %>% filter(variable == "WC_Change6_14"| variable == "WC_Change14_18")


library(dplyr)

WC_means_time_treatment <- WC_change_melted2 %>%
  # First, remove rows with missing Time values
  filter(!is.na(variable) & variable != "NA") %>%
  group_by(variable, Diet, Treatment) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value   = sd(value, na.rm = TRUE),
    n          = n(),
    .groups    = "drop"
  )

# Display the rows
WC_means_time_treatment



# Merge the individual-level and mean-level figures together 

library(ggplot2)
ggplot() +
    # Add a reference horizontal dashed line at y = 0
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 1) +
    # Connect the mean points for each Treatment group over Time
    geom_line(data = WC_means_time_treatment,
            aes(x = variable, y = mean_value, group = Treatment, color = Treatment),
            alpha = 1, size = 1.5) +
    # Plot the mean points
    geom_point(data = WC_means_time_treatment,
             aes(x = variable, y = mean_value, group = Treatment, color = Treatment),
             size = 3) +
    
    # Phase 2: Time 6 and 14 -- lines/points colored by Treatment
    geom_line(data = data_phase2,
              aes(x = variable, y = value, group = SubjectID, color = Treatment),
              alpha = 0.1, size = 1) +
    geom_point(data = data_phase2,
               aes(x = variable, y = value, group = SubjectID, color = Treatment),
               alpha = 0.1, size = 1) +
    # Phase 3: Time 14 and later -- lines/points colored by Treatment
    geom_line(data = data_phase3,
              aes(x = variable, y = value, group = SubjectID, color = Treatment ),
              alpha = 0.1,  size = 1) +
    geom_point(data = data_phase3,
               aes(x = variable, y = value, group = SubjectID, color = Treatment),
               alpha = 0.1, size = 1) +
    scale_x_discrete(labels = c("6", "14", "18") ) +
    scale_color_manual(values = c("HDG" = "#D0495B", "MED" = "#4DBBD6", "GreenMED" = "#049B86", 
                                  "Placebo" = "#7796A7", "aFMT" = "#B9772B" )) +
    labs(x = "Month",y = "Waist Circumference Change (cm)") +
    # Facet by Diet so that each diet group appears in its own panel
    facet_wrap(~ Diet, scales = "fixed") +
    
    theme(
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black", face = "bold", angle = 0),
        axis.text.y = element_text(size = 16, color = "black", face = "bold"),
        legend.title = element_blank(),
        legend.position = "right",
        strip.text = element_text(size = 16, face = "bold")
    )
    
# Output 8*5  


addWorksheet(wb, "Fig 1d Waist Circum individual")
writeData(wb, "Fig 1d Waist Circum individual", weight_change_melted2)

addWorksheet(wb, "Fig 1d Waist Circum group")
writeData(wb, "Fig 1d Waist Circum group", means_time_treatment)

#saveWorkbook(wb, "Source Data Fig 1d Weight.xlsx", overwrite = TRUE)




####################################################################################
# Plasma Fasting Insulin Change
####################################################################################
setwd("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/")

# Reading the metadata
library(readr)
library(dplyr)
library(reshape2)

metadata_biomarkers <- read_tsv(file = "../metadata_biomarkers.tsv", col_names = TRUE, show_col_types = FALSE)
head(metadata_biomarkers)
# A tibble: 6 × 39
#    SubjectID  Time Gender  Diet Treatment   Age Weight   BMI Pulse    WC   DBP   SBP   DEX Glucose Insulin Cholesterol  HDLc
#  <dbl> <dbl>  <dbl> <dbl> <chr>     <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>   <dbl>       <dbl> <dbl>
#1     4     0      0     3 aFMT       69.8   90.5  29.2    NA   116  70.5  158    105    99.1   11.7         170.  38.2
#2    13     0      1     2 Placebo    50.9   86.7  34.7    76   107  75.5  120.   110    91.0   11.3         193.  43.6
#3    15     0      0     3 Placebo    46.8   99.6  35.7    79   111  94.5  147     94    99.5   18.4         212.  40.5
#4    17     0      0     3 Placebo    45.9  106.   35.0    73   117  94    146.   115   114.    28.6         229.  36.3
#5    19     0      0     2 Placebo    55.9   93    28.1    61   111  71    123    120   109.     8.60        156.  40.2
#6    22     0      0     3 Placebo    50.4   80.2  29.5    75   105  77    112    101    94.1    9.47        191.  42.1
# ℹ 22 more variables: LDLc <dbl>, Triglycerides <dbl>, ALT <dbl>, AST <dbl>, BilT <dbl>, ALKP <dbl>, Iron <dbl>,
#   Ferritine <dbl>, Transferrin <dbl>, TransferrinSat <dbl>, FolicA <dbl>, VitB12 <dbl>, CRP <dbl>, IL6 <dbl>,
#   FFA <dbl>, TSH <dbl>, TnT <dbl>, VitD3 <dbl>, Leptin <dbl>, Adiponectin <dbl>, FetuinA <dbl>, Chemerin <dbl>

# Calculate insulin change (weight at time 6 - weight at time 0)
metadata_biomarkers$Insulin_Baseline <- 0
metadata_biomarkers <- metadata_biomarkers %>%
  group_by(SubjectID) %>%
  mutate(#Insulin_Change0_6 = Insulin[Time == 6] - Insulin[Time == 0], 
        Insulin_Change6_14 = Insulin[Time == 14] - Insulin[Time == 6],
        Insulin_Change14_18 = Insulin[Time == 18] - Insulin[Time == 14]) %>%
  ungroup()
  
# Ensure 'Diet' is treated as a factor for proper faceting 
metadata_biomarkers$Diet <- factor(metadata_biomarkers$Diet, levels = c(1, 2, 3), labels = c("HDG", "MED", "GreenMED")) 
# Convert 'Time' to a factor to treat time points as discrete categories 
#metadata$Time <- factor(metadata$Time, levels = c(0, 6, 14, 18), labels = c("0", "6", "14", "18"))

Insulin_change <- metadata_biomarkers[, c("SubjectID","Time", "Diet","Treatment", "Insulin_Baseline", "Insulin_Change6_14", "Insulin_Change14_18")]

Insulin_change_melted <- melt(Insulin_change, id=c("SubjectID","Time", "Diet","Treatment"))
head(Insulin_change_melted)
#  SubjectID Time     Diet Treatment         variable value
#1         4    0 GreenMED      aFMT Insulin_Baseline     0
#2        13    0      MED   Placebo Insulin_Baseline     0
#3        15    0 GreenMED   Placebo Insulin_Baseline     0
#4        17    0 GreenMED   Placebo Insulin_Baseline     0
#5        19    0      MED   Placebo Insulin_Baseline     0
#6        22    0 GreenMED   Placebo Insulin_Baseline     0


# Remove outliers SubjectID 225
#Insulin_change_melted2 <- Insulin_change_melted %>% filter(SubjectID != 225)

library(dplyr)
Insulin_means_time_treatment <- Insulin_change_melted %>%
  # First, remove rows with missing Time values
  filter(!is.na(variable) & variable != "NA") %>%
  group_by(variable, Diet, Treatment) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value   = sd(value, na.rm = TRUE),
    n          = n(),
    .groups    = "drop"
  )

# Display the rows
Insulin_means_time_treatment

#data_phase1 <- Insulin_change_melted %>% filter(variable == "Insulin_Baseline"| variable == "Insulin_Change0_6")
data_phase2 <- Insulin_change_melted %>% filter(variable == "Insulin_Baseline"| variable == "Insulin_Change6_14")
data_phase3 <- Insulin_change_melted %>% filter(variable == "Insulin_Change6_14"| variable == "Insulin_Change14_18")

# Merge the individual-level and mean-level figures together 
library(ggplot2)
ggplot() +
    # Add a reference horizontal dashed line at y = 0
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 1) +
    # Connect the mean points for each Treatment group over Time
    geom_line(data = Insulin_means_time_treatment,
            aes(x = variable, y = mean_value, group = Treatment, color = Treatment),
            alpha = 1, size = 1.5) +
    # Plot the mean points
    geom_point(data = Insulin_means_time_treatment,
             aes(x = variable, y = mean_value, group = Treatment, color = Treatment),
             size = 3) +
    
    # Phase 2: Time 6 and 14 -- lines/points colored by Treatment
    geom_line(data = data_phase2,
              aes(x = variable, y = value, group = SubjectID, color = Treatment),
              alpha = 0.1, size = 1) +
    geom_point(data = data_phase2,
               aes(x = variable, y = value, group = SubjectID, color = Treatment),
               alpha = 0.1, size = 1) +
    # Phase 3: Time 14 and later -- lines/points colored by Treatment
    geom_line(data = data_phase3,
              aes(x = variable, y = value, group = SubjectID, color = Treatment ),
              alpha = 0.1,  size = 1) +
    geom_point(data = data_phase3,
               aes(x = variable, y = value, group = SubjectID, color = Treatment),
               alpha = 0.1, size = 1) +
    scale_x_discrete(labels = c("6", "14", "18") ) +
    scale_color_manual(values = c("HDG" = "#D0495B", "MED" = "#4DBBD6", "GreenMED" = "#049B86", 
                                  "Placebo" = "#7796A7", "aFMT" = "#B9772B" )) +
    labs(x = "Month",y = "Insulin Change (pmol/l)") +
    # Facet by Diet so that each diet group appears in its own panel
    facet_wrap(~ Diet, scales = "fixed") +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black", face = "bold", angle = 0),
        axis.text.y = element_text(size = 16, color = "black", face = "bold"),
        legend.title = element_blank(),
        legend.position = "right",
        strip.text = element_text(size = 16, face = "bold")
    )
    
# Output 8*5  


addWorksheet(wb, "Fig 1d Insulin individual")
writeData(wb, "Fig 1d Insulin individual", weight_change_melted2)

addWorksheet(wb, "Fig 1d Insulin group")
writeData(wb, "Fig 1d Insulin group", means_time_treatment)

saveWorkbook(wb, "Source Data Fig 1d.xlsx", overwrite = TRUE)


