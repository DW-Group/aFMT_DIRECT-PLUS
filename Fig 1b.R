
##################################################################################
# Density plots for Weight, Age, Waist circumference, and Gender at Month 6
##################################################################################

library(ggplot2)
library(dplyr)

setwd("~/mydata/Project_aFMT/Metabolites/fecal/")
metadata <- read_tsv(file = "../metadata.tsv", col_names = TRUE, show_col_types = FALSE)


###############################################
# Density plot for Weight at M6
###############################################

metadata_weight_6 <- metadata[metadata$Time=="6",c("SubjectID", "Weight", "Treatment")]
ggplot(metadata_weight_6, aes(x = Weight, fill=Treatment)) +
    #geom_density(fill = "#F2AE40", alpha = 0.5) +  # Smooth density curve
    geom_density( alpha = 0.4) +
    labs(title = "Weight", x = "", y = "") +
    scale_x_continuous(limits = c(0, 150)) +  # Set x-axis range from 0 to 100
    scale_fill_manual(values = c("Placebo" = "#7796A7", "aFMT" = "#B9772B")) +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        #panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 22, face = "bold"),
        axis.title.x = element_text(size = 22, face = "bold"),
        axis.title.y = element_text(vjust = 0.5, size = 22, face = "bold"),
        axis.text.x = element_text(size = 20, color = "black", face = "bold"),
        axis.text.y = element_text(size = 0, color = "black", face = "bold"),
        legend.title = element_blank(),
        legend.position = "right",
        strip.text = element_text(size = 26, face = "bold")
    )
# Output 6*3



###############################################
# Density plot for Waist circumference at M6
###############################################

metadata_wc_6 <- metadata[metadata$Time=="6",c("SubjectID", "WC", "Treatment")]
ggplot(metadata_wc_6, aes(x = WC, fill=Treatment)) +
    #geom_density(fill = "#F2AE40", alpha = 0.5) +  # Smooth density curve
    geom_density(alpha = 0.4) +
    labs(title = "Waist circumference", x = "", y = "") +
    scale_x_continuous(limits = c(50, 150)) +  # Set x-axis range from 0 to 100
    scale_fill_manual(values = c("Placebo" = "#7796A7", "aFMT" = "#B9772B")) +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        #panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 22, face = "bold"),
        axis.title.x = element_text(size = 22, face = "bold"),
        axis.title.y = element_text(vjust = 0.5, size = 22, face = "bold"),
        axis.text.x = element_text(size = 20, color = "black", face = "bold"),
        axis.text.y = element_text(size = 0, color = "black", face = "bold"),
        legend.title = element_blank(),
        legend.position = "right",
        strip.text = element_text(size = 26, face = "bold")
    )
# Output 6*3


###############################################
# Density plot for Age distribution at M6
###############################################

metadata_age_6 <- metadata[metadata$Time=="6",c("SubjectID", "Age", "Treatment")]
ggplot(metadata_age_6, aes(x = Age, fill=Treatment)) +
    #geom_density(fill = "#69b3a2", alpha = 0.5) +  # Smooth density curve
    geom_density(alpha = 0.4) +  # Smooth density curve
    scale_fill_manual(values = c("Placebo" = "#7796A7", "aFMT" = "#B9772B")) +
    labs(title = "Age", 
            x = "", 
            y = ""
            ) +
    scale_x_continuous(limits = c(0, 100)) +  # Set x-axis range from 0 to 100
    theme(
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        #panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 22, face = "bold"),
        axis.title.x = element_text(size = 22, face = "bold"),
        axis.title.y = element_text(vjust = 0.5, size = 22, face = "bold"),
        axis.text.x = element_text(size = 20, color = "black", face = "bold"),
        axis.text.y = element_text(size = 0, color = "black", face = "bold"),
        legend.title = element_blank(),
        legend.position = "right",
        strip.text = element_text(size = 26, face = "bold")
    )
# Output 6*3




###############################################
# Density plot for Gender at M6
###############################################

metadata_gender <- metadata %>% select("SubjectID", "Gender", "Treatment") %>% distinct()
library(ggplot2)
library(dplyr)

# Count the number of each gender
gender_counts <- metadata_gender %>%
  group_by(Treatment, Gender) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = (Count / sum(Count)) * 100)  # Convert to percentage

library(dplyr)
gender_counts <- gender_counts %>%
  mutate(Gender = ifelse(Gender == 0, "Male", "Female"))

# Create a stacked bar plot
ggplot(gender_counts, aes(x = Percentage, y = Treatment, fill = Gender)) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.7) +
  #labs(title = "Gender Distribution", x = "", y = "Number of People") +
  #scale_y_continuous(labels = scales::percent_format()) +  # Converts y-axis to percentage
  scale_fill_manual(values = c( "Female" = "#D95F02", "Male" = "#46826d")) + # Blue for Male, Red for Female
  #theme_minimal() +
  #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    #ggplot(metadata_gender, aes(x = as.factor(Gender))) +
    #geom_bar(fill = "#D95F02", alpha = 0.7) +  # Bar chart since Gender is categorical
    labs(title = "Gender (%)", x = "", y = "") +
    #scale_x_discrete(labels = c("0" = "Male", "1" = "Female")) + 
    #theme_minimal()
    theme(
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        #panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 22, face = "bold"),
        axis.title.x = element_text(size = 22, face = "bold"),
        axis.title.y = element_text(vjust = 0.5, size = 22, face = "bold"),
        axis.text.x = element_text(size = 20, color = "black", face = "bold"),
        axis.text.y = element_text(size = 20, color = "black", face = "bold"),
        legend.title = element_blank(),
        legend.position = "right",
        strip.text = element_text(size = 26, face = "bold")
    )
# Output 6*2
 
