setwd("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/")

# Reading the metadata
library(readr)
library(dplyr)
library(reshape2)

md <- read_tsv(file = "metadata2.tsv", col_names = TRUE, show_col_types = FALSE)
metadata <- read_tsv(file = "metadata_directplus.tsv", col_names = TRUE, show_col_types = FALSE)

head(md)
# A tibble: 6 × 11
#  sno_time SubjectID  Time   Age Weight   BMI    WC Gender Richness  Diet Treatment
#  <chr>        <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl>  <dbl>    <dbl> <dbl> <chr>    
#1 4_0              4     0  69.8   90.5  29.2   116      0      214     3 aFMT     
#2 4_14             4    14  70.9   85.9  27.7   105      0      276     3 aFMT     
#3 4_18             4    18  71.2   85.9  27.7   108      0      242     3 aFMT     
#4 4_6              4     6  70.2   82.6  26.7   106      0      209     3 aFMT     
#5 13_0            13     0  50.9   86.7  34.7   107      1      128     2 Placebo  
#6 13_14           13    14  52.1   84.1  33.7   107      1      132     2 Placebo  
head(metadata)
# A tibble: 6 × 10
#    sno  time   age weight   BMI pulse    WC   DBP   SBP   DEX
#  <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#1   169     0  74.9   88.1  29.8    59   112  86.5  172     97
#2   169     6  75.4   84.1  28.4    54   104  76.5  154.   108
#3   169    18  76.4   83.3  28.2    53   101  83    155     NA
#4   263     0  48.9  108.   34.8    78   118 100.   161     98
#5   263     6  49.4   94.9  30.6    58   107  66.5  113    101
#6   263    18  50.4  105.   33.8    76   111  90.5  132    101

# Calculate weight change (weight at time 6 - weight at time 0)
md$weight_Baseline <- 0
md <- md %>%
  group_by(SubjectID) %>%
  mutate(WC_Change0_6 = Weight[Time == 6] - Weight[Time == 0]) %>%
  ungroup()
  
metadata$weight_Baseline <- 0
metadata <- metadata %>%
  group_by(sno) %>%
  mutate(WC_Change0_6 = weight[time == 6] - weight[time == 0]) %>%
  ungroup()
  
md_afmt <- md[, c("SubjectID","WC_Change0_6")] %>% distinct()
metadata_wc <- metadata[,c("sno","WC_Change0_6")] %>% distinct()
names(metadata_wc)[1] <- "SubjectID"

md_afmt$df <- "aFMT+Placebo\n(90 subjects)"
metadata_wc$df <- "DIRECT-PLUS\n(294 subjects)"
merged_data <- rbind(md_afmt, metadata_wc)

library(ggplot2)
ggplot(merged_data, aes(x = df, y = WC_Change0_6), alpha = 0.5) +
    # Add a reference horizontal dashed line at y = 0
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 1) +
    # Draw the violin plots
    geom_violin(aes(color = df), size=1, alpha = 0.5) +
    # Overlay jittered points to show individual data
    geom_point(aes(color = df), 
               position = position_jitter(width = 0.1), 
               size = 1, 
               alpha = 0.5) +
    scale_color_manual(values = c("DIRECT-PLUS\n(294 subjects)" = "#779", "aFMT+Placebo\n(90 subjects)" = "#F2AE40")) +
    # Add a horizontal crossbar at the median for each group
    stat_summary(fun.data = function(x) {
        data.frame(y = median(x), ymin = median(x), ymax = median(x))
    },
    geom = "crossbar", width = 0.3,  color = "black") +
    labs(title = "Weight Change by Diets (kg, month 0-6)", x = "", y = "Month 0-6") +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        #panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 16, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(hjust = 2, size = 16, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black", face = "bold"),
        axis.text.y = element_text(size = 20, color = "black", face = "bold"),
        legend.title = element_blank(),
        legend.position = "none"
        #strip.text = element_text(size = 16, face = "bold")
    )
    
# Output 5*5    

merged_data_out <- merged_data %>%
    mutate(
        WC_Change0_6 = round(WC_Change0_6, 2),
        df = gsub("\n", " ", df)
    )

write.table(
    merged_data_out,
    "Source Data Fig 1c.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)


# Output source data to files
library(openxlsx)

wb <- createWorkbook()

addWorksheet(wb, "Fig 1c aFMT vs DIRECT-PLUS")
writeData(wb, "Fig 1c aFMT vs DIRECT-PLUS", merged_data_out)

saveWorkbook(wb, "Source Data Fig 1c.xlsx", overwrite = TRUE)
