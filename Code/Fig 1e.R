
#############################################################################
# Part 1: PCoA analysis of Species
#############################################################################


setwd("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/")

# Reading the metadata
library(readr)
library(dplyr)
library(reshape2)
library(vegan) 
library(ggplot2)

abundance <- read_tsv(file = "merged_abundance_table_filtered.tsv", col_names = TRUE, show_col_types = FALSE)
#metadata <- read_tsv(file = "metadata_directplus.tsv", col_names = TRUE, show_col_types = FALSE)
abundance[1:5, 1:6]
# A tibble: 5 × 6
#  SampleID Treatment Time  SubjectID  Diet t__SGB10068
#  <chr>    <chr>     <chr>     <dbl> <dbl>       <dbl>
#1 A004     aFMT      A             4     3     0      
#2 A013     Placebo   A            13     2     0      
#3 A015     Placebo   A            15     3     0      
#4 A017     Placebo   A            17     3     0      
#5 A019     Placebo   A            19     2     0.00031
abundance <- abundance[abundance$Time!="A",]
abundance$Time <- factor(abundance$Time, levels = c("A", "B", "D", "C"), 
                    labels = c("0", "6", "14", "18")) 
                
## 1) Species matrix (rows = samples, cols = species)
species_data <- abundance %>%
  select(starts_with("t__")) %>%           # all species/features
  as.data.frame()
rownames(species_data) <- abundance$SampleID

# (Optional) drop all-zero samples to avoid NaNs in Bray–Curtis
keep <- rowSums(species_data) > 0
species_data <- species_data[keep, ]


## 2) Metadata aligned to species_data order
meta <- abundance %>%
    select(SampleID, SubjectID, Treatment, Time, Diet) %>%
    filter(SampleID %in% rownames(species_data)) %>%
    arrange(match(SampleID, rownames(species_data))) %>%
    mutate(
        SubjectID = factor(SubjectID),
        Treatment = factor(Treatment, levels = c("Placebo","aFMT")),
        Time      = factor(Time, levels = c("6","14","18")),       # keep A/B/D; relabel if desired
        Diet      = factor(Diet, levels = c(1,2,3), labels = c("HDG","MED","GreenMED"))
    )
stopifnot(all(meta$SampleID == rownames(species_data)))


## 3) Bray–Curtis dissimilarities
dist_mat <- vegdist(species_data, method = "bray")
# -------------------------------
# Perform PCoA (Classical MDS)
# -------------------------------
# k = 2 returns the first two axes; eig = TRUE returns eigenvalues.
pcoa_res <- cmdscale(dist_mat, k = 2, eig = TRUE)

# Calculate the percentage of explained variance.
# (If negative eigenvalues are present, sum only the positive ones.)
ev <- pcoa_res$eig
pos_ev_sum <- sum(ev[ev > 0])
perc1 <- round(100 * ev[1] / pos_ev_sum, 1)
perc2 <- round(100 * ev[2] / pos_ev_sum, 1)


meta2 <- meta %>% mutate(Treatment_Time = paste0(Treatment,Time))
## 4) PERMANOVA (marginal tests) with SubjectID as permutation stratum
set.seed(123)
adon_main <- adonis2(
  dist_mat ~ Treatment_Time,
  data         = meta2,
  permutations = 999,
  strata       = meta$SubjectID,   # respects repeated measures
  by           = "margin"          # report marginal (Type II-like) tests
)
print(adon_main)   # report R^2 and P for Treatment_Time
#Permutation test for adonis under reduced model
#Marginal effects of terms
#Blocks:  strata 
#Permutation: free
#Number of permutations: 999

#adonis2(formula = dist_mat ~ Treatment_Time, data = meta2, permutations = 999, by = "margin", strata = meta$SubjectID)
#                Df SumOfSqs      R2     F Pr(>F)    
#Treatment_Time   5    1.159 0.01977 0.944  0.001 ***
#Residual       234   57.445 0.98023                 
#Total          239   58.604 1.00000                 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##########################################################################################

library(dplyr)
library(ggplot2)
library(ggExtra)
library(cowplot)

pcoa_df <- abundance %>% 
  select(SampleID, Treatment, Time) %>%
  mutate(PC1 = pcoa_res$points[, 1],
         PC2 = pcoa_res$points[, 2],
         Treatment = as.factor(Treatment))  # Convert Diet to a factor for discrete coloring

pcoa_df2 <- pcoa_df[pcoa_df$Time != "0", ]

# 1. Reorder factor levels as needed
pcoa_df2 <- pcoa_df2 %>%
  mutate(Treatment_Time = paste(Treatment, Time, sep = "_M")) %>%
  mutate(Treatment_Time = factor(Treatment_Time, 
                                 levels = c("aFMT_M6", "Placebo_M6",
                                            "aFMT_M14", "Placebo_M14",
                                            "aFMT_M18", "Placebo_M18"))) %>% arrange(Treatment_Time)

# 2. Main scatter plot
p_main <- ggplot(pcoa_df2, aes(x = PC1, y = PC2, color = Treatment_Time)) +
  geom_point(size = 2, alpha = 1) + 
  labs(#title = "Species abundance",
       x = paste0("PCoA 1 (", perc1, "%)"),
       y = paste0("Species abundance\nPCoA 2 (", perc2, "%)")) +
  theme_minimal() +
  scale_color_manual(values = c("Placebo_M6" = "#A2B7C9", 
                                  "Placebo_M14" = "#7796A7", 
                                  "Placebo_M18" = "#4B6584", 
                                  "aFMT_M6" = "#D8B07F", 
                                  "aFMT_M14" = "#B9772B",
                                  "aFMT_M18" = "#A64D20")
                                  ) + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0, size = 18, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    legend.title = element_blank(),
    legend.position = "bottom"
  )

# 3. Top box plot (for PC1)
p_top <- ggplot(pcoa_df2, aes(x = Treatment_Time, y = PC1, color = Treatment_Time, fill = Treatment_Time)) +
  geom_boxplot(outliers=F, alpha=0.7) +
  coord_flip() +
  scale_x_discrete(limits = levels(pcoa_df2$Treatment_Time)) +
  theme_minimal() +
  scale_color_manual(values = c("Placebo_M6" = "#A2B7C9", 
                                  "Placebo_M14" = "#7796A7", 
                                  "Placebo_M18" = "#4B6584", 
                                  "aFMT_M6" = "#D8B07F", 
                                  "aFMT_M14" = "#B9772B",
                                  "aFMT_M18" = "#A64D20")
                                  ) + 
    scale_fill_manual(values = c("Placebo_M6" = "#A2B7C9", 
                                  "Placebo_M14" = "#7796A7", 
                                  "Placebo_M18" = "#4B6584", 
                                  "aFMT_M6" = "#D8B07F", 
                                  "aFMT_M14" = "#B9772B",
                                  "aFMT_M18" = "#A64D20")
                                  ) + 
  theme(panel.grid.major.y = element_blank(), 
        #panel.grid.minor.x = element_blank(),
        #axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        legend.position = "none")

# 4. Right box plot (for PC2)
p_right <- ggplot(pcoa_df2, aes(x = Treatment_Time, y = PC2, color = Treatment_Time, fill = Treatment_Time)) +
  geom_boxplot(outliers=F, alpha=0.7) +
  #coord_flip() +
  scale_x_discrete(limits = levels(pcoa_df2$Treatment_Time)) +
  theme_minimal() +
  scale_color_manual(values = c("Placebo_M6" = "#A2B7C9", 
                                  "Placebo_M14" = "#7796A7", 
                                  "Placebo_M18" = "#4B6584", 
                                  "aFMT_M6" = "#D8B07F", 
                                  "aFMT_M14" = "#B9772B",
                                  "aFMT_M18" = "#A64D20")
                                  ) + 
  scale_fill_manual(values = c("Placebo_M6" = "#A2B7C9", 
                                  "Placebo_M14" = "#7796A7", 
                                  "Placebo_M18" = "#4B6584", 
                                  "aFMT_M6" = "#D8B07F", 
                                  "aFMT_M14" = "#B9772B",
                                  "aFMT_M18" = "#A64D20")
                                  ) + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        #axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

# 5. Combine using cowplot
combined <- plot_grid(
  plot_grid(p_top, NULL, ncol = 1, rel_heights = c(1, 0.05)),
  plot_grid(p_main, p_right, ncol = 2, rel_widths = c(1, 0.3)),
  nrow = 2, rel_heights = c(0.3, 1)
)
combined

# Output 8*8

#write.table(pcoa_df2, file = "Source Data Fig 1e Species.tsv", sep = '\t')

# Output source data to files
library(openxlsx)

wb <- createWorkbook()

addWorksheet(wb, "Fig 1e Species_data")
writeData(wb, "Fig 1e Species_data", pcoa_df2)

#saveWorkbook(wb, "Source Data Fig 1e.xlsx", overwrite = TRUE)




#############################################################################
# Part 2: PCoA analysis of Metabolites
#############################################################################

setwd("~/mydata/Project_aFMT/Metabolites/fecal/")

# Reading the metadata
library(readr)
library(dplyr)
library(reshape2)
library(vegan) 
library(ggplot2)
library(ggExtra) 
library(cowplot)

#metabolites_abundance <- read.table(file = "Batch_normalizedData_90p.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = F)
metabolites_abundance <- read.table(file = "~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/8.phenotype_prediction/Fecal_metabolites_Batch-norm_Imputed_Data_90p2.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = F)

metabolites_abundance <- metabolites_abundance %>%
  mutate(sno_time = paste(SubjectID, Month, sep = "_"))

metadata <- read_tsv(file = "../metadata.tsv", col_names = TRUE, show_col_types = FALSE)
metadata <-  metadata[,c("sno_time","SubjectID","Diet","Treatment","Time")] %>% distinct()
metabolites_abundance2 <- merge(metadata, metabolites_abundance, id="sno_time")

rownames(metabolites_abundance2)<- metabolites_abundance2$CLIENT_SAMPLE_ID
# Remove data from Month 0
metabolites_abundance2 <- metabolites_abundance2[metabolites_abundance2$Month != "0", ]
metabolites_abundance2[1:5, 1:11]

metabolites_data <- metabolites_abundance2[,-c(1:11)]
# Use Bray–Curtis distance
dist_mat <- vegdist(metabolites_data, method = "bray", na.rm = T)

# 3. Perform PCoA (Classical MDS)
# -------------------------------
# k = 2 returns the first two axes; eig = TRUE returns eigenvalues.
pcoa_res <- cmdscale(dist_mat, k = 2, eig = TRUE)

# Calculate the percentage of explained variance.
# (If negative eigenvalues are present, sum only the positive ones.)
ev <- pcoa_res$eig
pos_ev_sum <- sum(ev[ev > 0])
perc1 <- round(100 * ev[1] / pos_ev_sum, 1)
perc2 <- round(100 * ev[2] / pos_ev_sum, 1)

# -------------------------------
# 4. Combine PCoA results with metadata
# -------------------------------
pcoa_df <- metabolites_abundance2 %>% 
  select(Treatment, Month) %>%
  mutate(PC1 = pcoa_res$points[, 1],
         PC2 = pcoa_res$points[, 2],
         Treatment = as.factor(Treatment))  # Convert Diet to a factor for discrete coloring

# 5. Reorder factor levels as needed
pcoa_df <- pcoa_df %>%
  mutate(Treatment_Time = paste(Treatment, Month, sep = "_M")) %>%
  mutate(Treatment_Time = factor(Treatment_Time, 
                                 levels = c("aFMT_M6", "Placebo_M6",
                                            "aFMT_M14", "Placebo_M14",
                                            "aFMT_M18", "Placebo_M18")))

ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Treatment_Time)) +
  geom_point(size = 2, alpha = 1) + 
  labs(#title = "Species abundance",
       x = paste0("PCoA 1 (", perc1, "%)"),
       y = paste0("Metabolite level\nPCoA 2 (", perc2, "%)")) +
  theme_minimal() +
  scale_color_manual(values = c("Placebo_M6" = "#A2B7C9", 
                                  "Placebo_M14" = "#7796A7", 
                                  "Placebo_M18" = "#4B6584", 
                                  "aFMT_M6" = "#D8B07F", 
                                  "aFMT_M14" = "#B9772B",
                                  "aFMT_M18" = "#A64D20")
                                  ) 

# Remove outliers for which PC2 is larger than 0.2
pcoa_df2 <- pcoa_df[pcoa_df$PC2 > -0.2, ]

# There are some outliers, check their name
pcoa_outliers <- pcoa_df[pcoa_df$PC2 > 0.2, ]
unique(pcoa_outliers$Sample)
# [1] "BRIG-05698" "BRIG-05282" "BRIG-05834" "BRIG-05824" "BRIG-05500" "BRIG-05392" "BRIG-05384"


# 2. Main scatter plot
p_main <- ggplot(pcoa_df2, aes(x = PC1, y = PC2, color = Treatment_Time)) +
  geom_point(size = 2, alpha = 1) + 
  labs(title = "Metabolite level\n",
       x = paste0("PCoA 1 (", perc1, "%)"),
       y = paste0("PCoA 2 (", perc2, "%)")) +
  theme_minimal() +
  scale_color_manual(values = c("Placebo_M6" = "#A2B7C9", 
                                  "Placebo_M14" = "#7796A7", 
                                  "Placebo_M18" = "#4B6584", 
                                  "aFMT_M6" = "#D8B07F", 
                                  "aFMT_M14" = "#B9772B",
                                  "aFMT_M18" = "#A64D20")
                                  ) + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0, size = 18, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    legend.title = element_blank(),
    legend.position = "bottom"
  )

# 3. Top box plot (for PC1)
p_top <- ggplot(pcoa_df2, aes(x = Treatment_Time, y = PC1, color = Treatment_Time, fill = Treatment_Time)) +
  geom_boxplot(outliers=F, alpha=0.8) +
  coord_flip() +
  scale_x_discrete(limits = levels(pcoa_df2$Treatment_Time)) +
  theme_minimal() +
  scale_color_manual(values = c("Placebo_M6" = "#A2B7C9", 
                                  "Placebo_M14" = "#7796A7", 
                                  "Placebo_M18" = "#4B6584", 
                                  "aFMT_M6" = "#D8B07F", 
                                  "aFMT_M14" = "#B9772B",
                                  "aFMT_M18" = "#A64D20")
                                  ) + 
    scale_fill_manual(values = c("Placebo_M6" = "#A2B7C9", 
                                  "Placebo_M14" = "#7796A7", 
                                  "Placebo_M18" = "#4B6584", 
                                  "aFMT_M6" = "#D8B07F", 
                                  "aFMT_M14" = "#B9772B",
                                  "aFMT_M18" = "#A64D20")
                                  ) + 
  theme(panel.grid.major.y = element_blank(), 
        #panel.grid.minor.x = element_blank(),
        #axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        legend.position = "none")

# 4. Right box plot (for PC2)
p_right <- ggplot(pcoa_df2, aes(x = Treatment_Time, y = PC2, color = Treatment_Time, fill = Treatment_Time)) +
  geom_boxplot(outliers=F, alpha=0.8) +
  #coord_flip() +
  scale_x_discrete(limits = levels(pcoa_df2$Treatment_Time)) +
  theme_minimal() +
  scale_color_manual(values = c("Placebo_M6" = "#A2B7C9", 
                                  "Placebo_M14" = "#7796A7", 
                                  "Placebo_M18" = "#4B6584", 
                                  "aFMT_M6" = "#D8B07F", 
                                  "aFMT_M14" = "#B9772B",
                                  "aFMT_M18" = "#A64D20")
                                  ) + 
  scale_fill_manual(values = c("Placebo_M6" = "#A2B7C9", 
                                  "Placebo_M14" = "#7796A7", 
                                  "Placebo_M18" = "#4B6584", 
                                  "aFMT_M6" = "#D8B07F", 
                                  "aFMT_M14" = "#B9772B",
                                  "aFMT_M18" = "#A64D20")
                                  ) + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        #axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

# 5. Combine using cowplot
combined <- plot_grid(
  plot_grid(p_top, NULL, ncol = 1, rel_heights = c(1, 0.05)),
  plot_grid(p_main, p_right, ncol = 2, rel_widths = c(1, 0.3)),
  nrow = 2, rel_heights = c(0.3, 1)
)
combined

# Output 8*8

#write.table(pcoa_df2, file = "Source Data Fig 1e Metabolites.tsv", sep = '\t')

addWorksheet(wb, "Fig 1e Metabolite_data")
writeData(wb, "Fig 1e Metabolite_data", pcoa_df2)

saveWorkbook(wb, "Source Data Fig 1e.xlsx", overwrite = TRUE)

