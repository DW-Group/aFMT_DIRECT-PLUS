###########################################################################
# Identify metabolites with higher Decrease-Persistence in aFMT
###########################################################################


setwd("~/mydata/Project_aFMT/Metabolites/fecal")
# Load necessary library 
library(dplyr) 
library(tidyr) # Ensure tidyr is loaded for pivot_longer
library(pheatmap)
library(oddsratio)
library(stringr)
library(tibble)
library(broom)

# Read the dataset (assuming it is in a TSV format with one row per person)

metadata <- read.table("../metadata.tsv", header=T, sep="\t", check.names = FALSE)
metadata <- metadata %>% filter(Time != 18)

######################### Co-variables Processing Start #####################
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

# Categorize Age into quartiles with only data at month 6
time_6_df <- metadata %>% filter(Time == 6)
time_6_df$Age_Quartile <- cut(
  time_6_df$Age,
  breaks = quantile(time_6_df$Age, probs = seq(0, 1, 0.25), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("1", "2", "3", "4"))
metadata <- merge(metadata, time_6_df[,c("SubjectID", "Age_Quartile")], by = "SubjectID")
metadata$Age_Quartile <- as.factor(metadata$Age_Quartile)

# Categorize Richness into quartiles with only data at month 6
time_6_df <- metadata %>% filter(Time == 6)
time_6_df$Richness_Quartile <- cut(
  time_6_df$Richness,
  breaks = quantile(time_6_df$Richness, probs = seq(0, 1, 0.25), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("1", "2", "3", "4"))
metadata <- merge(metadata, time_6_df[,c("SubjectID", "Richness_Quartile")], by = "SubjectID")
metadata$Richness_Quartile <- as.factor(metadata$Richness_Quartile)

head(metadata)
#metadata <- metadata %>% select("SubjectID", "Diet", "Treatment", "Age_Quartile", "Gender", "WC_change_Quartile", "Weight_change_Quartile", "Richness_Quartile") %>% distinct()
dim(metadata)
#[1] 270  17

metadata$Diet <- as.factor(metadata$Diet)
metadata$Gender <- as.factor(metadata$Gender)
metadata$Treatment <- as.factor(metadata$Treatment)

######################### Co-variables Processing END #####################


##### Prepare metadata: select relevant columns and ensure SubjectID is available
#metabolite_metadata2 <- metadata[, c("SubjectID", "Diet", "Treatment")] %>% distinct()
metabolite_data <- read.table("Batch_normalizedData_90p.tsv", header=T, sep="\t", check.names = FALSE)
metabolite_data[1:4, 1:15]
metabolite_data[is.na(metabolite_data)] <- 0

metabolite_data <- metabolite_data %>% mutate(sno_time = paste0(SubjectID, "_", Month))
merged_data <- merge(metadata, metabolite_data[, -c(2,4,5,7)], by = "sno_time")
merged_data$Diet <- as.factor(merged_data$Diet)
merged_data[1:14, 1:25]



# Function to determine quartile assignment
# ---- helpers ----
assign_quartile <- function(value, q1, q2, q3, q4) {
  if (is.na(value)) {
    NA_character_
  } else if (value <= q1) {
    "Q1"
  } else if (value <= q2) {
    "Q2"
  } else if (value <= q3) {
    "Q3"
  } else {
    "Q4"
  }
}

library(dplyr)
library(tidyr)
library(rlang)


library(dplyr)
library(tidyr)

perform_metabolite_transition_analysis <- function(data, metabolite_cols, metadata = NULL) {
  results <- data.frame()

  for (metabolite in metabolite_cols) {

    # Wide: A/B/D columns for this metabolite per subject
    temp_data <- data %>%
      select(SubjectID, Diet, Treatment, GROUP_NUMBER, !!sym(metabolite)) %>%
      distinct() %>%
      pivot_wider(
        names_from  = GROUP_NUMBER,
        values_from = !!sym(metabolite),
        values_fill = NA
      )

    # Ensure we have A, B, D
    if (!all(c("A","B","D") %in% names(temp_data))) {
      message("Skipping ", metabolite, " (missing A/B/D)")
      next
    }

    # Compute quartiles globally for this metabolite
    vals <- suppressWarnings(as.numeric(data[[metabolite]]))
    qs <- stats::quantile(vals, probs = c(0.25, 0.5, 0.75, 1), na.rm = TRUE)
    q1 <- qs[[1]]; q2 <- qs[[2]]; q3 <- qs[[3]]; q4 <- qs[[4]]

    temp_data <- temp_data %>%
      mutate(
        Quartile_A = vapply(A, assign_quartile, character(1), q1, q2, q3, q4),
        Quartile_B = vapply(B, assign_quartile, character(1), q1, q2, q3, q4),
        Quartile_D = vapply(D, assign_quartile, character(1), q1, q2, q3, q4),

        Persistently_Loss = ifelse(
          (Quartile_A == "Q4" & Quartile_B %in% c("Q1","Q2") & Quartile_D %in% c("Q1","Q2")) | 
          (Quartile_A == "Q3" & Quartile_B == "Q1" & Quartile_D == "Q1"), 
          1L, 0L)
      )

    # Clean table (one row per subject)
    temp_clean <- temp_data %>%
      select(SubjectID, Diet, Treatment, Persistently_Loss) %>%
      tidyr::drop_na() %>%
      distinct()

    ##########################################################################################
    # ---------- Unadjusted OR (Overall aFMT vs Placebo) ----------
    counts_overall <- temp_clean %>%
      group_by(Treatment) %>%
      summarise(
        loss    = sum(Persistently_Loss == 1L, na.rm = TRUE),
        no_loss = sum(Persistently_Loss == 0L, na.rm = TRUE),
        .groups = "drop"
      )

    OR_unadjusted <- NA_real_
    OR_unadj_ci_low <- NA_real_
    OR_unadj_ci_high <- NA_real_
    OR_unadj_p <- NA_real_

    if (all(c("aFMT","Placebo") %in% counts_overall$Treatment)) {
      a <- counts_overall$loss[counts_overall$Treatment == "aFMT"]
      b <- counts_overall$no_loss[counts_overall$Treatment == "aFMT"]
      c_ <- counts_overall$loss[counts_overall$Treatment == "Placebo"]
      d <- counts_overall$no_loss[counts_overall$Treatment == "Placebo"]

      mat <- matrix(c(a, b, c_, d), nrow = 2, byrow = TRUE,
                    dimnames = list(Treatment = c("aFMT","Placebo"),
                                    Outcome = c("Loss","NoLoss")))
      ft <- tryCatch(stats::fisher.test(mat), error = function(e) NULL)
      if (!is.null(ft)) {
        OR_unadjusted   <- unname(ft$estimate)
        OR_unadj_ci_low <- unname(ft$conf.int[1])
        OR_unadj_ci_high<- unname(ft$conf.int[2])
        OR_unadj_p      <- ft$p.value
      }
    }
    
    ##########################################################################################
    # ---------- Optional adjusted GLM (guarded) ----------
    odds_ratio <- NA_real_; ci_low <- NA_real_; ci_high <- NA_real_; p_val <- NA_real_
    
    if (!is.null(metadata)) {
      # de-duplicate on key
      md_use <- metadata %>% distinct(SubjectID, .keep_all = TRUE)
      td_use <- temp_clean %>% distinct(SubjectID, .keep_all = TRUE) %>% select(SubjectID,Persistently_Loss)

      temp_data2 <- left_join(md_use, td_use, by = "SubjectID")

      if (all(c("Persistently_Loss","Treatment","Diet") %in% names(temp_data2))) {
        temp_data2 <- temp_data2 %>%
          mutate(Treatment = factor(Treatment, levels = c("Placebo","aFMT")))

        vars <- c("Persistently_Loss","Treatment","Diet","Age_Quartile","Gender",
                  "WC_change_Quartile","Weight_change_Quartile","Richness_Quartile")

        temp_cc <- temp_data2 %>% filter(complete.cases(across(any_of(vars))))

        ok_arms <- length(unique(temp_cc$Treatment)) == 2
        ok_events <- sum(temp_cc$Persistently_Loss == 1L, na.rm = TRUE) >= 2

        if (ok_arms && ok_events) {
          model <- tryCatch(
            glm(Persistently_Loss ~ Treatment + Diet + Age_Quartile + Gender +
                  WC_change_Quartile + Weight_change_Quartile + Richness_Quartile,
                data = temp_cc, family = binomial()),
            error = function(e) NULL
          )
          if (!is.null(model) && "TreatmentaFMT" %in% names(coef(model))) {
            beta <- coef(model)["TreatmentaFMT"]
            odds_ratio <- unname(exp(beta))
            ci_vals <- tryCatch(
              unname(exp(confint.default(model, "TreatmentaFMT"))),
              error = function(e) c(NA_real_, NA_real_)
            )
            ci_low  <- ci_vals[1]
            ci_high <- ci_vals[2]
            p_val   <- summary(model)$coefficients["TreatmentaFMT","Pr(>|z|)"]
          }
        }
      }
    }

    # ---------- Summaries ----------
    summary_tables <- list(
      overall_all = temp_clean %>%
        summarise(n = n(), Persistently_Loss = sum(Persistently_Loss, na.rm = TRUE), .groups = "drop") %>%
        mutate(Diet = "Overall", Treatment = "Overall"),

      overall_diet = temp_clean %>%
        group_by(Diet) %>%
        summarise(n = n(), Persistently_Loss = sum(Persistently_Loss, na.rm = TRUE), .groups = "drop") %>%
        mutate(Treatment = "Overall"),

      overall_treatment = temp_clean %>%
        group_by(Treatment) %>%
        summarise(n = n(), Persistently_Loss = sum(Persistently_Loss, na.rm = TRUE), .groups = "drop") %>%
        mutate(Diet = "Overall"),

      combined = temp_clean %>%
        group_by(Diet, Treatment) %>%
        summarise(n = n(), Persistently_Loss = sum(Persistently_Loss, na.rm = TRUE), .groups = "drop")
    )

    for (nm in names(summary_tables)) {
      summary_tables[[nm]] <- summary_tables[[nm]] %>%
        mutate(
          Metabolite = metabolite,
          Persistently_Loss_p = ifelse(n > 0, round(Persistently_Loss * 100 / n, 2), NA_real_),
          # attach ORs (same for every row of this metabolite)
          odds_ratio   = odds_ratio,
          ci_low = ci_low,
          ci_high= ci_high,
          p_val = p_val
          #OR_adjusted     = odds_ratio,
          #OR_adj_ci_low   = ci_low,
          #OR_adj_ci_high  = ci_high,
          #OR_adj_p        = p_val
        ) %>%
        select(Metabolite, Diet, Treatment, n, Persistently_Loss, Persistently_Loss_p,
               odds_ratio, ci_low, ci_high, p_val)
               #,OR_adjusted, OR_adj_ci_low, OR_adj_ci_high, OR_adj_p)
    }

    final_summary <- bind_rows(summary_tables) %>%
      mutate(Diet = as.character(Diet), Treatment = as.character(Treatment))

    results <- bind_rows(results, final_summary)
    cat("Processed metabolite:", metabolite, "\n")
  }

  results
}


# Specify your metabolite columns (update with the actual names from your data)
metabolite_cols <- colnames(merged_data)[22:ncol(merged_data)]
metabolite_transition_results <- perform_metabolite_transition_analysis(merged_data, metabolite_cols, metadata = metadata)
print(metabolite_transition_results)


# View the resulting summary table
metabolite_transition_results2 <- metabolite_transition_results %>% select(Metabolite, odds_ratio) %>% distinct() %>% filter(odds_ratio > 2)
print(metabolite_transition_results2)
dim(metabolite_transition_results2)
#[1] 514   2

# There were 287 metabolites persistently decreased!


########################################################################################################
# assign fecal chemical real names
########################################################################################################

library(dplyr)
library(readr)
fecal_chem_annotation <- read_tsv("../Fecal_chem_annotation.tsv",
  na = c("", "NA", "N/A"),
  locale = locale(encoding = "UTF-8"),
  guess_max = 10000
)

fecal_chem_annotation_clean <- fecal_chem_annotation[, c("CHEM_ID", "SUPER_PATHWAY", "SUB_PATHWAY", "CHEMICAL_NAME")]
names(fecal_chem_annotation_clean)[1] <- "Metabolite"
head(fecal_chem_annotation_clean)
# A tibble: 6 × 4
#  Metabolite SUPER_PATHWAY          SUB_PATHWAY                            CHEMICAL_NAME              
#       <dbl> <chr>                  <chr>                                  <chr>                      
#1         30 Lipid                  Mevalonate Metabolism                  mevalonate                 
#2         35 Amino Acid             Glutamate Metabolism                   S-1-pyrroline-5-carboxylate
#3         49 Amino Acid             Polyamine Metabolism                   putrescine                 
#4         50 Amino Acid             Polyamine Metabolism                   spermidine                 
#5         54 Nucleotide             Purine Metabolism, Adenine containing  1-methyladenine            
#6         55 Cofactors and Vitamins Nicotinate and Nicotinamide Metabolism 1-methylnicotinamide       

metabolite_transition_results2 <- metabolite_transition_results %>% filter(!is.na(odds_ratio)) %>% distinct()
dim(metabolite_transition_results2)/12
#[1] 1569.0000000   0.8333333

fecal_chem_annotation_clean$Metabolite <- as.factor(fecal_chem_annotation_clean$Metabolite)
fecal_chem_annotation_merged <- left_join(metabolite_transition_results2, fecal_chem_annotation_clean, by="Metabolite")
dim(fecal_chem_annotation_merged)
#[1] 18828   13

head(fecal_chem_annotation_merged, 20)
#   Metabolite    Diet Treatment  n Persistently_Loss Persistently_Loss_p odds_ratio ci_low ci_high       p_val
#1          30 Overall   Overall 77                11               14.29  0.1647243   0.04    0.53 0.004551804
#2          30       1   Overall 14                 0                0.00  0.1647243   0.04    0.53 0.004551804
#3          30       2   Overall 29                 8               27.59  0.1647243   0.04    0.53 0.004551804
#4          30       3   Overall 34                 3                8.82  0.1647243   0.04    0.53 0.004551804
#5          30 Overall      aFMT 37                 3                8.11  0.1647243   0.04    0.53 0.004551804
#6          30 Overall   Placebo 40                 8               20.00  0.1647243   0.04    0.53 0.004551804
#7          30       1      aFMT  7                 0                0.00  0.1647243   0.04    0.53 0.004551804
#8          30       1   Placebo  7                 0                0.00  0.1647243   0.04    0.53 0.004551804

#   SUPER_PATHWAY           SUB_PATHWAY               CHEMICAL_NAME
#1          Lipid Mevalonate Metabolism                  mevalonate
#2          Lipid Mevalonate Metabolism                  mevalonate
#3          Lipid Mevalonate Metabolism                  mevalonate
#4          Lipid Mevalonate Metabolism                  mevalonate
#5          Lipid Mevalonate Metabolism                  mevalonate
#6          Lipid Mevalonate Metabolism                  mevalonate
#7          Lipid Mevalonate Metabolism                  mevalonate
#8          Lipid Mevalonate Metabolism                  mevalonate


# Write the results to an output file (tab-delimited)
per_metabolite_transition_results <- fecal_chem_annotation_merged
#output_file <- "~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/6.ER_per_species_results/per_metabolite_decrease_persistence_results.tsv"
#write.table(per_metabolite_transition_results, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)


library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "Decreased-Persist Metabolites")
writeData(wb, "Decreased-Persist Metabolites", per_metabolite_transition_results)

saveWorkbook(wb, "Supplementary Table 5.xlsx", overwrite = TRUE)


########################################################################################################
# Reorder your data frame by the column Persistently_Gain_p
library(dplyr)

# remove (unselect) metabolites where, in diet = 3, the aFMT group has a lower “Persistently_Gain_p” than the placebo group
fecal_chem_annotation_merged_filt <- fecal_chem_annotation_merged %>% filter(CHEMICAL_NAME %in% fecal_chem_annotation_merged$CHEMICAL_NAME) %>% 
  group_by(Metabolite) %>%
  filter(
    if_all(1, ~ {
      p_aFMT <- Persistently_Loss_p[Diet == 3 & Treatment == "aFMT"]
      p_placebo <- Persistently_Loss_p[Diet == 3 & Treatment == "Placebo"]
      length(p_aFMT) == 0 || length(p_placebo) == 0 || p_aFMT > p_placebo
    })
  ) %>%
  ungroup()

fecal_chem_annotation_merged2 <- fecal_chem_annotation_merged_filt  %>% filter(Diet=="Overall" & Treatment=="Overall")  %>% filter(Persistently_Loss > 2 & odds_ratio > 1 & odds_ratio <=15 & ci_high < 250  & SUPER_PATHWAY != "NA") %>% distinct()
dim(fecal_chem_annotation_merged2)
#[1] 151  13

# Reorder metabolites names
fecal_chem_annotation_merged2$CHEMICAL_NAME <- factor(fecal_chem_annotation_merged2$CHEMICAL_NAME,
  levels = fecal_chem_annotation_merged2$CHEMICAL_NAME[order(fecal_chem_annotation_merged2$odds_ratio, decreasing = F)])


############################################################################
# Plot metabolites with best decrease-persistence, ranked by odds ratio
############################################################################

# Select metabolites with bete coefficients > 0 (odds ratios > 1)
# Filter out rows whose CHEMICAL_NAME begins with "X-"

head(fecal_chem_annotation_merged, 12)

# Add a new column by linking Treatment and Diet
fecal_chem_annotation_merged3 <- fecal_chem_annotation_merged %>% mutate(Treatment_diet = paste0(Treatment, "_", Diet))

# Remove rows only for Diets
fecal_chem_annotation_merged_filtered <- fecal_chem_annotation_merged3 %>% filter(!Treatment_diet %in% c("Overall_1", "Overall_2", "Overall_3")) %>% filter(!grepl("^X-", CHEMICAL_NAME)) %>% filter(CHEMICAL_NAME %in% fecal_chem_annotation_merged2$CHEMICAL_NAME) 

dim(fecal_chem_annotation_merged_filtered)
#[1] 1341  14
dim(fecal_chem_annotation_merged_filtered)/9
#[1] 149.000000  1.555556


#To remove metabolites where Persistently_Loss_p in Diet = 2 & Treatment = Placebo is greater than that in Diet = 2 & Treatment = aFMT, you can do the following:
# Step 1: Filter rows to only those from Diet 2
diet2 <- fecal_chem_annotation_merged_filtered %>% filter(Diet == 2)
# Step 2: Get Persistently_Gain_p for each Metabolite in both Treatment groups
diet2_wide <- diet2 %>%
  select(Metabolite, Treatment, Persistently_Loss_p) %>%
  pivot_wider(names_from = Treatment, values_from = Persistently_Loss_p)
# Step 3: Identify metabolites to KEEP
metabolites_to_keep <- diet2_wide %>%
  filter(Placebo <= aFMT | is.na(Placebo) | is.na(aFMT)) %>%
  pull(Metabolite)
# Step 4: Filter the original dataset
filtered_data <- fecal_chem_annotation_merged_filtered %>%
  filter(Metabolite %in% metabolites_to_keep)
dim(filtered_data)/9
#[1] 213.000000  1.555556
head(filtered_data)


# Select the top 30 metabolites relevant to weight change
library(dplyr)
library(stringr)
library(forcats)   # just forcats
#filtered_data <- filtered_data %>% mutate(CHEMICAL_NAME = str_replace_all(CHEMICAL_NAME, c("alpha" = "α", "beta" = "β", "gamma" = "γ")) )

filtered_data2 <- filtered_data %>% select("Metabolite", "odds_ratio", "CHEMICAL_NAME",SUPER_PATHWAY,SUB_PATHWAY) %>% distinct() %>%
    filter(odds_ratio > 2)  %>% distinct()
# Reorder metabolites names
filtered_data2$CHEMICAL_NAME <- factor(filtered_data2$CHEMICAL_NAME,
  levels = filtered_data2$CHEMICAL_NAME[order(filtered_data2$odds_ratio, decreasing = F)])

## ensure the greek letters being outputted correctly by library them
library(showtext)
showtext_auto()

############################################################################
# Remove some metabolites
metabolites_rm <- c("100008932","100006051","100003588","100001026","100020204","799","100004509","100000700","100002082","100002113","100006430","361","100015609","100003926")
filtered_data2 <- filtered_data2 %>% filter(!Metabolite %in% metabolites_rm)

p <- ggplot(filtered_data2, aes(x = CHEMICAL_NAME, 
    # Originial odds ratio
    y = odds_ratio, ymin = odds_ratio, ymax = odds_ratio)) +
    # Log tranformation odds ratio
    #y = log(odds_ratio), ymin = log(ci_low), ymax = log(ci_high))) +
    geom_pointrange(size = 0.5, color = "#1F78B4") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    coord_flip() +
    labs(x = "",
         y = "Odds Ratio of participants\naFMT vs. Placebo",
         title = "Decrease-Persistence") +
    #scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, by = 2) ) +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 12, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 12, color = "black", face = "bold", angle=0),
        axis.text.y = element_text(size = 10, color = "black"),
        text = element_text(family = "Arial"),
        legend.title = element_blank(),  # Remove legend title
        legend.position = "right",
        strip.text = element_text(size = 16, face = "bold")  # Modify facet title style
    )
# Output portrait 5*6
p
#ggsave("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/Figure3f. Decrease−Persistence metabolites filtered odds ratio.pdf", plot = p, width = 5.5, height = 6)


library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "Decreased-Persist Metabolites")
writeData(wb, "Decreased-Persist Metabolites", filtered_data2)

saveWorkbook(wb, "Source Data Fig 3d.xlsx", overwrite = TRUE)




############################################################################
# Create a stacked bar chart for SUPER_PATHWAY
############################################################################

library(ggplot2)
library(cowplot)

# Create the stacked bar plot for SUPER_PATHWAY
 #"Cofactors and Vitamins" "Lipid"                  "Energy"                 "Nucleotide"               "Carbohydrate"           "Xenobiotics"      
#distinct_colors <- c("Amino Acid"="firebrick", "#377EB8", "#4DAF4A", "#fabed4", "#984EA3", "#FF7F00", "#ffd8b1", "#A65628")
distinct_colors <- c("Amino Acid"="firebrick", "Cofactors and Vitamins"="#377EB8", "Lipid"="#4DAF4A", "Peptide" = "orange", "Nucleotide"="#984EA3", "Xenobiotics"= "#ffd8b1", "Energy" = "#A65628", "Carbohydrate" = "skyblue")

# Here we assume each CHEMICAL_NAME may have multiple rows with different SUPER_PATHWAY values.

p <- ggplot(filtered_data2, aes(x = CHEMICAL_NAME, fill = SUPER_PATHWAY)) +
    geom_bar(alpha=0.8) +
    coord_flip() +
    labs(fill = "Super pathway\n(SP)") +
    #theme_minimal() +
    scale_fill_manual(values = distinct_colors) +
    labs(title="SP")+
    theme(
        panel.background = element_blank(),
        #axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "white", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 11, face = "bold"),
        axis.title.x = element_text(size = 0, face = "bold"),
        axis.title.y = element_text(size = 0, face = "bold"),
        axis.text.x = element_text(size = 0, color = "white"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.ticks.x = element_blank(),
        legend.title = element_text( face = "bold"),
        legend.position = "right",
        strip.text = element_text(size = 16, face = "bold")  # Modify facet title style
    )
    
# Output portrait 5.2*6 for top 30
p
#ggsave("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/Figure3f. Decrease−Persistence metabolites filtered super pathway.pdf", plot = p, width = 5.3, height = 6)




############################################################################
# Assign metabolites to more interpretable functional categories
############################################################################

library(readr)
library(dplyr)
library(ggplot2)

setwd('~/Downloads')
decrease_persistence_df <- read_tsv("decrease_persistence_functional_categories.tsv", show_col_types = FALSE)
head(decrease_persistence_df)
# A tibble: 6 × 2
#  metabolite                             functional_category                       
#  <fct>                                  <chr>                                     
#1 3-hydroxypropanoate                               Amino-acid and nucleotide metabolism      
#2 2,3-dihydroxy-5-methylthio-4-pentenoate (DMTPA)*  Amino-acid and nucleotide metabolism      
#3 1,2-dilinoleoyl-digalactosylglycerol (18:2/18:2)* Structural and membrane lipid metabolism  
#4 margaroylcarnitine (C17)*                         Fatty-acid oxidation and energy metabolism
#5 oleoyl-oleoyl-glycerol (18:1/18:1) [2]*           Structural and membrane lipid metabolism  
#6 uridine 5'-monophosphate (UMP)                    Amino-acid and nucleotide metabolism      

# Capitalize only metabolites that start with a letter, while leaving those starting with numbers/symbols unchanged.
decrease_persistence_df <- decrease_persistence_df %>%
  mutate(
    metabolite = if_else(
      str_detect(metabolite, "^[A-Za-z]"),
      str_replace(metabolite, "^[A-Za-z]", toupper(str_sub(metabolite, 1, 1))),
      metabolite
    ),
    # Keep the figure order from the TSV (top to bottom)
    metabolite = factor(metabolite, levels = rev(metabolite))
  )

#cat_colors <- c(
#  "Gut-brain and tryptophan-related signaling" = "#B27AC9",
#  "Plant-derived and antioxidant metabolites" = "#DDAA3B",
#  "Fatty-acid oxidation and energy metabolism" = "#74B66A",
#  "Structural and membrane lipid metabolism" = "#79B98B",
#  "Bile acid metabolism" = "#F1B77B",
#  "Amino-acid, peptide, and nucleotide metabolism" = "#C65A5A",
#  "Dietary xenobiotics and food additives" = "#7AA6D6"
#)

cat_colors <- c(
  "Gut-brain and tryptophan-related signaling" = "purple2",# "#6C5CE7",
  "Plant-derived and antioxidant metabolites" = "#00A896",
  "Fatty-acid oxidation and energy metabolism" = "#028090",
  "Structural and membrane lipid metabolism" = "#A8DADC",
  "Bile acid metabolism" = "#F4A261",
  "Amino-acid, peptide, and nucleotide metabolism" = "#E76F51",
  "Dietary xenobiotics and food additives" = "#457B9D"
)


## ensure the greek letters being outputted correctly by library them
library(showtext)
showtext_auto()

decrease_persistence_df2 <- decrease_persistence_df
decrease_persistence_df2$metabolite <- as.character(decrease_persistence_df2$metabolite)
decrease_persistence_df2[which(decrease_persistence_df2$metabolite=="2,3-dihydroxy-5-methylthio-4-pentenoate (DMTPA)*"),"metabolite"] = "2,3-dihy...-5-methy...-4-pent... (DMTPA)*"
decrease_persistence_df2[which(decrease_persistence_df2$metabolite=="1,2-dilinoleoyl-digalactosylglycerol (18:2/18:2)*"),"metabolite"] = "1,2-dilino...-digala...glycerol (18:2/18:2)*"
decrease_persistence_df2[which(decrease_persistence_df2$metabolite=="1-linoleoyl-2-linolenoyl-digalactosylglycerol (18:2/18:3)*"),"metabolite"] = "1-lino...-2-lino...-di...glycerol (18:2/18:3)*"


decrease_persistence_df2 <- decrease_persistence_df2 %>%
  mutate(metabolite = factor(metabolite, levels = rev(metabolite)))
  
p <- ggplot(decrease_persistence_df2, aes(x = metabolite, fill = functional_category)) +
    geom_bar(alpha=1) +
    coord_flip() +
    labs(fill = "Functional Categories (FC)") +
    #theme_minimal() +
    scale_fill_manual(values = cat_colors) +
    labs(title="FC") +
    theme(
        panel.background = element_blank(),
        text = element_text(family = "Arial"),
        #axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "white", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 12, face = "bold"),
        axis.title.x = element_text(size = 0, face = "bold"),
        axis.title.y = element_text(size = 0, face = "bold"),
        axis.text.x = element_text(size = 0, color = "white"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        legend.position = "right",
        strip.text = element_text(size = 16, face = "bold")  # Modify facet title style
    )
p
# Save to PDF
# ggsave("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/Figure3f. Decrease−Persistence metabolites filtered functional category.pdf", plot = p, width = 8.4, height = 6)




############################################################################
# Bubble plot: Rows = CHEMICAL_NAME, Columns = Treatment_diet Group, 
# Bubble size and color indicate Persistently_Gain_p 
############################################################################

fecal_chem_annotation_merged_bubble <- filtered_data
# Define the column order
fecal_chem_annotation_merged_bubble$Treatment_diet <- factor(fecal_chem_annotation_merged_bubble$Treatment_diet, 
    levels=c("Overall_Overall", "aFMT_Overall", "Placebo_Overall",
            "aFMT_1", "Placebo_1", "aFMT_2", "Placebo_2", "aFMT_3", "Placebo_3"))

# Select Metabolites in the same order of the forest plot
fecal_chem_annotation_merged_bubble$CHEMICAL_NAME <- factor(
  fecal_chem_annotation_merged_bubble$CHEMICAL_NAME,
  levels = unique(fecal_chem_annotation_merged_bubble$CHEMICAL_NAME[order(fecal_chem_annotation_merged_bubble$odds_ratio, decreasing = FALSE)])
)

# Subset the data frame to only include rows with CHEMICAL_NAME in the vector
persistent_decrease_subset_df <- fecal_chem_annotation_merged_bubble %>%
  filter(Metabolite %in% filtered_data2$Metabolite)
  
# Range of the persistent loss percentage
range(persistent_decrease_subset_df$Persistently_Loss_p)
#[1]  0.00 31.25

# Define a named vector of 6 new color codes
new_colors <- c(
    "Overall_Overall"      = "#8B0000",  # Deep DarkRed
    "aFMT_Overall"         = "#D95F02",  # Strong gold-orange
    "Placebo_Overall"      = "#4F7CAC",  # Bold blue-gray
    "aFMT_1"   = "#E6AB02",  # Vivid bright orange/gold
    "aFMT_2"   = "#008B8B",  # Strong dark cyan
    "aFMT_3"   = "#66A61E",  # Rich teal
    "Placebo_1"= "#E6AB02",  # (Same as aFMT_diet1)
    "Placebo_2"= "#008B8B",  # (Same as aFMT_diet2)
    "Placebo_3"= "#66A61E"   # (Same as aFMT_diet3)
)


# Bubble plot: Rows = CHEMICAL_NAME, Columns = Treatment_diet Group, Bubble size and color indicate Persistently_Gain_p 
# bubble_plot
p <- ggplot(persistent_decrease_subset_df, aes(x = Treatment_diet, y = CHEMICAL_NAME)) +
    geom_point(aes(size = Persistently_Loss_p, color = Treatment_diet, alpha=Persistently_Loss_p/31) ) +
    scale_color_manual(values = new_colors) +
    labs(title = "Within-group % of participants",
         x = "",
         y = "",
         size = "Within-group\n% of participants",
         color = "Engrafted groups") +
    guides(color = "none", alpha = "none") +
    scale_x_discrete(
        labels = c("Overall", "aFMT", "Placebo","aFMT_HDG","Placebo_HDG",
                   "aFMT_MED","Placebo_MED", "aFMT_GreenMED", "Placebo_GreenMED")) + # Set custom x-axis labels
    theme(
        panel.background = element_blank(),
        axis.line = element_line(color = "grey80"),
        panel.border = element_rect(color = "grey80", fill = NA, size = 1),
        plot.title = element_text(hjust = 0, size = 12, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 9, color = "black", vjust = 1,hjust = 0, angle=-20),
        axis.text.y = element_text(size = 8, color = "black", face = "bold"),
        legend.title = element_text( face = "bold"),
        legend.position = "right",
        strip.text = element_text(size = 16, face = "bold")  # Modify facet title style
    )
p
# Output portrait 7*5
# ggsave("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/Figure3f. Decrease−Persistence metabolites filtered bubble plot.pdf", plot = p, width = 7.5, height = 6)


