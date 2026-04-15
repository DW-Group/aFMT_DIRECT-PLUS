
# Figure da. Define metabolites change categories

# Schematic of the definition of metabolite abundance change
library(tibble)
library(dplyr)

make_row <- function(q0, q6, q14, cat)
  tibble(Month = c(0,6,14),
         Quartile = c(q0, q6, q14),
         Category = cat)

dat <- dplyr::bind_rows(
  make_row("Q1","Q3", NA,  "Increase"),
  make_row("Q1","Q4", NA,  "Increase"),
  make_row("Q2","Q4", NA,  "Increase"),

  make_row("Q4","Q1", NA,  "Decrease"),
  make_row("Q4","Q2", NA,  "Decrease"),
  make_row("Q3","Q1", NA,  "Decrease"),

  make_row("Q1","Q3","Q3","Increase-Persistence"),
  make_row("Q1","Q3","Q4","Increase-Persistence"),
  make_row("Q1","Q4","Q3","Increase-Persistence"),
  make_row("Q1","Q4","Q4","Increase-Persistence"),
  make_row("Q2","Q4","Q4","Increase-Persistence"),

  make_row("Q3","Q1","Q1","Decrease-Persistence"),
  make_row("Q4","Q1","Q1","Decrease-Persistence"),
  make_row("Q4","Q1","Q2","Decrease-Persistence"),
  make_row("Q4","Q2","Q1","Decrease-Persistence"),
  make_row("Q4","Q2","Q2","Decrease-Persistence")
) %>%
  mutate(
    Month    = factor(Month, levels = c(0,6,14), ordered = TRUE),
    Quartile = factor(Quartile, levels = c("Q1","Q2","Q3","Q4"), ordered = TRUE),
    Category = factor(Category,
      levels = c("Decrease","Increase","Decrease-Persistence","Increase-Persistence"))
  )

# quick sanity check
#dplyr::count(dat, Category, Month, Quartile)
dat
# A tibble: 48 × 3
#   Month Quartile Category
#   <ord> <ord>    <fct>   
# 1 0     Q1       Increase
# 2 6     Q3       Increase
# 3 14    NA       Increase

library(dplyr)
library(tidyr)
library(ggplot2)

# 1) add a line id per Category (starts a new id at Month==0)
dat2 <- dat %>%
  mutate(Month = factor(Month, levels = c(0,6,14), ordered = TRUE)) %>%
  group_by(Category) %>%
  mutate(line_id = cumsum(Month == levels(Month)[1])) %>%
  ungroup()

# 2) plot-ready: drop NA quartiles, map to numeric, add small per-line x-offset
dat_plot <- dat2 %>%
  filter(!is.na(Quartile)) %>%
  mutate(
    Month_num = as.numeric(as.character(Month)),
    q_num     = as.integer(Quartile)
  ) %>%
  group_by(Category) %>%
  mutate(
    n_lines = n_distinct(line_id),
    offset  = (line_id - (n_lines + 1)/2) * 0.5,   # tweak 0.12 if you want more/less spread
    x       = Month_num + offset
  ) %>%
  ungroup()

p <- ggplot(dat_plot, aes(x = x, y = q_num, group = interaction(Category, line_id),color=line_id)) +
    geom_line(alpha = 0.9, linewidth = 0.6) +
    geom_point(size = 1) +
    facet_wrap(~Category, nrow = 1) +
    geom_hline(yintercept = c(1,2,3,4), linetype = "dashed", linewidth = 0.3, color = "grey80") +
    geom_vline(xintercept = c(0,6,14), linetype = 3, linewidth = 0.3, color = "grey60") +
    scale_x_continuous(breaks = c(0,6,14), labels = c("0","6","14"), name = "Month") +
    scale_y_continuous(breaks = 1:4, labels = paste0("Q", 1:4), limits = c(1,4), name = "Quartile") +
    #theme_classic(base_size = 12) +
    theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "grey", fill = NA, size = 1),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold",vjust = 2),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        plot.title = element_text(hjust = 0, size = 14, color = "black", face = "bold"),
        strip.text = element_text(size = 12, face = "bold") ) +
    labs(title = "Patterns of metabolite abundance change" ) 

ggsave("~/mydata/Project_aFMT/metaphlan_v411/metaphlan_all_species/Figure3a. Category of metabolite abundance change.pdf", plot = p, width = 10, height = 2.5)
