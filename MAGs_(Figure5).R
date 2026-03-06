# MAGs analysis (Figure 5)


library(readxl)
library(dplyr)
library(tidyr)
library(rlang)
library(stringr)
library(purrr)
library(ggplot2)
library(forcats)
library(RColorBrewer)

infile <- "~/Documents/African_Microbiome_Group/Long_COVID_study/Masters_research/Results/GTDBTK/MAGs.xlsx"
top_n <- 30
pseudocount <- 1e-6

df <- readxl::read_excel(infile)
names(df) <- tolower(names(df))   

df <- df %>%
  rename(
    mag        = any_of("mag"),
    color      = any_of("color"),
    phylum_raw = any_of("phylum"),
    phylum2    = any_of("phylum2"),
    dataset    = any_of("group"),        
    coverage   = any_of("coverage"),
    domain     = any_of("domain"),
    class      = any_of("class"),
    order      = any_of("order"),
    family     = any_of("family"),
    genus      = any_of("genus"),
    species    = any_of("species")
  )

if ("phylum2" %in% names(df) && !all(is.na(df$phylum2))) {
  df <- df %>% mutate(phylum = phylum2)
} else {
  df <- df %>% mutate(phylum = phylum_raw)
}

if (!"mag" %in% names(df)) df$mag <- NA_character_

df <- df %>% mutate(coverage = as.numeric(coverage))

df <- df %>%
  mutate(
    species = if_else(is.na(species) | species == "", NA_character_, as.character(species)),
    genus   = if_else(is.na(genus)   | genus == "",   NA_character_, as.character(genus)),
    mag     = if_else(is.na(mag)     | mag == "",     NA_character_, as.character(mag)),
    species = case_when(
      !is.na(species) ~ species,
      is.na(species) & !is.na(genus) ~ genus,
      TRUE ~ mag
    )
  )

sample_candidates <- intersect(c("sample","sampleid","sample_id","library","run","sample_name"), names(df))
if (length(sample_candidates) > 0) {
  sample_col <- sample_candidates[1]
  message("Found sample column '", sample_col, "'. Normalizing per-sample to relative abundance.")
  df <- df %>%
    rename(sample = !!sym(sample_col)) %>%
    group_by(sample) %>%
    mutate(rel_abund = if_else(sum(coverage, na.rm = TRUE) > 0,
                               coverage / sum(coverage, na.rm = TRUE), 0)) %>%
    ungroup()
  value_col <- "rel_abund"
} else {
  message("No sample column found. Using coverage values directly for plotting (rel_abund := coverage).")
  df <- df %>% mutate(rel_abund = coverage)
  value_col <- "rel_abund"
}

req_cols <- c("mag","dataset","phylum","genus","species","coverage")
missing_cols <- setdiff(req_cols, names(df))
if (length(missing_cols) > 0) stop("Missing required columns after normalization: ", paste(missing_cols, collapse = ", "))

group_stats <- df %>%
  group_by(dataset, species) %>%
  summarise(
    mean_val = mean(.data[[value_col]], na.rm = TRUE),
    sd_val   = sd(.data[[value_col]], na.rm = TRUE),
    n        = sum(!is.na(.data[[value_col]])),
    se       = if_else(n > 0, sd_val / sqrt(n), NA_real_),
    .groups = "drop"
  )

is_placeholder <- function(s) {
  s <- as.character(s)
  pattern <- regex("(^CAG-|^UBA|^UMGS|^ER|\\ssp\\d+\\b|\\bsp\\d+\\b|\\bsp\\b\\s*\\d*)", ignore_case = TRUE)
  stringr::str_detect(s, pattern)
}
group_stats <- group_stats %>% mutate(is_placeholder = is_placeholder(species))

wide_stats <- group_stats %>%
  select(dataset, species, mean_val) %>%
  pivot_wider(names_from = dataset, values_from = mean_val, values_fill = 0)

ds_cols <- setdiff(names(wide_stats), "species")
if (length(ds_cols) < 2) stop("Less than two unique datasets found in 'dataset' column. Need two groups to compare.")
gA <- ds_cols[1]; gB <- ds_cols[2]
message("Comparing groups: '", gA, "' vs '", gB, "'")

wide_stats <- wide_stats %>%
  mutate(diff = .data[[gB]] - .data[[gA]],
         abs_diff = abs(diff),
         overall_mean = rowMeans(across(all_of(ds_cols)), na.rm = TRUE))

top_by_absdiff <- wide_stats %>% arrange(desc(abs_diff)) %>% slice_head(n = top_n) %>% pull(species)
top_by_overall  <- wide_stats %>% arrange(desc(overall_mean)) %>% slice_head(n = top_n) %>% pull(species)

plot1_df <- group_stats %>%
  filter(species %in% top_by_overall) %>%
  mutate(species = factor(species, levels = rev(top_by_overall)))   

species_list <- unique(df$species)
test_results <- map_dfr(species_list, function(sp) {
  sub <- df %>% filter(species == sp)
  valsA <- sub %>% filter(dataset == gA) %>% pull(.data[[value_col]])
  valsB <- sub %>% filter(dataset == gB) %>% pull(.data[[value_col]])
  nA <- sum(!is.na(valsA)); nB <- sum(!is.na(valsB))
  if (nA >= 2 && nB >= 2) {
    t <- wilcox.test(valsB, valsA, alternative = "two.sided", exact = FALSE)
    tibble(species = sp, nA = nA, nB = nB, p_value = t$p.value)
  } else {
    tibble(species = sp, nA = nA, nB = nB, p_value = NA_real_)
  }
})

test_results <- test_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         sig = case_when(
           is.na(p_adj) ~ NA_character_,
           p_adj < 0.001 ~ "***",
           p_adj < 0.01  ~ "**",
           p_adj < 0.05  ~ "*",
           TRUE ~ "ns"
         ))

plot1_ann <- group_stats %>% left_join(test_results, by = "species")


my_colors <- c("Control" = "#008f00", "Long COVID" = "#ff9300")
 

dodge <- position_dodge(width = 0.8)

p1 <- ggplot(plot1_df, aes(x = species, y = mean_val, fill = dataset)) +
  geom_col(position = dodge, width = 0.7) +
  coord_flip() +
  labs(x = "", y = ifelse(value_col == "rel_abund", "Mean relative abundance", "Mean Coverage"),
       title = paste0("")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.position = "top",
        plot.title = element_text(face = "bold")) +
  scale_fill_manual(values = my_colors, name = "Group")

annot_positions1 <- plot1_df %>%
  group_by(species) %>%
  summarise(ypos = max(mean_val, na.rm = TRUE), .groups = "drop") %>%
  left_join(test_results, by = "species") %>%
  mutate(ypos = if_else(is.na(ypos), 0, ypos * 1.05 + 1e-12)) %>%
  filter(!is.na(p_adj))

if (nrow(annot_positions1) > 0) {
  p1 <- p1 + geom_text(data = annot_positions1, aes(x = species, y = ypos, label = sig),
                       inherit.aes = FALSE, hjust = 0.5, size = 3)
}

print(p1)

plot2_df <- group_stats %>%
  filter(species %in% top_by_absdiff) %>%
  mutate(species = factor(species, levels = rev(top_by_absdiff)))

p2 <- ggplot(plot2_df, aes(x = species, y = mean_val, fill = dataset)) +
  geom_col(position = dodge, width = 0.7) +
  coord_flip() +
  labs(x = "", y = ifelse(value_col == "rel_abund", "Mean relative abundance", "Mean Coverage"),
       title = paste0("")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.position = "top",
        plot.title = element_text(face = "bold")) +
  scale_fill_manual(values = my_colors, name = "Group")

annot_positions2 <- plot2_df %>%
  group_by(species) %>%
  summarise(ypos = max(mean_val, na.rm = TRUE), .groups = "drop") %>%
  left_join(test_results, by = "species") %>%
  mutate(ypos = if_else(is.na(ypos), 0, ypos * 1.05 + 1e-12)) %>%
  filter(!is.na(p_adj))

if (nrow(annot_positions2) > 0) {
  p2 <- p2 + geom_text(data = annot_positions2, aes(x = species, y = ypos, label = sig),
                       inherit.aes = FALSE, hjust = 0.5, size = 3)
}

print(p2)

full_species <- group_stats %>% filter(!is_placeholder) %>% pull(species) %>% unique()
if (length(full_species) == 0) {
  message("No full (non-placeholder) species found for p3.")
  p3 <- NULL
} else {
  overall_rank_full <- group_stats %>%
    filter(!is_placeholder) %>%
    group_by(species) %>%
    summarise(overall_mean = mean(mean_val, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(overall_mean)) %>%
    slice_head(n = top_n) %>%
    pull(species)
  
  plot_p3_df <- group_stats %>%
    filter(species %in% overall_rank_full) %>%
    mutate(species = factor(species, levels = rev(overall_rank_full)))
  
  plot3_ann <- plot_p3_df %>% left_join(test_results, by = "species")
  
  p3 <- ggplot(plot_p3_df, aes(x = species, y = mean_val, fill = dataset)) +
    geom_col(position = dodge, width = 0.7) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_flip() +
    labs(x = "", y = ifelse(value_col == "rel_abund", "Mean relative abundance", "Mean Coverage"),
         title = paste0("")) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black"),
          legend.position = "top",
          plot.title = element_text(face = "bold")) +
    scale_fill_manual(values = my_colors, name = "Group")
  
  ann_pos3 <- plot3_ann %>%
    group_by(species) %>%
    summarise(ypos = max(mean_val, na.rm = TRUE), p_adj = first(p_adj), sig = first(sig), .groups = "drop") %>%
    mutate(ypos = if_else(is.na(ypos), 0, ypos * 1.05 + 1e-12)) %>%
    filter(!is.na(p_adj))
  
  if (nrow(ann_pos3) > 0) p3 <- p3 + geom_text(data = ann_pos3, aes(x = species, y = ypos, label = sig),
                                               inherit.aes = FALSE, size = 3)
  print(p3)
}

placeholder_species <- group_stats %>% filter(is_placeholder) %>% pull(species) %>% unique()
if (length(placeholder_species) == 0) {
  message("No placeholder species found; skipping p4.")
  p4 <- NULL
} else {
  overall_rank_placeholder <- group_stats %>%
    filter(is_placeholder) %>%
    group_by(species) %>%
    summarise(overall_mean = mean(mean_val, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(overall_mean)) %>%
    slice_head(n = top_n) %>%
    pull(species)
  
  plot_p4_df <- group_stats %>%
    filter(species %in% overall_rank_placeholder) %>%
    mutate(species = factor(species, levels = rev(overall_rank_placeholder)))
  
  plot4_ann <- plot_p4_df %>% left_join(test_results, by = "species")
  
  p4 <- ggplot(plot_p4_df, aes(x = species, y = mean_val, fill = dataset)) +
    geom_col(position = dodge, width = 0.7) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_flip() +
    labs(x = "", y = ifelse(value_col == "rel_abund", "Mean relative abundance", "Mean Coverage"),
         title = paste0("")) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black"),
          legend.position = "top",
          plot.title = element_text(face = "bold")) +
    scale_fill_manual(values = my_colors, name = "Group")
  
  ann_pos4 <- plot4_ann %>%
    group_by(species) %>%
    summarise(ypos = max(mean_val, na.rm = TRUE), p_adj = first(p_adj), sig = first(sig), .groups = "drop") %>%
    mutate(ypos = if_else(is.na(ypos), 0, ypos * 1.05 + 1e-12)) %>%
    filter(!is.na(p_adj))
  
  if (nrow(ann_pos4) > 0) p4 <- p4 + geom_text(data = ann_pos4, aes(x = species, y = ypos, label = sig),
                                               inherit.aes = FALSE, size = 3)
  print(p4)
}



# ---- libraries ----
library(readxl)
library(dplyr)
library(ggplot2)
library(forcats)
library(ggtext)   
library(tibble)

phylum_colors <- c(
  "Bacillota"         = "#ff7e79",
  "Bacteroidota"      = "#009e73",
  "Actinomycetota"    = "#332288",
  "Pseudomonadota"    = "#f0e442",
  "Bacillota_A"    = "#ffd579",
  "Bacillota_C"    = "#fffd78",
  "Verrucomicrobiota" = "#D55E00",
  "Fusobacteriota"    = "#0072B2",
  "Elusimicrobiota"   = "#999933",
  "Synergistota"      = "#117733",
  "Cyanobacteriota" = "#28d910",
  "Desulfobacterota"             = "#ffb8ff"
)

MAG_df <- read_excel("/Users/louisburger22/Documents/African_Microbiome_Group/Long_COVID_study/Masters_research/Results/GTDBTK/MAGs.xlsx")

MAG_df <- MAG_df %>%
  rename_with(~ .x, everything()) %>%  
  mutate(
    Phylum = as.character(Phylum),
    Group  = as.character(Group),
    Color  = as.character(Color)
  ) %>%
  filter(!is.na(Phylum) & Phylum != "")

counts_group <- MAG_df %>%
  group_by(Phylum, Group) %>%
  summarise(MAGs = n(), .groups = "drop")

counts_sum <- counts_group %>%
  group_by(Phylum) %>%
  summarise(total = sum(MAGs), .groups = "drop") %>%
  arrange(total)

phylum_levels <- counts_sum$Phylum

phylum_color_df <- tibble(Phylum = phylum_levels) %>%
  mutate(color = phylum_colors[Phylum]) %>%
  mutate(color = if_else(is.na(color) | color == "", "#BEBEBE", color))

counts_group <- counts_group %>% mutate(Phylum = factor(Phylum, levels = phylum_levels))
counts_sum   <- counts_sum    %>% mutate(Phylum = factor(Phylum, levels = phylum_levels))
phylum_color_df <- phylum_color_df %>% mutate(Phylum = factor(Phylum, levels = phylum_levels))

group_cols <- c("Long COVID" = "#ff9300", "Control" = "#008f00")

x_max <- max(counts_sum$total, na.rm = TRUE)
text_offset <- max(1, ceiling(x_max * 0.03))
left_space <- max(1, ceiling(x_max * 0.08))   

use_dash_instead_of_dot <- FALSE   

p_fixed <- ggplot(counts_group, aes(x = MAGs, y = Phylum, fill = Group)) +
  geom_col(position = "stack", color = "NA", size = 0.15) +
  geom_text(
    data = counts_sum,
    mapping = aes(x = total, y = Phylum, label = total),
    inherit.aes = FALSE,
    hjust = -0.05,
    size = 3.2
  ) +
  {
    if (!use_dash_instead_of_dot) {
      geom_point(
        data = phylum_color_df,
        mapping = aes(x = -left_space/2, y = Phylum),
        inherit.aes = FALSE,
        show.legend = FALSE,
        colour = phylum_color_df$color,
        size = 4.5,
        shape = 19
      )
    } else {
      dash_half_len <- left_space * 0.3
      geom_segment(
        data = phylum_color_df,
        mapping = aes(x = -left_space/2 - dash_half_len, xend = -left_space/2 + dash_half_len,
                      y = Phylum, yend = Phylum),
        inherit.aes = FALSE,
        show.legend = FALSE,
        colour = phylum_color_df$color,
        linewidth = 2
      )
    }
  } +
  labs(x = "Number of MAGs", y = "", title = "") +
  scale_fill_manual(values = group_cols, na.value = "grey70") +
  scale_x_continuous(expand = c(0, 0), limits = c(-left_space, x_max + text_offset * 4)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(margin = margin(r = 6)),   
    axis.title.y = element_text(margin = margin(r = 6), face = "bold"),
    axis.ticks.y = element_line(color = "black", linetype = "dashed", linewidth = 0.4),
    axis.ticks.x = element_line(color = "black", linetype = "dashed", linewidth = 0.4),
    panel.border = element_rect(fill = NA, color = "black", size = 0.6),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    legend.key = element_rect(fill = NA),
    plot.margin = margin(6, 40, 6, 6)
  )

print(p_fixed)

p1_legend <- p1 + theme(legend.position = "bottom", legend.direction = "horizontal")
shared_legend <- cowplot::get_legend(p1_legend)

p1_noleg <- p1 + theme(legend.position = "none")
p2_noleg <- p_fixed + theme(legend.position = "none")
p3_noleg <- p3 + theme(legend.position = "none")
p4_noleg <- p4 + theme(legend.position = "none")

top_pair <- plot_grid(p1_noleg, p2_noleg, labels = c("A", "B"), label_size = 14,
                      ncol = 2, rel_widths = c(0.48, 0.52), align = "hv")

legend_row <- ggdraw(shared_legend)

diversity_block <- plot_grid(
  top_pair,
  legend_row,
  ncol = 1,
  rel_heights = c(1, 0.12)   
)

abundance_block <- plot_grid(
  p3_noleg, p4_noleg,
  labels = c("C", "D"),
  label_size = 14,
  ncol = 2,
  rel_widths = c(0.48, 0.52),
  align = "hv"
)

final_plot2 <- plot_grid(
  diversity_block,
  abundance_block,
  ncol = 1,
  rel_heights = c(1, 1.2)    
)

print(final_plot2)

