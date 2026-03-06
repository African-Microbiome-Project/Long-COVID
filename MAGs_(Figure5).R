library(tidyverse)
library(ape)
library(readxl)
library(vegan)
library(cowplot)

lc_tree_file <- "/Users/louisburger22/Documents/African_Microbiome_Group/Long_COVID_study/Masters_research/Results/GTDBTK/Testing/LC_bac.tre"
control_tree_file <- "/Users/louisburger22/Documents/African_Microbiome_Group/Long_COVID_study/Masters_research/Results/GTDBTK/Testing/control_bac.tre"
combined_tree_file <- "/Users/louisburger22/Documents/African_Microbiome_Group/Long_COVID_study/Masters_research/Results/GTDBTK/Testing/combined_longcovid_control_tree.tre"

lc_summary_file <- path.expand("~/Documents/African_Microbiome_Group/Long_COVID_study/Masters_research/Results/GTDBTK/LC_gtdbtk.bac120.summary.tsv")
lc_arc_summary_file <- path.expand("~/Documents/African_Microbiome_Group/Long_COVID_study/Masters_research/Results/GTDBTK/Testing/gtdbtk.ar53.summary.tsv")
ctrl_summary_file <- path.expand("~/Documents/African_Microbiome_Group/Long_COVID_study/Masters_research/Results/GTDBTK/control_gtdbtk.bac120.summary.tsv")
out_meta_file <- path.expand("~/Documents/African_Microbiome_Group/Long_COVID_study/Masters_research/Results/GTDBTK/combined_MAG_metadata_for_iTOL.tsv")

mag_excel_file <- path.expand("~/Documents/African_Microbiome_Group/Long_COVID_study/Masters_research/Results/GTDBTK/MAGs.xlsx")

top_n <- 20
pseudocount <- 1e-6
normalize_per_sample <- TRUE

group_cols <- c("Control" = "#008f00", "Long_COVID" = "#ff9300")

base_theme <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black"),
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

alpha_theme <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black"),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

normalize_group <- function(x) {
  case_when(
    x %in% c("LongCOVID", "Long_COVID", "Long COVID") ~ "Long_COVID",
    x %in% c("Control", "control", "CTRL") ~ "Control",
    TRUE ~ as.character(x)
  )
}

clean_summary <- function(df, source_label) {
  if (!"user_genome" %in% names(df)) stop("No genome ID column found in summary file.")
  
  tax_candidates <- c("classification", "pplacer_taxonomy", "fastani_taxonomy", "closest_placement_taxonomy")
  tax_found <- intersect(tax_candidates, names(df))
  tax_col <- if (length(tax_found) > 0) tax_found[1] else NULL
  
  if (!"classification" %in% names(df)) {
    df$classification <- if (!is.null(tax_col)) df[[tax_col]] else NA_character_
  }
  
  df %>%
    mutate(across(everything(), as.character)) %>%
    select(any_of(c("user_genome", "classification", "fastani_ani", "fastani_af", "msa_percent", "translation_table", "warnings"))) %>%
    separate(
      classification,
      into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
      sep = ";",
      fill = "right",
      remove = FALSE
    ) %>%
    mutate(
      across(domain:species, ~ if_else(is.na(.) | . == "", NA_character_, str_remove(., "^[a-z]__"))),
      dataset = source_label
    )
}

is_placeholder <- function(x) {
  x <- as.character(x)
  str_detect(x, regex("(^CAG-|^UBA|^UMGS|^ER|\\ssp\\d+\\b|\\bsp\\d+\\b|\\bsp\\b\\s*\\d*)", ignore_case = TRUE))
}

make_annotation_df <- function(plot_df, test_results) {
  plot_df %>%
    group_by(Species) %>%
    summarise(ypos = max(mean_val, na.rm = TRUE), .groups = "drop") %>%
    left_join(test_results, by = "Species") %>%
    mutate(ypos = ypos * 1.05 + 1e-12) %>%
    filter(!is.na(p_adj))
}

tree_lc <- ape::read.tree(lc_tree_file)
tree_control <- ape::read.tree(control_tree_file)

prefixed_lc <- FALSE
prefixed_ctrl <- FALSE

if (any(tree_lc$tip.label %in% tree_control$tip.label)) {
  tree_lc$tip.label <- paste0("LC_", tree_lc$tip.label)
  tree_control$tip.label <- paste0("CTRL_", tree_control$tip.label)
  prefixed_lc <- TRUE
  prefixed_ctrl <- TRUE
}

root2 <- ape::read.tree(text = "(r1:0.0001,r2:0.0001);")
tmp_tree <- ape::bind.tree(root2, tree_lc, where = which(root2$tip.label == "r1"))
combined_tree <- ape::bind.tree(tmp_tree, tree_control, where = which(tmp_tree$tip.label == "r2"))
combined_tree <- ape::collapse.singles(combined_tree)
combined_tree <- ape::ladderize(combined_tree)

ape::write.tree(combined_tree, file = combined_tree_file)

lc_df <- readr::read_tsv(lc_summary_file, col_types = readr::cols(.default = "c"))
lc_arc_df <- readr::read_tsv(lc_arc_summary_file, col_types = readr::cols(.default = "c"))
ctrl_df <- readr::read_tsv(ctrl_summary_file, col_types = readr::cols(.default = "c"))

lc_clean <- clean_summary(lc_df, "Long_COVID")
lc_arc_clean <- clean_summary(lc_arc_df, "Long_COVID")
ctrl_clean <- clean_summary(ctrl_df, "Control")

if (prefixed_lc) {
  lc_clean <- lc_clean %>% mutate(user_genome = paste0("LC_", user_genome))
}

if (prefixed_ctrl) {
  ctrl_clean <- ctrl_clean %>% mutate(user_genome = paste0("CTRL_", user_genome))
}

combined_meta <- bind_rows(lc_clean, lc_arc_clean, ctrl_clean) %>%
  distinct(user_genome, .keep_all = TRUE)

itol_meta <- combined_meta %>%
  select(user_genome, dataset, phylum, genus, species, fastani_ani, fastani_af, msa_percent)

readr::write_tsv(itol_meta, out_meta_file)

mag_df <- readxl::read_excel(mag_excel_file)
names(mag_df) <- tolower(names(mag_df))

required_mag_cols <- c("mag", "coverage", "group")
missing_mag_cols <- setdiff(required_mag_cols, names(mag_df))
if (length(missing_mag_cols) > 0) {
  stop("Missing required columns in MAGs.xlsx: ", paste(missing_mag_cols, collapse = ", "))
}

if (!"phylum" %in% names(mag_df)) mag_df$phylum <- NA_character_
if (!"genus" %in% names(mag_df)) mag_df$genus <- NA_character_
if (!"species" %in% names(mag_df)) mag_df$species <- NA_character_
if (!"color" %in% names(mag_df)) mag_df$color <- NA_character_

mag_df <- mag_df %>%
  rename(
    MAG = mag,
    Coverage = coverage,
    Group = group,
    Phylum = phylum,
    Genus = genus,
    Species = species,
    Color = color
  ) %>%
  mutate(
    MAG = as.character(MAG),
    Group = normalize_group(Group),
    Phylum = as.character(Phylum),
    Genus = as.character(Genus),
    Species = as.character(Species),
    Color = as.character(Color),
    Coverage = as.numeric(Coverage),
    Species = if_else(is.na(Species) | Species == "", NA_character_, Species),
    Species = if_else(is.na(Species) & !is.na(Genus) & Genus != "", Genus, Species),
    Species = if_else(is.na(Species), MAG, Species),
    Sample = str_extract(MAG, "^[A-Za-z]+\\d+")
  ) %>%
  mutate(
    Sample = if_else(is.na(Sample), str_remove(MAG, "_.*$"), Sample)
  )

if (normalize_per_sample) {
  abundance_long <- mag_df %>%
    group_by(Sample, Group, Species) %>%
    summarise(Coverage = sum(Coverage, na.rm = TRUE), .groups = "drop") %>%
    group_by(Sample) %>%
    mutate(
      sample_total = sum(Coverage, na.rm = TRUE),
      value = dplyr::if_else(sample_total > 0, Coverage / sample_total, 0)
    ) %>%
    ungroup() %>%
    select(-sample_total)
  
  y_lab <- "Mean relative abundance"
} else {
  abundance_long <- mag_df %>%
    group_by(Sample, Group, Species) %>%
    summarise(value = sum(Coverage, na.rm = TRUE), .groups = "drop")
  
  y_lab <- "Mean Coverage"
}

group_stats <- abundance_long %>%
  group_by(Group, Species) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val = sd(value, na.rm = TRUE),
    n = sum(!is.na(value)),
    se = if_else(n > 0, sd_val / sqrt(n), NA_real_),
    .groups = "drop"
  ) %>%
  mutate(is_placeholder = is_placeholder(Species))

group_names <- unique(group_stats$Group)
if (length(group_names) < 2) stop("Need at least two groups in MAG data.")

gA <- "Control"
gB <- "Long_COVID"

if (!all(c(gA, gB) %in% group_names)) {
  gA <- group_names[1]
  gB <- group_names[2]
}

wide_stats <- group_stats %>%
  select(Group, Species, mean_val) %>%
  pivot_wider(names_from = Group, values_from = mean_val, values_fill = 0) %>%
  mutate(
    diff = .data[[gB]] - .data[[gA]],
    abs_diff = abs(diff),
    log2fc = log2((.data[[gB]] + pseudocount) / (.data[[gA]] + pseudocount))
  )

top_by_absdiff <- wide_stats %>%
  arrange(desc(abs_diff)) %>%
  slice_head(n = top_n) %>%
  pull(Species)

overall_ranking <- group_stats %>%
  group_by(Species) %>%
  summarise(overall_mean = mean(mean_val, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(overall_mean)) %>%
  pull(Species)

top_overall <- head(overall_ranking, top_n)

plot_top_overall_df <- group_stats %>%
  filter(Species %in% top_overall) %>%
  mutate(Species = factor(Species, levels = rev(top_overall)))

plot_top_change_df <- group_stats %>%
  filter(Species %in% top_by_absdiff) %>%
  mutate(Species = factor(Species, levels = rev(top_by_absdiff)))

full_species <- group_stats %>%
  filter(!is_placeholder) %>%
  group_by(Species) %>%
  summarise(overall_mean = mean(mean_val, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(overall_mean)) %>%
  slice_head(n = top_n) %>%
  pull(Species)

plot_full_df <- group_stats %>%
  filter(Species %in% full_species) %>%
  mutate(Species = factor(Species, levels = rev(full_species)))

placeholder_species <- group_stats %>%
  filter(is_placeholder) %>%
  group_by(Species) %>%
  summarise(overall_mean = mean(mean_val, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(overall_mean)) %>%
  slice_head(n = top_n) %>%
  pull(Species)

plot_placeholder_df <- group_stats %>%
  filter(Species %in% placeholder_species) %>%
  mutate(Species = factor(Species, levels = rev(placeholder_species)))

species_list <- unique(abundance_long$Species)

test_results <- purrr::map_dfr(species_list, function(sp) {
  sub <- abundance_long %>% filter(Species == sp)
  valsA <- sub %>% filter(Group == gA) %>% pull(value)
  valsB <- sub %>% filter(Group == gB) %>% pull(value)
  nA <- sum(!is.na(valsA))
  nB <- sum(!is.na(valsB))
  
  if (nA >= 2 && nB >= 2) {
    wt <- wilcox.test(valsB, valsA, alternative = "two.sided", exact = FALSE)
    tibble(Species = sp, p_value = wt$p.value, nA = nA, nB = nB)
  } else {
    tibble(Species = sp, p_value = NA_real_, nA = nA, nB = nB)
  }
}) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    sig = case_when(
      is.na(p_adj) ~ NA_character_,
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

dodge <- position_dodge(width = 0.8)

p1 <- ggplot(plot_top_overall_df, aes(x = Species, y = mean_val, fill = Group)) +
  geom_col(position = dodge, width = 0.7) +
  geom_errorbar(aes(ymin = mean_val - se, ymax = mean_val + se), width = 0.2, position = dodge) +
  coord_flip() +
  labs(x = "", y = y_lab, title = "") +
  scale_fill_manual(values = group_cols, name = "Group") +
  base_theme

ann1 <- make_annotation_df(plot_top_overall_df, test_results)
if (nrow(ann1) > 0) {
  p1 <- p1 + geom_text(
    data = ann1,
    aes(x = Species, y = ypos, label = sig),
    inherit.aes = FALSE,
    size = 3
  )
}

p2 <- ggplot(plot_top_change_df, aes(x = Species, y = mean_val, fill = Group)) +
  geom_col(position = dodge, width = 0.7) +
  geom_errorbar(aes(ymin = mean_val - se, ymax = mean_val + se), width = 0.2, position = dodge) +
  coord_flip() +
  labs(x = "", y = y_lab, title = "") +
  scale_fill_manual(values = group_cols, name = "Group") +
  base_theme

ann2 <- make_annotation_df(plot_top_change_df, test_results)
if (nrow(ann2) > 0) {
  p2 <- p2 + geom_text(
    data = ann2,
    aes(x = Species, y = ypos, label = sig),
    inherit.aes = FALSE,
    size = 3
  )
}

p3 <- ggplot(plot_full_df, aes(x = Species, y = mean_val, fill = Group)) +
  geom_col(position = dodge, width = 0.7) +
  geom_errorbar(aes(ymin = mean_val - se, ymax = mean_val + se), width = 0.2, position = dodge) +
  coord_flip() +
  labs(x = "", y = y_lab, title = "") +
  scale_fill_manual(values = group_cols, name = "Group") +
  base_theme

ann3 <- make_annotation_df(plot_full_df, test_results)
if (nrow(ann3) > 0) {
  p3 <- p3 + geom_text(
    data = ann3,
    aes(x = Species, y = ypos, label = sig),
    inherit.aes = FALSE,
    size = 3
  )
}

p4 <- NULL
if (nrow(plot_placeholder_df) > 0) {
  p4 <- ggplot(plot_placeholder_df, aes(x = Species, y = mean_val, fill = Group)) +
    geom_col(position = dodge, width = 0.7) +
    geom_errorbar(aes(ymin = mean_val - se, ymax = mean_val + se), width = 0.2, position = dodge) +
    coord_flip() +
    labs(x = "", y = y_lab, title = "") +
    scale_fill_manual(values = group_cols, name = "Group") +
    base_theme
  
  ann4 <- make_annotation_df(plot_placeholder_df, test_results)
  if (nrow(ann4) > 0) {
    p4 <- p4 + geom_text(
      data = ann4,
      aes(x = Species, y = ypos, label = sig),
      inherit.aes = FALSE,
      size = 3
    )
  }
}

counts_group <- mag_df %>%
  filter(!is.na(Phylum), Phylum != "") %>%
  count(Phylum, Group, name = "MAGs")

counts_sum <- counts_group %>%
  group_by(Phylum) %>%
  summarise(total = sum(MAGs), .groups = "drop") %>%
  arrange(total)

phylum_levels <- counts_sum$Phylum

counts_group <- counts_group %>%
  mutate(Phylum = factor(Phylum, levels = phylum_levels))

counts_sum <- counts_sum %>%
  mutate(Phylum = factor(Phylum, levels = phylum_levels))

x_max <- max(counts_sum$total, na.rm = TRUE)
text_offset <- max(1, ceiling(x_max * 0.03))

p_fixed <- ggplot(counts_group, aes(x = MAGs, y = Phylum, fill = Group)) +
  geom_col(position = "stack", color = "black", linewidth = 0.15) +
  geom_text(
    data = counts_sum,
    aes(x = total, y = Phylum, label = total),
    inherit.aes = FALSE,
    hjust = -0.05,
    size = 3.2
  ) +
  labs(x = "Number of MAGs", y = "", title = "") +
  scale_fill_manual(values = group_cols, na.value = "grey70") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, x_max + text_offset * 4)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.ticks.y = element_line(color = "black", linewidth = 0.4),
    axis.ticks.x = element_line(color = "black", linewidth = 0.4),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.6),
    legend.position = "right",
    legend.key = element_rect(fill = NA),
    plot.margin = margin(6, 20, 6, 6)
  )

otu_df <- abundance_long %>%
  group_by(Sample, Species) %>%
  summarise(val = sum(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Species, values_from = val, values_fill = 0) %>%
  arrange(Sample)

sample_meta <- abundance_long %>%
  distinct(Sample, Group) %>%
  arrange(Sample)

site_mat <- otu_df %>%
  column_to_rownames("Sample") %>%
  as.matrix()

mode(site_mat) <- "numeric"
site_mat[is.na(site_mat)] <- 0

alpha_df <- tibble(Sample = rownames(site_mat)) %>%
  mutate(
    richness = rowSums(site_mat > 0),
    shannon = vegan::diversity(site_mat, index = "shannon"),
    simpson = vegan::diversity(site_mat, index = "simpson")
  ) %>%
  left_join(sample_meta, by = "Sample")

if (length(unique(alpha_df$Group)) == 2) {
  print(wilcox.test(shannon ~ Group, data = alpha_df, exact = FALSE))
  print(wilcox.test(simpson ~ Group, data = alpha_df, exact = FALSE))
  print(wilcox.test(richness ~ Group, data = alpha_df, exact = FALSE))
} else {
  print(kruskal.test(shannon ~ Group, data = alpha_df))
  print(kruskal.test(simpson ~ Group, data = alpha_df))
  print(kruskal.test(richness ~ Group, data = alpha_df))
}

p_shannon <- ggplot(alpha_df, aes(x = Group, y = shannon, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = group_cols) +
  labs(x = "", y = "Shannon index", title = "Alpha diversity — Shannon") +
  alpha_theme

p_simpson <- ggplot(alpha_df, aes(x = Group, y = simpson, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = group_cols) +
  labs(x = "", y = "Simpson index", title = "Alpha diversity — Simpson") +
  alpha_theme

p_richness <- ggplot(alpha_df, aes(x = Group, y = richness, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = group_cols) +
  labs(x = "", y = "Species richness", title = "Alpha diversity — Richness") +
  alpha_theme

bc_dist <- vegan::vegdist(site_mat, method = "bray")
pcoa_res <- ape::pcoa(bc_dist)

pc_axes <- as.data.frame(pcoa_res$vectors[, 1:min(2, ncol(pcoa_res$vectors)), drop = FALSE])
names(pc_axes) <- c("Axis1", "Axis2")[1:ncol(pc_axes)]

pcoa_df <- tibble(Sample = rownames(site_mat)) %>%
  bind_cols(pc_axes) %>%
  left_join(sample_meta, by = "Sample")

if (!is.null(pcoa_res$values$Relative_eig) && length(pcoa_res$values$Relative_eig) >= 2) {
  pc1_var <- round(pcoa_res$values$Relative_eig[1] * 100, 1)
  pc2_var <- round(pcoa_res$values$Relative_eig[2] * 100, 1)
} else {
  pc1_var <- NA
  pc2_var <- NA
}

p_pcoa <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = Group)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.15, show.legend = FALSE) +
  scale_color_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  labs(
    x = paste0("PCoA1 (", pc1_var, "%)"),
    y = paste0("PCoA2 (", pc2_var, "%)"),
    title = "PCoA (Bray-Curtis)"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black")
  )

adonis_res <- vegan::adonis2(bc_dist ~ Group, data = sample_meta, permutations = 999, by = "margin")
bd <- vegan::betadisper(bc_dist, sample_meta$Group)
bd_test <- anova(bd)

p1_legend <- p1 + theme(legend.position = "bottom", legend.direction = "horizontal")
shared_legend <- cowplot::get_legend(p1_legend)

p1_noleg <- p1 + theme(legend.position = "none")
p2_noleg <- p_fixed + theme(legend.position = "none")
p3_noleg <- p3 + theme(legend.position = "none")
p4_noleg <- if (!is.null(p4)) p4 + theme(legend.position = "none") else ggplot() + theme_void()

top_pair <- cowplot::plot_grid(
  p1_noleg, p2_noleg,
  labels = c("A", "B"),
  label_size = 14,
  ncol = 2,
  rel_widths = c(0.48, 0.52),
  align = "hv"
)

legend_row <- cowplot::ggdraw(shared_legend)

diversity_block <- cowplot::plot_grid(
  top_pair,
  legend_row,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

abundance_block <- cowplot::plot_grid(
  p3_noleg, p4_noleg,
  labels = c("C", "D"),
  label_size = 14,
  ncol = 2,
  rel_widths = c(0.48, 0.52),
  align = "hv"
)

final_plot <- cowplot::plot_grid(
  diversity_block,
  abundance_block,
  ncol = 1,
  rel_heights = c(1, 1.2)
)

print(head(itol_meta, 10))
print(head(alpha_df, 12))
print(adonis_res)
print(bd_test)

print(p1)
print(p2)
print(p3)
if (!is.null(p4)) print(p4)
print(p_fixed)
print(final_plot)

print(p_shannon)
print(p_simpson)
print(p_richness)
print(p_pcoa)