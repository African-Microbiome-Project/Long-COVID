library(tidyverse)
library(data.table)
library(maaslin3)
library(pheatmap)
library(ggrepel)
library(vegan)

longcov_file_path  <- "/Users/louisburger22/Documents/African_Microbiome_Group/Long_COVID_study/Masters_research/Results/Functional_analysis/Long_covid/pathways_merged.tsv"
control_file_path  <- "/Users/louisburger22/Documents/African_Microbiome_Group/Long_COVID_study/Masters_research/Results/Functional_analysis/Control/pathways_merged.tsv"

longcov_raw_df <- fread(longcov_file_path, sep = "\t", header = TRUE, data.table = FALSE, stringsAsFactors = FALSE)
control_raw_df <- fread(control_file_path, sep = "\t", header = TRUE, data.table = FALSE, stringsAsFactors = FALSE)

colnames(longcov_raw_df)[1] <- "Feature"
colnames(control_raw_df)[1] <- "Feature"

longcov_sample_names <- setdiff(colnames(longcov_raw_df), "Feature")
control_sample_names <- setdiff(colnames(control_raw_df), "Feature")

merged_features_df <- full_join(longcov_raw_df, control_raw_df, by = "Feature")

unwanted_pattern <- "(?i)UNMAPPED|UNINTEG|UNINTEGRATED|UNASSIGNED|UNINTR"
merged_filtered_df <- merged_features_df %>%
  filter(!str_detect(Feature, regex(unwanted_pattern, ignore_case = TRUE)))

merged_stratified_df <- merged_filtered_df %>%
  filter(str_detect(Feature, "\\|"))

split_feature_df <- merged_stratified_df %>%
  mutate(tmp_pair = str_split_fixed(Feature, "\\|", 2),
         pathway_name_raw = str_trim(tmp_pair[, 1]),
         organism_raw = str_trim(tmp_pair[, 2])) %>%
  select(-tmp_pair)

clean_organism_fn <- function(x) {
  x %>%
    str_replace_all("([kpcofgst]__)", "") %>%
    str_replace_all("[\\[\\]\\(\\)]", "") %>%
    str_replace_all("[;|_]", " ") %>%
    str_squish()
}

split_feature_df <- split_feature_df %>%
  mutate(organism_name = clean_organism_fn(organism_raw),
         pathway_name = str_squish(pathway_name_raw))

sample_cols_vector <- setdiff(colnames(split_feature_df), c("Feature", "pathway_name_raw", "organism_raw", "organism_name", "pathway_name"))

stratified_long_df <- split_feature_df %>%
  pivot_longer(cols = all_of(sample_cols_vector), names_to = "sample_id", values_to = "abundance") %>%
  mutate(abundance = as.numeric(abundance)) %>%
  select(pathway_name, organism_name, sample_id, abundance)

pathway_sample_sum_df <- stratified_long_df %>%
  group_by(pathway_name, sample_id) %>%
  summarise(pathway_abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

pathway_matrix_wide_df <- pathway_sample_sum_df %>%
  pivot_wider(names_from = pathway_name, values_from = pathway_abundance, values_fill = 0) %>%
  column_to_rownames(var = "sample_id")

all_sample_ids_vector <- rownames(pathway_matrix_wide_df)

metadata_inferred_df <- tibble(sample = all_sample_ids_vector,
                               Group = case_when(
                                 sample %in% longcov_sample_names ~ "LongCOVID",
                                 sample %in% control_sample_names ~ "Control",
                                 TRUE ~ "Unknown"
                               ))

metadata_inferred_df <- metadata_inferred_df %>% filter(Group %in% c("LongCOVID", "Control"))
samples_to_keep <- metadata_inferred_df$sample

pathway_matrix_for_maaslin_df <- pathway_matrix_wide_df[samples_to_keep, , drop = FALSE]

metadata_for_maaslin_df <- metadata_inferred_df %>%
  mutate(Group = as.character(Group),
         Group = case_when(
           Group %in% c("LongCOVID", "Long_COVID", "Long COVID") ~ "Long_COVID",
           Group %in% c("Control", "control", "CTRL") ~ "Control",
           TRUE ~ Group
         )) %>%
  column_to_rownames(var = "sample")

metadata_for_maaslin_df$Group <- factor(metadata_for_maaslin_df$Group, levels = c("Control", "Long_COVID"))
metadata_for_maaslin_df <- metadata_for_maaslin_df[rownames(pathway_matrix_for_maaslin_df), , drop = FALSE]

maaslin_output_dir <- path.expand("~/Documents/African_Microbiome_Group/Long_COVID_study/Masters_research/Results/Functional_analysis/maaslin_pathways_run")
dir.create(maaslin_output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(42)
maaslin3::maaslin3(
  input_data = pathway_matrix_for_maaslin_df,
  input_metadata = metadata_for_maaslin_df,
  output = maaslin_output_dir,
  fixed_effects = c("Group"),
  normalization = "TSS",
  transform = "LOG",
  reference = "Group,Control",
  standardize = TRUE,
  cores = 1
)

maaslin_all_results_file <- list.files(maaslin_output_dir, pattern = "all_results.*\\.tsv$", full.names = TRUE)[1]
if (is.na(maaslin_all_results_file) || !file.exists(maaslin_all_results_file)) stop("Could not find all_results TSV in output directory: ", maaslin_output_dir)
maaslin_results_df <- data.table::fread(maaslin_all_results_file, data.table = FALSE)

pw_mat <- as.matrix(pathway_matrix_for_maaslin_df)
if (nrow(pw_mat) < 3) stop("Need at least 3 samples in pathway matrix for ordination.")
row_sums <- rowSums(pw_mat, na.rm = TRUE)
row_sums[row_sums == 0] <- 1
prop_mat <- sweep(pw_mat, 1, row_sums, FUN = "/")

meta_df <- metadata_for_maaslin_df %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::rename(Group = Group) %>%
  mutate(
    Group = as.character(Group),
    Group_display = case_when(
      Group %in% c("LongCOVID", "Long_COVID", "Long COVID") ~ "Long COVID",
      Group %in% c("Control", "control", "CTRL") ~ "Control",
      TRUE ~ NA_character_
    ),
    Group_display = factor(Group_display, levels = c("Control", "Long COVID"))
  ) %>%
  filter(!is.na(Group_display))

samples_keep <- intersect(rownames(prop_mat), meta_df$sample)
if (length(samples_keep) < 3) stop("Not enough overlapping samples between pathway matrix and metadata for NMDS.")
prop_mat <- prop_mat[samples_keep, , drop = FALSE]

meta_for_ord <- meta_df %>%
  filter(sample %in% samples_keep) %>%
  column_to_rownames("sample")
meta_for_ord <- meta_for_ord[rownames(prop_mat), , drop = FALSE]

dist_bc <- vegan::vegdist(prop_mat, method = "bray", na.rm = TRUE)

set.seed(42)
nmds_try <- try(vegan::metaMDS(prop_mat, distance = "bray", k = 2, trymax = 50, autotransform = FALSE), silent = TRUE)
if (inherits(nmds_try, "try-error")) nmds_try <- vegan::metaMDS(prop_mat, distance = "bray", k = 2, trymax = 200, autotransform = FALSE)
nmds_res <- nmds_try

nmds_df <- as.data.frame(nmds_res$points) %>%
  tibble::rownames_to_column("sample") %>%
  rename(NMDS1 = MDS1, NMDS2 = MDS2) %>%
  left_join(meta_for_ord %>% rownames_to_column("sample") %>% select(sample, Group_display), by = "sample")

centroids <- nmds_df %>%
  group_by(Group_display) %>%
  summarise(cent_NMDS1 = mean(NMDS1, na.rm = TRUE),
            cent_NMDS2 = mean(NMDS2, na.rm = TRUE),
            .groups = "drop")

nmds_joined <- left_join(nmds_df, centroids, by = "Group_display")

meta_for_perm <- meta_for_ord %>%
  tibble::rownames_to_column("sample") %>%
  filter(!is.na(Group_display)) %>%
  tibble::column_to_rownames("sample")

keep <- intersect(rownames(meta_for_perm), rownames(as.matrix(dist_bc)))
meta_for_perm <- meta_for_perm[keep, , drop = FALSE]
dist_bc_perm <- as.dist(as.matrix(dist_bc)[keep, keep])

adon_res <- vegan::adonis2(dist_bc_perm ~ Group_display, data = meta_for_perm, permutations = 999, by = "margin")
adon_tab <- as.data.frame(adon_res)

if ("Group_display" %in% rownames(adon_tab)) {
  ad_p <- signif(as.numeric(adon_tab["Group_display", "Pr(>F)"]), 3)
  ad_R2 <- signif(as.numeric(adon_tab["Group_display", "R2"]), 3)
} else {
  ad_p <- signif(as.numeric(adon_tab[1, "Pr(>F)"]), 3)
  ad_R2 <- signif(as.numeric(adon_tab[1, "R2"]), 3)
}

my_theme <- theme_classic() +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(color = "black"))

if (!exists("group_cols")) group_cols <- c("Control" = "#008f00", "Long COVID" = "#ff9300")
if (!all(c("Control", "Long COVID") %in% names(group_cols))) group_cols <- c("Control" = "#008f00", "Long COVID" = "#ff9300")

p_centroid <- ggplot() +
  geom_segment(data = nmds_joined,
               aes(x = NMDS1, y = NMDS2, xend = cent_NMDS1, yend = cent_NMDS2, color = Group_display),
               alpha = 0.25, size = 0.4, show.legend = FALSE) +
  geom_point(data = nmds_df, aes(x = NMDS1, y = NMDS2, color = Group_display), size = 2.4, alpha = 0.95) +
  geom_point(data = centroids, aes(x = cent_NMDS1, y = cent_NMDS2, fill = Group_display), shape = 23, size = 6, color = "black", show.legend = FALSE) +
  scale_color_manual(name = "Group", values = group_cols[c("Control", "Long COVID")]) +
  scale_fill_manual(values = group_cols[c("Control", "Long COVID")], guide = FALSE) +
  labs(title = paste0("NMDS (Bray-Curtis) — centroid plot    stress=", round(nmds_res$stress, 3)),
       subtitle = paste0("PERMANOVA: p = ", ad_p, "    R2 = ", ad_R2),
       x = "NMDS1", y = "NMDS2") +
  my_theme +
  theme(legend.position = "right")

print(p_centroid)

################################################################################
###################### Preparing data for pllotting ############################
################################################################################

# Libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(forcats)

group_cols <- c("Control" = "#008f00", "Long COVID" = "#ff9300")

get_meta_df <- function() {
  if (exists("metadata_for_maaslin_df")) {
    metadata_for_maaslin_df %>%
      as.data.frame() %>%
      tibble::rownames_to_column("sample") %>%
      dplyr::select(sample, Group)
  } else if (exists("metadata_inferred_df")) {
    metadata_inferred_df %>%
      dplyr::select(sample, Group)
  } else {
    stop("No metadata object found. Provide metadata_for_maaslin_df or metadata_inferred_df.")
  }
}

meta_df <- get_meta_df() %>%
  mutate(
    Group = as.character(Group),
    Group_display = case_when(
      Group %in% c("LongCOVID", "Long_COVID", "Long COVID") ~ "Long COVID",
      Group %in% c("Control", "control", "CTRL") ~ "Control",
      TRUE ~ NA_character_
    ),
    Group_display = factor(Group_display, levels = c("Control", "Long COVID"))
  ) %>%
  filter(!is.na(Group_display))

q_choice <- dplyr::case_when(
  "qval_individual" %in% colnames(maaslin_results_df) ~ "qval_individual",
  "qval_joint" %in% colnames(maaslin_results_df) ~ "qval_joint",
  TRUE ~ NA_character_
)
if (is.na(q_choice)) stop("No q-value column found (expected qval_individual or qval_joint).")

my_theme <- theme_classic() +
  theme(
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(color = "black"),
    legend.title = element_blank(),
    legend.position = "right"
  )

volcano_plot_df <- maaslin_results_df %>%
  mutate(
    coef = as.numeric(coef),
    qvalue = as.numeric(.data[[q_choice]])
  ) %>%
  filter(!is.na(coef) & !is.na(qvalue)) %>%
  mutate(
    regulation = case_when(
      qvalue <= 0.05 & coef > 0 ~ "Up-regulated",
      qvalue <= 0.05 & coef < 0 ~ "Down-regulated",
      TRUE ~ "Not significant"
    ),
    regulation = factor(regulation, levels = c("Up-regulated", "Down-regulated", "Not significant"))
  )

volcano_color_map <- c(
  "Up-regulated" = "#ff9300",
  "Down-regulated" = "#008f00",
  "Not significant" = "grey60"
)

top_n_labels <- 12
top_hits_df <- volcano_plot_df %>% arrange(qvalue) %>% slice_head(n = top_n_labels)

ggplot(volcano_plot_df, aes(x = coef, y = -log10(qvalue + 1e-300), color = regulation)) +
  geom_point(alpha = 0.9, size = 2) +
  geom_vline(xintercept = 0, color = "gray60") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  scale_color_manual(values = volcano_color_map, breaks = levels(volcano_plot_df$regulation)) +
  geom_text_repel(data = top_hits_df, aes(label = feature), size = 3, max.overlaps = 20, box.padding = 0.3) +
  labs(x = "Effect size (coef)", y = "-log10(q-value)", title = "Volcano: MaAsLin3 pathways") +
  my_theme

top_n_heat <- 50

feat_ranked <- maaslin_results_df %>%
  mutate(
    coef = as.numeric(coef),
    qvalue = as.numeric(.data[[q_choice]]),
    abs_coef = abs(coef)
  ) %>%
  filter(!is.na(feature) & !is.na(coef) & !is.na(qvalue)) %>%
  arrange(qvalue, desc(abs_coef)) %>%
  distinct(feature, .keep_all = TRUE)

sig_features <- feat_ranked %>%
  slice_head(n = top_n_heat) %>%
  pull(feature)

sig_features <- intersect(sig_features, colnames(pathway_matrix_for_maaslin_df))
if (length(sig_features) < 2) stop("Too few selected pathways found in pathway matrix for heatmap.")

heatmat <- t(pathway_matrix_for_maaslin_df[, sig_features, drop = FALSE])

sample_order <- intersect(meta_df$sample, colnames(heatmat))
if (length(sample_order) > 0) heatmat <- heatmat[, sample_order, drop = FALSE]

row_sds <- apply(heatmat, 1, sd, na.rm = TRUE)
keep_rows <- which(is.finite(row_sds) & row_sds > 0)
if (length(keep_rows) == 0) stop("All selected heatmap pathways have zero variance.")
heatmat <- heatmat[keep_rows, , drop = FALSE]

heatmat_z <- t(scale(t(heatmat)))
heatmat_z[!is.finite(heatmat_z)] <- 0

annot_df <- meta_df %>%
  dplyr::select(sample, Group_display) %>%
  tibble::column_to_rownames("sample")
annot_df <- annot_df[colnames(heatmat_z), , drop = FALSE]
annot_df$Group_display <- as.character(annot_df$Group_display)

anno_col_list <- list(Group_display = group_cols)

ha_top <- ComplexHeatmap::HeatmapAnnotation(
  df = annot_df,
  col = anno_col_list,
  annotation_height = grid::unit(6, "mm")
)

col_fun <- circlize::colorRamp2(c(min(heatmat_z), 0, max(heatmat_z)), c("#2166ac", "white", "#b2182b"))

legend_param <- list(title = "Z-score", title_position = "topcenter", legend_direction = "vertical")

ht <- ComplexHeatmap::Heatmap(
  heatmat_z,
  name = "zscore",
  col = col_fun,
  top_annotation = ha_top,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = grid::gpar(fontsize = 9),
  row_names_rot = 0,
  row_names_max_width = grid::unit(12, "cm"),
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  heatmap_legend_param = legend_param
)

grid::grid.newpage()
ComplexHeatmap::draw(ht, heatmap_legend_side = "right", padding = grid::unit(c(4, 4, 4, 4), "mm"))