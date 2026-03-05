
library(tidyverse)
library(stringr)
library(ggpubr)
library(ggplot2)
library(vegan)     
library(ggrepel)
library(RColorBrewer)
library(dplyr)
library(forcats)


LC_paths <- list.files(
  path = "~/Documents/African_Microbiome_Group/Long_COVID_study/Masters_research/Results/Mphlan4_tax/Long_covid_tax",
  pattern = "*_taxonomy_profile.tsv",
  full.names = TRUE
)
LC_names <- str_extract(basename(LC_paths), "L[1-7]")

control_paths <- list.files(
  path = "~/Documents/African_Microbiome_Group/Long_COVID_study/Masters_research/Results/Mphlan4_tax/Control_tax",
  pattern = "*taxonomy_profile.tsv",
  full.names = TRUE
)
control_names <- str_extract(basename(control_paths), "C[1-6]")

file_paths <- c(LC_paths, control_paths)
sample_names <- c(LC_names, control_names)

taxa_list <- lapply(seq_along(file_paths), function(i) {
  read_tsv(
    file_paths[i],
    comment   = "#",
    col_names = c("clade_name", "NCBI_tax_id", "relative_abundance"),
    show_col_types = FALSE
  ) %>%
    transmute(
      taxon         = clade_name,
      rel_abundance = relative_abundance,
      sample        = sample_names[i]
    )
})

taxa_df <- bind_rows(taxa_list)

abundance_matrix <- taxa_df %>%
  pivot_wider(names_from = sample, values_from = rel_abundance, values_fill = 0) %>%
  column_to_rownames("taxon")


abundance_t <- t(abundance_matrix)  

# ---------------------------
# Compute alpha diversity metrics
# ---------------------------

shannon_vals <- diversity(abundance_t, index = "shannon")

invsimpson_vals <- diversity(abundance_t, index = "invsimpson")

richness_vals <- specnumber(abundance_t)

pielou_vals <- shannon_vals / log(ifelse(richness_vals > 0, richness_vals, NA))

alpha_diversity_df <- tibble(
  sample = rownames(abundance_t),
  shannon = shannon_vals,
  invsimpson = invsimpson_vals,
  richness = richness_vals,
  pielou = pielou_vals
) %>%
  mutate(
    group = case_when(
      str_detect(sample, "^L[1-7]$") ~ "Long COVID",
      str_detect(sample, "^C[1-6]$") ~ "Control",
      TRUE ~ "other"
    )
  )


indices <- c("shannon", "invsimpson", "richness", "pielou")
pvals <- sapply(indices, function(idx) {
  x <- alpha_diversity_df[[idx]]
  df <- alpha_diversity_df %>% filter(!is.na(.data[[idx]]))
  tryCatch({
    wt <- wilcox.test(df[[idx]] ~ df$group, exact = FALSE) 
    wt$p.value
  }, error = function(e) NA_real_)
})

p_adj <- p.adjust(pvals, method = "BH")
stat_results <- tibble(
  index = indices,
  p_value = pvals,
  p_adj = p_adj
)
print(stat_results)


theme_LB <- function() {
  theme_bw(base_size = 14) +
    theme(
      text = element_text(family = "serif", face = "plain"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title = element_text(face = "plain"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    )
}

my_comparisons <- list(c("Long COVID", "Control"))
my_colors <- c("Control" = "#008f00", "Long COVID" = "#ff9300")

make_alpha_plot <- function(df, index_col, ylab, palette = my_colors) {
  lab_y <- max(df[[index_col]], na.rm = TRUE) + 0.2 * diff(range(df[[index_col]], na.rm = TRUE))
  ggplot(df, aes_string(x = "group", y = index_col, fill = "group")) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 2) +
    stat_compare_means(
      method = "wilcox.test",
      comparisons = my_comparisons,
      label = "p.signif",
      tip.length = 0.02,
      label.y = lab_y,
      hide.ns = FALSE
    ) +
    theme_LB() +
    labs(title = paste0(ylab, ""), y = ylab, x = NULL) +
    scale_fill_manual(values = palette) +
    theme(legend.position = "none")
}

alpha_plot_shannon   <- make_alpha_plot(alpha_diversity_df, "shannon", "Shannon Index")
alpha_plot_simpson   <- make_alpha_plot(alpha_diversity_df, "invsimpson", "Inverse Simpson")
alpha_plot_richness  <- make_alpha_plot(alpha_diversity_df, "richness", "Observed richness")
alpha_plot_pielou    <- make_alpha_plot(alpha_diversity_df, "pielou", "Pielou's evenness")

print(alpha_plot_shannon)
print(alpha_plot_simpson)
print(alpha_plot_richness)
print(alpha_plot_pielou)

# ---------------------------
# Beta diversity (Bray-Curtis) and PERMANOVA
# ---------------------------
bray_dist <- vegdist(abundance_t, method = "bray")
pcoa_res <- cmdscale(bray_dist, k = 2, eig = TRUE)

pcoa_df <- as.data.frame(pcoa_res$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$sample <- rownames(pcoa_df)
pcoa_df <- pcoa_df %>%
  mutate(group = case_when(
    str_detect(sample, "^L[1-7]$") ~ "Long COVID",
    str_detect(sample, "^C[1-6]$") ~ "Control",
    TRUE ~ "other"
  ))

# PERMANOVA
metadata <- data.frame(sample = rownames(abundance_t)) %>%
  mutate(group = case_when(
    str_detect(sample, "^L[1-7]$") ~ "Long COVID",
    str_detect(sample, "^C[1-6]$") ~ "Control",
    TRUE ~ "other"
  ))

adonis_result <- adonis2(bray_dist ~ group, data = metadata, permutations = 999)
print(adonis_result)

bd <- betadisper(bray_dist, metadata$group)
bd_test <- permutest(bd, permutations = 999)
print(bd)       
print(bd_test)  

eig_vals <- pcoa_res$eig
var_explained <- eig_vals / sum(eig_vals) * 100
x_lab <- paste0("PCoA1 (", round(var_explained[1], 1), "%)")
y_lab <- paste0("PCoA2 (", round(var_explained[2], 1), "%)")

pval <- adonis_result$`Pr(>F)`[1]
r2 <- round(adonis_result$R2[1], 3)

beta_plot <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = group, fill = group, label = sample)) +
  stat_ellipse(type = "t", geom = "polygon", alpha = 0.2, show.legend = FALSE) +
  geom_point(size = 4) +
  geom_text_repel(show.legend = FALSE, max.overlaps = Inf) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "",
    subtitle = sprintf("", r2, pval),
    x = x_lab, y = y_lab, color = "Group", fill = "Group"
  ) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme_LB() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12))

print(beta_plot)


## Relative abundance plots ##
# -------------------------
library(tidyverse)
library(forcats)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(readxl)
library(maaslin3)

long_covid_dir <- "/Users/louisburger22/Documents/African_Microbiome_Group/Long_COVID_study/Masters_research/Results/Mphlan4_tax/Long_covid_tax"
control_dir    <- "/Users/louisburger22/Documents/African_Microbiome_Group/Long_COVID_study/Masters_research/Results/Mphlan4_tax/Control_tax"
metadata_file  <- "~/Documents/African_Microbiome_Group/Long_COVID_study/LC_metadata.xlsx"
out_base       <- "maaslin3"
dir.create(out_base, showWarnings = FALSE)

TOP_N <- 20           
TOP_N_GENUS <- 10     
TOPN_select <- 50     
TOP_N_SPECIES <- 10

pseudocount <- 1e-6
prevalence_thresh_run <- 0   
abund_thresh_run     <- 0

read_norm_mp <- function(path) {
  fname <- basename(path)
  samp  <- str_extract(fname, "^[LC][0-9]+") %>% coalesce(str_extract(fname, "^[A-Za-z0-9\\-]+"))
  readr::read_tsv(path, comment = "#",
                  col_names = c("clade","taxid","abund","additional"),
                  show_col_types = FALSE) %>%
    mutate(sample = toupper(trimws(as.character(samp))),
           abund  = as.numeric(abund))
}

clean_taxonomy_cols <- function(df, ranks = c("k","p","c","o","f","g","s")) {
  df %>%
    mutate(across(all_of(ranks), ~ as.character(.x))) %>%
    mutate(across(all_of(ranks), ~ ifelse(.x == "" | is.na(.x), NA_character_, .x))) %>%
    mutate(across(all_of(ranks), ~ str_remove(.x, "^([a-zA-Z]__|[a-zA-Z]_+|_+)"))) %>%
    mutate(across(all_of(ranks), ~ ifelse(.x == "", NA_character_, .x)))
}

make_safe_name <- function(x) {
  x %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("[^A-Za-z0-9_\\-]", "_") %>%
    str_replace_all("^([0-9])", "_\\1")
}

make_display_name <- function(x, drop_make_unique_suffix = TRUE, first_only = FALSE) {
  res <- x
  if (drop_make_unique_suffix) res <- str_replace(res, "_\\d+$", "")
  if (first_only) res <- sub("_", " ", res)           
  res <- str_replace_all(res, "_", " ")
  str_trim(res)
}

long_files <- list.files(long_covid_dir, "*_taxonomy_profile.tsv", full.names = TRUE)
control_files <- list.files(control_dir, "*_taxonomy_profile.tsv", full.names = TRUE)
all_files <- c(long_files, control_files)
if (length(all_files) == 0) stop("No MetaPhlAn files found — check long_covid_dir / control_dir")

mp_raw <- purrr::map_dfr(all_files, read_norm_mp)

mp_parsed <- mp_raw %>%
  separate(clade, into = c("k","p","c","o","f","g","s"), sep = "\\|", fill = "right") %>%
  clean_taxonomy_cols(c("k","p","c","o","f","g","s")) %>%
  mutate(s = ifelse(is.na(s), NA_character_, s),
         sample = toupper(trimws(as.character(sample))))

summarise_rank_global <- function(df, rank_col) {
  df %>%
    filter(!is.na(.data[[rank_col]])) %>%
    group_by(Group = NA_character_, sample, taxon = .data[[rank_col]]) %>% 
    summarise(abund = sum(abund), .groups = "drop") %>%
    group_by(sample) %>%
    mutate(abund = 100 * abund / sum(abund)) %>%
    ungroup()
}

summarise_rank <- function(df, rank_col) {
  df %>%
    filter(!is.na(.data[[rank_col]])) %>%
    group_by(Group, sample, taxon = .data[[rank_col]]) %>%
    summarise(abund = sum(abund), .groups = "drop") %>%
    group_by(sample) %>%
    mutate(abund = 100 * abund / sum(abund)) %>%
    ungroup()
}

mp_parsed <- mp_parsed %>% mutate(Group = ifelse(str_detect(sample, "^L"), "Long COVID",
                                                 ifelse(str_detect(sample, "^C"), "Control", NA_character_)))

junk_genus_pattern <- regex("^GGB\\d+|^UBA|uncultured|metagenome|unclassified|Candidatus|^Ca\\.", ignore_case = TRUE)
junk_species_pattern <- regex("^GGB\\d+|^UBA|uncultured|metagenome|unclassified|candidatus|^Ca\\.\\bsp[\\._]|_sp_|\\bsp\\b", ignore_case = TRUE)

phylum_raw <- summarise_rank(mp_parsed, "p")

genus_raw <- mp_parsed %>%
  filter(!is.na(g)) %>%
  mutate(genus = str_remove_all(g, "^_+")) %>%
  mutate(genus = if_else(str_detect(genus, junk_genus_pattern), "Other", genus)) %>%
  group_by(Group, sample, genus) %>%
  summarise(abund = sum(abund), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(abund = 100 * abund / sum(abund)) %>%
  ungroup()

species_raw_all <- mp_parsed %>%
  filter(!is.na(s)) %>%                                  
  transmute(sample, species = s, phylum = p, abund) %>%  
  filter(!str_detect(species, junk_species_pattern)) %>% 
  group_by(sample, species, phylum) %>%
  summarise(abund = sum(as.numeric(abund)), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(abund = 100 * abund / sum(abund)) %>%           
  ungroup()

# quick archaeal report
archaeal_check <- mp_parsed %>%
  filter(!is.na(g), !is.na(k)) %>%
  mutate(genus_clean = str_remove_all(g, "^_+")) %>%
  filter(!str_detect(genus_clean, junk_genus_pattern)) %>%
  filter(str_to_lower(k) %in% c("archaea", "archaeal")) %>%
  group_by(genus_clean) %>%
  summarise(total_abund = sum(abund), .groups = "drop") %>%
  arrange(desc(total_abund))

if (nrow(archaeal_check) == 0) {
  message("No archaeal genera detected (after cleaning/filtering).")
} else {
  message("Archaeal genera detected (cleaned name and total % abundance across samples):")
  print(archaeal_check)
}

genus_by_domain <- mp_parsed %>%
  filter(!is.na(g), !is.na(k)) %>%
  mutate(genus_clean = str_remove_all(g, "^_+"),
         genus_clean = if_else(str_detect(genus_clean, junk_genus_pattern), "Other", genus_clean)) %>%
  group_by(k, genus_clean) %>%
  summarise(total = sum(abund), .groups = "drop")

top_bacterial <- genus_by_domain %>%
  filter(str_to_lower(k) %in% c("bacteria", "bacterial")) %>%
  arrange(desc(total)) %>%
  slice_head(n = TOP_N_GENUS) %>%
  pull(genus_clean)

archaeal_genera <- genus_by_domain %>%
  filter(str_to_lower(k) %in% c("archaea", "archaeal")) %>%
  arrange(desc(total)) %>%
  pull(genus_clean) %>%
  unique()

top_genus <- union(top_bacterial, archaeal_genera)

message("Top bacterial genera (top ", TOP_N_GENUS, "): ", paste(head(top_bacterial, 30), collapse = ", "))
if (length(archaeal_genera) > 0) message("Archaeal genera included: ", paste(archaeal_genera, collapse = ", "))

genus_df <- genus_raw %>%
  mutate(taxon = if_else(genus %in% top_genus, genus, "Other")) %>%
  group_by(Group, sample, taxon) %>%
  summarise(abund = sum(abund), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(abund = 100 * abund / sum(abund)) %>%
  ungroup()

phylum_df <- phylum_raw %>%
  mutate(taxon = recode(taxon,
                        "Firmicutes" = "Bacillota",
                        "Bacteroidetes" = "Bacteroidota",
                        "Actinobacteria" = "Actinomycetota",
                        "Proteobacteria" = "Pseudomonadota",
                        "Lentisphaerae" = "Lentisphaerota",
                        "Elusimicrobia" = "Elusimicrobiota",
                        "Fusobacteria" = "Fusobacteriota",
                        "Synergistetes" = "Synergistota",
                        "Tenericutes" = "Mycoplasmatota",
                        "Verrucomicrobia" = "Verrucomicrobiota",
                        "Euryarchaeota" = "Methanobacteriota",
                        .default = taxon)) %>%
  mutate(taxon = if_else(str_detect(taxon, regex("^Candidatus", ignore_case = TRUE)) |
                           str_detect(taxon, regex("unclassified$", ignore_case = TRUE)),
                         "Other", taxon)) %>%
  group_by(Group, sample, taxon) %>%
  summarise(abund = sum(abund), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(abund = 100 * abund / sum(abund)) %>%
  ungroup() %>%
  mutate(taxon = fct_lump_n(as.factor(taxon), n = TOP_N, other_level = "Other") %>% as.character())

order_taxa_flip <- function(df) {
  taxon_order <- df %>%
    group_by(taxon) %>%
    summarise(total = sum(abund), .groups = "drop") %>%
    arrange(total) %>%
    pull(taxon)
  if ("Other" %in% taxon_order) taxon_order <- c("Other", setdiff(taxon_order, "Other"))
  taxon_order <- rev(taxon_order)
  df %>% mutate(taxon = factor(taxon, levels = taxon_order))
}

phylum_df <- order_taxa_flip(phylum_df)
genus_df  <- order_taxa_flip(genus_df)

# ---------------- Color maps  ----------------
phylum_colors <- c(
  "Bacillota"        = "#E69F00",
  "Bacteroidota"     = "#009E73",
  "Actinomycetota"   = "#332288",
  "Pseudomonadota"   = "#F0E442",
  "Lentisphaerota"   = "#CC79A7",
  "Mycoplasmatota"   = "#56B4E9",
  "Verrucomicrobiota"= "#D55E00",
  "Fusobacteriota"   = "#0072B2",
  "Elusimicrobiota"  = "#999933",
  "Synergistota"     = "#117733",
  "Methanobacteriota"= "#882255",
  "Other"            = "grey80"
)

genus_colors <- c(
  "Prevotella"         = "#E6AB02",
  "Blautia"            = "#377EB8",
  "Faecalibacterium"   = "#4DAF4A",
  "Bifidobacterium"    = "#984EA3",
  "Bacteroides"        = "#FF7F00",
  "Escherichia"        = "#A65628",
  "Phocaeicola"        = "#F781BF",
  "Clostridium"        = "#6A3D9A",
  "Ruminococcus"       = "#66C2A5",
  "Alistipes"          = "#FFD92F",
  "Parabacteroides"    = "#8DA0CB",
  "Roseburia"          = "#FC8D62",
  "Coprococcus"        = "#B3DE69",
  "Dorea"              = "#A6CEE3",
  "Dysosmobacter"      = "#FB9A99",
  "Eubacterium"        = "#CAB2D6",
  "Mediterraneibacter" = "#FFED6F",
  "Oscillibacter"      = "#80B1D3",
  "Phascolarctobacterium" = "#FDBF6F",
  "Prevotellamassilia" = "#B2DF8A",
  "Streptococcus"      = "#33A02C",
  "Methanobrevibacter" = "#E41A1C",
  "Methanosphaera"     = "black",
  "Other"              = "grey90"
)

color_map <- c(phylum_colors, genus_colors)

# ----------------- Plots: phylum & genus stacked bars -----------------
bacillota_order <- phylum_df %>%
  filter(taxon == "Bacillota") %>%
  arrange(desc(abund)) %>%
  pull(sample)

phylum_df <- phylum_df %>%
  mutate(sample = factor(sample, levels = bacillota_order))

prevotella_order <- genus_df %>%
  group_by(sample) %>%
  summarise(prev_abund = sum(abund[taxon == "Prevotella"], na.rm = TRUE)) %>%
  arrange(desc(prev_abund), sample) %>%
  pull(sample)


genus_df <- genus_df %>%
  mutate(sample = factor(sample, levels = prevotella_order))

ggplot(phylum_df, aes(x = sample, y = abund, fill = taxon)) +
  geom_col() +
  facet_grid(~ Group, scales = "free_x", space = "free") +
  labs(x = "Sample", y = "Relative abundance (%)", fill = "Phylum") +
  scale_fill_manual(values = color_map) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5))

ggplot(genus_df, aes(x = sample, y = abund, fill = taxon)) +
  geom_col() +
  facet_grid(~ Group, scales = "free_x", space = "free") +
  labs(x = "Sample", y = "Relative abundance (%)", fill = "Genus") +
  scale_fill_manual(values = color_map) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5))

########### Bacillota/Bacteroidota ratio Suplementary figure 3 #############################
pseudocount <- ifelse(exists("pseudocount"), pseudocount, 1e-6)

bac_bact_ratio_df <- phylum_df %>%
  filter(taxon %in% c("Bacillota", "Bacteroidota")) %>%
  select(Group, sample, taxon, abund) %>%
  pivot_wider(names_from = taxon,
              values_from = abund,
              values_fill = list(abund = 0)) %>%
  mutate(
    Bacillota    = coalesce(Bacillota, 0),
    Bacteroidota = coalesce(Bacteroidota, 0),
    ratio        = (Bacillota + pseudocount) / (Bacteroidota + pseudocount),
    log2ratio    = log2(ratio)
  )

print(head(bac_bact_ratio_df, n = 20))

group_summary <- bac_bact_ratio_df %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    mean_ratio = mean(ratio, na.rm = TRUE),
    median_ratio = median(ratio, na.rm = TRUE),
    sd_ratio = sd(ratio, na.rm = TRUE),
    mean_log2 = mean(log2ratio, na.rm = TRUE),
    median_log2 = median(log2ratio, na.rm = TRUE),
    .groups = "drop"
  )

print(group_summary)

if (length(unique(na.omit(bac_bact_ratio_df$Group))) == 2) {
  wt <- tryCatch(
    wilcox.test(ratio ~ Group, data = bac_bact_ratio_df),
    error = function(e) NULL
  )
  if (!is.null(wt)) {
    message("Wilcoxon test p-value: ", signif(wt$p.value, 3))
  } else {
    message("Wilcoxon test unavailable (insufficient data or error).")
  }
}

group_colors <- c("Long COVID" = "#ff9300",   
                  "Control"    = "#008f00") 

ggplot(bac_bact_ratio_df, aes(x = Group, y = ratio, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Group), width = 0.15, height = 0, size = 2) +
  scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(y = "Bacillota / Bacteroidota (log10 scale)", x = "Group") +
  theme_bw() +
  theme(legend.position = "none")


########### Combine print ############

my_comparisons <- list(c("Control", "Long COVID"))

make_alpha_plot <- function(df, index_col, ylab, palette = my_colors) {
  rng <- range(df[[index_col]], na.rm = TRUE)
  lab_y <- max(df[[index_col]], na.rm = TRUE) + 0.15 * diff(rng)
  ggplot(df, aes(x = group, y = .data[[index_col]], fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 2) +
    stat_compare_means(
      method = "wilcox.test",
      comparisons = my_comparisons,
      label = "p.signif",
      tip.length = 0.02,
      label.y = lab_y,
      hide.ns = FALSE
    ) +
    theme_LB() +
    labs(title = NULL, y = ylab, x = NULL) +
    scale_fill_manual(values = palette, name = "Group") +
    theme(legend.position = "none")  
}

alpha_plot_shannon <- make_alpha_plot(alpha_diversity_df, "shannon", "Shannon Index") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())


centroids <- pcoa_df %>%
  group_by(group) %>%
  summarise(centroid1 = mean(PCoA1), centroid2 = mean(PCoA2), .groups = "drop")

pcoa_join <- left_join(pcoa_df, centroids, by = "group")

beta_plot_centroids <- ggplot(pcoa_join, aes(x = PCoA1, y = PCoA2, color = group)) +
  geom_segment(aes(x = centroid1, y = centroid2, xend = PCoA1, yend = PCoA2, color = group),
               alpha = 0.6, size = 0.6, show.legend = FALSE) +
  geom_point(size = 4, show.legend = FALSE) +
  geom_point(data = centroids, aes(x = centroid1, y = centroid2, fill = group),
             shape = 23, size = 5, stroke = 0.8, color = "black", show.legend = FALSE) +
  
  geom_text_repel(aes(label = sample), show.legend = FALSE, max.overlaps = Inf) +
  coord_fixed() +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme_LB() +
  labs(title = NULL,
       subtitle = sprintf("", adonis_result$R2[1], adonis_result$`Pr(>F)`[1]),
       x = x_lab, y = y_lab) +
  theme(legend.position = "none", plot.subtitle = element_text(size = 10))


p1 <- alpha_plot_shannon
p2 <- beta_plot_centroids

combined_plots <- plot_grid(
  p1, p2,
  labels = c("A", "B"),
  label_size = 14,
  ncol = 2,
  align = "hv",
  rel_widths = c(0.48, 0.52)
)

legend_plot <- ggplot(alpha_diversity_df, aes(x = group, y = shannon, fill = group, color = group)) +
  geom_point(aes(color = group, fill = group), size = 3) +
  scale_fill_manual(values = my_colors, name = "Group") +
  scale_color_manual(values = my_colors, name = "Group") +
  theme_LB() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.box.spacing = unit(0.4, "lines"),
        axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        panel.border = element_blank(), panel.grid = element_blank())

shared_legend <- get_legend(legend_plot)

final_plot <- plot_grid(
  combined_plots,
  ggdraw(shared_legend),
  ncol = 1,
  rel_heights = c(1, 0.12)  
)

print(final_plot)

bacillota_order <- phylum_df %>%
  filter(taxon == "Bacillota") %>%
  arrange(desc(abund)) %>%
  pull(sample)

phylum_df <- phylum_df %>%
  mutate(sample = factor(sample, levels = bacillota_order))

prevotella_order <- genus_df %>%
  group_by(sample) %>%
  summarise(prev_abund = sum(abund[taxon == "Prevotella"], na.rm = TRUE)) %>%
  arrange(desc(prev_abund), sample) %>%
  pull(sample)


genus_df <- genus_df %>%
  mutate(sample = factor(sample, levels = prevotella_order))


p3 <- ggplot(phylum_df, aes(x = sample, y = abund, fill = taxon)) +
  geom_col(width = 1) +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  labs(x = "Sample", y = "Relative abundance (%)", fill = "Phylum") +
  scale_fill_manual(values = color_map) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1))


p4 <- ggplot(genus_df, aes(x = sample, y = abund, fill = taxon)) +
  geom_col(width = 1) +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  labs(x = "Sample", y = "Relative abundance (%)", fill = "Genus") +
  scale_fill_manual(values = color_map) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5))

# Top block: diversity plots with shared legend underneath
legend_plot <- ggplot(alpha_diversity_df, aes(x = group, y = shannon, fill = group, color = group)) +
  geom_point(aes(color = group, fill = group), size = 3) +
  scale_fill_manual(values = my_colors, name = "Group") +
  scale_color_manual(values = my_colors, name = "Group") +
  theme_LB() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.box.spacing = unit(0.4, "lines"),
        axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        panel.border = element_blank(), panel.grid = element_blank())

shared_legend <- get_legend(legend_plot)
diversity_block <- plot_grid(
  plot_grid(p1, p2, labels = c("A", "B"), label_size = 14, ncol = 2, rel_widths = c(0.48, 0.52)),
  ggdraw(shared_legend),
  ncol = 1,
  rel_heights = c(1, 0.12)
)

# Bottom block: abundance plots
abundance_block <- plot_grid(
  p3, p4,
  labels = c("C", "D"),
  label_size = 14,
  ncol = 2,
  rel_widths = c(0.48, 0.52)
)

# Final combined layout
final_plot2 <- plot_grid(
  diversity_block,
  abundance_block,
  ncol = 1,
  rel_heights = c(1, 1.2)  
)

print(final_plot2)


##### Differential abundance plot
genus_raw <- genus_raw %>%
  filter(!str_detect(genus, regex("^GGB\\d+|^UBA|uncultured|metagenome|unclassified|candidatus|^Ca\\.", ignore_case = TRUE)))

genus_safe <- genus_raw %>%
  mutate(
    genus_safe = genus %>%
      str_replace_all("\\s+", "_") %>%                     
      str_remove_all("[^A-Za-z0-9_\\-]") %>%               
      str_remove("^([0-9])")                               
  )

genus_abund_wide <- genus_safe %>%
  select(sample, genus_safe, abund) %>%
  group_by(genus_safe, sample) %>%                         
  summarise(abund = sum(abund), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = abund, values_fill = 0) %>%
  column_to_rownames("genus_safe")

sample_metadata <- genus_raw %>%
  select(sample, Group) %>%
  distinct() %>%
  arrange(sample) %>%
  mutate(Group = factor(Group)) %>%                        
  column_to_rownames("sample")

if (!all(colnames(genus_abund_wide) %in% rownames(sample_metadata))) {
  stop("Mismatch between sample names in abundance matrix and sample metadata. Check sample IDs.")
}
sample_metadata <- sample_metadata[colnames(genus_abund_wide), , drop = FALSE]

sample_metadata$Group <- relevel(sample_metadata$Group, ref = "Control")

output_dir <- "maaslin2_output_genera"
fit <- Maaslin2(
  input_data     = genus_abund_wide,
  input_metadata = sample_metadata,
  output         = output_dir,
  fixed_effects  = c("Group"),
  normalization  = "TSS",
  transform      = "AST",
  reference      = "Control"
)

sig_file <- file.path(output_dir, "significant_results.tsv")
if (!file.exists(sig_file)) stop("MaAsLin2 output file not found: ", sig_file)
sig_results <- read.delim(sig_file, stringsAsFactors = FALSE)

plot_df <- sig_results %>%
  mutate(
    Genus = str_remove(feature, "^X_"),
    Direction = if_else(coef > 0, "Enriched in Long COVID", "Depleted in Long COVID"),
    label = case_when(
      qval < 0.001 ~ "***",
      qval < 0.01  ~ "**",
      qval < 0.05  ~ "*",
      TRUE ~ ""
    ),
    hjust = if_else(coef > 0, -0.1, 1.1)
  ) %>%
  arrange(coef)

DA_plot <- ggplot(plot_df, aes(x = reorder(Genus, coef), y = coef, fill = Direction)) +
  geom_col(width = 1) +  # bars touching
  coord_flip() +
  geom_hline(yintercept = 0, color = "black") +
  geom_text(aes(label = label, hjust = hjust), size = 5) +
  scale_fill_manual(values = c(
    "Enriched in Long COVID" = "#ff9300",
    "Depleted in Long COVID" = "#008f00"
  )) +
  labs(x = NULL, y = "Effect size (Coefficient)", title = "", fill = NULL) +
  theme_bw(base_size = 14) +
  theme(
    axis.title.x = element_text(color = "black"),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank()
  )


final_da_plot <- plot_grid(
  DA_plot, labels = "E",
  label_size = 14,
  ncol = 1
)

print(final_da_plot)

