# Figure 3 script
## This script won't run without "Alpha_Beta_Relabund_DAFigure2.R" script (This script is a follow up after Figure 2's script)
### Therefore, tun this script after running "Alpha_Beta_Relabund_DAFigure2.R"

species_map <- species_raw_all %>%
  distinct(species, phylum) %>%
  arrange(phylum, species) %>%
  mutate(species_safe = make_safe_name(species))

species_map$species_safe <- make.unique(species_map$species_safe, sep = "_")

species_long_safe <- species_raw_all %>%
  left_join(species_map, by = c("species", "phylum"))

species_abund_wide_df <- species_long_safe %>%
  group_by(sample, species_safe) %>%
  summarise(abund = sum(abund), .groups = "drop") %>%
  pivot_wider(names_from = species_safe, values_from = abund, values_fill = 0) %>%
  arrange(sample)

if ("sample" %in% colnames(species_abund_wide_df)) {
  species_abund_wide_df <- column_to_rownames(species_abund_wide_df, "sample")
}

rownames(species_abund_wide_df) <- toupper(rownames(species_abund_wide_df))

meta_raw <- if (file.exists(metadata_file)) {
  readxl::read_excel(metadata_file, sheet = 1) %>% as_tibble()
} else {
  tibble(sample = character(0))
}

if (nrow(meta_raw) > 0) {
  names(meta_raw) <- names(meta_raw) %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("[:punct:]", "") %>%
    tolower()
  
  sid_matches <- which(grepl("^sample$|^sampleid$|^sample_id$|^id$", names(meta_raw)))
  sid_col <- if (length(sid_matches) > 0) names(meta_raw)[sid_matches[1]] else names(meta_raw)[1]
  
  meta_raw <- meta_raw %>%
    rename(sample = !!sym(sid_col))
  
  meta_raw$sample <- toupper(trimws(as.character(meta_raw$sample)))
}

meta_from_excel <- meta_raw %>%
  filter(sample %in% rownames(species_abund_wide_df)) %>%
  distinct(sample, .keep_all = TRUE)

missing_samples <- setdiff(rownames(species_abund_wide_df), meta_from_excel$sample)

minimal_rows <- tibble(sample = missing_samples) %>%
  mutate(
    Group = case_when(
      str_detect(sample, "^L") ~ "Long COVID",
      str_detect(sample, "^C") ~ "Control",
      TRUE ~ NA_character_
    )
  )

sample_metadata_complete <- bind_rows(meta_from_excel, minimal_rows) %>%
  distinct(sample, .keep_all = TRUE) %>%
  arrange(sample) %>%
  column_to_rownames("sample")

sample_metadata_complete$Group <- as.character(sample_metadata_complete$Group)

sample_metadata_complete$Group <- case_when(
  !is.na(sample_metadata_complete$Group) & str_to_lower(sample_metadata_complete$Group) %in% c("control", "controls", "ctrl") ~ "Control",
  !is.na(sample_metadata_complete$Group) & str_to_lower(sample_metadata_complete$Group) %in% c("long covid", "long_covid", "long-covid", "longcovid") ~ "Long COVID",
  TRUE ~ sample_metadata_complete$Group
)

missing_idx <- which(is.na(sample_metadata_complete$Group) | sample_metadata_complete$Group == "")

if (length(missing_idx) > 0) {
  rn <- rownames(sample_metadata_complete)[missing_idx]
  inferred <- case_when(
    str_detect(rn, regex("^L", ignore_case = TRUE)) ~ "Long COVID",
    str_detect(rn, regex("^C", ignore_case = TRUE)) ~ "Control",
    TRUE ~ NA_character_
  )
  sample_metadata_complete$Group[missing_idx] <- inferred
}

sample_metadata_complete$Group <- factor(sample_metadata_complete$Group, levels = c("Control", "Long COVID"))

if (!("overweight_flag" %in% names(sample_metadata_complete))) {
  sample_metadata_complete$overweight_flag <- NA_character_
}

if ("comorbidities" %in% names(sample_metadata_complete)) {
  sample_metadata_complete$overweight_flag <- ifelse(
    !is.na(sample_metadata_complete$comorbidities) &
      str_detect(tolower(as.character(sample_metadata_complete$comorbidities)), "overweight|obes|obesity"),
    "Yes",
    sample_metadata_complete$overweight_flag
  )
}

sample_metadata_complete$overweight_flag <- case_when(
  tolower(as.character(sample_metadata_complete$overweight_flag)) %in% c("yes", "y", "true", "1", "overweight", "obese") ~ "Yes",
  tolower(as.character(sample_metadata_complete$overweight_flag)) %in% c("no", "n", "false", "0") ~ "No",
  TRUE ~ as.character(sample_metadata_complete$overweight_flag)
)

if (!("severe_flag" %in% names(sample_metadata_complete))) {
  sample_metadata_complete$severe_flag <- NA_character_
}

if ("covid19severity" %in% names(sample_metadata_complete)) {
  sample_metadata_complete$severe_flag <- ifelse(
    !is.na(sample_metadata_complete$covid19severity) &
      str_detect(tolower(as.character(sample_metadata_complete$covid19severity)), "sever"),
    "Severe",
    sample_metadata_complete$severe_flag
  )
}

controls <- rownames(sample_metadata_complete)[
  !is.na(sample_metadata_complete$Group) & sample_metadata_complete$Group == "Control"
]

longcovid_all <- rownames(sample_metadata_complete)[
  !is.na(sample_metadata_complete$Group) & sample_metadata_complete$Group == "Long COVID"
]

overweight_cases <- rownames(sample_metadata_complete)[
  !is.na(sample_metadata_complete$overweight_flag) &
    tolower(as.character(sample_metadata_complete$overweight_flag)) == "yes" &
    sample_metadata_complete$Group == "Long COVID"
]

severe_cases <- rownames(sample_metadata_complete)[
  !is.na(sample_metadata_complete$severe_flag) &
    sample_metadata_complete$severe_flag == "Severe" &
    sample_metadata_complete$Group == "Long COVID"
]

message(
  "After normalization -> Controls: ", length(controls),
  "; Long COVID: ", length(longcovid_all),
  "; Overweight cases: ", length(overweight_cases),
  "; Severe cases: ", length(severe_cases)
)

print(table(sample_metadata_complete$Group, useNA = "ifany"))

run_maaslin_allfeatures <- function(case_samples, case_name, abund_df, sample_meta, outdir_base,
                                    pseudocount = 1e-6, prevalence_thresh = 0, abund_thresh = 0) {
  outdir <- file.path(outdir_base, paste0("Control_vs_", case_name))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  if (length(case_samples) == 0) {
    message("No cases for ", case_name, " -> skipping")
    return(NULL)
  }
  
  controls_local <- rownames(sample_meta)[
    !is.na(sample_meta$Group) & sample_meta$Group == "Control"
  ]
  
  keep <- intersect(c(controls_local, case_samples), rownames(abund_df))
  
  if (length(keep) < 2) {
    message("Not enough samples for ", case_name, " -> skipping")
    return(NULL)
  }
  
  abund_sub <- abund_df[keep, , drop = FALSE]
  meta_sub <- sample_meta[keep, , drop = FALSE]
  
  meta_sub$ComparisonGroup <- factor(
    ifelse(rownames(meta_sub) %in% controls_local, "Control", case_name),
    levels = c("Control", case_name)
  )
  
  abund_mat <- as.data.frame(abund_sub)
  
  if (any(abund_mat == 0)) {
    abund_mat <- abund_mat + pseudocount
  }
  
  feat_prev <- colSums(abund_mat > 0) / nrow(abund_mat)
  feat_mean <- colMeans(abund_mat)
  
  keep_feats <- names(feat_prev)[feat_prev >= prevalence_thresh & feat_mean >= abund_thresh]
  
  if (length(keep_feats) < 1) {
    stop("No features to analyze after prefiltering")
  }
  
  abund_filtered <- abund_mat[, keep_feats, drop = FALSE]
  meta_filtered <- meta_sub[rownames(abund_filtered), , drop = FALSE]
  
  message(
    "Running MaAsLin3 for Control vs ", case_name,
    " (samples=", nrow(abund_filtered),
    ", features=", ncol(abund_filtered), ")"
  )
  
  fit <- tryCatch(
    {
      maaslin3(
        input_data = abund_filtered,
        input_metadata = meta_filtered,
        output = outdir,
        fixed_effects = c("ComparisonGroup"),
        normalization = "TSS",
        transform = "LOG",
        reference = "Control"
      )
    },
    error = function(e) {
      message("MaAsLin3 error for ", case_name, ": ", e$message)
      NULL
    }
  )
  
  if (is.null(fit)) {
    return(NULL)
  }
  
  message("MaAsLin3 finished for ", case_name, " -> outdir: ", outdir)
  outdir
}

out_long <- run_maaslin_allfeatures(
  longcovid_all,
  "Long_COVID",
  species_abund_wide_df,
  sample_metadata_complete,
  out_base,
  pseudocount = pseudocount,
  prevalence_thresh = prevalence_thresh_run,
  abund_thresh = abund_thresh_run
)

out_over <- run_maaslin_allfeatures(
  overweight_cases,
  "Overweight_LongCOVID",
  species_abund_wide_df,
  sample_metadata_complete,
  out_base,
  pseudocount = pseudocount,
  prevalence_thresh = prevalence_thresh_run,
  abund_thresh = abund_thresh_run
)

out_sev <- run_maaslin_allfeatures(
  severe_cases,
  "Severe_LongCOVID",
  species_abund_wide_df,
  sample_metadata_complete,
  out_base,
  pseudocount = pseudocount,
  prevalence_thresh = prevalence_thresh_run,
  abund_thresh = abund_thresh_run
)

read_maaslin_signed <- function(maaslin_outdir) {
  if (is.null(maaslin_outdir) || !dir.exists(maaslin_outdir)) {
    return(NULL)
  }
  
  files <- list.files(maaslin_outdir, pattern = "all_results|results", full.names = TRUE)
  if (length(files) == 0) {
    files <- list.files(maaslin_outdir, full.names = TRUE)
  }
  if (length(files) == 0) {
    return(NULL)
  }
  
  all_res <- read.delim(files[1], stringsAsFactors = FALSE, check.names = FALSE)
  
  if (!"feature" %in% names(all_res)) {
    stop("No 'feature' column found in MaAsLin results.")
  }
  
  if (!"coef" %in% names(all_res)) {
    stop("No 'coef' column found in MaAsLin results.")
  }
  
  q_candidates <- c("qval", "q_val", "qvalue", "q_value", "q")
  qcol <- intersect(q_candidates, names(all_res))[1]
  
  if (is.na(qcol)) {
    qcol <- names(all_res)[grepl("^q", names(all_res), ignore.case = TRUE)][1]
  }
  
  if (is.na(qcol)) {
    p_candidates <- c("pval", "p_val", "pvalue", "p_value", "p")
    pcol <- intersect(p_candidates, names(all_res))[1]
    if (is.na(pcol)) {
      pcol <- names(all_res)[grepl("^p", names(all_res), ignore.case = TRUE)][1]
    }
    
    if (!is.na(pcol)) {
      all_res$qval <- as.numeric(all_res[[pcol]])
    } else {
      all_res$qval <- 1
    }
  } else {
    all_res$qval <- as.numeric(all_res[[qcol]])
  }
  
  all_res %>%
    mutate(
      feature_safe = str_remove(as.character(feature), "^X_"),
      coef = as.numeric(coef),
      qval = ifelse(is.na(qval), 1, qval),
      qval = pmax(qval, 1e-300),
      signed = -log10(qval) * sign(coef),
      sig_label = case_when(
        qval < 0.001 ~ "***",
        qval < 0.01 ~ "**",
        qval < 0.05 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    select(feature_safe, signed, sig_label)
}

count_dups <- function(tbl) {
  if (is.null(tbl)) {
    return(tibble(feature_safe = character(), n = integer()))
  }
  
  tbl %>%
    as_tibble() %>%
    count(feature_safe, name = "n") %>%
    filter(n > 1)
}

dedup_maaslin <- function(tbl) {
  if (is.null(tbl)) {
    return(NULL)
  }
  
  tbl %>%
    as_tibble() %>%
    group_by(feature_safe) %>%
    slice_max(order_by = abs(signed), n = 1, with_ties = FALSE) %>%
    ungroup()
}

res_long_df <- read_maaslin_signed(out_long)
res_over_df <- read_maaslin_signed(out_over)
res_sev_df <- read_maaslin_signed(out_sev)

message("Duplicates in res_long_df:")
print(count_dups(res_long_df))

message("Duplicates in res_over_df:")
print(count_dups(res_over_df))

message("Duplicates in res_sev_df:")
print(count_dups(res_sev_df))

res_long_df <- dedup_maaslin(res_long_df)
res_over_df <- dedup_maaslin(res_over_df)
res_sev_df <- dedup_maaslin(res_sev_df)

message("After deduplication — rows:")
message("res_long_df: ", ifelse(is.null(res_long_df), 0, nrow(res_long_df)))
message("res_over_df: ", ifelse(is.null(res_over_df), 0, nrow(res_over_df)))
message("res_sev_df: ", ifelse(is.null(res_sev_df), 0, nrow(res_sev_df)))

map_tbl <- species_map %>%
  select(species, species_safe, phylum) %>%
  distinct(species_safe, .keep_all = TRUE)

topN_safe <- c(
  if (!is.null(res_long_df)) res_long_df$feature_safe else character(0),
  if (!is.null(res_over_df)) res_over_df$feature_safe else character(0),
  if (!is.null(res_sev_df)) res_sev_df$feature_safe else character(0)
) %>%
  unique()

combined <- tibble(species_safe = topN_safe) %>%
  left_join(map_tbl, by = "species_safe") %>%
  left_join(
    if (!is.null(res_long_df)) {
      res_long_df %>% rename(signed_Long = signed, sig_Long = sig_label)
    } else {
      tibble(feature_safe = character(), signed_Long = numeric(), sig_Long = character())
    },
    by = c("species_safe" = "feature_safe")
  ) %>%
  left_join(
    if (!is.null(res_over_df)) {
      res_over_df %>% rename(signed_Over = signed, sig_Over = sig_label)
    } else {
      tibble(feature_safe = character(), signed_Over = numeric(), sig_Over = character())
    },
    by = c("species_safe" = "feature_safe")
  ) %>%
  left_join(
    if (!is.null(res_sev_df)) {
      res_sev_df %>% rename(signed_Sev = signed, sig_Sev = sig_label)
    } else {
      tibble(feature_safe = character(), signed_Sev = numeric(), sig_Sev = character())
    },
    by = c("species_safe" = "feature_safe")
  ) %>%
  mutate(
    signed_Long = coalesce(signed_Long, 0),
    signed_Over = coalesce(signed_Over, 0),
    signed_Sev = coalesce(signed_Sev, 0),
    sig_Long = coalesce(sig_Long, ""),
    sig_Over = coalesce(sig_Over, ""),
    sig_Sev = coalesce(sig_Sev, "")
  ) %>%
  distinct(species_safe, .keep_all = TRUE)

phylum_recode <- c(
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
  "Euryarchaeota" = "Methanobacteriota"
)

combined <- combined %>%
  mutate(
    phylum = as.character(phylum),
    phylum = ifelse(phylum %in% names(phylum_recode), phylum_recode[phylum], phylum),
    phylum = if_else(
      str_detect(phylum, regex("^Candidatus", ignore_case = TRUE)) |
        str_detect(phylum, regex("unclassified$", ignore_case = TRUE)) |
        is.na(phylum) |
        phylum == "",
      "Other",
      phylum
    ),
    species_display = if_else(
      is.na(species) | species == "",
      make_display_name(species_safe, drop_make_unique_suffix = TRUE),
      make_display_name(make_safe_name(species), drop_make_unique_suffix = TRUE)
    ),
    Control = -signed_Long,
    Long_COVID = signed_Long,
    Overweight = signed_Over,
    Severe = signed_Sev
  )

hm_df <- combined %>%
  select(
    species_safe,
    species_display,
    phylum,
    Control,
    Long_COVID,
    Overweight,
    Severe,
    sig_Long,
    sig_Over,
    sig_Sev
  ) %>%
  arrange(phylum, species_display)

cols <- c("Control", "Long_COVID", "Overweight", "Severe")

mat <- as.matrix(hm_df[, cols, drop = FALSE])
rownames(mat) <- hm_df$species_display
mode(mat) <- "numeric"

plusminus_mat <- apply(
  mat,
  c(1, 2),
  function(x) ifelse(x > 0, "+", ifelse(x < 0, "-", ""))
)

sig_mat <- matrix(
  "",
  nrow = nrow(hm_df),
  ncol = length(cols),
  dimnames = list(hm_df$species_display, cols)
)

for (i in seq_len(nrow(hm_df))) {
  rn <- hm_df$species_display[i]
  sig_mat[rn, "Long_COVID"] <- hm_df$sig_Long[i]
  sig_mat[rn, "Overweight"] <- hm_df$sig_Over[i]
  sig_mat[rn, "Severe"] <- hm_df$sig_Sev[i]
  sig_mat[rn, "Control"] <- ""
}

message("Heatmap rows (species): ", nrow(hm_df))

layer_fun_text <- function(j, i, x, y, width, height, fill) {
  rn <- rownames(mat)[i]
  cn <- colnames(mat)[j]
  
  pm_texts <- plusminus_mat[cbind(rn, cn)]
  star_texts <- sig_mat[cbind(rn, cn)]
  
  grid.text(
    pm_texts,
    x = x,
    y = y,
    gp = gpar(fontsize = ifelse(nrow(hm_df) > 80, 8, 10), fontface = "bold")
  )
  
  idx_star <- which(star_texts != "")
  if (length(idx_star) > 0) {
    grid.text(
      star_texts[idx_star],
      x = x[idx_star] + unit(0.06, "npc"),
      y = y[idx_star],
      gp = gpar(fontsize = ifelse(nrow(hm_df) > 80, 7, 9), col = "black")
    )
  }
}

lim <- max(abs(mat), na.rm = TRUE)
if (is.na(lim) || lim == 0) {
  lim <- 1
}

col_fun <- colorRamp2(c(-lim, 0, lim), c("#4575b4", "white", "#d73027"))

present_phyla <- unique(hm_df$phylum)

phylum_color_map <- c(
  "Bacillota" = "#E69F00",
  "Bacteroidota" = "#009E73",
  "Actinomycetota" = "#332288",
  "Pseudomonadota" = "#F0E442",
  "Lentisphaerota" = "#CC79A7",
  "Elusimicrobiota" = "#999933",
  "Fusobacteriota" = "#0072B2",
  "Synergistota" = "#117733",
  "Mycoplasmatota" = "#56B4E9",
  "Verrucomicrobiota" = "#D55E00",
  "Methanobacteriota" = "#882255",
  "Other" = "grey80"
)

missing_phyla <- setdiff(present_phyla, names(phylum_color_map))
if (length(missing_phyla) > 0) {
  extra_cols <- setNames(rep("grey80", length(missing_phyla)), missing_phyla)
  phylum_color_map <- c(phylum_color_map, extra_cols)
}

phylum_color_map <- phylum_color_map[present_phyla]

row_ha <- rowAnnotation(
  Phylum = anno_simple(
    as.character(hm_df$phylum),
    col = phylum_color_map[as.character(hm_df$phylum)],
    border = FALSE
  ),
  width = unit(6, "mm")
)

coef_legend <- Legend(
  title = expression(-log[10](FDR) * sign(coef)),
  col_fun = col_fun,
  direction = "vertical"
)

phylum_legend <- Legend(
  title = "Phylum",
  labels = names(phylum_color_map),
  legend_gp = gpar(fill = phylum_color_map, col = "black"),
  direction = "vertical",
  ncol = 1
)

ht <- Heatmap(
  mat,
  name = "signed",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = factor(hm_df$phylum, levels = present_phyla),
  row_title = rep("", length(present_phyla)),
  row_title_gp = gpar(fontsize = 0),
  row_names_side = "right",
  row_names_gp = gpar(fontsize = ifelse(nrow(hm_df) > 80, 7, 9)),
  column_names_gp = gpar(fontsize = 11),
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  left_annotation = row_ha,
  show_heatmap_legend = FALSE,
  layer_fun = layer_fun_text
)

draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  heatmap_legend_list = list(coef_legend, phylum_legend),
  padding = unit(c(2, 8, 2, 2), "mm")
)



# Top 50 differentially abundant species
top_n_taxa <- 50

hm_df_top <- combined %>%
  mutate(rank_effect = abs(Long_COVID)) %>%
  arrange(desc(rank_effect), phylum, species_display) %>%
  slice_head(n = top_n_taxa) %>%
  select(
    species_safe,
    species_display,
    phylum,
    Control,
    Long_COVID,
    Overweight,
    Severe,
    sig_Long,
    sig_Over,
    sig_Sev
  ) %>%
  arrange(phylum, species_display)

cols <- c("Control", "Long_COVID", "Overweight", "Severe")

mat <- as.matrix(hm_df_top[, cols, drop = FALSE])
rownames(mat) <- hm_df_top$species_display
mode(mat) <- "numeric"

plusminus_mat <- apply(
  mat,
  c(1, 2),
  function(x) ifelse(x > 0, "+", ifelse(x < 0, "-", ""))
)

sig_mat <- matrix(
  "",
  nrow = nrow(hm_df_top),
  ncol = length(cols),
  dimnames = list(hm_df_top$species_display, cols)
)

for (i in seq_len(nrow(hm_df_top))) {
  rn <- hm_df_top$species_display[i]
  sig_mat[rn, "Long_COVID"] <- hm_df_top$sig_Long[i]
  sig_mat[rn, "Overweight"] <- hm_df_top$sig_Over[i]
  sig_mat[rn, "Severe"] <- hm_df_top$sig_Sev[i]
  sig_mat[rn, "Control"] <- ""
}

layer_fun_text <- function(j, i, x, y, width, height, fill) {
  rn <- rownames(mat)[i]
  cn <- colnames(mat)[j]
  
  pm_texts <- plusminus_mat[cbind(rn, cn)]
  star_texts <- sig_mat[cbind(rn, cn)]
  
  grid.text(
    pm_texts,
    x = x,
    y = y,
    gp = gpar(fontsize = 10, fontface = "bold")
  )
  
  idx_star <- which(star_texts != "")
  if (length(idx_star) > 0) {
    grid.text(
      star_texts[idx_star],
      x = x[idx_star] + unit(0.06, "npc"),
      y = y[idx_star],
      gp = gpar(fontsize = 9, col = "black")
    )
  }
}

lim <- max(abs(mat), na.rm = TRUE)
if (is.na(lim) || lim == 0) lim <- 1

col_fun <- colorRamp2(c(-lim, 0, lim), c("#4575b4", "white", "#d73027"))

present_phyla <- unique(hm_df_top$phylum)

missing_phyla <- setdiff(present_phyla, names(phylum_color_map))
if (length(missing_phyla) > 0) {
  phylum_color_map <- c(
    phylum_color_map,
    setNames(rep("grey80", length(missing_phyla)), missing_phyla)
  )
}

phylum_color_map_top <- phylum_color_map[present_phyla]

row_ha <- rowAnnotation(
  Phylum = anno_simple(
    as.character(hm_df_top$phylum),
    col = phylum_color_map_top[as.character(hm_df_top$phylum)],
    border = FALSE
  ),
  width = unit(6, "mm")
)

coef_legend <- Legend(
  title = expression(-log[10](FDR) * sign(coef)),
  col_fun = col_fun,
  direction = "vertical"
)

phylum_legend <- Legend(
  title = "Phylum",
  labels = names(phylum_color_map_top),
  legend_gp = gpar(fill = phylum_color_map_top, col = "black"),
  direction = "vertical",
  ncol = 1
)

ht <- Heatmap(
  mat,
  name = "signed",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = factor(hm_df_top$phylum, levels = present_phyla),
  row_title = rep("", length(present_phyla)),
  row_title_gp = gpar(fontsize = 0),
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 11),
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  left_annotation = row_ha,
  show_heatmap_legend = FALSE,
  layer_fun = layer_fun_text
)

draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  heatmap_legend_list = list(coef_legend, phylum_legend),
  padding = unit(c(2, 8, 2, 2), "mm")
)
