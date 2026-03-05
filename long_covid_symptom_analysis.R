# long_covid_symptom_analysis

required_pkgs <- c("readxl", "dplyr", "tidyr", "stringr", "ggplot2", "forcats", "readr", "janitor")
to_install <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if(length(to_install)) install.packages(to_install)

library(readxl); library(dplyr); library(tidyr); library(stringr)
library(ggplot2); library(forcats); library(readr); library(janitor)


metadata_file_path <- "/Users/louisburger22/Documents/African_Microbiome_Group/Long_COVID_study/Long_COVID_metadata.xlsx"
metadata_sheet_name <- "Sheet1"
plot_top_n <- 0                
ordering_method <- "overall"   


output_directory <- dirname(metadata_file_path)
output_summary_csv <- file.path(output_directory, "long_covid_symptom_summary_fixed.csv")
output_plot_png <- file.path(output_directory, "long_covid_symptoms_barplot_fixed.png")

message("Reading metadata from: ", metadata_file_path)
metadata_table <- readxl::read_excel(metadata_file_path, sheet = metadata_sheet_name)
message("Metadata shape: ", paste0(dim(metadata_table)[1], " rows x ", dim(metadata_table)[2], " cols"))

looks_like_1_to_10 <- function(column_values) {
  numeric_values <- suppressWarnings(as.numeric(as.character(column_values)))
  numeric_values <- numeric_values[!is.na(numeric_values)]
  if(length(numeric_values) < 3) return(FALSE)
  mean((numeric_values >= 1 & numeric_values <= 10), na.rm = TRUE) > 0.5
}
candidate_symptom_columns <- names(metadata_table)[sapply(metadata_table, looks_like_1_to_10)]
message("Detected ", length(candidate_symptom_columns), " candidate symptom/rating columns.")

parse_symptom_column <- function(column_name) {
  parts <- strsplit(column_name, ":")[[1]] %>% trimws()
  symptom_name <- NA_character_
  timepoint_label <- NA_character_
  if(length(parts) >= 3) {
    raw_symptom <- parts[2]
    symptom_name <- str_replace(raw_symptom, "^\\d+\\.?\\s*", "") %>% str_squish()
    raw_timepart <- parts[length(parts)]
    raw_timepart_clean <- raw_timepart %>%
      str_replace_all("\\(.*?\\)", "") %>%
      str_replace_all("1-10|1 10|1–10", "") %>%
      str_squish()
    low_time <- tolower(raw_timepart_clean)
    if(str_detect(low_time, "before")) timepoint_label <- "A"
    if(str_detect(low_time, "3 month") | str_detect(low_time, "3months") | str_detect(low_time, "3 mth")) timepoint_label <- "B"
    if(str_detect(low_time, "current|now|present|ongoing")) timepoint_label <- "C"
  } else {
    symptom_name <- column_name
  }
  list(symptom = symptom_name, timepoint = timepoint_label)
}

parsed_columns <- tibble(original_colname = candidate_symptom_columns) %>%
  rowwise() %>%
  mutate(parsed = list(parse_symptom_column(original_colname)),
         symptom_name = parsed$symptom,
         timepoint_label = parsed$timepoint) %>%
  dplyr::select(-parsed) %>% ungroup()

parsed_columns <- parsed_columns %>%
  mutate(timepoint_label = ifelse(
    is.na(timepoint_label),
    case_when(
      str_detect(tolower(original_colname), "before") ~ "A",
      str_detect(tolower(original_colname), "3 month|3 months|3_months|3months|3 mth|3m") ~ "B",
      str_detect(tolower(original_colname), "current|now|present|ongoing") ~ "C",
      TRUE ~ NA_character_
    ),
    timepoint_label
  ))

if(any(is.na(parsed_columns$timepoint_label))) {
  if(length(candidate_symptom_columns) %% 3 == 0) {
    norm_prefix <- function(name) {
      n <- tolower(name)
      n <- str_replace(n, "[_\\.-]*[123]$", "")
      n <- str_replace(n, "\\(1-10 scale\\)", "")
      n <- str_replace_all(n, "[:]", " ")
      str_squish(n)
    }
    parsed_columns <- parsed_columns %>% mutate(prefix = sapply(original_colname, norm_prefix))
    grouping_info <- parsed_columns %>% group_by(prefix) %>% mutate(count_in_group = n()) %>% ungroup()
    auto_assign <- grouping_info %>% filter(count_in_group == 3) %>%
      arrange(prefix, original_colname) %>%
      group_by(prefix) %>%
      mutate(timepoint_label = c("A","B","C")) %>%
      ungroup() %>%
      select(original_colname, timepoint_label)
    parsed_columns <- parsed_columns %>%
      select(-prefix) %>%
      left_join(auto_assign, by = "original_colname", suffix = c("", ".assigned")) %>%
      mutate(timepoint_label = coalesce(timepoint_label.assigned, timepoint_label)) %>%
      select(-timepoint_label.assigned)
  }
}

unmapped_columns <- parsed_columns %>% filter(is.na(timepoint_label))
if(nrow(unmapped_columns) > 0) {
  message("WARNING: The following candidate symptom columns could not be automatically assigned to a timepoint (they will be excluded from plotting):")
  print(unmapped_columns$original_colname)
  message("If any of these are symptom columns, either rename them in the spreadsheet or add a manual mapping in the script.")
}

symptom_columns_for_pivot <- parsed_columns %>% filter(!is.na(timepoint_label)) %>% pull(original_colname)

if(length(symptom_columns_for_pivot) == 0) {
  stop("No symptom columns were mapped to a timepoint. Check column names or add manual mappings.")
} else {
  message("Pivoting ", length(symptom_columns_for_pivot), " symptom columns (those mapped to a timepoint).")
}

metadata_table[symptom_columns_for_pivot] <- lapply(metadata_table[symptom_columns_for_pivot], function(x) as.character(x))

symptom_long_form <- metadata_table %>%
  mutate(.row_id = row_number()) %>%
  pivot_longer(cols = all_of(symptom_columns_for_pivot),
               names_to = "original_colname",
               values_to = "rating_value") %>%
  left_join(parsed_columns, by = "original_colname") %>%
  mutate(rating_value_numeric = suppressWarnings(as.numeric(rating_value))) %>%
  dplyr::select(.row_id, original_colname, symptom_name, timepoint_label, rating_value_numeric)

symptom_summary_table <- symptom_long_form %>%
  filter(!is.na(timepoint_label), !is.na(symptom_name)) %>%
  group_by(symptom_name, timepoint_label) %>%
  summarise(
    mean_severity = mean(rating_value_numeric, na.rm = TRUE),
    median_severity = median(rating_value_numeric, na.rm = TRUE),
    count_non_missing = sum(!is.na(rating_value_numeric)),
    .groups = "drop"
  )

symptom_overall_means <- symptom_summary_table %>%
  group_by(symptom_name) %>%
  summarise(
    overall_mean = mean(mean_severity, na.rm = TRUE),
    current_mean = mean_severity[timepoint_label == "Current"] %>% { if(length(.)==0) NA_real_ else .[1] },
    .groups = "drop"
  )

symptom_summary_full <- symptom_summary_table %>% left_join(symptom_overall_means, by = "symptom_name")

if(ordering_method == "A" && any(!is.na(symptom_overall_means$current_mean))) {
  ordered_symptoms <- symptom_overall_means %>% arrange(desc(current_mean)) %>% pull(symptom_name)
} else {
  ordered_symptoms <- symptom_overall_means %>% arrange(desc(overall_mean)) %>% pull(symptom_name)
}
if(plot_top_n > 0) ordered_symptoms <- head(ordered_symptoms, plot_top_n)

plotting_table <- symptom_summary_full %>%
  filter(symptom_name %in% ordered_symptoms) %>%
  mutate(symptom_name = factor(symptom_name, levels = rev(ordered_symptoms)),
         timepoint_label = factor(timepoint_label, levels = c("A", "B", "C")))

if(nrow(plotting_table) == 0) stop("No symptoms available for plotting after parsing.")

plot_height <- max(4, length(unique(plotting_table$symptom_name)) * 0.35)

ggplot(plotting_table, aes(x = symptom_name, y = mean_severity, fill = timepoint_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", size = 0.15) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
  scale_fill_manual(values = c(
    "A"    = "#6BAE60",  
    "B"  = "#C44E52",  
    "C"   = "#E6A157"   
  )) +
  labs(x = "Symptoms", y = "Mean symptom severity (1-10)", fill = "Timepoint",
       title = "") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.ticks.y = element_line(color = "black", linetype = "dashed", linewidth = 0.4),
    axis.ticks.x = element_line(color = "black", linetype = "dashed", linewidth = 0.4),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# Stats

alpha <- 0.05
dot_size <- 3.2
asterisk_size <- 6    
outline_linewidth <- 0.6
asterisk_color <- "#D62828"  

pairwise_results <- pairwise_results %>% mutate(p_value = coalesce(p_adj, p))

if(!exists("plotting_table")) stop("plotting_table must exist (mean severity per symptom/timepoint).")
symptom_order <- plotting_table %>%
  group_by(symptom_name) %>%
  summarise(overall_mean = mean(mean_severity, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(overall_mean)) %>%
  pull(symptom_name)

pairwise_norm <- pairwise_results %>%
  mutate(gmin = pmin(as.character(group1), as.character(group2)),
         gmax = pmax(as.character(group1), as.character(group2)),
         pair_key = paste(gmin, gmax, sep = "__")) %>%
  dplyr::select(symptom_name, pair_key, p_value)

pairs <- tibble::tribble(
  ~label,               ~g1,        ~g2,
  "A_vs_B",  "A",   "B",
  "A_vs_C",  "A",   "C",
  "B_vs_C", "B", "C"
) %>%
  mutate(pair_key = paste(pmin(g1,g2), pmax(g1,g2), sep = "__"))

p_df <- expand.grid(symptom_name = symptom_order, label = pairs$label, stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  left_join(pairs, by = "label") %>%
  mutate(symptom_name = factor(symptom_name, levels = rev(symptom_order))) %>%   
  left_join(pairwise_norm, by = c("symptom_name","pair_key")) %>%
  mutate(is_sig = !is.na(p_value) & p_value < alpha)

build_panel <- function(pair_label, panel_title = "") {
  df <- p_df %>% filter(label == pair_label)
  df_plot <- df %>% filter(!is.na(p_value)) 
  df_sig  <- df_plot %>% filter(is_sig)
  df_ns   <- df_plot %>% filter(!is_sig)
  
  ggplot() +
    geom_point(data = df_ns, aes(x = p_value, y = symptom_name), shape = 16, size = dot_size, color = "black", na.rm = TRUE) +
    geom_text(data = df_sig, aes(x = p_value, y = symptom_name, label = "*"), size = asterisk_size, color = asterisk_color, vjust = 0.45, na.rm = TRUE) +
    geom_vline(xintercept = alpha, linetype = "dashed", color = "grey50") +
    scale_x_continuous(limits = c(0,1), breaks = c(0.05, 1), labels = c("0.05", "1")) +
    labs(x = "p-value", y = NULL, title = panel_title) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = outline_linewidth),
      axis.text.y = element_blank(),    # HIDE symptom labels as requested
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(color = "black", linetype = "dashed", linewidth = 0.4),
      axis.text.x = element_text(size = 10),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )
}

p_before_3m   <- build_panel("A_vs_B",  panel_title = "A vs B")
p_before_curr <- build_panel("A_vs_C",  panel_title = "A vs C")
p_3m_curr     <- build_panel("B_vs_C", panel_title = "B vs C")


combined_three <- p_before_3m + p_before_curr + p_3m_curr + plot_layout(ncol = 3, widths = c(1,1,1))
print(combined_three)