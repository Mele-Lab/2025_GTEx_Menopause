###### executed in hpc
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
shhh(library(Seurat))
shhh(library(SeuratDisk))
shhh(library(SingleCellExperiment))
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(ggplot2))
shhh(library(dreamlet))
shhh(library(zenith))
shhh(library(scater))
shhh(library(RColorBrewer))

# ---- Read input data ----
sc_data <- readRDS("/X/Data/scRNAseq/FT2025_menopause/converted_data_28b2cc25-efd7-4232-93ef-0e3153c7b218.rds")
cell_metadata <- as.data.frame(colData(sc_data))

# ---- Create donor-level metadata ----
cell_metadata$donor <- as.character(cell_metadata$donor_id)
mdata_donor <- cell_metadata[!duplicated(cell_metadata$donor), ]

# Parse numeric age from development_stage
mdata_donor$Age <- as.numeric(str_extract(mdata_donor$development_stage, "\\d+"))

# Create contraceptive_usage variable and update menstrual_phase_at_collection (those only available for premeno)
mdata_donor_clean <- mdata_donor %>%
  mutate(
    contraceptive_usage = ifelse(
      donor_menopausal_status == "premenopausal" & menstrual_phase_at_collection == "not applicable",
      1,
      0
    ),
    menstrual_phase_at_collection = na_if(menstrual_phase_at_collection, "not applicable")
  )

# Optional: sanity check required fields
required_vars <- c("donor", "Age", "self_reported_ethnicity", "donor_menopausal_status", "menstrual_phase_at_collection", "contraceptive_usage", "donor_BMI_at_collection")
stopifnot(all(required_vars %in% colnames(mdata_donor_clean)))

# Part 1: function to calculate proportions by a given cell label
calculate_props_by_label <- function(cell_label_col, output_path) {
  # Rename for uniform processing
  props.df <- cell_metadata %>%
    rename(cell_label = !!sym(cell_label_col)) %>%
    group_by(donor, cell_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(donor) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup() %>%
    as.data.frame()

  # Add metadata and dataset label
  props.df$dataset <- "FT_2025_Lengyel"

  props.df <- props.df %>%
    left_join(mdata_donor_clean, by = "donor")

  # Select columns of interest
  final_cols <- c("donor", "cell_label", "n", "freq", "dataset",
                  "Age", "self_reported_ethnicity", "donor_menopausal_status")
  final_cols <- final_cols[final_cols %in% colnames(props.df)]
  props.df <- props.df %>% select(all_of(final_cols))

  # Save
  saveRDS(props.df, output_path)
}

calculate_props_by_label("cell_type", "/X/Data/scRNAseq/FT2025_menopause/temp_replication/proportions_allDonors_celltype.rds")

#### Part 2) select the proportion file to read
props.df <- readRDS("/X/Data/scRNAseq/FT2025_menopause/temp_replication/proportions_allDonors_celltype.rds")

# Rename 'cell_label' to 'cell_type' temporarily for downstream compatibility
props.df <- props.df %>%
  rename(cell_type = cell_label)

# ---- Summarize cell type stats ----
celltype.stats.df <- props.df %>%
  group_by(cell_type, .drop = FALSE) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq)) %>%
  mutate(celltype.label = paste0(cell_type, ' (n=', n, '; freq=', round(freq, 2), ')')) %>%
  as.data.frame()

# ---- CLR transform function ----
Geom_mean <- function(x) exp(mean(log(x)))
CLR <- function(D) log2(D / Geom_mean(D))

clr_func <- function(df) {
  wide_df <- reshape2::dcast(df, donor ~ cell_type, value.var = "freq")
  rownames(wide_df) <- wide_df$donor
  wide_df$donor <- NULL
  wide_df[is.na(wide_df)] <- 0
  
  positive_vals <- wide_df[wide_df > 0]
  if (length(positive_vals) == 0) stop("All values are zero — cannot compute pseudocount")
  pseudocount <- 1/3 * min(positive_vals)
  
  clr_matrix <- apply(wide_df + pseudocount, 2, CLR)
  return(t(clr_matrix))  # rows = cell types
}

# ---- Apply CLR transformation ----
props_clr.df <- clr_func(props.df)

# --- Model function per cell type: here only model 2 (age only, adjusted for ethnicity) ---

lm_by_ct_age <- function(cell_type, df = props_clr.df, md = mdata_donor) {
  message(paste0("Modeling cell type: ", cell_type))
  
  vec_i <- df[rownames(df) == cell_type, ]
  df_i <- data.frame(freq = as.numeric(vec_i), donor = colnames(df))
  df_i <- merge(df_i, md, by = "donor")
  
  # Format covariates
  df_i$Age <- as.numeric(df_i$Age)
  df_i$self_reported_ethnicity <- factor(df_i$self_reported_ethnicity)
  
  # --- Model 2: freq ~ Age + self_reported_ethnicity ---
  mod2 <- lm(freq ~ Age + self_reported_ethnicity, data = df_i)
  tidy2 <- broom::tidy(mod2, conf.int = TRUE) %>%
    mutate(celltype = cell_type, model = "age_only")
  
  return(tidy2)
}

# Run model 2 for all cell types
cell_types <- unique(props.df$cell_type[!is.na(props.df$cell_type)])
tidy_mod.list <- lapply(cell_types, lm_by_ct_age, df = props_clr.df, md = mdata_donor)
tidy_mod.df <- do.call(rbind, tidy_mod.list)

# Adjust p-values for each term (model = "age_only" only)
tidy_mod.df <- tidy_mod.df %>%
  group_by(term, model) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  ungroup()

# Extract age effect results
age_results <- tidy_mod.df %>%
  filter(model == "age_only", term == "Age") %>%
  mutate(direction = ifelse(estimate > 0, "pos", "neg"))

# Add celltype labels
age_results <- merge(age_results, celltype.stats.df[, c("cell_type", "celltype.label")], 
                     by.x = "celltype", by.y = "cell_type", all.x = TRUE)

# Save results
saveRDS(age_results, "/X/Data/scRNAseq/FT2025_menopause/temp_replication/FT_CoDA_celltype_results_age_only.rds")
