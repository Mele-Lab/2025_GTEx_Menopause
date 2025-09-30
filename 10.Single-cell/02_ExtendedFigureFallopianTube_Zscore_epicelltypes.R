library(tidyverse)
library(dplyr)
library(ggplot2)
library(writexl)

# --- Load data ---
props_all <- readRDS("/X/scRNAseq/FT2025_menopause/proportions_allDonors_celltype.rds")
cell_metadata <- readRDS("/X/scRNAseq/FT2025_menopause/cell_metadata.rds")

# --- Donor-level metadata ---
cell_metadata$donor <- as.character(cell_metadata$donor_id)
metadata <- cell_metadata[!duplicated(cell_metadata$donor), ]

# Parse numeric Age
metadata$Age <- as.numeric(stringr::str_extract(metadata$development_stage, "\\d+"))

# --- Join proportions with metadata ---
props_df <- props_all %>%
  left_join(metadata, by = "donor") %>%
  mutate(freq = as.numeric(freq))

# --- Create donor x cell_type matrix ---
wide_df <- props_df %>%
  select(donor, cell_label, freq) %>%
  pivot_wider(names_from = donor, values_from = freq, values_fill = 0) %>%
  column_to_rownames("cell_label")

wide_mat <- as.matrix(wide_df)

# --- Z-score across donors (row-wise scaling) ---
zscore_matrix <- t(scale(t(wide_mat)))

# --- Order donors by age ---
ordered_donors <- metadata %>%
  filter(donor %in% colnames(zscore_matrix)) %>%
  arrange(Age) %>%
  pull(donor)

zscore_matrix <- zscore_matrix[, ordered_donors]

# --- Convert to long format with age info ---
zscore_long <- as.data.frame(zscore_matrix) %>%
  rownames_to_column(var = "celltype") %>%
  pivot_longer(-celltype, names_to = "donor", values_to = "zscore") %>%
  left_join(metadata %>% select(donor, Age), by = "donor")

# --- Filter for selected epithelial cell types ---
selected_cells <- c("peg cell", 
                    "ciliated epithelial cell", 
                    "fallopian tube secretory epithelial cell")

zscore_filtered <- zscore_long %>%
  filter(celltype %in% selected_cells) %>%
  # Rename the cell type
  mutate(celltype = ifelse(celltype == "fallopian tube secretory epithelial cell", 
                           "secretory epithelial cell", celltype))


# Base orange colors for each cell type

cell_colors <- c(
  "peg cell" = "#E69F00",
  "secretory epithelial cell" = "#FC8D62",
  "ciliated epithelial cell" = "#FDBF6F"
)

# plot harmonized with the theme used
p_filtered_harmonized <- ggplot(zscore_filtered, aes(x = Age, y = zscore)) +
  geom_point(aes(color = celltype), size = 2) +
  geom_smooth(aes(color = celltype), method = "loess", se = FALSE, linewidth = 1) +
  scale_color_manual(values = cell_colors, name = "Cell Type") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Age",
    y = "Z-score"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.y = element_text(size = 15, color = "black"),
    axis.text.x = element_text(size = 15, color = "black", angle = 0, hjust = 0.5),  # horizontal
    axis.title = element_text(size = 15, color = "black"),
    plot.title = element_text(size = 15, hjust = 0.5, color = "black"),
    legend.text = element_text(size = 15, color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey80", size = 0.5, fill = NA),
    legend.position = "right"
  )

p_filtered_harmonized


# --- Save plot ---
ggsave("/X/SC_ANALYSES/replication_final_analyses/02_ExtendedFigureFallopianTube_Zscore_epicelltypes.png",
       p_filtered_harmonized, bg = "white", width = 10, height = 6, dpi = 300)

# --- Save plot (PDF, vector format) ---
ggsave("/X/SC_ANALYSES/replication_final_analyses/02_ExtendedFigureFallopianTube_Zscore_epicelltypes.pdf", p_filtered_harmonized, bg = "white", width = 10, height = 6, device = cairo_pdf)


#below i copied the input data as *.xlsx in case you need to replicate/modify the plot:

write_xlsx(zscore_filtered, 
           "/X/SC_ANALYSES/replication_final_analyses/input_data_for_Fig2E_Zscore_selected_celltypes_age.xlsx")



