library(dplyr)
library(vroom)
library(stringr)


tissue_colors <- c(
  "Uterus" = "#6DB6FFFF",
  "Ovary" = "#009292FF",
  "Vagina" = "#B66DFFFF",
  "BreastMammaryTissue" = "#490092FF",
  "Endocervix" = "#FF6DB6FF",
  "Ectocervix" = "#FFB6DBFF",
  "FallopianTube" = "#DB6D00FF"
)



# Load donor IDs that only have myometrium
only_myometrium_id <- readRDS('/X/Laura/03.Image_processing/Second_filtering_images/Uterus_filtered_myometrium_images.rds')$Subject.ID

# Define tissues excluding Uterus
tissues <- c("Vagina", "Ovary", "Endocervix", "Ectocervix", "FallopianTube", "BreastMammaryTissue")

# Process non-Uterus tissues
results <- bind_rows(lapply(tissues, function(tissue) {
  tile_features <- vroom(paste0("/X/Allal/RNAPath/tile_features_balanced_classes_1000/balanced_tile_features_", tissue, "_1000.csv"))

  tile_features %>%
    count(label) %>%
    mutate(tissue = tissue) %>%
    arrange(desc(n))
}))

# Process Uterus separately with donor filtering
tile_features_uterus <- vroom("/X/Allal/RNAPath/tile_features_balanced_classes_1000/balanced_tile_features_Uterus_1000.csv")

# Remove endometrium rows for donors that only have myometrium
tile_features_uterus_filtered <- tile_features_uterus %>%
  filter(!(str_extract(slide_name, "GTEX-\\w+") %in% only_myometrium_id & label == "endometrium")) %>%
  count(label) %>%
  mutate(tissue = "Uterus") %>%
  arrange(desc(n))

# Concatenate Uterus results with other tissues
final_results <- bind_rows(results, tile_features_uterus_filtered)

print(final_results)



tissue_order <- c("BreastMammaryTissue", "Vagina", "Ovary", "Uterus", 
             "FallopianTube", "Endocervix", "Ectocervix")



# Ensure tissue is a factor with the specified order
final_results <- final_results %>%
  mutate(tissue = factor(tissue, levels = tissue_order))

# Ensure label is treated as a character for proper alphabetical sorting
final_results <- final_results %>%
  mutate(label = as.character(label))

# Arrange first by tissue (using the factor order), then by label alphabetically
final_results <- final_results %>%
  arrange(tissue, label)

# Print the final_results dataframe
print(final_results)


final_results <- final_results %>%
  filter(!(label == "gynecomastoid_hyperplasia" & tissue == "BreastMammaryTissue"))

tissue_order <- readRDS('/X/Allal/Differential_expression/old_varpar/tissue_order_violin.rds')




library(ggplot2)
library(dplyr)

final_results <- final_results %>%
  mutate(label_tissue = paste(tissue, label, sep = ".")) %>%  # Match the format of tissue_order
  mutate(label_tissue = factor(label_tissue, levels = tissue_order)) %>%  # Set factor levels to match the order
  arrange(label_tissue)  # Arrange the dataframe by the new order
# Create the plot
notiles_all_tissues_plot <- ggplot(final_results, aes(y = label_tissue, x = n, fill = tissue)) +
  geom_bar(stat = "identity") +  # Create a barplot with the unique label-tissue pairs
  scale_fill_manual(values = tissue_colors) +  # Apply custom colors for each tissue
  theme_classic() +
  labs(
    title = "",
    x = "No. tiles",
    y = ""
  ) +
scale_y_discrete(labels = function(x) sub("^.*\\.", "", x)) +  # Remove everything before and including the first dot
  theme(
    text = element_text(family = "Arial", size = 18, color = "black"),  
    axis.text.y = element_text(angle = -45, hjust = 1, vjust = 1, size = 20, color = "black", family = "Arial"),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18, color = "black", family = "Arial"),
    legend.position = "none",  # Remove the legend
    plot.margin = margin(t = 50, r = 50, b = 50, l = 10)  # Increase right and bottom margins
  )

# Save the plot
ggsave("no_tiles_all_tissues.svg", plot = notiles_all_tissues_plot, device = "svg", width = 4.5, height = 11)


write.csv(final_results, "no_tiles_all_subtissues.csv", row.names = FALSE)






small_tissues <- c("Ectocervix", "Endocervix", "FallopianTube")

final_results_big <- final_results %>%
  filter(!tissue %in% small_tissues)



# Create the plot
notiles_big_tissues_plot <- ggplot(final_results_big, aes(y = label_tissue, x = n, fill = tissue)) +
  geom_bar(stat = "identity") +  # Create a barplot with the unique label-tissue pairs
  scale_fill_manual(values = tissue_colors) +  # Apply custom colors for each tissue
  theme_classic() +
  labs(
    title = "",
    x = "No. tiles",
    y = ""
  ) +
  scale_y_discrete(labels = function(x) sub(" - .*", "", x)) +  # Remove the tissue part from y-axis labels
  theme(
    text = element_text(family = "Arial", size = 18, color = "black"),  
    axis.text.y = element_text(angle = -45, hjust = 1, vjust = 1, size = 20, color = "black", family = "Arial"),  
    axis.text.x = element_text(angle = -90, hjust = 1, size = 18, color = "black", family = "Arial"),
    legend.position = "none",  # Remove the legend
    plot.margin = margin(t = 50, r = 50, b = 50, l = 10)  # Increase right and bottom margins
  )

# Save the plot
ggsave("no_tiles_big_tissues.svg", plot = notiles_big_tissues_plot, device = "svg", width = 4, height = 8)
