.libPaths(c("/gpfs/home/bsc/bsc110937/Rlibs", .libPaths()))
library(ggplot2)
library(dplyr)
library(reshape)
library(yaml)
library(tidyr)

tissues <- c("Vagina", "Uterus", "Ovary", "Endocervix", "Ectocervix", "FallopianTube", "BreastMammaryTissue")
input<- "Desktop/TFM/bsc83671/GTEx_v8/"


covariate_colors <- c(
  "HardyScale"   = "#6DB6FFFF",  # HS
  "IschemicTime" = "#FFB6DBFF",  # IT
  "age"          = "#FF6DB6FF",  # AGE
  "ancestry"     = "#009292FF",  # AN
  "bmi"          = "#DB6D00FF"   # BMI
)




yaml_file <- paste0(input, "Allal/RNAPath/resources/clusters.yaml")
yaml_data <- yaml.load_file(yaml_file)

# View the structure of the data to ensure it is read correctly
str(yaml_data)
results_folder <- paste0(input, "Allal/Differential_expression/variance_partition/results/mean_tiles")

color_palette <- function(tissue) {
  # Check if the tissue exists in the yaml_data
  if (!tissue %in% names(yaml_data)) {
    stop("Tissue not found: ", tissue)
  }
  
  # Access the color definitions for the tissue
  color_defs <- yaml_data[[tissue]]$colors
  
  # Convert RGBA vectors to hexadecimal colors
  hex_colors <- lapply(color_defs, function(rgba) {
    if (length(rgba) != 4) {
      warning("Invalid color format for tissue component. Expected 4 values (RGBA).")
      return(NA)
    }
    rgb(rgba[1], rgba[2], rgba[3], rgba[4], maxColorValue = 255)
  })
  
  # Return as named vector (works better with ggplot)
  return(unlist(hex_colors))
}


load_results <- function(tissue,subtissue) {

    results_path <- paste0(results_folder, "/",tissue, "/", subtissue, ".rds")
    results <- readRDS(results_path)
    return(results)
}


# Initialize an empty list to store the results
feature_list <- list()

# Define the covariates of interest
covariates_of_interest <- c("bmi", "age", "ancestry")

for (tissue in tissues) {
  # List all .rds files in the results folder for the given tissue
  files <- list.files(paste0(results_folder, "/", tissue), pattern = ".rds$", full.names = FALSE)
  
  # Extract subtissue names by removing the suffix
  subtissues <- gsub("_results_glm.rds", "", files)
  subtissues <- subtissues[!grepl("_summary", subtissues)]
  subtissues <- gsub(".rds", "", subtissues)
  subtissues <- subtissues[subtissues != "gynecomastoid_hyperplasia"]

  
  # Initialize an empty data frame to store the proportion explained by each covariate
  covariate_proportions <- data.frame(matrix(ncol = 0, nrow = length(subtissues)))
  rownames(covariate_proportions) <- subtissues
  
  # Loop through subtissues
  for (subtissue in subtissues) {
    # Load the results for the given tissue and subtissue
    results <- load_results(tissue, subtissue)
    results$Residuals <- NULL
    
    # Identify the covariates by excluding feature rows
    covariates <- colnames(results)[!colnames(results) %in% rownames(results)]
    
    # Filter covariates to include only BMI, age, and ancestry
    covariates <- covariates[covariates %in% covariates_of_interest]
    
    # Initialize a temporary data frame for the current subtissue
    current_subtissue_proportions <- numeric(length(covariates))
    names(current_subtissue_proportions) <- covariates
    
    # Calculate the total variance explained by the covariates of interest
    total_variance_of_interest <- 0
    for (covariate in covariates) {
      covariate_values <- results[[covariate]]
      total_variance_of_interest <- total_variance_of_interest + sum(covariate_values, na.rm = TRUE)
    }
    
    # Loop through each covariate and calculate the proportion of variance explained
    for (covariate in covariates) {
      # Extract the current covariate's values (column)
      covariate_values <- results[[covariate]]
      
      # Calculate the proportion of variance explained by the current covariate
      covariate_value <- sum(covariate_values, na.rm = TRUE)  # Summing the covariate values
      covariate_proportion <- (covariate_value / total_variance_of_interest) * 100
      
      # Store the proportion for this covariate in the current subtissue
      current_subtissue_proportions[covariate] <- covariate_proportion
    }
    
    # Store the results for the current subtissue in the main data frame
    covariate_proportions[subtissue, covariates] <- current_subtissue_proportions
  }
  
  # Store the results in the feature_list for each tissue
  feature_list[[tissue]] <- covariate_proportions
}

# Combine all results into a single data frame
all_results <- do.call(rbind, feature_list)

# Print or save the final result
print(all_results)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# Define the desired tissue order
tissue_order <- rev(c("BreastMammaryTissue", "Vagina", "Ovary", "Uterus", 
                  "FallopianTube", "Endocervix", "Ectocervix"))

# First, order the Tissue factor levels based on the age column in all_results
ordered_tissues <- all_results %>%
  rownames_to_column(var = "Tissue") %>%
  arrange(age) %>%
  pull(Tissue)

# Reshape the data into long format, excluding the "Residuals" column
long_results <- all_results %>%
  rownames_to_column(var = "Tissue") %>%
  gather(key = "Covariate", value = "Proportion", -Tissue) %>%
  filter(Covariate != "Residuals")  # Exclude the "Residuals" column

# Extract TissueGroup from Tissue (assuming 'TissueGroup.Subtype' format)
long_results <- long_results %>%
  mutate(TissueGroup = sub("\\..*", "", Tissue)) %>%
  arrange(desc(TissueGroup), desc(Tissue))  # Extracts TissueGroup before the dot

# Convert TissueGroup to an ordered factor based on tissue_order
long_results$TissueGroup <- factor(long_results$TissueGroup, levels = tissue_order, ordered = TRUE)

# Sort long_results by TissueGroup order first, then Tissue alphabetically
long_results <- long_results %>%
  arrange(TissueGroup, Tissue)

# Convert Tissue to a factor with levels in the new sorted order
long_results$Tissue <- factor(long_results$Tissue, levels = unique(long_results$Tissue))

# Create horizontal stacked bar plots
p <- ggplot(long_results, aes(x = Tissue, y = Proportion, fill = Covariate)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bars
  coord_flip() +  # Flip coordinates for horizontal bars
  theme_classic() +  # Classic theme
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    legend.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 11, color = "black"),
    plot.title = element_text(size = 14, face = "bold", color = "black")
  ) +
  labs(title = "",
       x = "",
       y = "% of Variance Explained")

# Save the plot
ggsave(paste0(input,"/Laura/Figures_plots/varpar_stacked_barplot_by_subtissue.svg"), plot = p, device = "svg")


long_results <- long_results %>%
  mutate(TissueGroup = sub("\\..*", "", Tissue))  # Keep only the part before "."

# Aggregate proportions by averaging within each tissue group
aggregated_results <- long_results %>%
  group_by(TissueGroup, Covariate) %>%
  summarise(Proportion = mean(Proportion, na.rm = TRUE), .groups = "drop")

# Extract the proportion of variance explained by age
age_order <- aggregated_results %>%
  filter(Covariate == "age") %>%
  arrange(Proportion) %>%  # Sort ascending
  pull(TissueGroup)  # Extract ordered tissue names

# Convert TissueGroup to a factor with levels sorted by age proportion
aggregated_results <- aggregated_results %>%
  mutate(TissueGroup = factor(TissueGroup, levels = age_order))


p<-ggplot(aggregated_results, aes(x = TissueGroup, y = Proportion, fill = Covariate)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bars
  coord_flip() +  # Flip coordinates for horizontal bars
  theme_classic() +  # Classic theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(
    text = element_text(family = "Arial" , color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    legend.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 11, color = "black"),
    plot.title = element_text(size = 14, face = "bold", color = "black")
  ) +  
    scale_fill_manual(values = covariate_colors) +  # Apply custom colors
   # Angle X-axis labels for readability
  labs(title = "",
       x = "",
       y = "% of Variance Explained")
ggsave(paste0(input,"/Laura/Figures_plots/varpar_stacked_barplot_by_tissue.svg"),plot = p, device = "svg")

small_Tissues <- c("Ectocervix", "Endocervix", "FallopianTube")

long_results_filtered <- long_results %>%
  filter(!TissueGroup %in% small_Tissues)

# View the first few rows of the filtered dataframe
head(long_results_filtered)



covariate_colors <- c(
  "HardyScale"   = "#6DB6FFFF",  # HS
  "IschemicTime" = "#FFB6DBFF",  # IT
  "age"          = "#FF6DB6FF",  # AGE
  "ancestry"     = "#009292FF",  # AN
  "bmi"          = "#DB6D00FF"   # BMI
)




p2 <- ggplot(long_results_filtered, aes(x = Tissue, y = Proportion, fill = Covariate)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.7) +  # Stacked bars
  coord_flip() +  # Flip coordinates for horizontal bars
  scale_fill_manual(values = covariate_colors) +  # Apply custom colors
  theme_classic() +  # Classic theme
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.text.x = element_text(angle = -90, hjust = 1, size = 14, color = "black"),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 14, color = "black"),
    legend.text = element_blank(),
    legend.title = element_blank(),
    plot.title = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "none"
  ) +
  labs(title = "",
       x = "",
       y = "Mean variance explained 
excluding residuals (%)")

# Save the plot
ggsave(paste0(input,"/Laura/Figures_plots/varpar_stacked_barplot_big_subtissue.svg"), plot = p2, device = "svg", width = 3)

# Save the plot




p3 <- ggplot(long_results, aes(x = Tissue, y = Proportion, fill = Covariate)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.7) +  # Stacked bars
  coord_flip() +  # Flip coordinates for horizontal bars
  scale_fill_manual(values = covariate_colors) +  # Apply custom colors
  theme_classic() +  # Classic theme
  theme(
    text = element_text( color = "black"),
    axis.text.x = element_text(angle = -90, hjust = 1, size = 14, color = "black"),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 14, color = "black"),
    legend.text = element_blank(),
    legend.title = element_blank(),
    plot.title = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "none"
  ) +
  labs(title = "",
       x = "",
       y = "Mean variance explained 
excluding residuals (%)")

# Save the plot
ggsave(paste0(input,"/Laura/Figures_plots/varpar_stacked_barplot_all_subtissue.svg"), plot = p3, device = "svg", width = 2.5)
ggsave(paste0(input,"/Laura/Figures_plots/varpar_stacked_barplot_all_subtissue.svg"), plot = p3, device = "svg", width = 2.5)
pdf("Desktop/TFM/bsc83671/GTEx_v8/Laura/Figure_plots/figS2_I.pdf", width = 2.5) 
p3
dev.off()

svg("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/Figure_plots/figS2_I.svg", width = 2.5, pointsize = 12)
p3
dev.off()

pdf("Desktop/Figures/figS2_I.pdf", width = 2.5) 
p3
dev.off()

svg("Desktop/Figures/figS2_I.svg", width = 2.5)
p3
dev.off()




