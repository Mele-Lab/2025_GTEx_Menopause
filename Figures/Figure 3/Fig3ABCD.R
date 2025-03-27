
##### Functions for calculating the slope changes and plottin
###################
###libraries
library(segmented)  # For davies.test
library(ggplot2)
library(dplyr)
library(tidyr)
library(stats)


# Option 1: Maximum slope using LOESS
find_max_slope_point_loess <- function(age, measurement, span = 0.75) {
  # Ensure data is sorted by age
  ord <- order(age)
  age <- age[ord]
  measurement <- measurement[ord]
  
  # Fit LOESS model
  loess_fit <- loess(measurement ~ age, span = span)
  
  # Create a fine grid for evaluating the derivative
  grid_x <- seq(min(age), max(age), length.out = 200)
  grid_y <- predict(loess_fit, newdata = grid_x)
  
  # Calculate numerical derivative
  slopes <- diff(grid_y) / diff(grid_x)
  
  # Find the point with maximum absolute slope
  max_slope_idx <- which.max(abs(slopes))
  max_slope_point <- grid_x[max_slope_idx]
  
  return(max_slope_point)
}


analyze_organ_trajectories_gen <- function(long_data, 
                                           p_threshold = 0.05, 
                                           slope_threshold = 0.1
) {
  organs = unique(long_data$organ)
  organ_level = ifelse(("Uterus" %in% organs), TRUE, FALSE) # Check if we are at the organ level
  
  # Making sure it is a factor
  long_data$organ = factor(long_data$organ)
  # Initialize an empty named list
  # Finding the age of maximum slope change for the trajectories with significant davies test 
  
  slope_changes <- list()
  for(i in seq_along(unique(long_data$organ))) {
    organ_name <- as.character(unique(long_data$organ)[i])
    
    # Get data for this organ
    organ_data <- long_data %>% 
      filter(organ == organ_name)
    
    # Ensure enough data points
    if(nrow(organ_data) < 5) next
    
    # Create data frame for model fitting
    x_values <- organ_data$age
    y_values <- organ_data$measurement
    model_data <- data.frame(y = y_values, x = x_values)
    
    # Fit a linear model first (required for davies.test)
    linear_model <- lm(y ~ x, data = model_data)
    
    # Perform Davies test for a non-zero change in slope
    davies_result <- davies.test(linear_model)
    print(organ_name)
    print(davies_result)
    # Only find maximum slope point if Davies test is significant
    if (davies_result$p.value < p_threshold) {
      segmented_model <- segmented(linear_model, seg.Z = ~x)
   
      
      # Find the point with maximum slope
      max_slope_point <- find_max_slope_point_loess(x_values, y_values)
  
      print(max_slope_point)
      # Only add to list if point was found
      if (!is.null(max_slope_point)) {
        slope_changes[[organ_name]] <- max_slope_point
      }
    }
  }
  
  # Create data frame for slope change lines (only for organs with changes)
  slope_change_df <- if(length(slope_changes) > 0) {
    do.call(rbind, lapply(names(slope_changes), function(org) {
      data.frame(
        organ = org,
        age = slope_changes[[org]]
      )
    }))
  } else {
    NULL
  }
  
  
  # Create plot with custom colors for reproductive organs
  organ_colors <- c("#490092FF","#6DB6FFFF","#009292FF","#B66DFFFF")
  if (organ_level){
    organ_colors <- c("Uterus"="#6DB6FFFF", "Ovary" = "#009292FF","Vagina" = "#B66DFFFF")
  }
  p <- ggplot() +
    # Add individual points
    geom_point(data = long_data, 
               aes(x = age, y = measurement, color = organ), 
               alpha = 0.2, size = 0.3) +
    # Add smoothed trajectories
    geom_smooth(data = long_data, 
                aes(x = age, y = measurement, color = organ),
                linewidth = 1)
  
  # Add vertical lines only if slope changes were detected
  if(!is.null(slope_change_df) > 0) {
    p <- p + 
      geom_vline(data = slope_change_df,
                 aes(xintercept = age, color = organ),
                 linetype = "dashed",
                 alpha = 0.7)
  }
  ###if it's organ trajectories, write "Organ", if tissue - "Tissue"
  ylabel = ifelse(organ_level, "Probability (being old organ)", "Probability (being old tissue)")
  color_label = ifelse(organ_level, "Organ","Tissue")
  # Customize appearance
  p <- p + ylim(0,1)+
    scale_color_manual(values = organ_colors) +
    labs(x = "Age (years)", 
         y = ylabel,
         color = color_label) +
    #Change font sizes on axes and labels
    theme_minimal(base_size = 15) +
    theme(axis.text.y = element_text(size = 15, color = "black"),
          axis.text.x = element_text(size = 15, color = "black"),
          axis.title = element_text(size = 15, color = "black"),
          plot.title = element_text(size = 15, hjust = 0.5, color = "black"),
          legend.text = element_text(size = 15, color = "black"),
          panel.grid = element_blank(),   
          legend.position = 'bottom',
          legend.direction = "vertical",# Remove all grid lines
          # axis.line.x = element_line(size = 0.5, color = "grey80"), # Add x-axis line
          # axis.line.y = element_line(size = 0.5, color = "grey80"),
          panel.border = element_rect(color = "grey80", size = 0.5, fill = NA)) 
  
  return(list(
    plot = p,
    slope_changes = slope_changes
    #abs_slope_changes = changes_abs
  ))
}

#### 1. Load the data
prfx = paste0(input, "/Ole/paper/for_plotting/")
### CNN data
merged_df = readRDS(paste0(prfx, "CNN_classification.rds"))
### Gene expression data
merged_df = readRDS(paste0(prfx, "GE_classification.rds"))
### Subtissue data
tissues = c("Vagina", "Uterus", "Ovary")
merged_df = readRDS(paste0(prfx, tissues[2],"_subtissues_classification.rds"))

### 2. Analyze the data
# Analyze the CNN data
results = analyze_organ_trajectories_gen(merged_df, p_threshold = 0.05)

# Plot the plot
p<-plot(results$plot)
p <- results$plot + 
  theme(
    legend.title = element_blank(),
    legend.direction = "horizontal",
    # legend.position = c(0.77, 0.15),  # Move legend inside the bottom-right
    # legend.background = element_rect(fill = alpha("white", 0.6), color = NA),  # Semi-transparent background
    # legend.key = element_blank()  # Remove borders around legend items
  )
p


ut_traj<-p
va_traj<-p
ov_traj<-p

trajectories<- ut_traj+va_traj+ov_traj

##Epithelium measurements
library(jsonlite)
library(ggplot2)
library(ggsignif)
library(gridExtra)

input<-"Desktop/TFM/bsc83671/GTEx_v8/"
measures <- read.csv(paste0(input, "/Allal/CellProfiler/results/Vagina/newSecondaryFilteredVagina_Image.csv"))
median_value_columns <- c("FileName_Vagina","Median_FilterObjects_AreaShape_MinorAxisLength", "Median_FilterObjects_AreaShape_MinorAxisLength", "Median_FilterObjects_AreaShape_MinorAxisLength")
image_median_values <- measures[,median_value_columns]
#all_measures <- read.csv("/gpfs/projects/bsc83/Projects/GTEx_v8/Allal/CellProfiler/results/Vagina/newSecondaryFilteredVagina_FilterObjects.csv")
#measures <- read.csv("thickness_measures/Vagina_epithelium_measures.csv")
#metadata <- read.csv("/gpfs/projects/bsc83/Projects/GTEx_v8/Laura/00.Data/GTEx_Subject_Phenotypes.GRU.csv", sep = '\t')
# Replace the suffix with an empty string
image_median_values$FileName_Vagina <- sub("_binary_mask_class_epithelium_compressed\\.jpg$", "", measures$FileName_Vagina)
image_median_values$FileName_Vagina <- sub("^([^-]+-[^-]+)-.*$", "\\1", image_median_values$FileName_Vagina)

tissue<- "Vagina"
metadata<-read.csv("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/00.Data/GTEx_Subject_Phenotypes.GRU.csv", sep="\t", header=TRUE)

filter<-readRDS(paste0(input,"/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))

metadata<- metadata[metadata$SUBJID %in% filter$Subject.ID,]

ancestry_file<- read.table(paste0(input,"/Laura/00.Data/admixture_inferred_ancestry.txt"))
ances<-ancestry_file[,c(1,3)]
colnames(ances)<- c("Donor", "Ancestry")


metadata <- merge(metadata,ancestry_file, by.x = "SUBJID", by.y = "V1", all.x = FALSE)



# Read the GTEx filtered IDs from the text file
#filtered_ids <- readLines("../filtered_ids/Vagina_filtered_ids.txt")

#filtered_ids <- fromJSON(filtered_ids)
# Substitute the last character '6' with '5' for each GTEx ID


# Filter the measures data to include only rows where GTEx_ID is in the modified list
#measures_filtered <- measures[measures$GTEx_ID %in% filtered_ids_modified, ]


# Merge AGE column from metadata into measures_filtered
measures_filtered <- merge(
  image_median_values,
  metadata[, c("SUBJID","AGE","BMI","DTHHRDY","TRISCH","COHORT","V2","V3","V4","V5","V6")],  # Keep only necessary columns from metadata
  by.x = "FileName_Vagina",
  by.y = "SUBJID",
  all.x = FALSE  # Keep all rows from measures_filtered
)

measures_filtered <- na.omit(measures_filtered)


#Check data distribution

# check normal distribution
## If the p-value > 0.05: Fail to reject the null hypothesis, meaning the data is likely normal.
## If the p-value ≤ 0.05: Reject the null hypothesis, meaning the data is not normal.

x_norm <- measures_filtered$Median_FilterObjects_AreaShape_MinorAxisLength

# Shapiro-Wilk test
shapiro_result <- shapiro.test(x_norm)
print(shapiro_result)
if (shapiro_result$p.value > 0.05) {
  print("Shapiro-Wilk Test: The data is likely normally distributed (p > 0.05).")
} else {
  print("Shapiro-Wilk Test: The data is NOT normally distributed (p <= 0.05).")
}

# Kolmogorov-Smirnov test
ks_result <- ks.test(x_norm, "pnorm", mean = mean(x_norm), sd = sd(x_norm))
print(ks_result)
if (ks_result$p.value > 0.05) {
  print("Kolmogorov-Smirnov Test: The data is likely normally distributed (p > 0.05).")
} else {
  print("Kolmogorov-Smirnov Test: The data is NOT normally distributed (p <= 0.05).")
}


# Load required libraries
library(ggplot2)

# Convert to data frame for ggplot
data <- data.frame(x_norm = x_norm)


measures_filtered$Median_FilterObjects_AreaShape_MinorAxisLength <- measures_filtered$Median_FilterObjects_AreaShape_MinorAxisLength



# Convert 'COHORT' to a factor
measures_filtered$COHORT <- factor(measures_filtered$COHORT)

# Create a mapping of factor levels to their corresponding original text values
map_COHORT <- setNames(levels(measures_filtered$COHORT), seq_along(levels(measures_filtered$COHORT)) - 1)


# Convert 'COHORT' to numeric (discretize it) and assign it to the column
measures_filtered$COHORT <- as.integer(measures_filtered$COHORT) - 1

# Extract hours and minutes, then convert to total minutes

measures_filtered$TRISCH <- as.integer(sub("^(\\d+) hour\\(s\\), (\\d+) minute\\(s\\)$", "\\1", measures_filtered$TRISCH)) * 60 + 
  as.integer(sub("^(\\d+) hour\\(s\\), (\\d+) minute\\(s\\)$", "\\2", measures_filtered$TRISCH))

# View the result
head(measures_filtered$TRISCH)

lm_model <- lm(
  Median_FilterObjects_AreaShape_MinorAxisLength ~ AGE + BMI + DTHHRDY + TRISCH + COHORT + 
    V2 ,
  data = measures_filtered
)

# View the model summary
summary(lm_model)
''




# Finer bins
measures_filtered$quinquennial_bins <- cut(
  measures_filtered$AGE,
  breaks = c(20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70),
  labels = c("20-25", "26-30", "31-35", "36-40", "41-45", "46-50", "51-55", "56-60", "61-65", "66-70"),
  right = TRUE
)

# Decades
measures_filtered$decade_bins <- cut(
  measures_filtered$AGE,
  breaks = c(20, 30, 40, 50, 60, 70),
  labels = c("20-30", "31-40", "41-50", "51-60", "61-70"),
  right = TRUE
)





compute_pairwise_comparisons <- function(continuous_variable, levels) {
  
  ttest <- pairwise.wilcox.test(
    continuous_variable,
    levels,
    p.adjust.method = "bonferroni"  # Adjust for multiple comparisons
  )
  return(ttest)
  
}

quinquennial_ttest <- compute_pairwise_comparisons(measures_filtered$Median_FilterObjects_AreaShape_MinorAxisLength, measures_filtered$quinquennial_bins)
decade_ttest <- compute_pairwise_comparisons(measures_filtered$Median_FilterObjects_AreaShape_MinorAxisLength, measures_filtered$decade_bins)




extract_significant_comparisons <- function(ttest_result, alpha = 0.05) {
  # Convert the matrix of p-values into a data frame
  significant_pairs <- as.data.frame(as.table(ttest_result$p.value))
  
  # Filter for significant results (p-value < alpha)
  significant_pairs <- significant_pairs[!is.na(significant_pairs$Freq) & significant_pairs$Freq < alpha, ]
  
  # Ensure the comparisons are returned as a list of character vectors
  comparisons_list <- lapply(1:nrow(significant_pairs), function(i) {
    c(as.character(significant_pairs[i, "Var1"]), as.character(significant_pairs[i, "Var2"]))
  })
  
  return(comparisons_list)
}

# Example usage:
quinquennial_comparisons <- extract_significant_comparisons(quinquennial_ttest)
decade_comparisons <- extract_significant_comparisons(decade_ttest)
# Print the resulting list of significant comparisons


# 1. Initialize a vector to store p-values for each comparison
comparison_pvalues <- vector("numeric", length = length(quinquennial_comparisons))



quinquennial_comparisons <- Filter(
  function(x) !any(is.na(x)) && length(x) == 2,
  quinquennial_comparisons)

# Check if quinquennial_comparisons is empty
if (length(quinquennial_comparisons) == 0) {
  message("No valid comparisons found. Skipping processing.")
} else {
  # Initialize a vector for p-values
  comparison_pvalues <- numeric(length(quinquennial_comparisons))
  
  # Loop through each comparison
  for (i in seq_along(quinquennial_comparisons)) {
    comparison <- quinquennial_comparisons[[i]]
    
    # Check if both elements exist in the p-value matrix
    if (comparison[1] %in% rownames(quinquennial_ttest$p.value) &&
        comparison[2] %in% colnames(quinquennial_ttest$p.value)) {
      # Extract the p-value
      comparison_pvalues[i] <- quinquennial_ttest$p.value[comparison[1], comparison[2]]
    } else {
      # Assign NA for missing comparisons
      comparison_pvalues[i] <- NA
    }
  }
  
  # Create a data frame of comparisons and their p-values
  significant_comparisons <- data.frame(
    comparison = sapply(quinquennial_comparisons, function(x) paste(x[1], x[2], sep = " vs ")),
    p_value = comparison_pvalues
  )
  
  # Filter for significant comparisons (p < 0.05)
  significant_comparisons <- significant_comparisons[!is.na(significant_comparisons$p_value) & significant_comparisons$p_value < 0.05, ]
  
  # Handle the case where no significant results are found
  if (nrow(significant_comparisons) == 0) {
    message("No significant comparisons found.")
  } else {
    print(significant_comparisons)
  }
}


is_consecutive_quinquennial <- function(age_range1, age_range2) {
  age_bounds1 <- as.numeric(unlist(regmatches(age_range1, gregexpr("\\d+", age_range1))))
  age_bounds2 <- as.numeric(unlist(regmatches(age_range2, gregexpr("\\d+", age_range2))))
  
  print(paste("Comparing:", age_bounds1[1], "vs", age_bounds2[1]))
  
  return(abs(age_bounds1[1] - age_bounds2[1]) == 5)  # Check if the difference is 5
}
# Filter the consecutive comparisons (age difference of 5 years)
quinquennial_comparisons <- Filter(function(pair) {
  is_consecutive_quinquennial(pair[1], pair[2])  # Correctly pass the arguments
}, quinquennial_comparisons)

# Print the result
print(quinquennial_comparisons)

is_consecutive_decade <- function(age_range1, age_range2) {
  # Extract the numeric lower bounds of the age ranges
  age_bounds1 <- as.numeric(unlist(regmatches(age_range1, gregexpr("\\d+", age_range1))))
  age_bounds2 <- as.numeric(unlist(regmatches(age_range2, gregexpr("\\d+", age_range2))))
  
  # Check if the two ranges are consecutive (difference of 10)
  return(abs(age_bounds1[1] - age_bounds2[1]) == 10)
}

# Filter the consecutive comparisons for decade-based comparisons
decade_comparisons <- Filter(function(pair) {
  is_consecutive_decade(pair[1], pair[2])
}, decade_comparisons)

# Print the result
print(decade_comparisons)

# Filter for significant comparisons

# 1. Initialize a vector to store p-values for each comparison
comparison_pvalues <- vector("numeric", length = length(decade_comparisons))

# 2. Loop through each comparison in the list


decade_comparisons <- Filter(
  function(x) !any(is.na(x)) && length(x) == 2,
  decade_comparisons)

# Check if decade_comparisons is empty
if (length(decade_comparisons) == 0) {
  message("No valid comparisons found. Skipping processing.")
} else {
  # Initialize a vector for p-values
  comparison_pvalues <- numeric(length(decade_comparisons))
  
  # Loop through each comparison
  for (i in seq_along(decade_comparisons)) {
    comparison <- decade_comparisons[[i]]
    
    # Check if both elements exist in the p-value matrix
    if (comparison[1] %in% rownames(decade_ttest$p.value) &&
        comparison[2] %in% colnames(decade_ttest$p.value)) {
      # Extract the p-value
      comparison_pvalues[i] <- decade_ttest$p.value[comparison[1], comparison[2]]
    } else {
      # Assign NA for missing comparisons
      comparison_pvalues[i] <- NA
    }
  }
  
  # Create a data frame of comparisons and their p-values
  significant_comparisons <- data.frame(
    comparison = sapply(decade_comparisons, function(x) paste(x[1], x[2], sep = " vs ")),
    p_value = comparison_pvalues
  )
  
  # Filter for significant comparisons (p < 0.05)
  significant_comparisons <- significant_comparisons[!is.na(significant_comparisons$p_value) & significant_comparisons$p_value < 0.05, ]
  
  # Handle the case where no significant results are found
  if (nrow(significant_comparisons) == 0) {
    message("No significant comparisons found.")
  } else {
    print(significant_comparisons)
  }
}


# Filter for significant comparisons

# 4. Plot the boxplot with geom_signif using significant comparisons
p2 <- ggplot(measures_filtered, aes(x = decade_bins, y = Median_FilterObjects_AreaShape_MinorAxisLength, fill = tissue, color = tissue)) +
  geom_jitter(width = 0.2, alpha = 1, size = 1, color = "#B66DFFFF",) +  # Jittered points with full alpha
  geom_boxplot(alpha = 0.7, fill = "#B66DFFFF", color = "black") +  # Violin plot with fill color
  labs(
    title = "",
    x = "Age (years)",
    y = "Thickness (px)"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
        legend.text = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        panel.border = element_rect(color = "grey80", size = 0.5, fill = NA)
  ) + 
  
  # Add geom_signif with significant comparisons
  geom_signif(
    comparisons = decade_comparisons,  # Only significant comparisons
    map_signif_level = TRUE,  # Show asterisks for significance
    step_increase = 0.1,  # Adjust vertical spacing between comparisons
    color = "black"
  )


fig3<-trajectories+p2

pdf("Desktop/TFM/bsc83671/GTEx_v8/Laura/Figure_plots/Fig3ABCD.pdf", width = 12.2, height = 5.3) 
fig3
dev.off()

