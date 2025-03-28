#Fig 1CDE
library(patchwork)
library(jsonlite)
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(segmented)  # For davies.test
library(ggplot2)
library(dplyr)
library(tidyr)
library(stats)


#C
#Read CNN's probabilities per-donor and per-tile
result_cnn<- readRDS("X/Laura/05.CNN/metrics_new_validation_filtered_perdonor_ACC.rds")
result_cnn_pertile<- readRDS("X/Laura/05.CNN/metrics_new_validation_filteredpertile_ACC.rds")
result_cnn$Metric<- "Per-donor"
result_cnn_pertile$Metric <- "Per-tile"
pal <- c("Uterus" = "#6DB6FFFF", "Ovary" = "#009292FF", "Vagina" = "#B66DFFFF", "Breast" = "#490092FF", "CervixEndocervix"= "#FF6DB6FF", "CervixEctocervix"="#FFB6DBFF", "FallopianTube"="#DB6D00FF")

global<- rbind(result_cnn, result_cnn_pertile)#
rownames(global)[4]<-"Breast"
global$Tissue<- c("Uterus", "Ovary", "Vagina", "Breast", "Uterus", "Ovary", "Vagina", "Breast")
global$Tissue<- factor(global$Tissue, levels= c("Uterus", "Ovary", "Vagina", "Breast"))
global$Organ<-global$Tissue

global$Tissue_offset <- ifelse(global$Metric == "Per-tile", 
                               as.numeric(global$Tissue) + 0.1,  # Shift per-tile metrics left
                               as.numeric(global$Tissue) - 0.1)  # Shift per-donor metrics right

cnn<-ggplot(global, aes(x = Tissue_offset, y = acc, color = Organ, alpha = Metric)) +
  geom_segment(aes(x = Tissue_offset, xend = Tissue_offset, y = 0, yend = acc), 
               size = 1) +  # Use tissue color for the segment
  geom_point(size = 4) +
  geom_text(aes(label = sprintf("%.2f", acc)), 
            vjust = -0.5, size = 5, color = "black") +  # Accuracy labels in black
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "grey50", size = 0.5) +  # Add horizontal dashed line
  scale_color_manual(values = pal[1:4]) +  # Use your custom color palette
  scale_alpha_manual(values = c("Per-tile" = 0.5, "Per-donor" = 1)) +  # Different alpha for metrics
  scale_x_continuous(breaks = 1:4, labels = levels(global$Tissue), expand = c(0.2, 0.2)) +
  scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
  labs(
    title = "",
    x = "",
    y = "Accuracy"
  ) +
  theme_minimal(base_size = 15) +
  theme(axis.text.y = element_text(size = 15, color = "black"), 
        axis.text.x = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 15, color = "black"), 
        plot.title = element_text(size = 15, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        panel.grid = element_blank(),                    # Remove all grid lines
        # axis.line.x = element_line(size = 0.5, color = "grey80"), # Add x-axis line
        # axis.line.y = element_line(size = 0.5, color = "grey80"),
        panel.border = element_rect(color = "grey80", size = 0.5, fill = NA),
        legend.position = "bottom" # Add full border
  )+
  guides(
    alpha = guide_legend(title = "", override.aes = list(size = 4, alpha = c(1, 0.5), shape = 16)),
    # shape = guide_legend(title = "", override.aes = list(size = 4)),
    color = "none",  # This removes the color legend
    # alpha = guide_legend(override.aes = list(shape = 14, size = 4)) # Ensures alpha doesn't show 'A'
  )
cnn 


#####Trajectories
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
tiss<-tissues[3]
merged_df = readRDS(paste0(prfx,tiss,"_subtissues_classification.rds"))

if(tiss == "Uterus"){
  merged_df<-merged_df[merged_df$organ %in% c("myometrium", "vessels"),]
}else if(tiss == "Vagina"){
  merged_df<-merged_df[merged_df$organ %in% c("stroma", "epithelium", "lamina_propria"),]
}else{
  merged_df<-merged_df[merged_df$organ %in% c("corpora", "vessels", "cortex"),]
}

### 2. Analyze the data
# Analyze the CNN data
results = analyze_organ_trajectories_gen(merged_df, p_threshold = 0.05)

# Plot the plot
p<-plot(results$plot)
p <- results$plot + 
  theme(
    legend.title = element_blank(),
    legend.direction = "vertical",
    legend.position = c(0.8, 0.15),  # Move legend inside the bottom-right
    legend.background = element_rect(fill = alpha("white", 0.6), color = NA),  # Semi-transparent background
    legend.key = element_blank()  # Remove borders around legend items
  )
p
# Save the plot
ggsave(paste0(prfx, "CNN_classification_plot.png"), plot = results$plot, width = 10, height = 6, dpi = 300)


cnn<-p
ge<-p
acc<-p
all_tog<- acc+cnn+ge

pdf("Desktop/Figures/Fig1_CDE.pdf", width = 12.2, height = 5.3) 
all_tog
dev.off()



