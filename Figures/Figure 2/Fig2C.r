.libPaths(c("/x/Rlibs", .libPaths()))
library(ggplot2)
library(dplyr)
library(reshape)
library(yaml)
library(tidyr)
input<-"X/"
tissues <- c("Vagina", "Uterus", "Ovary", "Endocervix", "Ectocervix", "FallopianTube", "BreastMammaryTissue")

yaml_file <-paste0(input, "/Allal/RNAPath/resources/clusters.yaml")
yaml_data <- yaml.load_file(yaml_file)

tissue_colors <- c(
  "Uterus" = "#6DB6FFFF",
  "Ovary" = "#009292FF",
  "Vagina" = "#B66DFFFF",
  "BreastMammaryTissue" = "#490092FF",
  "Endocervix" = "#FF6DB6FF",
  "Ectocervix" = "#FFB6DBFF",
  "FallopianTube" = "#DB6D00FF"
)


# View the structure of the data to ensure it is read correctly
str(yaml_data)
results_folder <- paste0(input, "/Allal/Differential_expression/variance_partition/results/mean_tiles")

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


violin_plot <- function(final_df) {
    # Reshape the data into long format
    df_long <- melt(final_df, varnames = c("Feature", "Subtissue"), value.name = "Age_Proportion")
    df_long$variable <- gsub("age_proportion_", "", df_long$variable)

    
    # Create the violin plot
    ggplot(df_long, aes(x = variable, y = value, fill = variable)) +
        geom_violin(trim = FALSE) +
        geom_boxplot(width = 0.05, outlier.shape = NA, color = "black") +  # Boxplot inside violin plot
        theme_classic() +
        labs(title = "Age Proportion Distribution by Subtissue", 
             x = "Subtissue", 
             y = "% of variance explained by age") +
        scale_fill_manual(values = color_palette(tissue)) # Customize colors if needed
}




# Initialize an empty list to store results for all tissues
all_tissues_list <- list()

for (tissue in tissues) {
    # Initialize an empty list to store the results for the current tissue
    feature_list <- list()

    # Get the list of files for the given tissue
    files <- list.files(paste0(results_folder, "/", tissue), pattern = ".rds$", full.names = FALSE)

    # Extract subtissue names
    subtissues <- gsub(".rds", "", files)
    subtissues <- subtissues[!grepl("_summary", subtissues)]

    # Loop through subtissues
    for(subtissue in subtissues) {
        # Load the results for the given tissue and subtissue
        results <- load_results(tissue, subtissue)

        # Check if results have at least one row
        if (nrow(results) == 0) next 

        # Initialize an empty dataframe for feature proportions
        feature_proportions <- data.frame(matrix(nrow = nrow(results), ncol = 1))

        for (i in 1:nrow(results)) {
            # Extract the current feature (row)
            feature_values <- results[i, ]

            # Get the value of 'age' for this feature
            age_value <- feature_values["age"]

            # Calculate the proportion explained by 'age'
            age_proportion <- age_value * 100

            # Store the result
            feature_proportions[i, 1] <- age_proportion
        }

        # Name the column appropriately
        colnames(feature_proportions) <- paste0("", subtissue)

        # Store the dataframe in the list
        feature_list[[subtissue]] <- feature_proportions
    }

    # Combine all subtissue results into one dataframe for this tissue
    if (length(feature_list) > 0) {
        tissue_df <- do.call(cbind, feature_list)
        all_tissues_list[[tissue]] <- tissue_df
    }
}

# Combine results from all tissues into one dataframe
final_df <- do.call(cbind, all_tissues_list)

# Save the final dataframe
  # Convert to long format for statistical testing
  final_df_long <- final_df %>%
    pivot_longer(cols = everything(),
                 names_to = "subtissue",
                 values_to = "proportion") %>%
    mutate(subtissue = gsub("age_proportion_", "", subtissue))



final_df_long$main_tissue <- sapply(strsplit(final_df_long$subtissue, "\\."), `[`, 1)

# Assign colors based on the main tissue name
final_df_long$color <- tissue_colors[final_df_long$main_tissue]


#plot <- violin_plot(final_df)


#tissue_order <- rev(c("BreastMammaryTissue", "Vagina", "Ovary", "Uterus", 
#                       "FallopianTube", "Endocervix", "Ectocervix"))

tissue_order <- rev(c("Uterus", "FallopianTube", "Ectocervix", "Ovary", 
                       "Endocervix", "Vagina", "BreastMammaryTissue"))


final_df_long$main_tissue <- factor(final_df_long$main_tissue, levels = tissue_order)

final_df_long <- final_df_long %>%
  arrange(main_tissue, subtissue)  # First order by main tissue, then by subtissue

final_df_long$subtissue <- factor(final_df_long$subtissue, levels = unique(final_df_long$subtissue))

final_df_long <- final_df_long[final_df_long$subtissue != "BreastMammaryTissue.gynecomastoid_hyperplasia", ]

unique(final_df_long$subtissue)


# Reorder subtissue by decreasing mean proportion
#final_df_long <- final_df_long %>%
#  group_by(subtissue) %>%
#  mutate(mean_proportion = mean(proportion, na.rm = TRUE)) %>%
#  ungroup() %>%
#  mutate(subtissue = reorder(subtissue, mean_proportion))



# Create a new column to specify alpha values
final_df_long <- final_df_long %>%
  mutate(alpha_value = ifelse(subtissue %in% c(
    "Ectocervix.epithelium", 
    "Endocervix.glandular_epithelium", 
    "FallopianTube.epithelium", 
    "Vagina.epithelium", 
    "Uterus.myometrium", 
    "Ovary.vessels", 
    "BreastMammaryTissue.lobule"
  ), 1, 0.8))





library(ggplot2)
library(dplyr)
library(grid)

# Define the subtissues to make bold
bold_subtissues <- c(
  "Ectocervix.epithelium", 
  "Endocervix.glandular_epithelium", 
  "FallopianTube.epithelium", 
  "Vagina.epithelium", 
  "Uterus.myometrium", 
  "Ovary.vessels", 
  "BreastMammaryTissue.lobule"
)


library(ggplot2)
library(dplyr)

# Create the plot
p <- ggplot(final_df_long, aes(x = proportion, y = subtissue, fill = main_tissue, alpha = alpha_value)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = NA, show.legend = FALSE) + # Hide boxplot in legend
  theme_classic() + 
  theme(
    legend.title = element_text(size = 14, family = "Arial", color = "black"),
    legend.text = element_text(size = 14, family = "Arial", color = "black"),
    text = element_text(family = "Arial", size = 14, color = "black"),
    axis.text = element_text(family = "Arial", size = 14, color = "black"),
    axis.title = element_text(family = "Arial", size = 14, color = "black"), 
    axis.title.y = element_text(margin = margin(r = 20), family = "Arial", size = 1, color = "black")
  ) +
  xlab("Variance explained by age for all features (%)") + 
  scale_fill_manual(values = tissue_colors, name = "Main Tissue") + # Add legend title for colors
  scale_alpha_identity() +  
  scale_y_discrete(labels = function(x) gsub("_", " ", sub("^.*\\.", "", x))) +  
  guides(
    fill = guide_legend(override.aes = list(alpha = 1, shape = NA), reverse = TRUE), # Reverse legend order
    alpha = "none" # Hide alpha legend
  )

# Modify y-axis labels to make selected subtissues bold
p <- p + theme(axis.text.y = element_text(face = ifelse(levels(final_df_long$subtissue) %in% bold_subtissues, "bold", "plain")))

# Save plot
svg("~/X/Laura/Figure_plots/Violin_plot_vp.svg",  width = 6.5, height = 6)
p
dev.off()



# Sort by main_tissue first, then by proportion in descending order
final_df_long <- final_df_long %>%
  arrange(main_tissue, desc(proportion))

# Reverse the subtissue order **within each main_tissue group** without affecting the main_tissue order
final_df_long <- final_df_long %>%
  group_by(main_tissue) %>%
  mutate(subtissue = factor(subtissue, levels = rev(unique(subtissue)))) %>%
  ungroup()

# Create the plot
p <- ggplot(final_df_long, aes(x = proportion, y = subtissue, fill = main_tissue, alpha = alpha_value)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = NA, show.legend = FALSE) +
  theme_classic() +
  theme(
    legend.title = element_text(size = 15, family = "Arial", color = "black"),
    legend.text = element_text(size = 15, family = "Arial", color = "black"),
    text = element_text(family = "Arial", size = 15, color = "black"),
    axis.text = element_text(family = "Arial", size = 15, color = "black"),
    axis.title = element_text(family = "Arial", size = 15, color = "black"),
    axis.title.y = element_blank()
  ) +
  xlab("Variance explained by age for all features (%)") +
  scale_fill_manual(values = tissue_colors, name = "Main Tissue") +
  scale_alpha_identity() +
  scale_y_discrete(labels = function(x) gsub("_", " ", sub("^.*\\.", "", x))) +  
  guides(
    fill = guide_legend(override.aes = list(alpha = 1, shape = NA), reverse = TRUE),
    alpha = "none"
  )

# Bold selected subtissues
p <- p + theme(axis.text.y = element_text(face = ifelse(levels(final_df_long$subtissue) %in% bold_subtissues, "bold", "plain")))

# Save plot
ggsave(p, file = "~/X/Laura/Figure_plots/Violin_plot_vp.svg", device = "svg", width = 8, height = 7.5)

# Create the plot without legend
p2 <- ggplot(final_df_long, aes(x = proportion, y = subtissue, fill = main_tissue, alpha = alpha_value)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = NA, show.legend = FALSE) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 15, color = "black"),
    axis.text = element_text( size = 15, color = "black"),
    axis.title = element_text(size = 15, color = "black"),
    axis.title.y = element_blank()
  ) +
  xlab("Variance explained by age 
  for all features (%)") +
  scale_fill_manual(values = tissue_colors) +
  scale_alpha_identity() +
  scale_y_discrete(labels = function(x) gsub("_", " ", sub("^.*\\.", "", x))) 

# Bold selected subtissues
p2 <- p2 + theme(axis.text.y = element_text(face = ifelse(levels(final_df_long$subtissue) %in% bold_subtissues, "bold", "plain")))

# Save plot
ggsave(p2, file = "X/Laura/Figure_plots/Violin_plot_vp.svg", device = "svg", width = 6, height = 7)


pdf("~/X/Figures/Fig2C.pdf", width = 6, height = 8)  # Adjust width and height as needed
#grid.arrange(grobTree(heatmaptot), grobTree(heat), ncol = 2, widths = c(0.7, 0.3))
p2
dev.off()

svg("~/X/Figures/Fig2C.svg",  width = 6, height = 8)
p2
dev.off()


#######ADD FEATURES TRAJECTORIES
input<-"X/"

library("MOFA2")
library("MOFAdata")
library(dplyr)
library(tidyr)
library(psych)
library(limma)
library(ggplot2)
library(edgeR)
library(readxl)
tissue<-"Uterus"

m<-readRDS("X/Allal/Differential_expression/variance_partition/results/mean_tiles/Uterus/myometrium.rds")
max_feat_myo<-rownames(m[order(m$age, decreasing = TRUE),])[1]

MOFAobject_trained<-readRDS(paste0(input, "/Laura/MOFA/MOFA_", tissue, "_",option,"_","10_ensg_endocorrected.rds"))
metadata_tot <- readRDS(paste0(input, "/Laura/00.Data/v10/", tissue, "/metadata.rds"))
filter<-readRDS(paste0(input, "/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))


#For using continuous ancestry as a covariate (otherwise we use it as a categorical variable)
metadata_tot$Ancestry<- NULL
ancestry_file<- read.table(paste0(input,"/Laura/00.Data/admixture_inferred_ancestry.txt"))
ances<-ancestry_file[,c(1,3)]
colnames(ances)<- c("Donor", "Ancestry")
metadata_tot<-merge(metadata_tot, ances, by= "Donor")
metadata<- metadata_tot[metadata_tot$Donor %in% filter$Subject.ID,]
metadata$Age_bins <- ifelse(metadata$Age<=35, "young",
                            ifelse(metadata$Age>35 & metadata$Age<60, "middle", "old"))
metadata$Age_bins<-as.factor(metadata$Age_bins)



metadata_MOFA<-metadata
colnames(metadata_MOFA)[1]<-"sample"
metadata_MOFA<- metadata_MOFA[metadata_MOFA$sample %in% colnames(MOFAobject_trained@data[[1]]$group1),]
metadata_MOFA$Sample<-NULL
ordered_metadata <- metadata_MOFA[order(metadata_MOFA$Age), ]
samples_metadata(MOFAobject_trained) <- metadata_MOFA
data <- get_data(MOFAobject_trained, 
                 views = "myometrium", 
                 as.data.frame = TRUE
)

data_with_age<-merge(data, metadata_MOFA[,c("sample","Age", "Ancestry", "BMI")], by= "sample")
d<-data_with_age[data_with_age$feature=="346_myometrium",]
filter<-readRDS(paste0("X/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
d<-d[d$sample %in% filter$Subject.ID,]
d<-d[order(d$value, decreasing=TRUE),]
d_max<-d[d$value==max(d$value), ]
d_min<-d[d$value==min(d$value), ]
d_max
d_min
myometrium<-d
head(d)
tail(d)

tissue<-"Vagina"
m<-readRDS("X/Allal/Differential_expression/variance_partition/results/mean_tiles/Vagina/epithelium.rds")
max_feat_epi<-rownames(m[order(m$age, decreasing = TRUE),])[1]

MOFAobject_trained<-readRDS(paste0(input, "/Laura/MOFA/MOFA_", tissue, "_",option,"_","10_ensg_peer2.rds"))
metadata_tot <- readRDS(paste0(input, "/Laura/00.Data/v10/", tissue, "/metadata.rds"))
filter<-readRDS(paste0(input, "/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))


#For using continuous ancestry as a covariate (otherwise we use it as a categorical variable)
metadata_tot$Ancestry<- NULL
ancestry_file<- read.table(paste0(input,"/Laura/00.Data/admixture_inferred_ancestry.txt"))
ances<-ancestry_file[,c(1,3)]
colnames(ances)<- c("Donor", "Ancestry")
metadata_tot<-merge(metadata_tot, ances, by= "Donor")
metadata<- metadata_tot[metadata_tot$Donor %in% filter$Subject.ID,]
metadata$Age_bins <- ifelse(metadata$Age<=35, "young",
                            ifelse(metadata$Age>35 & metadata$Age<60, "middle", "old"))
metadata$Age_bins<-as.factor(metadata$Age_bins)
###Ancestry bins
metadata$Ancestry_bins <-cut(metadata$Ancestry, 
                             breaks = 3, 
                             labels = c("Low", "Medium","High"),
                             include.lowest = TRUE)
metadata$Ancestry_bins<-as.factor(metadata$Ancestry_bins)

metadata_MOFA<-metadata
colnames(metadata_MOFA)[1]<-"sample"
metadata_MOFA<- metadata_MOFA[metadata_MOFA$sample %in% colnames(MOFAobject_trained@data[[1]]$group1),]
metadata_MOFA$Sample<-NULL
ordered_metadata <- metadata_MOFA[order(metadata_MOFA$Age), ]
samples_metadata(MOFAobject_trained) <- metadata_MOFA
data <- get_data(MOFAobject_trained, 
                 views = "epithelium", 
                 as.data.frame = TRUE
)

data_with_age<-merge(data, metadata_MOFA[,c("sample","Age", "Ancestry", "BMI")], by= "sample")
d<-data_with_age[data_with_age$feature=="313_epithelium",]#For age
filter<-readRDS(paste0("X/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
d<-d[d$sample %in% filter$Subject.ID,]
d<-d[order(d$value, decreasing = TRUE),]
d_max<-d[d$value==max(d$value), ]
d_min<-d[d$value==min(d$value), ]
epithelium<-d

total<-rbind(myometrium, epithelium)
colors<-c(tissue_colors[[1]], tissue_colors[[3]])
f<-ggplot(total, aes(x = Age, y = value, color = view, group = "view")) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "loess", se = TRUE, color = "black", fill = "grey80", alpha = 0.5)  +
  scale_color_manual(values = colors) +
  labs(title= " ", x = "Age (years)", y = "Feature value") +  # Label axes
  theme_minimal()+
  theme(axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        plot.title = element_text(size = 15, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 15),
        panel.border = element_rect(color = "grey80", size = 0.5, fill = NA),
        strip.text = element_text(size = 15, color = "black"))+  # Increase facet label size)+
  facet_wrap(~view, scales = "free_y", nrow=4)
# 
dev.off()

pdf("~/X/Figures/Fig2C_2.pdf", width = 3.5, height = 8)  # Adjust width and height as needed
#grid.arrange(grobTree(heatmaptot), grobTree(heat), ncol = 2, widths = c(0.7, 0.3))
f
dev.off()

svg("~/X/Figures/Fig2C_2.svg",  width = 3.5, height = 8)
f
dev.off()
