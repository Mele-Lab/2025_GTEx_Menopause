####Supp 1: A-B-C-D
library(dplyr)
library(ggplot2)

common_dir<- "~/Desktop/TFM/bsc83671/GTEx_v8/Laura/05.CNN/35yo_separation_FEMALE/"
tissues<- c("Uterus", "Ovary", "Vagina", "BreastMammaryTissue")
results<-list()
for (tissue in tissues){
  if (tissue =="Vagina"){
    mid_val<-paste0(common_dir, "/", tissue, "/selected_tile_ACC/hybridtrainval_and_mid_256_tile_Alinsaif_age30_vgg19_predictions_filtered", tissue, ".csv")
    traintest<-paste0(common_dir, "/", tissue, "/selected_tile_ACC/hybridtraintest_256_tile_Alinsaif_age30_vgg19_predictions_filtered", tissue, ".csv")
    yp<-paste0(common_dir, "/", tissue, "/selected_tile_ACC/hybridtrainyoung_post_256_tile_Alinsaif_age30_vgg19_predictions_filtered", tissue, ".csv")
    
  }else if (tissue =="BreastMammaryTissue"){
    mid_val<-paste0(common_dir, "/", tissue, "/selected_tile_ACC/hybridtrainval_and_mid_512_tile_Alinsaif_age30_vgg19_predictions_filtered", tissue, ".csv")
    traintest<-paste0(common_dir, "/", tissue, "/selected_tile_ACC/hybridtraintest_512_tile_Alinsaif_age30_vgg19_predictions_filtered", tissue, ".csv")
    yp<-paste0(common_dir, "/", tissue, "/selected_tile_ACC/hybridtrainyoung_post_512_tile_Alinsaif_age30_vgg19_predictions_filtered", tissue, ".csv")
    
  }else{
    mid_val<-paste0(common_dir, "/", tissue, "/selected_tile_ACC/hybridtrainval_and_mid_512_tile_Alinsaif_age30_vgg19_predictions_filtered", tissue, ".csv")
    traintest<-paste0(common_dir, "/", tissue, "/selected_tile_ACC/hybridtraintest_512_tile_Alinsaif_age30_vgg19_predictions_filtered", tissue, ".csv")
    yp<-paste0(common_dir, "/", tissue, "/selected_tile_ACC/hybridtrainyoung_post_512_tile_Alinsaif_age30_vgg19_predictions_filtered", tissue, ".csv")
    
  }
  mid_val<- read.csv(mid_val)
  mid_val$Donor <- gsub("^(\\w+-\\w+).*", "\\1", mid_val$tile)
  mid_val$Donor<- as.factor(mid_val$Donor)
  
  traintest<- read.csv(traintest)
  traintest$Donor <- gsub("^(\\w+-\\w+).*", "\\1", traintest$tile)
  traintest$Donor<- as.factor(traintest$Donor)
  
  yp<- read.csv(yp)
  yp$Donor <- gsub("^(\\w+-\\w+).*", "\\1", yp$tile)
  yp$Donor<- as.factor(yp$Donor)
  
  
  traintest2<- traintest[!(traintest$Donor %in% mid_val$Donor),]
  mid_val<- mid_val[!(mid_val$Donor %in% traintest$Donor),]
  yp<- yp[!(yp$Donor %in% traintest$Donor),]
  
  yp<-c()
  prob_file<-rbind(mid_val, traintest, yp)
  #prob_file<-mid_val
  prob_file<-rbind(mid_val, traintest2, yp)
  
  result <- prob_file %>%
    group_by(Donor) %>%
    summarise("old" = mean(prob_class0), "young" = mean(prob_class1), "SD_young"=sd(prob_class1))
  result$Tissue<- tissue
  results[[tissue]]<-result
  # write.csv(result, paste0("Desktop/TFM/bsc83671/GTEx_v8/Laura/13.Hybrid_probs/", tissue, "_hybrid_tile_probabilities_filtered_all_TILE.csv"))
  
}
##QUADRANTS
ovary_probs<-as.data.frame(results$Ovary)[,c(1,2)]
uterus_probs<-as.data.frame(results$Uterus)[,c(1,2)]

age_file<-read.csv("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/00.Data/GTEx_Subject_Phenotypes.GRU.csv", sep="\t", header=TRUE)
colnames(age_file)[2]<-"Donor"
df_combined<-merge(ovary_probs, uterus_probs, by="Donor")
colnames(df_combined)<-c("Donor", "Ovary_probs", "Uterus_probs")
df_combined<- merge(df_combined,age_file[,c(2,5)], by="Donor") 
ej<-df_combined[df_combined$Ovary_probs<0.25 & df_combined$Uterus_probs>0.5,]
custom_palette <- c("#FFFF6DFF", "#490092FF") # Example gradient colors

df_combined<-df_combined[order(df_combined$AGE),]
ej<-df_combined[df_combined$AGE>35 & df_combined$AGE<60,]
old_ovary_yong_ut<- nrow(ej[ej$Ovary_probs>0.5 & ej$Uterus_probs<0.5,])
old_ut_yong_ov<- nrow(ej[ej$Ovary_probs<0.5 & ej$Uterus_probs>0.5,])

quadrants <- ggplot(df_combined, aes(x = Uterus_probs, y = Ovary_probs, color = AGE)) +
  # Scatter plot
  geom_point(alpha = 0.7, size = 3) +
  # Quadrant lines
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50", size = 0.8) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", size = 0.8) +
  # Custom color gradient
  scale_color_gradientn(colors = custom_palette) +
  # Axis labels and title
  labs(
    title = " ",
    x = "P(post-menopausal) uterus",
    y = "P(post-menopausal) ovary",
    color = "Age"
  ) +
  # Minimal theme
  theme_minimal() +
  theme(axis.text.y = element_text(size = 14, color = "black"), 
        axis.text.x = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, color = "black"), 
        plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 14, color = "black"),
        panel.grid = element_blank(),    
        legend.title = element_text(size = 14),
        panel.border = element_rect(color = "grey80", size = 0.5, fill = NA),
        legend.position = "bottom"
  )
quadrants


######Distribution boxplot
library(ggplot2)
library(tidyr)

# Reshape data to long format
df_long <- df_combined %>%
  pivot_longer(cols = c(Uterus_probs, Ovary_probs), 
               names_to = "Tissue", 
               values_to = "Probability")
pal <- c("Uterus" = "#6DB6FFFF", "Ovary" = "#009292FF", "Vagina" = "#B66DFFFF", "BreastMammaryTissue" = "#490092FF", "CervixEndocervix"= "#FF6DB6FF", "CervixEctocervix"="#FFB6DBFF", "FallopianTube"="#DB6D00FF")
pal <- c("Uterus" = "#6DB6FFFF", "Ovary" = "#009292FF", "Vagina" = "#B66DFFFF", "BreastMammaryTissue" = "#490092FF", "CervixEndocervix"= "#FF6DB6FF", "CervixEctocervix"="#FFB6DBFF", "FallopianTube"="#DB6D00FF")

# Create the boxplot
boxplot_plot <- ggplot(df_long, aes(x = Tissue, y = Probability, fill = Tissue)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("Uterus_probs" = pal[[1]] , "Ovary_probs" = pal[[2]])) + # Customize colors
  labs(
    title = "",
    x = "",
    y = "Probability"
  ) +
    theme_minimal() +
      theme(axis.text.y = element_text(size = 14, color = "black"), 
            axis.text.x = element_text(size = 14, color = "black"),
            axis.title = element_text(size = 14, color = "black"), 
            plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
            legend.text = element_text(size = 14, color = "black"),
            panel.grid = element_blank(),    
            legend.title = element_text(size = 14),
            panel.border = element_rect(color = "grey80", size = 0.5, fill = NA),
            legend.position = "none"
      
  )

# Print the boxplot
boxplot_plot




###B
ge_acc<-data.frame(
  Organ=c("Uterus", "Ovary", "Vagina"),
  Accuracy=c(0.902, 0.935, 0.696)
)
ggplot(global, aes(x = Tissue_offset, y = acc, color = Organ, alpha = Metric)) +
  geom_point(size = 4) +
  geom_text(aes(label = sprintf("%.2f", acc)), 
            vjust = -0.5, size = 5, color = "black") +  # Accuracy labels in black
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "grey50", size = 0.5) +  # Add horizontal dashed line
  scale_color_manual(values = pal[1:4]) +  # Use your custom color palette
  scale_alpha_manual(values = c("Per-tile" = 0.5, "Per-donor" = 1)) +  # Different alpha for metrics
  scale_x_continuous(breaks = 1:4, labels = levels(global$Tissue), expand = c(0.2, 0.2)) +
  scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0)) 

ge_acc<-ggplot(ge_acc, aes(x = Organ, y = Accuracy, color = Organ)) +
  geom_segment(aes(x = Organ, xend = Organ, y = 0, yend = Accuracy), 
               size = 1) +  # Use the tissue color for the segment
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "grey50", size = 0.5) +  # Add horizontal dashed line
  
  geom_point(size = 4) +
  geom_text(aes(label = sprintf("%.3f", Accuracy)), 
            vjust = -1, size = 5, color = "black") +  # Accuracy labels in black
  scale_color_manual(values = pal) +
  scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
  labs(
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

###C
suppressMessages(library(edgeR))
suppressMessages(library(limma))
suppressMessages(library(dplyr))
suppressMessages(library(DEswan))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(parallel))
suppressMessages(library(data.table))
suppressMessages(library(gplots))

window.center = seq(21,70,1)
buckets.size = 10
interval.size = 0

tissue<-"Vagina"
outpath<- paste0("Desktop/TFM/bsc83671/GTEx_v8/Laura/DEswan/All_tissues/", tissue, "/")
# buckets.size = 10
# interval.size = 0
r<-readRDS(paste0("Desktop/TFM/bsc83671/GTEx_v8/Laura/DEswan/All_tissues/", tissue, "/DEswan_", tissue, "_female_results.rds"))

pval<-"adj.P.Val"
metadata_tot <- readRDS(paste0(input, "Laura/00.Data/v10/", tissue, "/metadata.rds"))
filter<-readRDS(paste0(input, "/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
metadata<- metadata_tot[(metadata_tot$Donor %in% filter$Subject.ID),]   

sample_counts <- sapply(window.center, function(center) {
  sum(metadata$Age >= (center - buckets.size - interval.size) & 
        (metadata$Age <= (center + buckets.size + interval.size)) &
        !(metadata$Age > (center - interval.size) & 
            metadata$Age < (center + interval.size)))  
})


# Convert to DataFrame for plotting
df_counts <- data.frame(Window = window.center, Samples = sample_counts)


####1. Extract the significant genes in at least one comparison
num_cores <- detectCores() - 1  # Use one less than max cores
### 1. Extract significant genes in parallel
extract_significant_genes <- function(sublist) {
  if (!is.null(sublist) && pval %in% colnames(sublist) && "gene.name.x" %in% colnames(sublist)) {
    return(sublist[sublist[[pval]] < 0.05, "gene.name.x", drop = FALSE])
  }
  return(NULL)
}

# Use parallel execution instead of lapply
significant_genes_list <- mclapply(r, extract_significant_genes, mc.cores = num_cores)
significant_genes_list <- significant_genes_list[!sapply(significant_genes_list, is.null)]
final_significant_genes <- unique(do.call(rbind, significant_genes_list))
### 2. Optimize P-value matrix generation
all_windows <- seq_along(r)
gene_list <- unique(unlist(lapply(r, rownames)))

# Create an empty data.table instead of a matrix
p_value_matrix <- data.table(matrix(NA_real_, nrow = length(gene_list), ncol = length(all_windows)))
setnames(p_value_matrix, as.character(all_windows))
setattr(p_value_matrix, "rownames", gene_list)

# Fill p-value matrix in parallel
p_value_matrix <- do.call(cbind, mclapply(all_windows, function(win) {
  if (!is.null(r[[win]]) && is.data.frame(r[[win]])) {
    p_values <- r[[win]][[pval]]
    names(p_values) <- rownames(r[[win]])
    return(p_values)
  }
  return(rep(NA_real_, length(gene_list)))  # Return NA if data is missing
}, mc.cores = num_cores))

# Convert to data frame
p_value_df <- as.data.frame(p_value_matrix)

### 3. Optimize Counting Significant Genes
count_significant <- function(p_df, thresholds = c(0.1, 0.05, 0.01, 0.001)) {
  return(sapply(thresholds, function(thresh) colSums(p_df < thresh, na.rm = TRUE)))
}

res.signif <- count_significant(p_value_df)
rownames(res.signif) <- window.center
colnames(res.signif) <- c(0.1, 0.05, 0.01, 0.001)

threshold_interest<- "0.05"
df_counts <- data.table(Window = window.center, Samples = sample_counts)
df_counts$Genes <- res.signif[, threshold_interest]
df_counts$Genes_normalized<- df_counts$Genes/ df_counts$Samples


uterus_counts<-df_counts
ovary_counts<-df_counts
vagina_counts<-df_counts

ovary_counts$Organ<-"Ovary"
uterus_counts$Organ<-"Uterus"
vagina_counts$Organ<-"Vagina"

t_t<-as.data.frame(rbind(ovary_counts, uterus_counts, vagina_counts))
max(t_t[t_t$Organ=="Uterus", "Genes_normalized"])
pal <- c("Uterus" = "#6DB6FFFF", "Ovary" = "#009292FF", "Vagina" = "#B66DFFFF", "BreastMammaryTissue" = "#490092FF",  "Endocervix"= "#FF6DB6FF", "Ectocervix"="#FFB6DBFF", "Fallopian Tube"="#DB6D00FF")

not_normalized_degs <- ggplot(t_t, aes(x = Window, y = Genes, group = Organ, color = Organ)) +
  geom_point()+
  # geom_line(aes(group = Organ), size = 1) +  # Add lines to connect the dots
  # geom_smooth(method = "loess", se = FALSE, span = 0.2, size = 1) +  # Add smoothing line
  # geom_smooth(method = "loess", se = FALSE, span = 0.5, size = 1) + 
  scale_color_manual(values = pal) +  # Adjust colors as needed
  scale_size(range = c(2, 6)) +
  scale_y_continuous(limits = c(0, NA))+#, expand = c(0, 0)) +
  
  labs(title = "", 
       x = "Age (years)", 
       y = paste0("#DEGs"), 
       color = "Organ", size = "Number of Samples") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 14, color = "black"), 
        axis.text.x = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, color = "black"), 
        plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank(),                    # Remove all grid lines
        # axis.line.x = element_line(size = 0.5, color = "grey80"), # Add x-axis line
        # axis.line.y = element_line(size = 0.5, color = "grey80"),
        panel.border = element_rect(color = "grey80", size = 0.5, fill = NA) # Add full border
  )



#####ALL
library(patchwork)
sup1_1<-quadrants+ ge_acc+not_normalized_degs+plot_layout(nrow = 1, heights = unit(1, "null"))
# Print the plot
pdf("Desktop/TFM/bsc83671/GTEx_v8/Laura/Figure_plots/figS1_ABCD_2.pdf", width = 11, height = 5) 
sup1_1
dev.off()

svg("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/Figure_plots/figS1_ABCD_2.svg", width = 11, height = 5, pointsize = 12)
sup1_1
dev.off()

pdf("Desktop/Figures/figS1_ABCD_2.pdf", width = 11, height = 5) 
sup1_1
dev.off()

svg("Desktop/Figures/figS1_ABCD_2.svg", width = 11, height = 5, pointsize = 12)
sup1_1
dev.off()
