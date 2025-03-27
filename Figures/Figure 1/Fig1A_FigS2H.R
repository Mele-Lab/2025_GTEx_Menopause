#Fig A
library(ggplot2)
library(reshape2) 
library(gridExtra)  
library(colorBlindness)
library(tidyr)
library(ComplexHeatmap)
library(cowplot)
library(patchwork)  # For combining the plots
library(RColorBrewer)
library(dplyr)
library(paletteer)
palette <- paletteer_d("colorBlindness::paletteMartin")


#Samples sizes for images and RNA data
pal <- c("Uterus" = "#6DB6FFFF", "Ovary" = "#009292FF", "Vagina" = "#B66DFFFF", "BreastMammaryTissue" = "#490092FF", "CervixEndocervix"= "#FF6DB6FF", "CervixEctocervix"="#FFB6DBFF", "FallopianTube"="#DB6D00FF")
tissues<-c("Uterus", "Ovary", "Vagina", "BreastMammaryTissue", "CervixEndocervix", "CervixEctocervix", "FallopianTube")
samples_list<-list()
probs_list<-list()
age_file<-read.csv("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/00.Data/GTEx_Subject_Phenotypes.GRU.csv", sep="\t", header=TRUE)

for (tissue in tissues){
  
  metadata_tot <- readRDS(paste0("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/00.Data/v10/", tissue, "/metadata.rds"))
  colnames(metadata_tot)<-ifelse(colnames(metadata_tot) =="Donor","SUBJID",
                                 ifelse(colnames(metadata_tot)=="Age", "AGE", colnames(metadata_tot)))
  metadata_tot<-metadata_tot[,c("SUBJID", "AGE", "BMI")]
  ancestry_file<- read.table("Desktop/TFM/bsc83671/GTEx_v8/Laura/00.Data/admixture_inferred_ancestry.txt")
  ances<-ancestry_file[,c(1,3,4)]
  colnames(ances)<- c("SUBJID", "Ancestry", "Ancestry_cat")
  metadata_tot<-merge(metadata_tot, ances, by= "SUBJID")
  
  metadata_tot$Samples<-"RNA-Seq"
  filter<-readRDS(paste0("Desktop/TFM/bsc83671/GTEx_v8/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
  colnames(filter)[4]<-"SUBJID"
  filter$Samples<-"Images"
  filter<- merge(filter[,c("SUBJID", "Samples")], age_file[,c(2,5, 12)], by="SUBJID")
  filter<- merge(filter, ances, by="SUBJID")
  filter$Sample_size<- length(unique(filter$SUBJID))
  metadata_tot<-metadata_tot[metadata_tot$SUBJID %in% filter$SUBJID,]
  metadata_tot$Sample_size<- length(unique(metadata_tot$SUBJID))
  
  #combined_df<-filter
  combined_df <- bind_rows(metadata_tot, filter)
  combined_df$Tissue<-tissue
  samples_list[[tissue]]<- combined_df
  
  
}

combined_df <- bind_rows(samples_list)

sample_df <- combined_df %>%
  group_by(Tissue, Samples) %>%
  summarise(Sample_size = n(), .groups = "drop")  # Use n() to count the number of entries

sample_df$Tissue<- factor(sample_df$Tissue, levels= tissues)
sample_df$Color_with_alpha <- mapply(function(tissue, sample_type) {
  color <- pal[tissue]  # Get the base color for the tissue
  if (sample_type == "RNA-Seq") {
    # Add transparency to the color for RNA-Seq (using proper hex alpha format)
    rgb_val <- col2rgb(color)
    transparent_color <- rgb(rgb_val[1], rgb_val[2], rgb_val[3], maxColorValue = 255, alpha = 64)  # Alpha 128 corresponds to 50% transparency
    return(transparent_color)
  } else {
    return(color)  # Use the solid color for Images
  }
}, sample_df$Tissue, sample_df$Samples)
sample_df$Organ<- factor(sample_df$Tissue, levels= rev(tissues))

combined_df$Organ <- factor(combined_df$Tissue, levels = rev(tissues))
combined_df$Sample_size <- factor(combined_df$Sample_size)

combined_df$Ancestry_num<-ifelse(combined_df$Ancestry_cat=="AFR", 1, 
                                 ifelse(combined_df$Ancestry_cat=="AMR", 2, 3))

endo<-combined_df[combined_df$Tissue =="CervixEndocervix" & combined_df$Samples =="Images",]
ecto<-combined_df[combined_df$Tissue =="CervixEctocervix" & combined_df$Samples =="Images",]
fall<-combined_df[combined_df$Tissue =="FallopianTube" & combined_df$Samples =="Images",]
ecto<-combined_df[combined_df$Tissue =="CervixEctocervix" & combined_df$Samples =="Images",]

table(endo$Ancestry_cat)
table(ecto$Ancestry_cat)
table(fall$Ancestry_cat)

# combined_df %>%
#   group_by(Organ) %>%
#   summarise(p_value = shapiro.test(Ancestry_num)$p.value)
# kruskal.test(Ancestry_num ~ Organ, data = combined_df)
# library(FSA)
# dunnTest(Ancestry.y ~ Organ, data = combined_df, method = "bonferroni")

# Ancestry distro is the same across organs

# Plot the data using ggplot2
bar_plot <- ggplot(sample_df, aes(y = Organ, x = Sample_size, fill = Color_with_alpha)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # Do not sum, just use single values  scale_fill_identity() +  # Use the custom colors with alpha transparency
  scale_fill_identity()+
  labs(x = "Sample Size", y = "Organ") +
  scale_x_reverse() +  # Reverse the x-axis so that 0 is on the right
  theme_minimal() +
  
  theme(axis.text.y = element_text(size = 15, color = "black"), 
        axis.text.x = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 15, color = "black"), 
        axis.title.y = element_blank(),
        
        plot.title = element_text(size = 15, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        panel.grid = element_blank(),    
        legend.title = element_text(size = 15),
        panel.border = element_rect(color = "grey80", size = 0.5, fill = NA),
        legend.position = "left"
  )

# Display the plot
bar_plot

# Create the violin plot for age distribution
violin_plot_age <- ggplot(combined_df[combined_df$Samples=="Images",], aes(x = AGE, y = Organ, fill = Organ, color=Organ)) +
  geom_boxplot(alpha = 0.3, width = 0.5) +  # Adjust boxplot width to make them narrower
  geom_violin(alpha = 0.5, scale = "width") +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  
  labs(x = "Age", y = "Organ") +
  theme_minimal() +
  theme(
    legend.position = "none",  # Legend at the bottom
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 15, color = "black"),
    
    axis.text.x = element_text(size = 15, color = "black"),
    plot.title = element_text(size = 15, hjust = 0.5, color = "black"),
    legend.text = element_text(size = 15, color = "black"),
    panel.grid = element_blank(),    
    legend.title = element_text(size = 15),
    panel.border = element_rect(color = "grey80", size = 0.5, fill = NA),
  )
violin_plot_bmi <- ggplot(combined_df[combined_df$Samples=="Images",], aes(x = BMI, y = Organ, fill = Organ, color=Organ)) +
  geom_boxplot(alpha = 0.3, width = 0.5) +  # Adjust boxplot width to make them narrower
  geom_violin(alpha = 0.5, scale = "width") +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  
  labs(x = "BMI", y = "Organ") +
  theme_minimal() +
  theme(
    legend.position = "none",  # Legend at the bottom
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 15, color = "black"),
    
    axis.text.x = element_text(size = 15, color = "black"),
    plot.title = element_text(size = 15, hjust = 0.5, color = "black"),
    legend.text = element_text(size = 15, color = "black"),
    panel.grid = element_blank(),    
    legend.title = element_text(size = 15),
    panel.border = element_rect(color = "grey80", size = 0.5, fill = NA),
  )

violin_plot_ancestry <- ggplot(combined_df[combined_df$Samples=="Images",], aes(x = Ancestry, y = Organ, fill = Organ, color=Organ)) +
  geom_boxplot(alpha = 0.3, width = 0.5) +  # Adjust boxplot width to make them narrower
  geom_violin(alpha = 0.5, scale = "width") +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  
  labs(x = "Ancestry", y = "Organ") +
  theme_minimal() +
  theme(
    legend.position = "none",  # Legend at the bottom
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 15, color = "black"),
    
    axis.text.x = element_text(size = 15, color = "black"),
    plot.title = element_text(size = 15, hjust = 0.5, color = "black"),
    legend.text = element_text(size = 15, color = "black"),
    panel.grid = element_blank(),    
    legend.title = element_text(size = 15),
    panel.border = element_rect(color = "grey80", size = 0.5, fill = NA),
  )
combined_plot <- grid.arrange(
  bar_plot, violin_plot_age,
  ncol = 2, 
  widths = c(1, 0.9)
  # Adjust relative width to make plots closer
)
combined_plot <- grid.arrange(
  bar_plot, violin_plot_bmi, violin_plot_ancestry,
  ncol = 3, 
  widths = c(1, 0.7, 0.7)
  
  # top = "Age distribution per tissue"
  # Adjust relative width to make plots closer
)

pdf("Desktop/TFM/bsc83671/GTEx_v8/Laura/Figure_plots/fig1_A.pdf", width = 8, height = 5) 
grid.arrange(grobTree(combined_plot))
dev.off()

svg("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/Figure_plots/fig1_A.svg", width = 8, height = 5, pointsize = 12)
grid.arrange(grobTree(combined_plot))
dev.off()

pdf("Desktop/Figures//fig1_A.pdf", width = 8, height = 5) 
grid.arrange(grobTree(combined_plot))
dev.off()

svg("Desktop/Figures/fig1_A.svg", width = 8, height = 5, pointsize = 12)
grid.arrange(grobTree(combined_plot))
dev.off()


###### 2H supp
pdf("Desktop/TFM/bsc83671/GTEx_v8/Laura/Figure_plots/figS2_H.pdf", width = 8.5, height = 5) 
grid.arrange(grobTree(combined_plot))
dev.off()

svg("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/Figure_plots/figS2_H.svg", width = 8.5, height = 5, pointsize = 12)
grid.arrange(grobTree(combined_plot))
dev.off()

pdf("Desktop/Figures//figS2_H.pdf", width = 8.5, height = 5) 
grid.arrange(grobTree(combined_plot))
dev.off()

svg("Desktop/Figures/figS2_H.svg", width = 8.5, height = 5, pointsize = 12)
grid.arrange(grobTree(combined_plot))
dev.off()

