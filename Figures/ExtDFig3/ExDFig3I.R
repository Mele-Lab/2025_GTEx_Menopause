##Plotting follicle count

library(ggplot2)
library(dplyr)
library(gridExtra)
tissue_colors <- c(
  "Uterus" = "#6DB6FFFF",
  "Ovary" = "#009292FF",
  "Vagina" = "#B66DFFFF",
  "BreastMammaryTissue" = "#490092FF",
  "Endocervix" = "#FF6DB6FF",
  "Ectocervix" = "#FFB6DBFF",
  "FallopianTube" = "#DB6D00FF"
  
)

#1. Read the follicle count file
fol_count<- read.table("/home/laura/BSC/Daniel/RNAPath/results_follicles/interior_follicles_count_Cyst_MatFoll_merged.txt")
fol_count$V1<- gsub("^(\\w+-\\w+).*", "\\1", fol_count$V1)
colnames(fol_count)<- c("Donor", "Count")
GTEx_Subject_Phenotypes.GRU <- read.csv('/X/00.Data/GTEx_phenotypes/GTEx_Subject_Phenotypes.GRU.csv', sep = '\t', header = TRUE)
filter<-readRDS("/X/03.Image_processing/Second_filtering_images/Ovary_final_filtered_images.rds")
GTEx_Subject_Phenotypes.GRU<- GTEx_Subject_Phenotypes.GRU[GTEx_Subject_Phenotypes.GRU$SUBJID %in% filter$Subject.ID,]
fol_count_df<- merge(fol_count, GTEx_Subject_Phenotypes.GRU[,c("SUBJID", "AGE")], by.x = "Donor", by.y ="SUBJID")
fol_count_df$Group<-ifelse(fol_count_df$AGE<30, "[20-29]", 
                           ifelse(fol_count_df$AGE>29 & fol_count_df$AGE<40, "[30-39]",
                                  ifelse(fol_count_df$AGE>39 & fol_count_df$AGE<50, "[40-49]",
                                         ifelse(fol_count_df$AGE>49 & fol_count_df$AGE<60, "[50-59]","[60-70]"))))

fol_mean<-fol_count_df %>% 
  group_by(Group)%>%
  summarise(
    mean_count = mean(Count, na.rm =TRUE),
    n=n()
  )

        
p2 <- ggplot(fol_count_df, aes(x = Group, y = Count)) +
  geom_violin(fill= tissue_colors[2], alpha = 0.7, width = 1)+
  geom_boxplot(fill= tissue_colors[2], width =0.2)+
  # geom_col(fill= tissue_colors[2])+
  
  theme_classic() +  # Classic theme
  theme(
    text = element_text(color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 14, face = "bold", color = "black")
  ) +
  labs(title = "",
       x = "",
       y = "Follicle count")
p2
# Save the plot
pdf("/X/Paper_material/Post_review_analyses/variance_partition/follicle_count_boxplots.pdf", width = 6, height = 6) # Adjust width and height as needed
p2
dev.off()


# Plot the data using ggplot2
bar_plot <- ggplot(fol_mean, aes(y = Group, x = n)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, fill = tissue_colors[2]) +  # Do not sum, just use single values  scale_fill_identity() +  # Use the custom colors with alpha transparency
  scale_fill_identity()+
  labs(x = "Sample Size", y = "Organ") +
  scale_x_reverse() +  # Reverse the x-axis so that 0 is on the right
  theme_minimal() +
  scale_y_discrete(limits = rev) +   # 🔄 reverse the order
  
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
# Create the violin plot for age distribution
violin_plot_count <- ggplot(fol_count_df, aes(x = Count, y = Group)) +
  geom_boxplot(alpha = 0.3, width = 0.5, fill = tissue_colors[2]) +  # Adjust boxplot width to make them narrower
  geom_violin(alpha = 0.5, scale = "width", fill =tissue_colors[2]) +
  scale_y_discrete(limits = rev) +   # 🔄 reverse the order
  
  labs(x = "Follicle count", y = "") +
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
  bar_plot, violin_plot_count,
  ncol = 2, 
  widths = c(0.5, 0.9)
  # Adjust relative width to make plots closer
)


pdf("/X/Paper_material/Post_review_analyses/variance_partition/follicle_count_with_sample_size.pdf", width = 6, height = 6) # Adjust width and height as needed
grid.arrange(grobTree(combined_plot))
dev.off()



