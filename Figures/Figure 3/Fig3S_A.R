library(ggplot2)

balanced_accuracy_all = c()
prfx<-"X/"

#Read accuracies
tissue_list = c("Ovary", "Uterus", "Vagina", "BreastMammaryTissue")
for (tissue in tissue_list[1:4]){
  balanced_accuracy = read.csv(paste0(prfx, "Ole/RNApath/NewestVal/", tissue,"/BalSampling_new_balanced_accuracy_1000_intuned10cv10_realmean_enetCVdonor.csv"))
  balanced_accuracy$tissue = tissue
  balanced_accuracy_all = rbind(balanced_accuracy_all, balanced_accuracy)
}
balanced_accuracy_all = balanced_accuracy_all[!balanced_accuracy_all$subtissue %in%c("nerve", "gynecomastoid_hyperplasia"),]

# Transform the dataframe to long format
df_long <- balanced_accuracy_all %>%
  pivot_longer(cols = c(balacc_donor, balacc_tile), 
               names_to = "accuracy_type", 
               values_to = "accuracy")
## add a column to df_long that combines tissue and subtissue columns
df_long$tissue[df_long$tissue == "BreastMammaryTissue"] = "Breast"
df_long$tissue_subtissue <- paste(df_long$tissue, df_long$subtissue, sep = "_")

# Adjust the position of the lollipops slightly
df_long <- df_long %>%
  mutate(position = as.numeric(factor(tissue_subtissue)) + ifelse(accuracy_type == "balacc_donor", -0.2, 0.2))

# Create a position label for x-axis
df_long <- df_long %>%
  mutate(subtissue_label = as.numeric(factor(tissue_subtissue)))
df_long$accuracy_type_name = ifelse(df_long$accuracy_type=="balacc_donor", "per-donor", "per-tile")

custom_colors_tissues <- c("Breast"="#490092FF", "Uterus" ="#6DB6FFFF", "Ovary" = "#009292FF","Vagina" = "#B66DFFFF")


df_long_th<-df_long[df_long$accuracy_type=="balacc_donor",]
df_long_th<-df_long_th[df_long_th$accuracy>0.75,]
df_long_th$subtissue
# Create the lollipop plot

p<-ggplot(df_long, aes(x = position, y = accuracy, color = tissue, alpha = accuracy_type_name, group = tissue)) +
  geom_segment(aes(x = position, xend = position, y = 0, yend = accuracy), size = 1) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "gray50", size = 0.5) +
  
  scale_x_continuous(breaks = unique(df_long$subtissue_label), labels = unique(df_long$tissue_subtissue)) +
  labs(title = "",
       x = "",
       y = "Balanced Accuracy",
       color = "") +
  # theme_minimal() + 
  ylim(0,1)+
  
  scale_color_manual(values = custom_colors_tissues) +
  scale_alpha_manual(values =  c("per-donor" = 1., "per-tile" = 0.65), name = "") +
  #theme(axis.text.x = element_text(, vjust = 0.5, hjust = 1))+
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black", angle = 45, hjust = 1),
        axis.title = element_text(size = 15, color = "black"),
        plot.title = element_text(size = 15, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        panel.grid = element_blank(),                    # Remove all grid lines
        axis.line.x = element_line(size = 0.5, color = "grey80"), # Add x-axis line
        axis.line.y = element_line(size = 0.5, color = "grey80"),
        panel.border = element_rect(color = "grey80", size = 0.5, fill = NA) )

pdf("X/Figures/Fig3S_A.pdf", width = 6.5, height = 4.5) 
p
dev.off()


