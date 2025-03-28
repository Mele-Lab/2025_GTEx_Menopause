## Fig 2S A
#Accuracies Top-1
acc_top1<-data.frame(
  Organ= c("Uterus", "Uterus", "Ovary", "Ovary", "Vagina", "Vagina",
           "Breast", "Breast", "Endocervix", "Endocervix", "Ectocervix", "Ectocervix",
           "Fallopian tube", "Fallopian tube"),
  Acc= c(0.96, 0.98, 0.95, 0.97, 0.95, 0.98, 0.96, 0.96, 0.95, 0.98, 0.95, 0.99, 0.92, 0.98),
  Metric= c( "Per-tile", "Per-donor","Per-tile", "Per-donor","Per-tile", "Per-donor","Per-tile", "Per-donor",
             "Per-tile", "Per-donor","Per-tile", "Per-donor","Per-tile", "Per-donor")
)
pal <- c("Uterus" = "#6DB6FFFF", "Ovary" = "#009292FF", "Vagina" = "#B66DFFFF", "Breast" = "#490092FF", "Endocervix"= "#FF6DB6FF", "Ectocervix"="#FFB6DBFF", "Fallopian tube"="#DB6D00FF")
names(pal)
acc_top1$Organ<-factor(acc_top1$Organ, levels=names(pal))
acc_top1$Tissue_offset <- ifelse(acc_top1$Metric == "Per-tile", 
                               as.numeric(acc_top1$Organ) + 0.25,  # Shift per-tile metrics left
                               as.numeric(acc_top1$Organ) - 0.25)  # Shift per-donor metrics right

cnn<-ggplot(acc_top1, aes(x = Tissue_offset, y = Acc, fill = "black", alpha = Metric)) +
  geom_segment(aes(x = Tissue_offset, xend = Tissue_offset, y = 0, yend = Acc), 
               size = 1) +  # Use tissue color for the segment
  geom_point(size = 4) +
  geom_text(aes(label = sprintf("%.2f", Acc)), 
            vjust = 0.5, hjust=-0.2, size = 5, color = "black", angle=90) +  # Accuracy labels in black
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "grey50", size = 0.5) +  # Add horizontal dashed line
  # scale_color_manual(values = "black") +  # Use your custom color palette
  scale_alpha_manual(values = c("Per-tile" = 0.2, "Per-donor" = 1)) +  # Different alpha for metrics
  scale_x_continuous(breaks = 1:7, labels = levels(acc_top1$Organ), expand = c(0.05, 0.05)) +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, by = 0.25), expand = c(0, 0))+
  labs(
    title = "",
    x = "",
    y = "Accuracy"
  ) +
  theme_minimal(base_size = 15) +
  theme(axis.text.y = element_text(size = 15, color = "black"), 
        axis.text.x = element_text(size = 15, color = "black", angle = 25, vjust= 0.8, hjust=0.5),
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

pdf("X/Figures//figS2_A.pdf", width = 5.5, height = 5) 
cnn
dev.off()

svg("X/Figures/figS2_A.svg", width = 5.5, height = 5, pointsize = 12)
cnn
dev.off()


