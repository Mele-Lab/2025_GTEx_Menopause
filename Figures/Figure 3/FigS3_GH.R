#Fig3S_F
##Trajectory of MOFA features
input<-"Desktop/TFM/bsc83671/GTEx_v8/"

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

max<-plot_top_weights(MOFAobject_trained,
                 view = "myometrium",
                 #view = "medulla",
                 factor = 4,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)
max<-as.character(max$data$feature[2])
data <- get_data(MOFAobject_trained, 
                 views = "myometrium", 
                 as.data.frame = TRUE
)

data_with_age<-merge(data, metadata_MOFA[,c("sample","Age", "Ancestry", "BMI")], by= "sample")
d<-data_with_age[data_with_age$feature==max,]
filter<-readRDS(paste0("Desktop/TFM/bsc83671/GTEx_v8/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
d<-d[d$sample %in% filter$Subject.ID,]
d_max<-d[d$value==max(d$value), ]
d_min<-d[d$value==min(d$value), ]
d_max
d_min
myometrium<-d

tissue<-"Vagina"
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
max<-plot_top_weights(MOFAobject_trained,
                      view = "epithelium",
                      #view = "medulla",
                      factor = 5,
                      nfeatures = 10,     # Top number of features to highlight
                      scale = T           # Scale weights from -1 to 1
)
max<-as.character(max$data$feature[1])

data_with_age<-merge(data, metadata_MOFA[,c("sample","Age", "Ancestry", "BMI")], by= "sample")
d<-data_with_age[data_with_age$feature==max,]
filter<-readRDS(paste0("Desktop/TFM/bsc83671/GTEx_v8/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
d<-d[d$sample %in% filter$Subject.ID,]
d<-d[order(d$value, decreasing = TRUE),]
d_max<-d[d$value==max(d$value), ]
d_min<-d[d$value==min(d$value), ]
epithelium<-d


tissue<-"Ovary"
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
                 views = "cortex", 
                 as.data.frame = TRUE
)
max<-plot_top_weights(MOFAobject_trained,
                      view = "cortex",
                      factor = 7,
                      nfeatures = 10,     # Top number of features to highlight
                      scale = T           # Scale weights from -1 to 1
)
max<-as.character(max$data$feature[1])

data_with_age<-merge(data, metadata_MOFA[,c("sample","Age", "Ancestry", "BMI")], by= "sample")
d<-data_with_age[data_with_age$feature==max,]
filter<-readRDS(paste0("Desktop/TFM/bsc83671/GTEx_v8/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
d<-d[d$sample %in% filter$Subject.ID,]
d<-d[order(d$value, decreasing = TRUE),]
d_max<-d[d$value==max(d$value), ]
d_min<-d[d$value==min(d$value), ]
cortex<-d
t_colors <- c(
  "myometrium" = "#6DB6FFFF",
  "cortex" = "#009292FF",
  "epithelium" = "#B66DFFFF"
)
#Gradient for vagina
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
total<- rbind(myometrium, epithelium, cortex)
total$view<-factor(total$view, levels=c("myometrium", "epithelium", "cortex"))

f<-ggplot(total, aes(x = Age, y = value, color = view, group = "view")) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "loess", se = TRUE, color = "black", fill = "grey80", alpha = 0.5)  +
  scale_color_manual(values=t_colors)+
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
  facet_wrap(~view, scales = "free_y", nrow=1)
# facet_wrap(~view, scales = "free_y", nrow=4)


pdf("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/Figure_plots/MOFA_features_m_e_c.pdf", width = 8, height = 4)  # Adjust width and height as needed
#grid.arrange(grobTree(heatmaptot), grobTree(heat), ncol = 2, widths = c(0.7, 0.3))
f
dev.off()

svg("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/Figure_plots/MOFA_features_m_e_c.svg",  width = 8, height = 4)
f
dev.off()

pdf("~/Desktop/Figures/Fig3S_GH.pdf", width = 8, height = 4)  # Adjust width and height as needed
#grid.arrange(grobTree(heatmaptot), grobTree(heat), ncol = 2, widths = c(0.7, 0.3))
f
dev.off()

svg("~/Desktop/Figures/Fig3S_GH.svg",  width = 8, height = 4)
f
dev.off()
