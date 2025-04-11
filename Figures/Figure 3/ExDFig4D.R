#Extended Data Fig 4D

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
input<- "/X/"
MOFAobject_trained<-readRDS(paste0(input, "/Laura/13.MOFA/MOFA_", tissue, "_",option,"_","10_ensg_CORRECTED.rds"))
metadata_tot <- readRDS(paste0(input, "/Laura/00.Data/v10/", tissue, "/metadata.rds"))
filter<-readRDS(paste0(input, "/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
metadata<- metadata_tot[metadata_tot$Donor %in% filter$Subject.ID,]
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

data_with_age<-merge(data, metadata_MOFA[,c("sample","Age", "BMI")], by= "sample")
d<-data_with_age[data_with_age$feature=="222_myometrium",]#For age

d_max<-d[d$value==max(d$value), ]
d_min<-d[d$value==min(d$value), ]

colors<-"#6DB6FFFF"
total<-d

f<-ggplot(total, aes(x = Age, y = value, color = view))+#, group = "view")) +
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
        strip.text = element_text(size = 15, color = "black"))  # Increase facet label size)+
f
