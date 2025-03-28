#Fig 2E
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library("MOFA2")
library("MOFAdata")

tissue_colors <- c(
  "Uterus" = "#6DB6FFFF",
  "Ovary" = "#009292FF",
  "Vagina" = "#B66DFFFF",
  "BreastMammaryTissue" = "#490092FF",
  "Endocervix" = "#FF6DB6FF",
  "Ectocervix" = "#FFB6DBFF",
  "FallopianTube" = "#DB6D00FF"
)
input<- "X/"
tissue<-"Vagina"

#Read variance partition results
m<-readRDS("X/Allal/Differential_expression/old_varpar/variance_partition/results/mean_tiles/Vagina/epithelium.rds")
max_feat_epi<-rownames(m[order(m$age, decreasing = TRUE),])[1]

#Read MOFA Object only to take the feature values needed
MOFAobject_trained<-readRDS(paste0(input, "/Laura/MOFA/MOFA_", tissue, "_only_matching_donors_10_ensg_CORRECTED.rds"))
metadata_tot <- readRDS(paste0(input, "/Laura/00.Data/v10/", tissue, "/metadata.rds"))
filter<-readRDS(paste0(input, "/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))

#For using continuous ancestry as a covariate (otherwise we use it as a categorical variable)
metadata_tot$Ancestry<- NULL
ancestry_file<- read.table(paste0(input,"/Laura/00.Data/admixture_inferred_ancestry.txt"))
ances<-ancestry_file[,c(1,3)]
colnames(ances)<- c("Donor", "Ancestry")
metadata_tot<-merge(metadata_tot, ances, by= "Donor")
metadata<- metadata_tot[metadata_tot$Donor %in% filter$Subject.ID,]
metadata_MOFA<-metadata
colnames(metadata_MOFA)[1]<-"sample"
metadata_MOFA<- metadata_MOFA[metadata_MOFA$sample %in% colnames(MOFAobject_trained@data[[1]]$group1),]
metadata_MOFA$Sample<-NULL
ordered_metadata <- metadata_MOFA[order(metadata_MOFA$Age), ]
samples_metadata(MOFAobject_trained) <- metadata_MOFA

#Select the view we need
data <- get_data(MOFAobject_trained, 
                 views = "epithelium", 
                 as.data.frame = TRUE
)

data_with_age<-merge(data, metadata_MOFA[,c("sample","Age", "Ancestry", "BMI")], by= "sample")

feature_ancestry<-"213_epithelium"
feature_age<-"313_epithelium"

d<-data_with_age[data_with_age$feature=="213_epithelium",]
filter<-readRDS(paste0("X/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
d<-d[d$sample %in% filter$Subject.ID,]

#Order by feature value and select max and mins
d<-d[order(d$value, decreasing = TRUE),]
d_max<-d[d$value==max(d$value), ]
d_min<-d[d$value==min(d$value), ]
epithelium<-d

tail(epithelium)


tissue<-"BreastMammaryTissue"

#Read variance partition results
m<-readRDS("X/Allal/Differential_expression/variance_partition/results/mean_tiles/BreastMammaryTissue/adipocyte.rds")
max_feat_cortex<-rownames(m[order(m$bmi, decreasing = TRUE),])[1]

#Read MOFA Object only to take the feature values needed
MOFAobject_trained<-readRDS(paste0(input, "/Laura/MOFA/MOFA_", tissue, "_only_matching_donors_10_ensg_CORRECTED.rds"))
metadata_tot <- readRDS(paste0(input, "/Laura/00.Data/v10/", tissue, "/metadata.rds"))
filter<-readRDS(paste0(input, "/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))


#For using continuous ancestry as a covariate (otherwise we use it as a categorical variable)
metadata_tot$Ancestry<- NULL
ancestry_file<- read.table(paste0(input,"/Laura/00.Data/admixture_inferred_ancestry.txt"))
ances<-ancestry_file[,c(1,3)]
colnames(ances)<- c("Donor", "Ancestry")
metadata_tot<-merge(metadata_tot, ances, by= "Donor")
metadata<- metadata_tot[metadata_tot$Donor %in% filter$Subject.ID,]

###Age bins
metadata$Age_bins <- ifelse(metadata$Age<=35, "young",
                            ifelse(metadata$Age>35 & metadata$Age<60, "middle", "old"))
metadata$Age_bins<-as.factor(metadata$Age_bins)

###Ancestry bins
metadata$Ancestry_bins <-cut(metadata$Ancestry, 
                             breaks = 3, 
                             labels = c("Low", "Medium","High"),
                             include.lowest = TRUE)
metadata$Ancestry_bins<-as.factor(metadata$Ancestry_bins)

###BMI

metadata_MOFA<-metadata
colnames(metadata_MOFA)[1]<-"sample"
metadata_MOFA<- metadata_MOFA[metadata_MOFA$sample %in% colnames(MOFAobject_trained@data[[1]]$group1),]
metadata_MOFA$Sample<-NULL
ordered_metadata <- metadata_MOFA[order(metadata_MOFA$Age), ]
samples_metadata(MOFAobject_trained) <- metadata_MOFA
data <- get_data(MOFAobject_trained, 
                 views = "adipocyte", 
                 as.data.frame = TRUE
)
nrow(data)
data_with_age<-merge(data, metadata_MOFA[,c("sample","Age", "Ancestry", "BMI")], by= "sample")
d<-data_with_age[data_with_age$feature=="90_adipocyte",]
filter<-readRDS(paste0("X/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
d<-d[d$sample %in% filter$Subject.ID,]

#Order by feature value and select max and mins
d_max<-d[d$value==max(d$value), ]
d_min<-d[d$value==min(d$value), ]
d<-d[order(d$value, decreasing= FALSE),]

tail(d)
adipo<-d

c <- c("#FFFF6DFF","#A47F7F", "#490092FF") 
epi<-ggplot(epithelium, aes(x = Ancestry, y = value, color = Ancestry, fill = Ancestry)) +
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.2,alpha=0.2, size = 1) +
  # geom_point(alpha = 0.5, size = 1.5) +
  # geom_smooth(method = "loess", se = TRUE, color = "black", fill = "grey80", alpha = 0.5)  +
  scale_fill_manual(values = c) +
  scale_color_manual(values = c) +
  # scale_x_discrete(labels = c("Low [0,0.33] " = "Low", 
  #                             "Medium (0.33,0.66]" = "Medium", 
  #                             "(0.66, 1]" = "High"))+
  labs(title= " ", x = "Ancestry", y = "Feature value") +  # Label axes
  theme_minimal()+
  theme(axis.text.y = element_text(size = 15, color = "black"),
        # axis.text.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        
        plot.title = element_text(size = 15, hjust = 0.5, color = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        legend.title = element_text(size = 15),
        panel.border = element_rect(color = "grey80", size = 0.5, fill = NA),
        strip.text = element_text(size = 15, color = "black"))  # Increase facet label size)+
# facet_wrap(~view, scales = "free_y")
# facet_wrap(~view, scales = "free_y", nrow=1)

br_bmi<- ggplot(adipo, aes(x = BMI,color="view", y = value)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "loess", se = TRUE, color = "black", fill = "grey80", alpha = 0.5)  +
  scale_color_manual(values= tissue_colors[[4]]) +
  labs(title= " ", x = "BMI", y = "Feature value") +  # Label axes
  theme_minimal()+
  theme(axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        plot.title = element_text(size = 15, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        legend.position = "bottom",
        panel.grid = element_blank(),
        legend.title = element_text(size = 15),
        panel.border = element_rect(color = "grey80", size = 0.5, fill = NA),
        strip.text = element_text(size = 15, color = "black"))  # Increase facet label size)+
# facet_wrap(~view, scales = "free_y", nrow=4)
# facet_wrap(~view, scales = "free_y", nrow=4)
epi <- ggplot(epithelium, aes(x = Ancestry, y = value, fill = Ancestry)) +
  geom_violin(color = "black", alpha = 0.6, size = 0.5) +  # Black border for violin
  geom_boxplot(width = 0.2, alpha = 0.2, size = 0.5) +  # Boxplot retains its default border color
  
  scale_fill_manual(values = c) +  # Colors for ancestry
  labs(title= " ", x = "Ancestry", y = "Feature value") +  
  theme_minimal() +
  theme(axis.text.y = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, color = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        legend.title = element_text(size = 15),
        panel.border = element_rect(color = "grey80", size = 0.5, fill = NA),
        strip.text = element_text(size = 15, color = "black"))  

ej_features<-epi+br_bmi+plot_layout(nrow = 2, heights = unit(1, "null"))
ej_features <- (epi + theme(plot.margin = margin(2, 2, 2, 2))) + 
  (br_bmi + theme(plot.margin = margin(2, 2, 2, 2))) +
  plot_layout(nrow = 2, heights = c(1, 1))
ej_features
# Print the plot


pdf("~/Desktop/Figures/Fig2E.pdf", width = 3.5, height = 7.5)  # Adjust width and height as needed
ej_features
dev.off()

svg("~/Desktop/Figures/Fig2E.svg",  width = 3.5, height = 7.5)
ej_features
dev.off()


###Tiles examples
input<- "X/"
tissue<-"FallopianTube"

#Read variance partition results
m<-readRDS("X/Allal/Differential_expression/old_varpar/variance_partition/results/mean_tiles/FallopianTube/epithelium.rds")
max_feat_epi<-rownames(m[order(m$age, decreasing = TRUE),])[1]

#Read MOFA Object only to take the feature values needed
MOFAobject_trained<-readRDS(paste0(input, "/Laura/MOFA/MOFA_", tissue, "_only_matching_donors_4_ensg_CORRECTED.rds"))
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

#Select the view we need
data <- get_data(MOFAobject_trained, 
                 views = "epithelium", 
                 as.data.frame = TRUE
)

data_with_age<-merge(data, metadata_MOFA[,c("sample","Age", "Ancestry", "BMI")], by= "sample")
d<-data_with_age[data_with_age$feature=="313_epithelium",]
filter<-readRDS(paste0("X/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
d<-d[d$sample %in% filter$Subject.ID,]

#Order by feature value and select max and mins
d<-d[order(d$value, decreasing = TRUE),]
d_max<-d[d$value==max(d$value), ]
d_min<-d[d$value==min(d$value), ]

#Endocervix tiles
input<- "X/"
tissue<-"CervixEndocervix"

#Read variance partition results
m<-readRDS("X/Allal/Differential_expression/old_varpar/variance_partition/results/mean_tiles/Endocervix/glandular_epithelium.rds")
max_feat_epi<-rownames(m[order(m$age, decreasing = TRUE),])[1]

#Read MOFA Object only to take the feature values needed
MOFAobject_trained<-readRDS(paste0(input, "/Laura/MOFA/MOFA_", tissue, "_only_matching_donors_4_ensg_CORRECTED.rds"))
metadata_tot <- readRDS(paste0(input, "/Laura/00.Data/v10/", tissue, "/metadata.rds"))
filter<-readRDS(paste0(input, "/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
metadata_tot$Ancestry<- NULL
ancestry_file<- read.table(paste0(input,"/Laura/00.Data/admixture_inferred_ancestry.txt"))
ances<-ancestry_file[,c(1,3)]
colnames(ances)<- c("Donor", "Ancestry")
metadata_tot<-merge(metadata_tot, ances, by= "Donor")
metadata<- metadata_tot[metadata_tot$Donor %in% filter$Subject.ID,]

#Add MOFA metadata
metadata_MOFA<-metadata
colnames(metadata_MOFA)[1]<-"sample"
metadata_MOFA<- metadata_MOFA[metadata_MOFA$sample %in% colnames(MOFAobject_trained@data[[1]]$group1),]
metadata_MOFA$Sample<-NULL
ordered_metadata <- metadata_MOFA[order(metadata_MOFA$Age), ]
samples_metadata(MOFAobject_trained) <- metadata_MOFA

#Select the view we need
data <- get_data(MOFAobject_trained, 
                 views = "glandular_epithelium", 
                 as.data.frame = TRUE
)

data_with_age<-merge(data, metadata_MOFA[,c("sample","Age", "Ancestry", "BMI")], by= "sample")
d<-data_with_age[data_with_age$feature=="313_glandular_epithelium",]
filter<-readRDS(paste0("X/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
d<-d[d$sample %in% filter$Subject.ID,]

#Order by feature value and select max and mins
d<-d[order(d$value, decreasing = TRUE),]
d_max<-d[d$value==max(d$value), ]
d_min<-d[d$value==min(d$value), ]



#Ectocervix tiles
tissue<-"CervixEctocervix"

#Read variance partition results
m<-readRDS("X/Allal/Differential_expression/old_varpar/variance_partition/results/mean_tiles/Ectocervix/epithelium.rds")
max_feat_epi<-rownames(m[order(m$age, decreasing = TRUE),])[1]

#Read MOFA Object only to take the feature values needed
MOFAobject_trained<-readRDS(paste0(input, "/Laura/MOFA/MOFA_", tissue, "_only_matching_donors_4_ensg_CORRECTED.rds"))
metadata_tot <- readRDS(paste0(input, "/Laura/00.Data/v10/", tissue, "/metadata.rds"))
filter<-readRDS(paste0(input, "/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))


#For using continuous ancestry as a covariate (otherwise we use it as a categorical variable)
metadata_tot$Ancestry<- NULL
ancestry_file<- read.table(paste0(input,"/Laura/00.Data/admixture_inferred_ancestry.txt"))
ances<-ancestry_file[,c(1,3)]
colnames(ances)<- c("Donor", "Ancestry")
metadata_tot<-merge(metadata_tot, ances, by= "Donor")
metadata<- metadata_tot[metadata_tot$Donor %in% filter$Subject.ID,]
metadata_MOFA<-metadata
colnames(metadata_MOFA)[1]<-"sample"
metadata_MOFA<- metadata_MOFA[metadata_MOFA$sample %in% colnames(MOFAobject_trained@data[[1]]$group1),]
metadata_MOFA$Sample<-NULL
ordered_metadata <- metadata_MOFA[order(metadata_MOFA$Age), ]
samples_metadata(MOFAobject_trained) <- metadata_MOFA

#Select the view we need
data <- get_data(MOFAobject_trained, 
                 views = "epithelium", 
                 as.data.frame = TRUE
)
data_with_age<-merge(data, metadata_MOFA[,c("sample","Age", "Ancestry", "BMI")], by= "sample")
d<-data_with_age[data_with_age$feature=="313_epithelium",]
filter<-readRDS(paste0("X/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
d<-d[d$sample %in% filter$Subject.ID,]

#Order by feature value and select max and mins
d<-d[order(d$value, decreasing = TRUE),]
d_max<-d[d$value==max(d$value), ]
d_min<-d[d$value==min(d$value), ]


