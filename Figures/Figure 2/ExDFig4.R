#Supplementary figures
#A --- Elastic-Net accs
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

pdf("X/Figures/Fig4S_A.pdf", width = 6.5, height = 4.5) 
p
dev.off()




#B --substructures' proportions
library(dplyr)
library(tidyverse)
tis<-c("Uterus", "Ovary", "Vagina", "BreastMammaryTissue", "CervixEndocervix", "CervixEctocervix", "FallopianTube")

organ_list<-list()
for (tissue in tis){
  
  if(tissue =="BreastMammaryTissue"){
    cell_prop<-read.csv(paste0("~/X/Laura/derived_proportions_Craig/", tissue, "/", tissue, "_compositional.csv"))
    substructures<- c("adipocyte", "lobule", "duct", "stroma", "nerve", "gynecomastoid_hyperplasia")
    
  }else if(tissue == "Ovary"){
    cell_prop<-read.csv(paste0("~/X/Laura/12.Tissues_substructures/", tissue, "_no_follicles_compositional.csv"))
    substructures<- c("cortex", "corpora", "medulla", "vessels")
    
  }else if (tissue =="Vagina"){
    cell_prop<-read.csv(paste0("~/X/Laura/12.Tissues_substructures/", tissue, "_compositional.csv"))
    substructures<- c("epithelium", "lamina_propria", "vessels", "stroma")
    
  }else if (tissue =="Uterus"){
    cell_prop<-read.csv(paste0("~/X/Laura/12.Tissues_substructures/", tissue, "_compositional.csv"))
    substructures<- c("myometrium", "endometrium", "vessels")
    
    
  }else if (tissue =="CervixEndocervix"){
    cell_prop<-read.csv(paste0("~/X/Laura/12.Tissues_substructures/Endocervix_compositional.csv"))
    substructures<- c("vessels", "glandular_epithelium", "stroma")
    
  }else if (tissue =="CervixEctocervix"){
    cell_prop<-read.csv(paste0("~/X/Laura/12.Tissues_substructures/Ectocervix_compositional.csv"))
    substructures<- c("epithelium", "stroma", "vessels", "glands")
    
    
  }else if (tissue == "FallopianTube"){
    cell_prop<-read.csv(paste0("~/X/Laura/12.Tissues_substructures/FallopianTube_compositional.csv"))
    substructures<- c("vessels", "lumen", "smooth_muscle", "epithelium", "stroma") ##always control for adipose, since it is contamination
    
  }
  
  metadata_tot <- readRDS(paste0("~/X/Laura/00.Data/v10/", tissue, "/metadata.rds"))
  filter<-readRDS(paste0("X/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
  metadata<- metadata_tot[metadata_tot$Donor %in% filter$Subject.ID,]
  cell_prop$Donor<-  gsub("^(\\w+-\\w+).*", "\\1", cell_prop$slide_id)
  metadata<- merge(metadata, cell_prop, by="Donor")
  m<-metadata[, c("Donor", "Age", substructures)]
  df_long <- m %>%
    pivot_longer(cols = substructures, 
                 names_to = "Subtissue", 
                 values_to = "Proportion")
  df_long$Tissue<- tissue
  organ_list<-rbind(organ_list, df_long)
  
}
pal <- c( "#6DB6FFFF", "#009292FF",  "#B66DFFFF", "#490092FF",   "#FF6DB6FF",
          "#FFB6DBFF", "#DB6D00FF")

organ_list$Organ<-factor(organ_list$Tissue, levels= tis)
organ_list$Tissue<-NULL
colnames(organ_list)[3]<-"Structure"
library(ggplot2)
library(patchwork)
d<-organ_list

wb_s<-createWorkbook()
d$Donor <- paste0("Donor", as.numeric(factor(d$Donor)))
d$Age<- ifelse(d$Age >=20 & d$Age, "20-29",
               ifelse(d$Age >=30 & d$Age<40, "30-39",
                      ifelse(d$Age >=40 & d$Age<50, "40-49",
                             ifelse(d$Age >=50 & d$Age<60, "50-59",
                                    ifelse(d$Age >=60 & d$Age<70, "60-69","70-79")))))

sheet_name<-"ExDFig4A"
addWorksheet(wb_s, sheet_name)
writeData(wb_s, sheet_name, d, startCol = 1, startRow = 3)  # Title

# Create individual plots for each tissue
plots <- lapply(unique(organ_list$Organ), function(tissue) {
  ggplot(subset(organ_list, Tissue == tissue), aes(x = Age, y = Proportion, color = Structure, group = Structure)) +
    geom_point(alpha = 0.2, size = 1.5) +
    geom_smooth(method = "loess", se = TRUE, aes(color = Structure, fill = Structure)) +
    scale_color_manual(values = pal) +  # Use custom colors but let ggplot assign them
    scale_fill_manual(values = pal) +   # Same for fill
    labs(title = paste(tissue), x = "Age (years)", y = "Proportion") +
    theme_minimal() +
    theme(    #  legend.position = "right",  # Keep legend on the right
      legend.justification = "left",  # Move legend closer to the plot
      legend.margin = margin(r = -20, l = 0, t = 0, b = 0),  # Reduce right margin
      legend.title= element_text(size=12, color= "black"),
      axis.text.y = element_text(size = 14, color = "black"), 
      axis.text.x = element_text(size = 14, color = "black"),
      axis.title = element_text(size = 14, color = "black"), 
      plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
      legend.text = element_text(size = 12, color = "black"),
      legend.position = "bottom",
      legend.direction = "vertical",  # Make the legend vertical
      panel.grid = element_blank(),  
      panel.border = element_rect(color = "grey80", size = 0.5, fill = NA))
})

# Combine all plots horizontally or in a grid
final_plot <- wrap_plots(plots, ncol = 7)  # Adjust ncol to change layout
###Test the changes with age

#Chose the tissue
cell_prop<- read.csv(paste0("~/X/Laura/12.Tissues_substructures/", tissue, "_pivot.csv"))  decon_results<-cell_prop
decon_results$Donor<- as.factor(gsub("^(\\w+-\\w+).*", "\\1", decon_results$slide_id))
decon_results$slide_id<-NULL
decon_results <- aggregate(. ~ Donor, data = decon_results, FUN = mean, na.rm = TRUE)decon_results<- as.data.frame(decon_results)
rownames(decon_results)<-decon_results$Donor#To input into the function: ROWS-CELL TYPES, COLUMNS-DONORS
df<-decon_results
decon_results$Donor<-NULL
DF_proportions<-(t(decon_results))sex <- "female"
metadata<-readRDS(paste0("~/X/Laura/00.Data/v10/", tissue, "/metadata.rds"))
filter<-readRDS(paste0("X/Laura/03.Image_processing/Second_filtering_images/", tissue, "_final_filtered_images.rds"))
metadata<- metadata[metadata$Donor %in% filter$Subject.ID,]
metadata<- merge(metadata, df, by= "Donor")
meta<- metadata

rownames(meta)<-meta$Donor
meta$donor<- meta$Donor
meta$Donor<-NULL
meta<- meta[, !colnames(meta) %in% c("PEER1", "PEER2", "ExonicRate", "RIN")]

lm_by_ct <- function(cell_type, df = DF_proportions, md = meta){
  print(cell_type)
  vec_i <- df[rownames(df)==cell_type, ]  df_i <- as.data.frame(vec_i)
  colnames(df_i) <- 'freq'
  df_i$donor <- rownames(df_i)
  df_i <- merge(df_i, md, by = 'donor')
  rownames(df_i) <- df_i$donor
  # df_i <- df_i[,-1]  df_i$Donor<-NULL
  df_i$donor<-NULL
  covariates <- colnames(md)[colnames(md) %in% c("HardyScale", "IschemicTime","RIN", "Cohort", "NucAcIsoBatch", "ExonicRate")]
  individual_traits<- c("Age", "Ancestry", "BMI")  # Formula  fmla <- paste(c(covariates, individual_traits), collapse = " + ")
  fmla <- paste0('freq ~ ',fmla)
  form <- as.formula(fmla)
  print(paste0('Fitting lmer: ',fmla))  # fit model
  print('lmer...')
  #mod <-  lmerTest::lmer(form, data = df_i)
  mod <- lm(formula = form, data = df_i)
  # tidy model
  tidy_mod <- broom.mixed::tidy(mod, conf.int = TRUE, effects = "fixed")  # tidy to dataframe
  cnames <- c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high")
  tidy_mod <- tidy_mod[,which(colnames(tidy_mod)%in%cnames)]
  tidy_mod.df <- as.data.frame(tidy_mod)
  tidy_mod.df$celltype <- cell_type
  return(tidy_mod.df)
}
## Apply function
tidy_mod.list <- sapply(rownames(DF_proportions), function(i) lm_by_ct(i), simplify = FALSE)
tidy_mod.df <- do.call("rbind", tidy_mod.list)# Split by phenotype and compute FDR
tidy_mod.by_phe <- split(tidy_mod.df, tidy_mod.df$term)
tidy_mod.Age <- tidy_mod.by_phe$Age
tidy_mod.Age$fdr <- p.adjust(tidy_mod.Age$p.value, 'BH')
tidy_mod.Age <- tidy_mod.Age[order(tidy_mod.Age$fdr),]table(tidy_mod.Age$fdr < 0.05)


# Print combined figure
pdf("X/Laura/Figure_plots/figS_Notes_A.pdf", width = 10, height = 4) 
final_plot
dev.off()

pdf("X/figS_Notes_A.pdf", width = 18, height = 5) 
final_plot
dev.off()




#C---Uterus samples
tissue<-"Uterus"
myo<- readRDS(paste0("X/Laura/03.Image_processing/Second_filtering_images/Uterus_filtered_myometrium_images.rds"))
tot<-readRDS(paste0("X/Laura/03.Image_processing/Second_filtering_images/", tissue, "_final_filtered_images.rds"))
metadata<-readRDS(paste0("~/X/Laura/00.Data/v10/", tissue, "/metadata.rds"))
metadata_myo<- metadata[metadata$Donor %in% myo$Subject.ID,]
metadata_tot<-metadata[metadata$Donor %in% tot$Subject.ID,]
ut<-data.frame(
  Number=c(nrow(myo),nrow(tot),nrow(metadata_myo), nrow(metadata_tot)),
  Samples=c("Only myometrium", "Both", "Only myometrium", "Both"),
  Data=c("Images","Images", "RNA-Seq", "RNA-Seq")
)
ut$Data<-factor(ut$Data)
ut$Samples<-factor(ut$Samples)
cols<-c("#B66DFFFF", "#490092FF")
sheet_name<-"SNFig1B"
addWorksheet(wb_s, sheet_name)
writeData(wb_s, sheet_name, ut, startCol = 1, startRow = 3)  # Title

bar_plot <- ggplot(ut, aes(y = Number, x = Data, fill = Samples)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # Do not sum, just use single values  scale_fill_identity() +  # Use the custom colors with alpha transparency
  geom_text(aes(label = Number), 
            position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 5, color = "black") +  # Adjust vertical position
  
  theme_minimal() +
  scale_fill_manual(values=cols)+
  theme(axis.text.y = element_text(size = 15, color = "black"), 
        axis.text.x = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 15, color = "black"), 
        axis.title.y = element_blank(),
        
        plot.title = element_text(size = 15, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        panel.grid = element_blank(),    
        legend.title = element_text(size = 15),
        panel.border = element_rect(color = "grey80", size = 0.5, fill = NA),
        legend.position = "bottom"
  )

# Display the plot
bar_plot
pdf("X/ExDFig4C.pdf", width =5, height = 5 )
bar_plot
dev.off()


#D -- endometrium variability in all age groups
tissue<-"Uterus"
cell_prop<-read.csv(paste0("~/X/Laura/12.Tissues_substructures/", tissue, "_compositional.csv"))
substructures<- c("myometrium", "endometrium", "vessels")

metadata_tot <- readRDS(paste0("~/X/Laura/00.Data/v10/", tissue, "/metadata.rds"))
filter<-readRDS(paste0("X/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
metadata<- metadata_tot[metadata_tot$Donor %in% filter$Subject.ID,]
cell_prop$Donor<-  gsub("^(\\w+-\\w+).*", "\\1", cell_prop$slide_id)
metadata<- merge(metadata, cell_prop, by="Donor")
m<-metadata[, c("Donor", "Age", substructures)]
m<-m[,c("Donor", "Age", "endometrium")]

m$Agebin<-ifelse(m$Age<=35, "Young", 
                 ifelse(m$Age>35 & m$Age<60, "Middle", "Old") )
m$Tissue<-tissue
m$Agebin<-factor(m$Agebin, levels=c("Young", "Middle", "Old"))
c <- c("#FFFF6DFF","#A47F7F", "#490092FF") 
sheet_name<-"SNFig1C"
addWorksheet(wb_s, sheet_name)
d<-m
d$Donor <- paste0("Donor", as.numeric(factor(d$Donor)))
d$Age<- ifelse(d$Age >=20 & d$Age, "20-29",
               ifelse(d$Age >=30 & d$Age<40, "30-39",
                      ifelse(d$Age >=40 & d$Age<50, "40-49",
                             ifelse(d$Age >=50 & d$Age<60, "50-59",
                                    ifelse(d$Age >=60 & d$Age<70, "60-69","70-79")))))

write.csv(d,"Descargas/endom.csv")
writeData(wb_s, sheet_name, d, startCol = 1, startRow = 3)  # Title

p <- ggplot(m, aes(x = Agebin, y = endometrium, fill = Agebin, color = Agebin)) +
  geom_violin(alpha = 0.4) +
  geom_jitter(width = 0.2, alpha = 1, size = 1) +
  #facet_wrap(~ Substructure, scales = "free_y") +
  scale_fill_manual(values = c) +
  scale_color_manual(values = c) +
  labs(
    title = " ",
    x = "Age",
    y = "Endometrium proportion"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 13, color = "black"), 
        axis.text.x = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13, color = "black"), 
        plot.title = element_text(size = 13, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        panel.grid = element_blank(),                    # Remove all grid lines
        panel.border = element_rect(color = "grey80", size = 0.5, fill = NA) # Add full border
  )
p

pdf("X/ExtDFigD.pdf", width = 5, height = 4.5) 
p
dev.off()


#EFG - Cell type deconvolution

shhh(library(plyr))
shhh(library(dplyr))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(ggplot2))
shhh(library(scater))
shhh(library(RColorBrewer))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(ggplot2))


# read in files 
tissue<- "Uterus"
if(tissue == "Ovary"){
  cell_prop<-read.csv("X/Laura/gtm-decon/Deconvolution_ovary_pivot_v10.csv", row.names = "X")
  decon_results<-cell_prop
  colnames(decon_results)<-c('Endothelial cells', 'Epithelial cells','Granulosa cells', 'Immune cells', 'Stromal cells',  'Smooth muscle cells','Theca cells')
  
}else{
  cell_prop<-read.csv("X/Laura/gtm-decon/Deconvolution_uterus_pivot_v10.csv", row.names = "X")
  decon_results<-cell_prop
  colnames(decon_results)<-c('B cells', 'Endothelial cells','Epithelial cells', 'Fibroblasts', 'Macrophages',  'NK cells','Pericytes', 'Smooth muscle cells', 'T cells') 
  
}
decon_results$Donor<- as.factor(gsub("^(\\w+-\\w+).*", "\\1", rownames(decon_results)))
rownames(decon_results)<-decon_results$Donor

decon_results<- as.data.frame(decon_results)
rownames(decon_results)<-decon_results$Donor

#TO INOUT INTO THE FUNCTION: ROWS-CELL TYPES, COLUMNS-DONORS
df<-decon_results
decon_results$Donor<-NULL
DF_proportions<-(t(decon_results))

sex <- "female"
metadata<-readRDS(paste0("~/X/Laura/00.Data/v10/", tissue, "/metadata.rds"))
if (tissue =="Uterus"){
  filter<- readRDS(paste0("X/Laura/03.Image_processing/Second_filtering_images/Uterus_filtered_myometrium_images.rds"))
  
}else{
  filter<-readRDS(paste0("X/Laura/03.Image_processing/Second_filtering_images/", tissue, "_final_filtered_images.rds"))
  
}
nrow(metadata)

metadata<- metadata[metadata$Donor %in% filter$Subject.ID,]
metadata<- merge(metadata, df, by= "Donor")
meta<- metadata
rownames(meta)<-meta$Donor
meta$donor<- meta$Donor
meta$Donor<-NULL
meta<- meta[, !colnames(meta) %in% c("PEER1", "PEER2")]

lm_by_ct <- function(cell_type, df = DF_proportions, md = meta){
  
  print(cell_type)
  vec_i <- df[rownames(df)==cell_type, ]
  df_i <- as.data.frame(vec_i)
  colnames(df_i) <- 'freq'
  df_i$donor <- rownames(df_i)
  df_i <- merge(df_i, md, by = 'donor')
  rownames(df_i) <- df_i$donor
  
  df_i$Donor<-NULL
  df_i$donor<-NULL  
  covariates <- colnames(md)[colnames(md) %in% c("HardyScale", "IschemicTime","RIN", "Cohort", "NucAcIsoBatch", "ExonicRate")]
  individual_traits<- c("Age", "Ancestry", "BMI")
  
  # Formula
  
  fmla <- paste(c(covariates, individual_traits), collapse = " + ")
  fmla <- paste0('freq ~ ',fmla)
  form <- as.formula(fmla)
  print(paste0('Fitting lmer: ',fmla))
  
  # fit model
  print('lm...')
  mod <- lm(formula = form, data = df_i)
  # tidy model
  tidy_mod <- broom.mixed::tidy(mod, conf.int = TRUE, effects = "fixed")
  
  # tidy to dataframe
  cnames <- c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high")
  tidy_mod <- tidy_mod[,which(colnames(tidy_mod)%in%cnames)]
  tidy_mod.df <- as.data.frame(tidy_mod)
  tidy_mod.df$celltype <- cell_type
  return(tidy_mod.df)
}
## Apply function
tidy_mod.list <- sapply(rownames(DF_proportions), function(i) lm_by_ct(i), simplify = FALSE)
tidy_mod.df <- do.call("rbind", tidy_mod.list)

# Split by phenotype and compute FDR
tidy_mod.by_phe <- split(tidy_mod.df, tidy_mod.df$term)
tidy_mod.Age <- tidy_mod.by_phe$Age
tidy_mod.Age$fdr <- p.adjust(tidy_mod.Age$p.value, 'BH')
tidy_mod.Age <- tidy_mod.Age[order(tidy_mod.Age$fdr),]

table(tidy_mod.Age$fdr < 0.05)


tidy_mod.Age$direction <- ifelse(tidy_mod.Age$estimate>0, 'pos', 'neg')
###plot
coda<- tidy_mod.Age
coda$significance <- ifelse(coda$p.value< 0.01, "ss","ns")
coda$celltype <- gsub("_", " ", coda$celltype)
coda$direction <- ifelse(coda$estimate < 0, "down", "up")

alpha_vec = c(ns = 0.5, ss = 1)
library(scales)


sheet_name<-"SNFig1K"
addWorksheet(wb_s, sheet_name)
writeData(wb_s, sheet_name, coda, startCol = 1, startRow = 3)  # Title

p_uter_all <- ggplot(coda, aes(x = estimate, y = celltype)) + 
  geom_point(aes(alpha = significance, color = direction), size = 4) + 
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high, alpha = significance, color = direction), fatten = 0.1, size=0.6) +
  geom_vline(xintercept = 0, linetype = "dashed",color = "grey50", size = 0.7) +
  scale_color_manual(values = c("up" = "#63B8FF", "down" = "#8B0000")) + 
  
  scale_alpha_manual(values = alpha_vec) +
  #facet_grid(celltype ~ ., space = "free", scale = "free") +
  
  scale_x_continuous(
    limits = c(-0.025, 0.025), 
    breaks = seq(-0.025, 0.025, by = 0.05),
    labels = label_number(accuracy = 0.02)  # Set the number of decimal places
  )  +
  xlab("CoDA estimate") +
  ylab(NULL)+
  # labs(title = "Ovary (n=169)")+
  
  labs(title = "Uterus (n=61)")+
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 13, color = "black"), 
        axis.text.x = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13, color = "black"), 
        plot.title = element_text(size = 13, hjust = 1, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.position = "top", 
        legend.key.size = unit(0.5, "lines"),
        legend.title = element_blank(),
        panel.grid = element_blank(),                    # Remove all grid lines
        axis.line.x = element_line(size = 0.5, color = "grey80"), # Add x-axis line
        axis.line.y = element_line(size = 0.5, color = "grey80"),
        panel.border = element_rect(color = "grey80", size = 0.5, fill = NA) # Add full border
  )
p_uter_all
decon_plots<-p_uter_all+ p_uter +p_ova

pdf("X/uteruscoda.pdf", width =5, height = 5 )
p_uter_all
dev.off()

pdf("X/Laura/Figure_plots/figS_Notes_C.pdf", width = 5, height = 5) 
p
dev.off()
