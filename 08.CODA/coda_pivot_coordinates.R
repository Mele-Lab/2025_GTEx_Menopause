#Compositional Analysis (CODA)
#For running this analysis, you need cell-type or structures proportions that have been previously transformed
#We transform proportion percentages into pivot coordinates and use them for this analysis

shhh <- suppressPackageStartupMessages
shhh(library(optparse))
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(ggplot2))
shhh(library(RColorBrewer))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(scales))
shhhlibrary(ggplot2))

# read in files 
tissue<- "Ovary"
experiments <- c("Deconvol", "Substructures")
experiment<- experiments[1]


if (experiment == "Deconvol"){
  cell_prop<-read.csv("~/X/Laura/gtm-decon/Deconvolution_ovary_pivot_v10.csv", row.names = "X")
  decon_results<-cell_prop
  #Uterus
  # colnames(decon_results)<-c('B cells', 'Endothelial cells','Epithelial cells', 'Fibroblasts', 'Macrophages',  'NK cells','Pericytes', 'Smooth muscle cells', 'T cells') 
  #Ovary
   colnames(decon_results)<-c('Endothelial cells', 'Epithelial cells','Granulosa cells', 'Immune cells', 'Stromal cells',  'Smooth muscle cells','Theca cells')
  decon_results$Donor<- as.factor(gsub("^(\\w+-\\w+).*", "\\1", rownames(decon_results)))
  rownames(decon_results)<-decon_results$Donor
  

}else{
  #Breast
  #cell_prop<-read.csv(paste0("~/X/Laura/derived_proportions_Craig/", tissue, "/", tissue, "_pivot.csv"))
  #Ovary
  # cell_prop<- read.csv(paste0("~/X/Laura/12.Tissues_substructures/", tissue, "_no_follicles_pivot.csv"))
  #Uterus/vagina
  cell_prop<- read.csv(paste0("~/X/Laura/12.Tissues_substructures/", tissue, "_pivot.csv"))
  
  decon_results<-cell_prop
  decon_results$Donor<- as.factor(gsub("^(\\w+-\\w+).*", "\\1", decon_results$slide_id))
  decon_results$slide_id<-NULL
  decon_results <- aggregate(. ~ Donor, data = decon_results, FUN = mean, na.rm = TRUE)
}

decon_results<- as.data.frame(decon_results)
rownames(decon_results)<-decon_results$Donor

#TO INOUT INTO THE FUNCTION: ROWS-CELL TYPES, COLUMNS-DONORS
df<-decon_results
decon_results$Donor<-NULL
DF_proportions<-(t(decon_results))
                              
sex <- "female"
metadata_tot<-readRDS(paste0("~/X/Laura/00.Data/v10/", tissue, "/metadata.rds"))

#Add continuous ancestry
metadata_tot$Ancestry<- NULL
ancestry_file<- read.table(paste0(input, "/Laura/00.Data/admixture_inferred_ancestry.txt"))
ances<-ancestry_file[,c(1,3)]
colnames(ances)<- c("Donor", "Ancestry")
metadata_tot<-merge(metadata_tot, ances, by= "Donor")
filter<-readRDS(paste0("/X/Laura/03.Image_processing/Second_filtering_images/", tissue, "_final_filtered_images.rds"))
metadata<- metadata_tot[metadata_tot$Donor %in% filter$Subject.ID,]
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
  print(paste0('Fitting lm: ',fmla))
  
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
tidy_mod.Age$direction <- ifelse(tidy_mod.Age$estimate>0, 'pos', 'neg')

###plot
coda<- tidy_mod.Age
coda$significance <- ifelse(coda$p.value< 0.01, "ss","ns")
coda$celltype <- gsub("_", " ", coda$celltype)
coda$direction <- ifelse(coda$estimate < 0, "down", "up")

alpha_vec = c(ns = 0.5, ss = 1)

p <- ggplot(coda, aes(x = estimate, y = celltype)) + 
  geom_point(aes(alpha = significance, color = direction), size = 4) + 
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high, alpha = significance, color = direction), fatten = 0.1, size=0.6) +
  geom_vline(xintercept = 0, linetype = "dashed",color = "grey50", size = 0.7) +
  scale_color_manual(values = c("up" = "#63B8FF", "down" = "#8B0000")) + 
  
  scale_alpha_manual(values = alpha_vec) +
  scale_x_continuous(
    limits = c(-0.025, 0.025), 
    breaks = seq(-0.025, 0.025, by = 0.05),
    labels = label_number(accuracy = 0.02)  # Set the number of decimal places
  )  +
  xlab("CoDA estimate") +
  ylab(NULL)+
  labs(title = "Cell-type proportion changes with age in ovary (n=169)")+
theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 14, color = "black"), 
        axis.text.x = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, color = "black"), 
        plot.title = element_text(size = 14, hjust = 1, color = "black"),
        legend.text = element_text(size = 14, color = "black"),
        legend.position = "top", 
        legend.key.size = unit(0.5, "lines"),
        legend.title = element_blank(),
        panel.grid = element_blank(),                    # Remove all grid lines
        axis.line.x = element_line(size = 0.5, color = "grey80"), # Add x-axis line
        axis.line.y = element_line(size = 0.5, color = "grey80"),
        panel.border = element_rect(color = "grey80", size = 0.5, fill = NA) # Add full border
  )

print(p)

pdf(paste0("~/X/Laura/12.Tissues_substructures/Figures/Final_v10/", tissue, "_deconvolution_ovary_coda_AGE_pivot.pdf"), width = 5, height = 5)  # Adjust width and height as needed
p
dev.off()

