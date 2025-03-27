#FigS1E
#Correlations
##TEST THE CORRELATIONS
library(car)
library(ComplexHeatmap)
library(brew)
library(ggplot2)
library(RColorBrewer)
library(corrplot)
plots_corr<-list()
tissues<-c("Uterus", "Ovary", "Vagina", "BreastMammaryTissue")

for (tissue in tissues){
  metadata_tot <- readRDS(paste0("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/00.Data/v10/", tissue, "/metadata.rds"))
  filter<-readRDS(paste0("Desktop/TFM/bsc83671/GTEx_v8/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
  #filter<- readRDS(paste0("Desktop/TFM/bsc83671/GTEx_v8/Laura/03.Image_processing/", tissue, "_filtered_MYO_AND_ENDO_images.rds"))
  
  metadata<- metadata_tot[metadata_tot$Donor %in% filter$Subject.ID,]
  metadata$Ancestry <- relevel(metadata$Ancestry, ref="EUR")  
  
  #Read substructures' proportions
  if(tissue %in% both_tissues){
    tissue<- "BreastMammaryTissue"
    cell_prop<-read.csv(paste0("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/derived_proportions_Craig/", tissue, "/", tissue, "_pivot.csv"))
    substructures<- c("adipocyte", "lobule")
    
  }else if(tissue == "Ovary"){
    cell_prop<-read.csv(paste0("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/12.Tissues_substructures/", tissue, "_no_follicles_pivot.csv"))
    substructures<- c("cortex", "corpora", "medulla")
    
  }else if (tissue =="Vagina"){
    cell_prop<-read.csv(paste0("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/12.Tissues_substructures/", tissue, "_pivot.csv"))
    substructures<- c("epithelium", "lamina_propria")
    
  }else if (tissue =="Uterus"){
    cell_prop<-read.csv(paste0("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/12.Tissues_substructures/", tissue, "_pivot.csv"))
    #cell_prop<-read.csv(paste0("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/12.Tissues_substructures/Uterus_atrophy_pivot.csv"))
    
    substructures<- c("myometrium", "endometrium")
    #substructures<- c("myometrium", "endometrium","atrophied_myometrium", "atrophied_endometrium")
    
  }
  cell_prop$Donor<-  gsub("^(\\w+-\\w+).*", "\\1", cell_prop$slide_id)
  metadata<- merge(metadata, cell_prop, by="Donor")
  peers<- c("PEER1", "PEER2")
  # Calculate Pearson correlation
  
  correlations <- data.frame(matrix(nrow = length(substructures), ncol = length(peers)))
  rownames(correlations) <- substructures
  colnames(correlations) <- peers
  
  # Loop through substructures and peers to calculate Pearson correlation
  for (columns in substructures) {
    for (peer in peers) {
      # Calculate Pearson correlation and store it in the dataframe
      correlations[columns, peer] <- cor(metadata[[peer]], metadata[[columns]], method = "spearman")
    }
  }
  
  normality_results <- list()
  
  # Check normality for PEERs
  for (peer in peers) {
    normality_results[[peer]] <- shapiro.test(metadata[[peer]])$p.value
  }
  
  # Check normality for substructures
  for (columns in substructures) {
    normality_results[[columns]] <- shapiro.test(metadata[[columns]])$p.value
  }
  
  all_corr<- list()
  all_corr[[tissue]]<- correlations
  
  peers<- c("PEER1", "PEER2")
  tec<- colnames(metadata)[colnames(metadata) %in% c("HardyScale","IschemicTime","RIN","ExonicRate","Cohort", "NucIsoBatch")]
  covariates <- c(substructures, peers)#, tec)
  individual_traits <- c("Age", "Ancestry", "BMI")
  
  fml_args_mod <- paste(c(covariates, individual_traits), collapse = " + ")
  
  # Create the full formula for your model
  full_formula <- as.formula(paste("dummy_response ~  ", paste(fml_args_mod,collapse = " ")))
  metadata$dummy_response <- rnorm(nrow(metadata))
  # This can be a constant or any placeholder
  
  # Fit the model with the dummy response variable
  mod <- lm(full_formula, data = metadata)
  
  # Calculate VIF
  vif_values <- vif(mod)
  vif_values<- as.data.frame(vif_values)
  
  
  # predictor_matrix <- sapply(metadata[, c(covariates,interaction_terms)], as.numeric)
  predictor_matrix <- sapply(metadata[, c(covariates)], as.numeric)
  
  # Calculate the correlation matrix
  cor_matrix <- cor(predictor_matrix, use = "complete.obs",  method = "spearman")
  plots_corr[[tissue]]<-corrplot(cor_matrix, type = "lower", method = "color", tl.col = "black", tl.srt = 20)
  
  pdf(paste0("Desktop/TFM/bsc83671/GTEx_v8/Laura/Figure_plots/", tissue,"_corr.pdf"))
  corrplot(plots_corr[[tissue]]$corr, type = "lower", method = "color", tl.col = "black", tl.srt = 20, 
           tl.cex = 1.5)
  dev.off()
  
  pdf(paste0("Desktop/Figures/", tissue,"_corr.pdf"))
      corrplot(plots_corr[[tissue]]$corr, type = "lower", method = "color", tl.col = "black", tl.srt = 20, 
               tl.cex = 1.5)
  dev.off()
}
