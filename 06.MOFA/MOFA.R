#####MOFA

library("MOFA2")
library("MOFAdata")
library(dplyr)
library(tidyr)
library(psych)
library(limma)
library(ggplot2)
library(edgeR)
library(readxl)

tissues<- c("Uterus", "Ovary", "Vagina", "BreastMammaryTissue", "FallopianTube", "CervixEctocervix", "CervixEndocervix")

if (getwd()== "/X/Laura" ){
    input<-"/X/"
}else{
  input<-"X/GTEx_v8/"
  
}

log_file<-paste0(input, "Laura/MOFA/log_file_training2.txt")
log_message <- function(message) {
  cat(paste(Sys.time(), "-", message, "\n"), file = log_file, append = TRUE)
}

log_message("Starting the process...\n")

options<- c("all_donors", "only_matching_donors")
option<-options[2]

log_message(paste0("Option: ", option, "\n"))

for (tissue in tissues){
    log_message(paste0("----Model for ", tissue, "\n"))
    
      if (getwd()== "/X/Laura" ){
        log_message("----Reading feature data----\n")
        
          ###Preprocess feature data
          if (tissue== "CervixEctocervix"){
            features<- read.csv(paste0(input, "/Allal/RNAPath/tile_features_balanced_classes_1000/balanced_tile_features_Ectocervix_1000.csv"))
          }else if  (tissue== "CervixEndocervix"){
            features<- read.csv(paste0(input, "/Allal/RNAPath/tile_features_balanced_classes_1000/balanced_tile_features_Endocervix_1000.csv"))
          }else{
            features<- read.csv(paste0(input, "/Allal/RNAPath/tile_features_balanced_classes_1000/balanced_tile_features_", tissue, "_1000.csv"))
            
          }
        
          features$Donor<-gsub("^(\\w+-\\w+).*", "\\1",features$slide_name)
          features$Donor<-as.factor(features$Donor)
          
          #Leave only the donors that passed the filter 
          filter<-readRDS(paste0(input, "/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
          features<-features[levels(features$Donor) %in% filter$Subject.ID,]
          features$slide_name<-NULL
          features$tile_number<-NULL
          
          ##Aggregate features per substructure and per donor
          df_long <- features %>%
            pivot_longer(cols = starts_with("feature"),  # Stack all feature columns
                         names_to = "Feature",
                         values_to = "Value") %>%
            group_by(Donor, label, Feature)%>%
            summarise(across(where(is.numeric), mean, na.rm = FALSE))%>%
            mutate(NewFeature = paste(Feature,label, sep = "_"))
          df_long$label<-NULL
          df_long$Feature<-NULL
          df_transformed <- df_long %>%
            pivot_wider(names_from = Donor, values_from = Value) %>%
            as.data.frame()
  
          #Normalize image features
          feature_data<- scale(df_transformed)
  
          ##Create one matrix per substructure
          list_of_substructures<-list()
          
          df_long <- features %>%
            pivot_longer(cols = starts_with("feature"),  # Stack all feature columns
                         names_to = "Feature",
                         values_to = "Value") %>% 
            group_by(Donor, label, Feature)%>%
            summarise(
              across(where(is.numeric), \(x) mean(x, na.rm = FALSE)),  # Updated function syntax
              .groups = "drop"  # Drop grouping after summarization
            )%>%
            mutate(
              Feature = gsub("feature_", "", Feature),  # Remove "feature_" prefix
              Feature = as.numeric(Feature)  # Convert to numeric for ordering
            ) %>%
            arrange(Feature) 
          
          # df_long$Feature<-as.numeric(gsub("^feature_", "",df_long$Feature))
          
          for (element in unique(df_long$label)){
            
            list_of_substructures[[element]]<-df_long[df_long[["label"]]==element,]
            list_of_substructures[[element]]$label<-NULL
            list_of_substructures[[element]]<-list_of_substructures[[element]]%>%
                  pivot_wider(names_from = Donor, values_from = Value) %>%
                  as.data.frame()
            rownames(list_of_substructures[[element]])<-list_of_substructures[[element]]$Feature
            list_of_substructures[[element]]$Feature<-NULL
            
            #Normalize image features
            list_of_substructures[[element]]<- scale(list_of_substructures[[element]])
            list_of_substructures[[element]]<-as.matrix(list_of_substructures[[element]])
            log_message(paste0("----Matrix for ", element, " done----\n"))
  
            
          }
          
          common_columns <- Reduce(intersect, lapply(list_of_substructures, colnames))
          
          # Step 2: Subset and reorder columns in each dataframe
          list_of_substructures <- lapply(list_of_substructures, function(df) df[, common_columns, drop = FALSE])
          
          log_message("----Reading gene expression data----\n")
  
          ###Read and normalize gene expression data
          counts_tot <- readRDS(paste0(input,"/Laura/00.Data/v10/", tissue, "/counts.rds"))
          tpm_tot <- readRDS(paste0(input,"/Laura/00.Data/v10/", tissue, "/tpm.rds"))
          metadata_tot <- readRDS(paste0(input, "/Laura/00.Data/v10/", tissue, "/metadata.rds"))
          
          #For using continuous ancestry as a covariate (otherwise we use it as a categorical variable)
          metadata_tot$Ancestry<- NULL
          ancestry_file<- read.table(paste0(input,"/Laura/00.Data/admixture_inferred_ancestry.txt"))
          ances<-ancestry_file[,c(1,3)]
          colnames(ances)<- c("Donor", "Ancestry")
          metadata_tot<-merge(metadata_tot, ances, by= "Donor")
          metadata<- metadata_tot[metadata_tot$Donor %in% filter$Subject.ID,]
          counts <- counts_tot[, metadata$Sample]
          tpm<-tpm_tot[, metadata$Sample]
          
    
          # 1.2.1 TPM>=0.2 in at least 20% of the tissue samples
          
          exprs_genes.tpm <- rownames(tpm)[apply(tpm, 1, function(x) sum(x>=1) ) >= 0.2*ncol(tpm)  ]  
          
          # 1.2.2 Count >=10 in at least 20% of the tissue samples
          exprs_genes.counts <- rownames(counts)[ apply(counts, 1, function(x) sum(x>=10) ) >= 0.2*ncol(counts)  ]  
          
          # 1.2.3. Intersect gene lists
          exprs_genes <- intersect(exprs_genes.tpm,
                                   exprs_genes.counts) 
          
          # Exclude chrY genes in female-only tissues
          gene_annotation  <- read.csv(paste0(input, "/Laura/00.Data/gencode.v39.annotation.bed"))
          Y_genes <- gene_annotation[gene_annotation$chr=="chrY",]$ensembl.id
          exprs_genes <- exprs_genes[!exprs_genes %in% Y_genes]
          
          dge <- DGEList(counts[exprs_genes,])
          dge <- calcNormFactors(dge)
          v <- voom(dge, design = NULL, normalize="quantile", save.plot=F, plot =F ) # samples are treated as replicates
          
          
          if(option == options[1]){
              donors_to_add<-as.vector(colnames(list_of_substructures[[element]])[!colnames(list_of_substructures[[element]]) %in% colnames(gene_expression)])
          
              for (col in donors_to_add) {
                gene_expression[[col]] <- NA
              }
              
              rownames(gene_expression)<- gsub("\\.\\d+$", "",rownames(gene_expression))
              gene_data<-gene_expression
              gene_data <- gene_data[, colnames(list_of_substructures[[1]])]
              
          }else{
            for (element in names(list_of_substructures)){
              list_of_substructures[[element]]<-list_of_substructures[[element]][,colnames(list_of_substructures[[element]]) %in% metadata$Donor]
            }
              
            #Remove batch effect
              tec<- colnames(metadata)[colnames(metadata) %in% c("HardyScale", "IschemicTime", "RIN", "Cohort", "NucIsoBatch", "ExonicRate", "PEER2")]
            # Create the design matrix
              design <- model.matrix(as.formula(paste("~", paste(tec, collapse= " + "))), data = metadata)
              
            # Fit the linear model
              fit <- lmFit(v, design)
              exprs_residuals <- resid(fit, v)
              
              if(g == genes[1]){
                  log_message("----Selecting highly variable genes----\n")
                  
                  # Calculate variance
                  gene_variances <- apply(exprs_residuals, 1, var)
                  
                  # Select top 5,000 most variable genes
                  selected_genes <- order(gene_variances, decreasing = TRUE)[1:5000]
              }else{
                log_message("----Selecting ANM genes----\n")
                selected_genes<-gene_annotation[gene_annotation$gene.name.x %in% gwas_anm,]
                selected_genes<-selected_genes[,"ensembl.id"]
                selected_genes<-selected_genes[selected_genes %in% rownames(exprs_residuals)]
              }
  
              filtered_gene_expression <- exprs_residuals[selected_genes, ]
              
              #gene_expression<- v$E
              gene_expression<-as.data.frame(filtered_gene_expression)
              
              colnames(gene_expression)<-gsub("^(\\w+-\\w+).*", "\\1", colnames(gene_expression))
          
          }
          
          gene_data<-gene_expression[, colnames(gene_expression) %in% colnames(list_of_substructures[[element]])]
          rownames(gene_data)<- gsub("\\.\\d+$","", rownames(gene_data))
  
          log_message("----Creating MOFA Object----\n")
          
          total_data<-list()
  
          for (element in names(list_of_substructures)){
            total_data[[element]]<-list_of_substructures[[element]]
          }
          total_data$gene_expression<-gene_data
          total_data$gene_expression<-as.matrix(total_data$gene_expression)
          
          #Create MOFA Object
          MOFAobject <- create_mofa(total_data)
          plot_data_overview(MOFAobject)
          
          data_opts <- get_default_data_options(MOFAobject)
          model_opts <- get_default_model_options(MOFAobject)
          if (tissue %in% c("FallopianTube", "CervixEndocervix", "CervixEctocervix")){
            model_opts$num_factors <- 4
  
          }else{
            model_opts$num_factors <- 10
            
          }
          train_opts <- get_default_training_options(MOFAobject)
          train_opts$convergence_mode <- "slow"
          train_opts$seed <- 42
          log_message("----Training MOFA----\n")
          
          MOFAobject <- prepare_mofa(MOFAobject,
                                     data_options = data_opts,
                                     model_options = model_opts,
                                     training_options = train_opts
          )
          log_message("----Saving the model----\n")
          
          MOFAobject_trained<-run_mofa(MOFAobject, outfile=paste0(input,"/Laura/MOFA/MOFA_", tissue, "_",option,"_",model_opts$num_factors,"_ensg.hdf5"), use_basilisk = TRUE)
          saveRDS(MOFAobject_trained,paste0(input, "/Laura/MOFA/MOFA_", tissue, "_", option,"_",model_opts$num_factors,"_ensg.rds"))
  
      }
}


