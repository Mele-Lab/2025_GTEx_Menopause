###############################################Supp table 1 - CNNs
library(writexl)
library(openxlsx)

combined_wb <- createWorkbook()
save_with_rownames <- function(df) {
  df <- cbind(Feature = rownames(df), df)  # Add row names as a column called "Feature"
  rownames(df) <- NULL  # Remove original row names
  return(df)
}

add_combined_sheet <- function(wb, sheet_name, file_paths, file_titles = NULL) {
  addWorksheet(wb, sheet_name)
  
  if (is.null(file_titles)) {
    file_titles <- sub("\\.rds$", "", basename(file_paths))# Default: Use filenames as titles
  }
  
  start_col <- 1
  for (i in seq_along(file_paths)) {
    table <- save_with_rownames(readRDS(file_paths[i]))
    
    # Write title (e.g., "Myometrium", "Duct", etc.)
    writeData(wb, sheet_name, file_titles[i], startCol = start_col, startRow = 2)
    
    # Write the actual table
    writeData(wb, sheet_name, table, startCol = start_col, startRow = 3)
    
    # Move to the next block (leave 2 columns gap for visual clarity)
    start_col <- start_col + ncol(table) + 2
  }
}


# Breast files
breast_files <- c(
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/BreastMammaryTissue/adipocyte.rds",
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/BreastMammaryTissue/duct.rds",
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/BreastMammaryTissue/lobule.rds",
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/BreastMammaryTissue/nerve.rds",
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/BreastMammaryTissue/stroma.rds",
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/BreastMammaryTissue/gynecomastoid_hyperplasia.rds"
)

# Uterus files
uterus_files <- c(
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/Uterus/myometrium.rds",
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/Uterus/endometrium.rds",
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/Uterus/vessels.rds"
)
ovary_files <- c(
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/Ovary/corpora.rds",
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/Ovary/medulla.rds",
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/Ovary/cortex.rds",
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/Ovary/vessels.rds"
)

vagina_files <- c(
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/Vagina/epithelium.rds",
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/Vagina/lamina_propria.rds",
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/Vagina/stroma.rds",
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/Vagina/vessels.rds"
)

endocervix_files <- c(
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/Endocervix/glandular_epithelium.rds",
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/Endocervix/stroma.rds",
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/Endocervix/vessels.rds"
)

ectocervix_files <- c(
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/Ectocervix/epithelium.rds",
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/Ectocervix/glands.rds",
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/Ectocervix/stroma.rds",
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/Ectocervix/vessels.rds"
)

fallopian_files <- c(
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/FallopianTube/epithelium.rds",
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/FallopianTube/lumen.rds",
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/FallopianTube/smooth_muscle.rds",
  "X/Allal/Differential_expression/variance_partition/results/mean_tiles/FallopianTube/stroma.rds"
)
########MEDIAN
# Breast files
breast_files_med <- c(
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/BreastMammaryTissue/adipocyte.rds",
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/BreastMammaryTissue/duct.rds",
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/BreastMammaryTissue/lobule.rds",
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/BreastMammaryTissue/nerve.rds",
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/BreastMammaryTissue/stroma.rds",
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/BreastMammaryTissue/gynecomastoid_hyperplasia.rds"
)

# Uterus files
uterus_files_med <- c(
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/Uterus/myometrium.rds",
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/Uterus/endometrium.rds",
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/Uterus/vessels.rds"
)
ovary_files_med <- c(
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/Ovary/corpora.rds",
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/Ovary/medulla.rds",
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/Ovary/cortex.rds",
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/Ovary/vessels.rds"
)

vagina_files_med <- c(
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/Vagina/epithelium.rds",
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/Vagina/lamina_propria.rds",
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/Vagina/stroma.rds",
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/Vagina/vessels.rds"
)

endocervix_files_med <- c(
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/Endocervix/glandular_epithelium.rds",
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/Endocervix/stroma.rds",
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/Endocervix/vessels.rds"
)

ectocervix_files_med <- c(
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/Ectocervix/epithelium.rds",
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/Ectocervix/glands.rds",
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/Ectocervix/stroma.rds",
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/Ectocervix/vessels.rds"
)

fallopian_files_med <- c(
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/FallopianTube/epithelium.rds",
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/FallopianTube/lumen.rds",
  "X/Allal/Differential_expression/variance_partition/results/median_tiles/FallopianTube/smooth_muscle.rds",
  "X/Allal/Differential_expression/variance_partition/results/median_tiles//FallopianTube/stroma.rds"
)

# Add both sheets to the same workbook
add_combined_sheet(combined_wb, "Uterus-mean", uterus_files)
add_combined_sheet(combined_wb, "Ovary-mean", ovary_files)
add_combined_sheet(combined_wb, "Vagina-mean", vagina_files)
add_combined_sheet(combined_wb, "Breast-mean", breast_files)
add_combined_sheet(combined_wb, "Endocervix-mean", endocervix_files)
add_combined_sheet(combined_wb, "Ectocervix-mean", ectocervix_files)
add_combined_sheet(combined_wb, "FallopianTube-mean", fallopian_files)
add_combined_sheet(combined_wb, "Uterus-median", uterus_files_med)
add_combined_sheet(combined_wb, "Ovary-median", ovary_files_med)
add_combined_sheet(combined_wb, "Vagina-median", vagina_files_med)
add_combined_sheet(combined_wb, "Breast-median", breast_files_med)
add_combined_sheet(combined_wb, "Endocervix-median", endocervix_files_med)
add_combined_sheet(combined_wb, "Ectocervix-median", ectocervix_files_med)
add_combined_sheet(combined_wb, "FallopianTube-median", fallopian_files_med)


# Save the final workbook with both sheets
saveWorkbook(combined_wb, "X/Laura/Supp_tables/Supplementary_table_1.xlsx", overwrite = TRUE)


###############################################Supp table 3 - DEGs

tissues<-c("Uterus", "Ovary", "Vagina", "BreastMammaryTissue", "CervixEndocervix", "CervixEctocervix", "FallopianTube")

combined_wb <- createWorkbook()
save_with_rownames <- function(df) {
  df <- cbind(Gene = rownames(df), df)  # Add row names as a column called "Feature"
  rownames(df) <- NULL  # Remove original row names
  return(df)
}

add_age_and_substructure_sheet <- function(wb, sheet_name, age_degs, sub_degs) {
  addWorksheet(wb, sheet_name)
  
  # Start with age_degs
  start_col <- 1
  writeData(wb, sheet_name, "Age_DEGs", startCol = start_col, startRow = 2)  # Title
  writeData(wb, sheet_name, save_with_rownames(age_degs), startCol = start_col, startRow = 3)
  
  # Move to the next column block
  start_col <- start_col + ncol(age_degs) + 3  # leave some gap
  
  # Add each substructure table
  for (sub in names(sub_degs)) {
    sub_table <- sub_degs[[sub]]
    writeData(wb, sheet_name, sub, startCol = start_col, startRow = 2)  # Title
    writeData(wb, sheet_name, save_with_rownames(sub_table), startCol = start_col, startRow = 3)
    start_col <- start_col + ncol(sub_table) + 3  # leave some gap between tables
  }
}

for (tissue in tissues){
  
  if(tissue =="BreastMammaryTissue"){
    substructures<- c("adipocyte", "lobule", "duct", "stroma", "nerve", "gynecomastoid_hyperplasia")
    
  }else if(tissue == "Ovary"){
    substructures<- c("cortex", "corpora", "medulla", "vessels")
    
  }else if (tissue =="Vagina"){
    substructures<- c("epithelium", "lamina_propria", "vessels", "stroma")
    
  }else if (tissue =="Uterus"){
    substructures<- c("myometrium", "endometrium", "vessels")
    
  }else if (tissue =="CervixEndocervix"){
    substructures<- c("vessels", "glandular_epithelium", "stroma")
    
    
  }else if (tissue =="CervixEctocervix"){
    substructures<- c("epithelium", "stroma", "vessels", "glands")
    
  }else if (tissue == "FallopianTube"){
    substructures<- c("vessels", "lumen", "smooth_muscle", "epithelium", "stroma") ##always control for adipose, since it is contamination
  }
  inpath<-paste0("~/X/Laura/12.Tissues_substructures/DEA/DEA_v10/", tissue, "/")
  age_f<-readRDS( paste0(inpath, "all_cov_no_inter/", tissue, "_all_cov_no_inter_AGE_covariates_and_traits.results.rds"))
  age_f<-age_f$Age
  if (!tissue %in% c("CervixEndocervix","CervixEctocervix", "FallopianTube")){
      sub_degs<- sub_f[sub_f$adj.P.Val<0.05,]
  }else{
    age_degs<- age_f[age_f$P.Value<0.01,]
  }
  sub_degs_t <- list()
  for (sub in substructures){
    sub_f<-readRDS( paste0(inpath,"one_cov_no_inter/", tissue, "_",sub, "_one_cov_no_inter_covariates_and_traits.results.rds"))
    sub_f<-sub_f[[sub]]
    if (!tissue %in% c("CervixEndocervix","CervixEctocervix", "FallopianTube")){
      sub_degs<- sub_f[sub_f$adj.P.Val<0.05,]
    }else{
      sub_degs<- sub_f[sub_f$P.Value<0.01,]
    }
    if(nrow(sub_degs)>0){
      sub_degs_t[[sub]]<-sub_degs
    }
  }
  if(length(sub_degs_t)>0){
    add_age_and_substructure_sheet(combined_wb, sheet_name = tissue, age_degs, sub_degs_t)
    
  }else{
    add_age_and_substructure_sheet(combined_wb, sheet_name = tissue, age_degs)
    
  }
}

saveWorkbook(combined_wb, "~/X/Laura/Supp_tables/Supplementary_table_DEGs.xlsx", overwrite = TRUE)

tissue<- "Uterus"
## myometrium-age degs
m<-readRDS("X/Laura/12.Tissues_substructures/DEA/DEA_v10/Uterus/all_cov_no_inter/Myometrium/Uterus_all_cov_no_inter_AGE_covariates_and_traits.results.rds")
all<-m$Age 
degs<-all[all$adj.P.Val<0.05,]
nrow(degs) 

write.xlsx(save_with_rownames(degs), "~/X/Laura/Supp_tables/myometrium.xlsx")
all<- readRDS( paste0(inpath, test, "/", tissue, "_", test, "_AGE_covariates_and_traits.results.rds"))


#############################################Supp table - CNNs
tissues<-c("Uterus", "Ovary", "Vagina", "BreastMammaryTissue")

file_counts <- data.frame(
  Tissue = character(),
  Folder = character(),
  Subfolder = character(),
  FileCount = integer(),
  stringsAsFactors = FALSE
)

# Base path
base_path <- "~/X/Laura/05.CNN/35yo_separation_FEMALE/"

# Folder-to-subfolder mapping
folder_mapping <- list(
  train = c("old", "young"),
  test = c("old", "young"),
  middle_age = "mid",
  validationset = "val",
  young_post = "young_post"
)

# Loop through each tissue
for (tissue in tissues) {
  
  # Construct path to tissue folder
  tissue_path <- file.path(base_path, tissue)
  
  # List all folders in the tissue directory
  all_folders <- list.files(tissue_path, full.names = TRUE)
  
  # Find folder that starts with "Tiles"
  tiles_folder <- all_folders[grepl("^Tiles", basename(all_folders))]
  
  if (length(tiles_folder) == 0) {
    warning(paste("No Tiles folder found for", tissue))
    next
  }
  
  # Loop through each main folder (train, test, etc.)
  for (main_folder in names(folder_mapping)) {
    subfolders <- folder_mapping[[main_folder]]
    
    # Make sure subfolders are always treated as a vector (even if length 1)
    if (!is.vector(subfolders)) subfolders <- c(subfolders)
    
    for (subfolder in subfolders) {
      # Full path to the subfolder
      subfolder_path <- file.path(tiles_folder, main_folder, subfolder)
      
      if (dir.exists(subfolder_path)) {
        # Count files
        file_count <- length(list.files(subfolder_path, full.names = TRUE))
        
        # Set folder name based on rules
        if (main_folder %in% c("train", "test")) {
          combined_folder <- paste0(main_folder, "-", subfolder)
        } else {
          combined_folder <- main_folder
        }
        
        # Add to dataframe
        file_counts <- rbind(
          file_counts,
          data.frame(
            Tissue = tissue,
            Folder = combined_folder,
            FileCount = file_count,
            stringsAsFactors = FALSE
          )
        )
      } else {
        warning(paste("Missing folder:", subfolder_path))
      }
    }
  }
}

# Preview results
print(file_counts)
# Create workbook
wb <- createWorkbook()
# Add the file counts dataframe as a sheet
addWorksheet(wb, "FileCounts")
writeData(wb, "FileCounts", file_counts,  startRow = 1)

#Add porbabilities
probs<-readRDS("~/X/Ole/paper/for_plotting/CNN_classification.rds")
addWorksheet(wb, "Probabilities")
# writeData(wb, "OtherSheet", other_df)
writeData(wb, "Probabilities", probs, startCol = 1, startRow = 1)

saveWorkbook(wb, "~/X/Laura/Supp_tables/Supplementary_table_CNNs.xlsx", overwrite = TRUE)
##Now add parameters
#Manually to parameters.xlsx

##Add accuracies
result_cnn<- readRDS("X/Laura/05.CNN/metrics_new_validation_filtered_perdonor_ACC.rds")
result_cnn$Organ<-c("Uterus", "Ovary", "Vagina", "Breast")
result_cnn_pertile<- readRDS("X/Laura/05.CNN/metrics_new_validation_filteredpertile_ACC.rds")
result_cnn_pertile$Organ<-c("Uterus", "Ovary", "Vagina", "Breast")
result_cnn$Metric<- "Per-donor"
result_cnn_pertile$Metric <- "Per-tile"
global<- rbind(result_cnn, result_cnn_pertile)#[result_cnn_pertile$Tissue %in% c("Uterus", "Ovary"),])
global$Tissue<-NULL

write.xlsx(global, "~/X/Laura/Supp_tables/accuracies.xlsx", rowNames = FALSE)


#####Supp- MOFA
library(MOFA2)
library(MOFAdata)
tissues<- c("Uterus", "Ovary", "Vagina", "BreastMammaryTissue", "FallopianTube", "CervixEctocervix", "CervixEndocervix")



wb <- createWorkbook()

for (tissue in tissues){
  
  if (tissue %in% c("FallopianTube", "CervixEctocervix", "CervixEndocervix")){
    MOFAobject_trained<-readRDS(paste0(input, "/Laura/MOFA/MOFA_", tissue, "_only_matching_donors_4_ensg_CORRECTED.rds"))
    
  }else{
    MOFAobject_trained<-readRDS(paste0(input, "/Laura/MOFA/MOFA_", tissue, "_only_matching_donors_10_ensg_CORRECTED.rds"))
    
  }
  ###Data conforming the figures
  
  metadata_tot <- readRDS(paste0(input, "/Laura/00.Data/v10/", tissue, "/metadata.rds"))
  filter<-readRDS(paste0(input, "/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
  metadata<- metadata_tot[metadata_tot$Donor %in% filter$Subject.ID,]
  metadata_MOFA<-metadata
  colnames(metadata_MOFA)[1]<-"sample"
  metadata_MOFA<- metadata_MOFA[metadata_MOFA$sample %in% colnames(MOFAobject_trained@data[[1]]$group1),]
  metadata_MOFA$Sample<-NULL
  ordered_metadata <- metadata_MOFA[order(metadata_MOFA$Age), ]
  samples_metadata(MOFAobject_trained) <- metadata_MOFA
  variance<-get_variance_explained(MOFAobject_trained)
  
  p<-correlate_factors_with_covariates(MOFAobject_trained, 
                                       covariates = c("Age"), 
                                       #covariates = colnames(metadata_MOFA_nonas), 
                                       plot="r",
                                       return_data = TRUE)
  p<-as.data.frame(p)
  
  p[2]<-p$Age
  p[1]<-rownames(p)
  colnames(p)<-c("Factor", "r")

  p_log <- tryCatch(
    correlate_factors_with_covariates(MOFAobject_trained, 
                                      covariates = c("Age"), 
                                      plot = "log_pval", 
                                      return_data = TRUE),
    error = function(e) {
      message("Error caught: ", e$message)
      # Ensure same structure, fill with 1s
      return(data.frame(Age = rep(1, nrow(p)), row.names = rownames(p)))
    }
  )
  p_log <- as.data.frame(p_log)
  
  # Compute adjusted p-values
  p$adj_p_val<-10^(-p_log$Age)
  p2<-plot_variance_explained(MOFAobject_trained, max_r2=max(variance$r2_per_factor$group1))
  # p3<-plot_variance_explained(MOFAobject_trained, plot_total = T)[[2]]
  variance_df<-p2$data
  addWorksheet(wb, paste0(tissue))
  
  writeData(wb, tissue, "Correlations", startCol = 1, startRow = 3)  # Title
  writeData(wb, tissue, p, startCol = 1, startRow = 4)
  writeData(wb, tissue, "Variance explained (%)", startCol = ncol(p)+2, startRow = 3)  # Title
  writeData(wb, tissue, variance_df, startCol = ncol(p)+2, startRow = 4)

  
  # factors <- get_factors(MOFAobject_trained, 
  #                        # factors = all, 
  #                        as.data.frame = TRUE
  # )
  # data <- get_data(MOFAobject_trained, 
  #                  # views = "myometrium", 
  #                  as.data.frame = TRUE
  # )
  # 
  # addWorksheet(wb, tissue)
  # 
  # # Write title (e.g., "Myometrium", "Duct", etc.)
  # writeData(wb, tissue, "Factors", startCol = 1, startRow = 3)
  # writeData(wb, tissue, factors, startCol = 1, startRow = 4)
  # writeData(wb, tissue, "Features", startCol = ncol(factors) +2, startRow = 2)
  # writeData(wb, tissue, data, startCol = ncol(factors) +2, startRow = 3)

}
saveWorkbook(wb, "~/X/Laura/Supp_tables/Supplementary_table_MOFA.xlsx", overwrite = TRUE)


data[data$view =="gene_expression",]

###########################Supp enrichments
#Enrichments for uterus, myometrium, ovary
wb <- createWorkbook()
#Enrichments MOFA

uterus_4<- readRDS("X/Laura/12.Tissues_substructures/Enrichment/Enrichments_GSEA/gsea_Uterus_goterms_MOFA4.rds")
uterus_4 <- data.frame(GO = rownames(uterus_4), uterus_4, row.names = NULL)
uterus_6<- readRDS("X/Laura/12.Tissues_substructures/Enrichment/Enrichments_GSEA/gsea_Uterus_goterms_MOFA6.rds")
uterus_6 <- data.frame(GO = rownames(uterus_6), uterus_6, row.names = NULL)

vagina<- readRDS("X/Laura/12.Tissues_substructures/Enrichment/Enrichments_GSEA/gsea_Vagina_goterms_MOFA.rds")
vagina <- data.frame(GO = rownames(vagina), vagina, row.names = NULL)

ovary<- readRDS("X/Laura/12.Tissues_substructures/Enrichment/Enrichments_GSEA/gsea_Ovary_goterms_MOFA.rds")
ovary <- data.frame(GO = rownames(ovary), ovary, row.names = NULL)

addWorksheet(wb, "Uterus - MOFA(F4)")
writeData(wb, "Uterus - MOFA(F4)", uterus_4, startCol = 1, startRow = 3)

addWorksheet(wb, "Uterus - MOFA(F6)")
writeData(wb, "Uterus - MOFA(F6)", uterus_6, startCol = 1, startRow = 3)

addWorksheet(wb, "Vagina - MOFA")
writeData(wb, "Vagina - MOFA", vagina, startCol = 1, startRow = 3)

addWorksheet(wb, "Ovary - MOFA")
writeData(wb, "Ovary - MOFA", ovary, startCol = 1, startRow = 3)

saveWorkbook(wb, "~/X/Laura/Supp_tables/Table_S5.xlsx", overwrite = TRUE)


##Github tables
ge<-readRDS("X/Laura/Supp_tables/GE_classification.rds")
ge$age<-NULL
cnn<-readRDS("X/Laura/Supp_tables/CNN_classification.rds")
cnn$age<-NULL
write.csv(ge, "X/Laura/Supp_tables/gene_expression_classification.csv", row.names = FALSE)
write.csv(cnn, "X/Laura/Supp_tables/cnn_classification.csv", row.names=FALSE)




