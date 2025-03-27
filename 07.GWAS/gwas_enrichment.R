####################################################################################

###Code for the heatmaps (data shown in Main Figure 4F) and Supplementary Figure XX)#

### date: 18.03.2025                                                                #

### script version v1.0                                                             #

### R version R4.1.2.                                                               #

### code written by: Natàlia Pujol Gualdo                                           #

#####################################################################################



# Load required libraries

library(dplyr)

library(tidyr)

library(pheatmap)

library(ggplot2)

library(tibble)  





# Define input path

input_path <- "/X/GTEx_v8/"



#Load candidate genes from Age at natural menopause coming from GWAS (PMID: 34349265) and WES studies (https://www.nature.com/articles/s41586-024-07931-x)



#For GWAS: Read and process the ANM table, which is the supplementary table found in PMID: 34349265



temp_gwas_anm <- read_excel(file.path(input_path, "Natalia/data/ANM_ruth_tables.xlsx"),
                            
                            sheet = "ST2", skip = 5) %>%
  
  separate_rows(Consensus, sep = " / ") %>%
  
  distinct(Consensus, .keep_all = TRUE) %>%
  
  mutate(Trait = "ANM", Genes = Consensus) 



gwas_anm <- dplyr::select(temp_gwas_anm, Trait, Genes)



#For WES; the 9 genes defined in the WES study (https://www.nature.com/articles/s41586-024-07931-x)

EWAS_anm <- data.frame(
  
  Trait = "ANM",
  
  Genes = c("BRCA2", "CHEK2", "ETAA1", "HELB", "HROB", "PALB2", "PNPLA8", "SAMHD1", "ZNF518A"),  
  
  stringsAsFactors = FALSE
  
)



# Combine and check if there are duplicates

ALL_ANM <- bind_rows(EWAS_anm, table_anm) %>%
  
  distinct(Genes, .keep_all = TRUE)



ALL_ANM<-dplyr::select(ALL_ANM,Trait,Genes)





###################### Loading DEA with age results adjusted by standard covariates and structures

# Define tissues and file paths

tissues <- c("Uterus", "Ovary", "Myometrium", "Vagina", "Breast")



# Load DEA results datasets for the 5 tissues that showed any DEG with age 

tissue_data <- list(
  
  Uterus = readRDS(paste0(input_path, "Laura/12.Tissues_substructures/DEA/DEA_v10/Uterus/all_cov_no_inter/Uterus_all_cov_no_inter_AGE_covariates_and_traits.results.rds")),
  
  Ovary = readRDS(paste0(input_path, "Laura/12.Tissues_substructures/DEA/DEA_v10/Ovary/all_cov_no_inter/Ovary_all_cov_no_inter_AGE_covariates_and_traits.results.rds")),
  
  Myometrium = readRDS(paste0(input_path, "Laura/12.Tissues_substructures/DEA/DEA_v10/Uterus/all_cov_no_inter/Myometrium/Uterus_all_cov_no_inter_AGE_covariates_and_traits.results.rds")),
  
  Vagina = readRDS(paste0(input_path, "Laura/12.Tissues_substructures/DEA/DEA_v10/Vagina/all_cov_no_inter/Vagina_all_cov_no_inter_AGE_covariates_and_traits.results.rds")),
  
  Breast = readRDS(paste0(input_path, "Laura/12.Tissues_substructures/DEA/DEA_v10/BreastMammaryTissue/all_cov_no_inter/BreastMammaryTissue_all_cov_no_inter_AGE_covariates_and_traits.results.rds")))



# Extract 'Age' subset and filter differentially expressed genes (p-adj < 0.05)

degs_list <- lapply(tissue_data, function(x) {
  
  age_subset <- x$Age
  
  degs <- age_subset %>%
    
    filter(adj.P.Val < 0.05) %>%
    
    select(gene.name.x) %>%
    
    pull()  # Extracts the gene names as a vector
  
  return(degs)
  
})



# Assign DEGs dynamically as variables

total_degs <- setNames(degs_list, paste0("degs_", names(tissue_data)))

list2env(total_degs, envir = .GlobalEnv)  # Create individual variables



# Print DEGs summary

cat("\nTotal DEGs count in each tissue:\n")

degs_counts <- sapply(total_degs, length)

print(degs_counts)



# Compute intersections of DEGs with GWAS genes

ANM_intersections <- setNames(lapply(names(total_degs), function(tissue) {
  
  intersect(total_degs[[tissue]], ALL_ANM$Genes)
  
}), names(total_degs))



# Combine all unique genes across tissues,t

ANM_combined_genes <- unique(unlist(ANM_intersections))

cat("\nTotal unique intersected genes across tissues:", length(ANM_combined_genes), "\n")



# Extract and combine logFC values for all intersected genes across tissues

combined_data <- bind_rows(lapply(names(tissue_data), function(tissue) {
  
  age_data <- tissue_data[[tissue]]$Age  # Extract age-specific data
  
  
  
  if (!"gene.name.x" %in% colnames(age_data)) 
    
    stop(paste("Dataset for", tissue, "does not contain a 'gene.name.x' column. Check column names."))
  
  
  
  age_data %>%
    
    filter(gene.name.x %in% ANM_combined_genes) %>%
    
    select(gene.name.x, logFC) %>%
    
    mutate(Tissue = tissue)
  
}))





# Reshape into matrix format with Tissue as row names

logFC_matrix <- combined_data %>%
  
  pivot_wider(names_from = gene.name.x, values_from = logFC, values_fill = NA) %>%
  
  column_to_rownames(var = "Tissue")  



##############################################

### Large heatmap which goes to the supp Fig #

##############################################





# Remove genes that come from WES evidence from the main heatmap

main_logFC_matrix <- logFC_matrix[,!colnames(logFC_matrix) %in% c("SAMHD1", "HROB", "PALB2", "ETAA1") ]

main_logFC_matrix<-t(main_logFC_matrix)



# Create an annotation matrix filled with NA to highlight which are DEGs

annot_matrix <- matrix(NA, nrow = nrow(main_logFC_matrix), ncol = ncol(main_logFC_matrix))



# Set row and column names to match main_logFC_matrix

rownames(annot_matrix) <- rownames(main_logFC_matrix)

colnames(annot_matrix) <- colnames(main_logFC_matrix)



# Add stars where genes overlap with degs_Breast, degs_Ovary, degs_Uterus, degs_Myometrium

annot_matrix[rownames(annot_matrix) %in% degs_Breast, "Breast"] <- "*"

annot_matrix[rownames(annot_matrix) %in% degs_Ovary, "Ovary"] <- "*"

annot_matrix[rownames(annot_matrix) %in% degs_Uterus, "Uterus"] <- "*"

annot_matrix[rownames(annot_matrix) %in% degs_Myometrium, "Myometrium"] <- "*"



# Replace NA values with empty strings (to avoid showing NA in the heatmap)

annot_matrix[is.na(annot_matrix)] <- ""





# Define high (max) and low (min) values for scaling

min_val <- min(main_logFC_matrix, na.rm = TRUE)

max_val <- max(main_logFC_matrix, na.rm = TRUE)



# Define color scale with fixed low/high values

heatmap_colors <- colorRampPalette(c("#63B8FF", "white", "#8B0000"))(100)



# Define breaks to ensure min -> blue, 0 -> white, max -> red

breaks <- seq(min_val, max_val, length.out = 101)



# Update the heatmap call

heatmap_plot1 <- pheatmap(t(main_logFC_matrix), 
                          
                          cluster_rows = FALSE,  # Flip clustering accordingly
                          
                          cluster_cols = TRUE, # Ensure correct orientation
                          
                          scale = "none",  
                          
                          color = heatmap_colors,
                          
                          breaks = breaks,  # Ensures white is at 0
                          
                          na_col = "grey70",  
                          
                          fontsize = 12, 
                          
                          fontsize_number = 15,  
                          
                          display_numbers = t(annot_matrix),  
                          
                          number_color = "black",  
                          
                          border_color = "grey50",
                          
                          main = "",
                          
                          angle_col = 45,
                          
                          cellwidth = 15, 
                          
                          cellheight = 15,
                          
                          treeheight_col = 0,
                          
                          treeheight_row = 0)



ggsave(filename = "/home/npujol/Desktop/FINAL_SUPP_big_heatmap_ANM_65GENES.pdf", 
       
       plot = heatmap_plot1,  
       
       width = 17, height = 5, dpi = 300)



##############################################

### Small heatmap which goes to the Fig.4F ###

##############################################



# Select only the four genes correctly 

small_logFC_matrix <- logFC_matrix[, c("SAMHD1", "HROB", "PALB2", "ETAA1"), drop = FALSE]



# Ensure small_logFC_matrix is not empty

if (nrow(small_logFC_matrix) == 0) stop("No matching genes found in logFC_matrix!")



# Create annotation matrix with the same dimensions

annot_matrix <- matrix(NA, nrow = nrow(small_logFC_matrix), ncol = ncol(small_logFC_matrix),
                       
                       dimnames = list(rownames(small_logFC_matrix), colnames(small_logFC_matrix)))



# Add stars to specific positions (Check if row names exist before assigning)

if ("Myometrium" %in% rownames(annot_matrix)) annot_matrix["Myometrium", "SAMHD1"] <- "*"

if ("Uterus" %in% rownames(annot_matrix)) annot_matrix["Uterus", "SAMHD1"] <- "*"

if ("Ovary" %in% rownames(annot_matrix)) annot_matrix["Ovary", "HROB"] <- "*"

if ("Ovary" %in% rownames(annot_matrix)) annot_matrix["Ovary", "PALB2"] <- "*"

if ("Ovary" %in% rownames(annot_matrix)) annot_matrix["Ovary", "ETAA1"] <- "*"



# Convert NA values in annot_matrix to blank spaces

annot_matrix[is.na(annot_matrix)] <- ""



# Define desired row order

desired_order <- c("Ovary", "Myometrium", "Uterus")



# Ensure rownames exist before reordering

existing_rows <- intersect(desired_order, rownames(small_logFC_matrix))

small_logFC_matrix <- small_logFC_matrix[existing_rows, , drop = FALSE]

annot_matrix <- annot_matrix[existing_rows, , drop = FALSE]



# Define high (max) and low (min) values for scaling

min_val <- min(small_logFC_matrix, na.rm = TRUE)

max_val <- max(small_logFC_matrix, na.rm = TRUE)



# Define color scale with fixed low/high values

heatmap_colors <- colorRampPalette(c("#63B8FF", "white", "#8B0000"))(100)



# Define breaks to ensure min -> blue, 0 -> white, max -> red

breaks <- seq(min_val, max_val, length.out = 101)



# Generate heatmap with fixed high/low colors and correct star placement

heatmap_plot2 <- pheatmap(small_logFC_matrix, 
                          
                          cluster_rows = FALSE,  
                          
                          cluster_cols = TRUE, 
                          
                          scale = "none",  
                          
                          color = heatmap_colors,
                          
                          breaks = breaks,  
                          
                          na_col = "grey70",  
                          
                          fontsize = 12, 
                          
                          fontsize_number = 20,  # Large stars
                          
                          display_numbers = annot_matrix,  # Overlay stars
                          
                          number_color = "black",  # Stars in black
                          
                          border_color = "grey50",
                          
                          angle_col = 45,
                          
                          cellwidth = 25,  
                          
                          cellheight = 25,  
                          
                          treeheight_col = 0,
                          
                          treeheight_row = 0)



heatmap_plot2





ggsave(filename = "/X/FINAL_small_heatmap_WES_ANM.pdf", 
       
       plot = heatmap_plot2,  
       
       width = 5, height = 5, dpi = 300)

