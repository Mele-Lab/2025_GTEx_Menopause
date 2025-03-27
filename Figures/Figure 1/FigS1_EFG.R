##Fig S1 EFG
heatmaps<-list()
for (tissue in c("Uterus", "Ovary", "Vagina")){
  tissue<- "Uterus"
  
  s<-"female"
  buckets.size<- 20
  interval.size<-0
  r<-readRDS(paste0("Desktop/TFM/bsc83671/GTEx_v8/Laura/DEswan/All_tissues/", tissue,"/DEswan_", tissue, "_", s,"_", buckets.size,"_", interval.size, "_results.rds"))
  
  window.center = seq(21,70,1)
  # window.center = seq(0,1,0.05)
  
  pvals<-c("P.Value", "adj.P.Val")
  #We choose depending on the results
  pval<-pvals[2]
  
  ####1. Extract the significant genes in at least one comparison
  
  library(parallel)
  library(data.table)
  
  # Detect the number of cores available
  num_cores <- detectCores() - 1  # Use one less than max cores
  
  ### 1. Extract significant genes in parallel
  extract_significant_genes <- function(sublist) {
    if (!is.null(sublist) && pval %in% colnames(sublist) && "gene.name.x" %in% colnames(sublist)) {
      return(sublist[sublist[[pval]] < 0.05, "gene.name.x", drop = FALSE])
    }
    return(NULL)
  }
  # Use parallel execution instead of lapply
  significant_genes_list <- mclapply(r, extract_significant_genes, mc.cores = num_cores)
  significant_genes_list <- significant_genes_list[!sapply(significant_genes_list, is.null)]
  final_significant_genes <- unique(do.call(rbind, significant_genes_list))

    ### 2. Optimize P-value matrix generation
  all_windows <- seq_along(r)
  gene_list <- unique(unlist(lapply(r, rownames)))
  
  # Create an empty data.table instead of a matrix
  p_value_matrix <- data.table(matrix(NA_real_, nrow = length(gene_list), ncol = length(all_windows)))
  setnames(p_value_matrix, as.character(all_windows))
  setattr(p_value_matrix, "rownames", gene_list)
  
  # Fill p-value matrix in parallel
  p_value_matrix <- do.call(cbind, mclapply(all_windows, function(win) {
    if (!is.null(r[[win]]) && is.data.frame(r[[win]])) {
      p_values <- r[[win]][[pval]]
      names(p_values) <- rownames(r[[win]])
      return(p_values)
    }
    return(rep(NA_real_, length(gene_list)))  # Return NA if data is missing
  }, mc.cores = num_cores))
  
  # Convert to data frame
  p_value_df <- as.data.frame(p_value_matrix)
  
  ### 3. Optimize Counting Significant Genes
  count_significant <- function(p_df, thresholds = c(0.05, 0.01, 0.001)) {
    return(sapply(thresholds, function(thresh) colSums(p_df < thresh, na.rm = TRUE)))
  }
  
  res.signif <- count_significant(p_value_df)
  rownames(res.signif) <- window.center
  colnames(res.signif) <- c(0.05, 0.01, 0.001)
  
  ### 5. Parallelize Sample Count Calculation
  sample_counts <- unlist(mclapply(window.center, function(center) {
    sum(metadata$Age >= (center - buckets.size) & metadata$Age < (center + buckets.size))
  }, mc.cores = num_cores))
  
  df_counts <- data.table(Window = window.center, Samples = sample_counts)
  df_counts$Genes <- res.signif[, "0.05"]
  
  ###4. Heatmap
  ##4.1. With z-score
  metadata_tot <- readRDS(paste0(input, "Laura/00.Data/v10/", tissue, "/metadata.rds"))
  counts_tot <- readRDS(paste0(input, "Laura/00.Data/v10/", tissue, "/counts.rds"))
  tpm_tot <- readRDS(paste0(input, "Laura/00.Data/v10/", tissue, "/tpm.rds"))
  metadata_tot$Ancestry<- NULL
  ancestry_file<- read.table(paste0(input, "/Laura/00.Data/admixture_inferred_ancestry.txt"))
  ances<-ancestry_file[,c(1,3)]
  colnames(ances)<- c("Donor", "Ancestry")
  metadata_tot<-merge(metadata_tot, ances, by= "Donor")
  
  filter<-readRDS(paste0(input, "/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
  metadata<- metadata_tot[(metadata_tot$Donor %in% filter$Subject.ID),]   
  
  z.score <- function(ft.matrix){
    t(scale(t(ft.matrix), center=TRUE, scale=TRUE))
  }
  
  counts_tot <- readRDS(paste0("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/00.Data/v10/", tissue, "/counts.rds"))
  tpm_tot <- readRDS(paste0("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/00.Data/v10/", tissue, "/tpm.rds"))
  counts <- counts_tot[, metadata$Sample]
  tpm<-tpm_tot[, metadata$Sample]
  
  ######RESIDUALS
  # library(dplyr)
  # residuals<-readRDS(paste0("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/12.Tissues_substructures/DEA/DEA_v10/",tissue,"/no_subs_no_age/", tissue,"_no_subs_no_age_residuals.results.rds"))
  # residuals <- as.data.frame(residuals)
  # residuals[1:5, 1:5]
  # residuals <- residuals %>%
  #   mutate(`GTEX-R55G-1626` = rowMeans(dplyr::select(., `GTEX-R55G-1626`, `GTEX-R55G-1626.1`), na.rm = TRUE)) %>%
  #   select(-`GTEX-R55G-1626.1`)
  # 
  # residuals<- as.data.frame(t(residuals))
  # residuals$Age<- metadata$Age
  # final_significant_genes<-final_significant_genes[rownames(final_significant_genes) %in% colnames(residuals), , drop = FALSE]
  # 
  # selected_cols <- intersect(c(rownames(final_significant_genes), "Age"), colnames(residuals))
  # genes_selected <- residuals[, selected_cols, drop = FALSE]
  # colnames(genes_selected)<-c(final_significant_genes$gene.name.x, "Age")
  # rownames(genes_selected)<-as.factor(gsub("^(\\w+-\\w+).*", "\\1", rownames(genes_selected)))
  # 
  
  ####TPMS
  exprs_genes.tpm <- rownames(tpm)[apply(tpm, 1, function(x) sum(x>=1) ) >= 0.2*ncol(tpm)  ]  
  
  # 1.2.2 Count >=6 in at least 20% of the tissue samples
  exprs_genes.counts <- rownames(counts)[ apply(counts, 1, function(x) sum(x>=10) ) >= 0.2*ncol(counts)  ]  
  
  # 1.2.3. Intersect gene lists
  exprs_genes <- intersect(exprs_genes.tpm,
                           exprs_genes.counts) 
  tpm<-tpm[exprs_genes,]
  
  tpm_matrix<- as.matrix(log2(tpm +1))
  tpm_matrix_zscored<-as.data.frame(t(z.score(tpm_matrix)))
  tpm_matrix_zscored$Age<- metadata$Age
  genes_selected<-tpm_matrix_zscored[,c(rownames(final_significant_genes), "Age")]
  colnames(genes_selected)<-c(final_significant_genes$gene.name.x, "Age")
  rownames(genes_selected)<-as.factor(gsub("^(\\w+-\\w+).*", "\\1", rownames(genes_selected)))
  genes_selected$Sample<-as.factor(rownames(genes_selected))
  
  z_scores <- genes_selected[, -which(colnames(genes_selected) %in% c("Sample"))]
  z_scores_avg <- aggregate(. ~ Age, data = z_scores, FUN = mean)
  z_scores_avg <- z_scores_avg[order(z_scores_avg$Age), ]
  
  # Convert to matrix format for the heatmap
  rownames(z_scores_avg) <- z_scores_avg$Age
  z_scores_matrix <- t(as.matrix(z_scores_avg[, -1]))  # Remove Age column
  
  # Define breaks and color palette
  
  selected_ages <- c(20,30, 40,50, 60, 70)
  
  # Replace all other column names with empty strings
  colnames(z_scores_matrix) <- ifelse(as.numeric(colnames(z_scores_matrix)) %in% selected_ages, 
                                      colnames(z_scores_matrix), "")

  min<--max(z_scores_matrix)

  max<- max(z_scores_matrix)

  
  pairs.breaks <- seq(min, max, by = 0.01)
  
  # Generate the color palette
  library(gplots)
  
  mycol <- gplots::colorpanel(n = length(pairs.breaks) - 1, low = "cyan", mid = "black", high = "yellow")
  # Plot the heatmap
  dev.off()
  pdf(paste0("Desktop/TFM/bsc83671/GTEx_v8/Laura/Figure_plots/FigS1_EFG/", tissue, "_ht.pdf"), width = 3, height = 2)

  heatmap.2(
    z_scores_matrix,
    cexRow = 0.7,                    # Adjust font sizes
    cexCol = 2,                      # Adjust font sizes for columns
    trace = "none",                  # No trace lines
    dendrogram = "row",              # Cluster rows (genes)
    breaks = pairs.breaks,           # Define the breaks for color scale
    col = mycol,                     # Use the color palette
    Rowv = TRUE,                     # Cluster rows
    Colv = FALSE,                    # Do not cluster columns (ages stay in order)
    lhei = c(0.2, 10),               # Adjust layout heights
    lwid = c(0.2, 3),                # Adjust layout widths
    key = TRUE,                      # Display legend
    key.title = "Z-score",           # Title of the legend
    key.xlab = "Z-score Range",      # Optional: Label for the x-axis of the legend
    labRow = NA,                     # Removes row names
    notecex = 0,                     # Ensure no annotations appear
    sepcolor = NA,                   # Disable border color
    sepwidth = c(0, 0),              # No separators
    srtCol = 0
    # Rotate column labels horizontally
  )
  dev.off()

  svg(paste0("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/Figure_plots/FigS1_EFG/", tissue, "_ht.svg"), width = 3, height = 2, pointsize = 12)
  heatmap.2(
    z_scores_matrix,
    cexRow = 0.7,                    # Adjust font sizes
    cexCol = 2,                      # Adjust font sizes for columns
    trace = "none",                  # No trace lines
    dendrogram = "row",              # Cluster rows (genes)
    breaks = pairs.breaks,           # Define the breaks for color scale
    col = mycol,                     # Use the color palette
    Rowv = TRUE,                     # Cluster rows
    Colv = FALSE,                    # Do not cluster columns (ages stay in order)
    lhei = c(0.2, 10),               # Adjust layout heights
    lwid = c(0.2, 3),                # Adjust layout widths
    key = TRUE,                      # Display legend
    key.title = "Z-score",           # Title of the legend
    key.xlab = "Z-score Range",      # Optional: Label for the x-axis of the legend
    labRow = NA,                     # Removes row names
    notecex = 0,                     # Ensure no annotations appear
    sepcolor = NA,                   # Disable border color
    sepwidth = c(0, 0),              # No separators
    srtCol = 0
    # Rotate column labels horizontally
  )
  dev.off()
  pdf(paste0("Desktop/Figures/FigS1_EFG/", tissue, "_ht.pdf"), width = 30, height = 30)
  heatmap.2(
    z_scores_matrix,
    cexRow = 0.7,                    # Adjust font sizes
    cexCol = 2,                      # Adjust font sizes for columns
    trace = "none",                  # No trace lines
    dendrogram = "row",              # Cluster rows (genes)
    breaks = pairs.breaks,           # Define the breaks for color scale
    col = mycol,                     # Use the color palette
    Rowv = TRUE,                     # Cluster rows
    Colv = FALSE,                    # Do not cluster columns (ages stay in order)
    lhei = c(0.2, 10),               # Adjust layout heights
    lwid = c(0.2, 3),                # Adjust layout widths
    key = TRUE,                      # Display legend
    key.title = "Z-score",           # Title of the legend
    key.xlab = "Z-score Range",      # Optional: Label for the x-axis of the legend
    labRow = NA,                     # Removes row names
    notecex = 0,                     # Ensure no annotations appear
    sepcolor = NA,                   # Disable border color
    sepwidth = c(0, 0),              # No separators
    srtCol = 0
    # Rotate column labels horizontally
  )
  dev.off()
  svg(paste0("Desktop/Figures/FigS1_EFG/", tissue, "_ht.svg"), width = 25, height = 25, pointsize = 12)
  heatmap.2(
    z_scores_matrix,
    cexRow = 0.7,                    # Adjust font sizes
    cexCol = 2,                      # Adjust font sizes for columns
    trace = "none",                  # No trace lines
    dendrogram = "row",              # Cluster rows (genes)
    breaks = pairs.breaks,           # Define the breaks for color scale
    col = mycol,                     # Use the color palette
    Rowv = TRUE,                     # Cluster rows
    Colv = FALSE,                    # Do not cluster columns (ages stay in order)
    lhei = c(0.2, 10),               # Adjust layout heights
    lwid = c(0.2, 3),                # Adjust layout widths
    key = TRUE,                      # Display legend
    key.title = "Z-score",           # Title of the legend
    key.xlab = "Z-score Range",      # Optional: Label for the x-axis of the legend
    labRow = NA,                     # Removes row names
    notecex = 0,                     # Ensure no annotations appear
    sepcolor = NA,                   # Disable border color
    sepwidth = c(0, 0),              # No separators
    srtCol = 0
    # Rotate column labels horizontally
  )
  dev.off()

  svg(paste0("Desktop/TFM/bsc83671/GTEx_v8/Laura/Figure_plots/FigS1_EFG/", tissue, "_legend.svg"), width = 3, height = 2)
  mycol <- colorpanel(n = 100, low = "cyan", mid = "black", high = "yellow")
  breaks <- seq(min, max, length.out = 100)
  labels <- c(min, 0, max)  # Only these three labels
  plot(1, type = "n", xlab = "", ylab = "", xlim = c(-6, 6), ylim = c(0, 1), axes = FALSE)
  image(x = breaks, y = 0, z = matrix(breaks, ncol = 1), col = mycol, axes = FALSE, xlab = "", ylab = "")
  box()
  axis(side = 1, at = labels, labels = labels, las = 1)
  title(main = "z-score", line = 1)
  
  dev.off()

  svg(paste0("Desktop/Figures/FigS1_EFG/", tissue, "_legend.svg"), width = 3, height = 2)
  mycol <- colorpanel(n = 100, low = "cyan", mid = "black", high = "yellow")
  breaks <- seq(min, max, length.out = 100)
  labels <- c(min, 0, max)  # Only these three labels
  plot(1, type = "n", xlab = "", ylab = "", xlim = c(-6, 6), ylim = c(0, 1), axes = FALSE)
  image(x = breaks, y = 0, z = matrix(breaks, ncol = 1), col = mycol, axes = FALSE, xlab = "", ylab = "")
  box()
  axis(side = 1, at = labels, labels = labels, las = 1)
  title(main = "z-score", line = 1)
  
  dev.off()
  
}



