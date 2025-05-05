## Functions for clustering and plotting
library(tidyverse)
library(splines)
library(cluster)
library(dendextend)
library(NbClust)
library(patchwork)

gene_trajectory_analysis <- function(
    counts_matrix,  # Genes x Donors matrix
    metadata, 
    age_column = 'Age', 
    output_dir = './', 
    span_t = 0.75,
    clustering_method = 'complete',
    cluster_selection_method = 'silhouette'  # Options: 'silhouette', 'elbow', 'manual'
) {
  # Z-score normalization
 
  z_scored_counts <- t(scale(t(counts_matrix)))
  ages <- metadata[[age_column]]
  
  # Fit LOESS curves
  loess_curves <- apply(z_scored_counts, 1, function(row) {
    fit <- loess(row ~ ages, span = span_t)
    predict(fit, newdata = seq(min(ages), max(ages), 1))
  }) %>% t()
  
  # Compute distance matrix
  distance_matrix <- dist(loess_curves)
  
  # Cluster selection
  cluster_selection <- switch(cluster_selection_method,
                              'silhouette' = {
                                nb_result <- NbClust(loess_curves, distance = "euclidean", 
                                                     min.nc = 2, max.nc = 10, 
                                                     method = clustering_method, 
                                                     index = "silhouette")
                                nb_result$Best.nc[1]
                              },
                              'elbow' = {
                                wss <- sapply(2:20, function(k) {
                                  kmeans(loess_curves, k, nstart = 10)$tot.withinss
                                })
                                plot(2:20, wss, type = "b")
                                readline("Select number of clusters based on elbow plot: ") %>% as.numeric()
                              },
                              'manual' = {
                                # Hierarchical clustering dendrogram for manual selection
                                hc <- hclust(distance_matrix, method = clustering_method)
                                plot(as.dendrogram(hc), main = "Hierarchical Clustering Dendrogram")
                                readline("Select number of clusters: ") %>% as.numeric()
                              },
                              'manual_done' = {
                                readline("Select number of clusters: ") %>% as.numeric()
                              },
                              ,
                              'decision_done' = {
                                ###final decisions on numbers of clusters per tissue
                                switch(tissue,
                                       "Uterus"= {14},
                                       "Vagina" = {2},
                                       "Ovary" = {10},
                                       "Myometrium"= {6}
                                       
                                       
                                       )
                              }
                              
  )
  
  # Perform clustering
  clusters <- cutree(hclust(distance_matrix, method = clustering_method), 
                     k = cluster_selection)
  
  # Plotting
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  pdf(file.path(output_dir, paste0('gene_trajectory_clusters_dea',cluster_selection,'res.pdf')), 
      width = 15, height = 5 * ceiling(cluster_selection / 2))
  par(mfrow = c(ceiling(cluster_selection / 2), 2))
  
  for (cluster in 1:cluster_selection) {
    cluster_genes <- loess_curves[clusters == cluster, ]
    
    matplot(t(cluster_genes), type = 'l', 
            col = rgb(0.5, 0.5, 0.5, 0.1), 
            lty = 1, 
            main = paste0('Cluster ', cluster, 
                          ' (n=', nrow(cluster_genes), ' genes)'),
            xlab = 'Age Progression', 
            ylab = 'Z-scored Expression')
    
    lines(colMeans(cluster_genes), col = 'red', lwd = 3)
  }
  dev.off()
  
  # Return results
  list(
    z_scored_counts = z_scored_counts,
    loess_curves = loess_curves,
    distance_matrix = as.matrix(distance_matrix),
    clusters = clusters,
    num_clusters = cluster_selection
  )
}



# Load required libraries

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Function to create cluster plots
plot_gene_clusters <- function(loess_curves, clusters, cluster_selection, tissue, single_row = TRUE) {
  # color_tissue depending on tissue
  get_cluster_color <- function(tissue) {case_when(
    tissue=="Uterus" ~ "#6DB6FFFF",
    tissue=="Ovary" ~ "#009292FF",
    tissue=="Vagina" ~ "#B66DFFFF",
    TRUE ~ "black"
  )
  }
  # Convert data to long format for ggplot
  plot_list <- list()
  
  for (cluster in 1:cluster_selection) {
    # Get genes for current cluster
    cluster_genes <- loess_curves[clusters == cluster, ]
    
    # Convert to long format
    cluster_df <- as.data.frame(t(cluster_genes)) %>%
      mutate(Age = 21:70) %>%  # Adjust x-axis range from 20 to 70
      pivot_longer(cols = -Age, 
                   names_to = "Gene",
                   values_to = "Expression")
    
    # Calculate mean expression
    mean_expr <- cluster_df %>%
      group_by(Age) %>%
      summarize(Mean_Expression = mean(Expression))
    
    # Create plot
    p <- ggplot() +
      geom_hline(yintercept = 0.0, color = "gray")+
      
      # Individual gene trajectories
      geom_line(data = cluster_df,
                aes(x = Age, y = Expression, group = Gene),
                color = get_cluster_color(tissue),
                alpha = 0.3) +
      # Mean trajectory
      # Replace the mean trajectory geom_line and theme sections with:
      geom_line(data = mean_expr,
                aes(x = Age, y = Mean_Expression),
                color = get_cluster_color(tissue),
                size = 1.5) +
      labs(
        title = sprintf("Cluster %d (n=%d genes)", 
                        cluster, nrow(cluster_genes)),
        x = "Age (years)",
        y = "Z-scored Expression"
      ) +
      theme_minimal(base_size = 15) +
      theme(axis.text.y = element_text(size = 15, color = "black"),
            axis.text.x = element_text(size = 15, color = "black"),
            axis.title = element_text(size = 15, color = "black"),
            plot.title = element_text(size = 15, hjust = 0.5, color = "black"),
            legend.text = element_text(size = 15, color = "black"),
            panel.grid = element_blank(),                    # Remove all grid lines
            # axis.line.x = element_line(size = 0.5, color = "grey80"), # Add x-axis line
            # axis.line.y = element_line(size = 0.5, color = "grey80"),
            panel.border = element_rect(color = "grey80", size = 0.5, fill = NA))+
      # Labels
      
      # Set axis limits
      scale_x_continuous(limits = c(20, 70)) +
      # Add some padding to y-axis
      scale_y_continuous(
        limits = function(x) c(min(x) - 0.1 * diff(range(x)),
                               max(x) + 0.1 * diff(range(x))),
        labels = function(x) sprintf("%.1f", x)  # This formats all numbers to one decimal place
      ) 
    
    plot_list[[cluster]] <- p
  }
  
  # Calculate layout dimensions
  n_rows <- if(single_row) 1 else ceiling(cluster_selection / 2)
  n_cols <- if(single_row) cluster_selection else 2
  
  # Combine plots
  # Save as PDF
  pdf(file.path(output_dir, 
                paste0('pretty_gene_trajectory_clusters_dea',
                       cluster_selection,'res.pdf')),
      width = if(single_row) 20 else 4, 
      height = if(single_row) 5 else 4 * n_rows)
  
  do.call(grid.arrange, 
          c(plot_list, 
            ncol = n_cols, 
            nrow = n_rows))
  
  dev.off()
  # Save as PNG
  png(file.path(output_dir, 
                paste0('pretty_gene_trajectory_clusters_dea',
                       cluster_selection,'res.png')),
      width = if(single_row) 2000 else 1000, 
      height = if(single_row) 500 else 500 * n_rows,
      res = 100)
  
  do.call(grid.arrange, 
          c(plot_list, 
            ncol = n_cols, 
            nrow = n_rows))
  dev.off()
}



# Function to create cluster plots
plot_gene_clusters_one_by_one <- function(loess_curves, clusters, cluster_selection, tissue, single_row = TRUE) {
  # color_tissue depending on tissue
  get_cluster_color <- function(tissue) {case_when(
    tissue=="Uterus" ~ "#6DB6FFFF",
    tissue=="Ovary" ~ "#009292FF",
    tissue=="Vagina" ~ "#B66DFFFF",
    TRUE ~ "black"
  )
  }
  # Convert data to long format for ggplot
  plot_list <- list()
  
  for (cluster in 1:cluster_selection) {
    # Get genes for current cluster
    cluster_genes <- loess_curves[clusters == cluster, ]
    
    # Convert to long format
    cluster_df <- as.data.frame(t(cluster_genes)) %>%
      mutate(Age = 21:70) %>%  # Adjust x-axis range from 20 to 70
      pivot_longer(cols = -Age, 
                   names_to = "Gene",
                   values_to = "Expression")
    
    # Calculate mean expression
    mean_expr <- cluster_df %>%
      group_by(Age) %>%
      summarize(Mean_Expression = mean(Expression))
    
    # Create plot
    p <- ggplot() +
      geom_hline(yintercept = 0.0, color = "gray", linetype = "dashed")+
      
      # Individual gene trajectories
      geom_line(data = cluster_df,
                aes(x = Age, y = Expression, group = Gene),
                color = get_cluster_color(tissue),
                alpha = 0.3) +
      # Mean trajectory
      # Replace the mean trajectory geom_line and theme sections with:
      geom_line(data = mean_expr,
                aes(x = Age, y = Mean_Expression),
                color = get_cluster_color(tissue),
                size = 1.5) +
      labs(
        title = sprintf("Cluster %d (n=%d genes)", 
                        cluster, nrow(cluster_genes)),
        x = "Age (years)",
        y = "Z-scored Expression"
      ) +
      theme_minimal(base_size = 15) +
      theme(axis.text.y = element_text(size = 15, color = "black"),
            axis.text.x = element_text(size = 15, color = "black"),
            axis.title = element_text(size = 15, color = "black"),
            plot.title = element_text(size = 15, hjust = 0.5, color = "black"),
            legend.text = element_text(size = 15, color = "black"),
            panel.grid = element_blank(),                    # Remove all grid lines
            # axis.line.x = element_line(size = 0.5, color = "grey80"), # Add x-axis line
            # axis.line.y = element_line(size = 0.5, color = "grey80"),
            panel.border = element_rect(color = "grey80", size = 0.5, fill = NA))+
      # Labels
      
      # Set axis limits
      scale_x_continuous(limits = c(20, 70)) +
      # Add some padding to y-axis
      scale_y_continuous(
        limits = function(x) c(min(x) - 0.1 * diff(range(x)),
                               max(x) + 0.1 * diff(range(x))),
        labels = function(x) sprintf("%.1f", x)  # This formats all numbers to one decimal place
      ) 
    
    ### highlight AMH genes
    if (length(intersect("ENSG00000104899.8", rownames(cluster_genes)))==1){
      p <- p+ geom_line(aes(y = loess_curves["ENSG00000104899.8",], x=mean_expr$Age), color = "#490092FF")
    }
    plot_list[[cluster]] <- p
    pdf(file.path(output_dir, 
                  paste0('pretty_4gene_clusters',
                         cluster,'res.pdf')),
        width = 4, 
        height = 4
    )
    print(p)
    
    
    dev.off()
  }
  
  
  # Combine plots
  
  # Save as PNG
  
  
}

# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(ggplot2)
library(enrichplot)
library(stringr)
library(cowplot)

dotplot_custom <- function(ego, showCategory = 20) {

  # Create and save dotplot with custom theme
  ego <- ego[order(ego$p.adjust, decreasing = FALSE), ]  # Order by adjusted p value
  
  # ego <- ego[order(ego$logOdds, decreasing = TRUE), ]  # Order by logOdds
  ego <- head(ego, showCategory)  # Select top_n terms
  
  dotplot_custom <- ggplot(ego, aes(x = logOdds, y = reorder(Description, logOdds), 
                                   color = p.adjust, size = Count)) +
    geom_point() +
    scale_color_gradient(low = "#8B0000", high = "#63B8FF") +
    labs(x = "Log odds ratio", title = "", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(size = 14, color = "black"), 
          axis.text.x = element_text(size = 14, color = "black"),
          axis.title = element_text(size = 14, color = "black"), 
          plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
          legend.text = element_text(size = 14, color = "black"),
          panel.grid = element_blank(),                    
          panel.border = element_rect(color = "grey80", size = 0.5, fill = NA) 
    )+
    scale_y_discrete(labels = function(x) {
      stringr::str_wrap(x, width = 30)  # Wraps text without regex issues
    })
  return(dotplot_custom)
}

perform_enrichment_analysis <- function(gene_clusters, background_genes, output_dir) {
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Convert Ensembl IDs to Entrez IDs for compatibility
  ensembl_to_entrez <- select(org.Hs.eg.db, 
                              keys = background_genes,
                              columns = c("ENSEMBL", "ENTREZID"),
                              keytype = "ENSEMBL")
  
  # Remove duplicates and NAs
  ensembl_to_entrez <- ensembl_to_entrez[!duplicated(ensembl_to_entrez$ENSEMBL), ]
  ensembl_to_entrez <- na.omit(ensembl_to_entrez)
  
  # Convert background genes
  background_entrez <- ensembl_to_entrez$ENTREZID
  
  # Initialize results list
  enrichment_results <- list()
  
  # Process each cluster
  for (cluster_id in names(gene_clusters)) {
    # Convert cluster genes to Entrez IDs
    cluster_genes <- gene_clusters[[cluster_id]]
    cluster_entrez <- ensembl_to_entrez$ENTREZID[ensembl_to_entrez$ENSEMBL %in% cluster_genes]
    
    # Create cluster directory
    cluster_dir <- file.path(output_dir, paste0("cluster_", cluster_id))
    dir.create(cluster_dir, showWarnings = FALSE)
    
    # Initialize cluster results
    cluster_results <- list()
    
    # 1. GO Analysis (BP, MF, CC)
    ont <- c("BP")  # Can be expanded to include "MF", "CC"
    
      ego <- enrichGO(gene = cluster_entrez,
                      universe = background_entrez,
                      OrgDb = org.Hs.eg.db,
                      ont = ont,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)
      
      ego_sum <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
      
      # Calculate additional metrics
      if (!is.null(ego_sum) && nrow(ego_sum@result) > 0) {
        ego_sum@result$GeneRatio_num <- as.numeric(sub("/.*", "", ego_sum@result$GeneRatio)) / 
          as.numeric(sub(".*/", "", ego_sum@result$GeneRatio))
        
        ego_sum@result$BgRatio_num <- as.numeric(sub("/.*", "", ego_sum@result$BgRatio)) / 
          as.numeric(sub(".*/", "", ego_sum@result$BgRatio))
        
        ego_sum@result$logOdds = log((ego_sum@result$GeneRatio_num/ego_sum@result$BgRatio_num)/
                                       (1-ego_sum@result$GeneRatio_num)/(1-ego_sum@result$BgRatio_num))
        
        # Save results to file
        write.csv(ego@result, 
                  file = file.path(cluster_dir, paste0("GO_", ont, "_results.csv")))
        
        # Store results in list
        cluster_results[[ont]] <- list(
          ego = ego,
          ego_sum = ego_sum,
          output_path = cluster_dir
        )
      }
    
    
    # Add cluster results to main results list
    enrichment_results[[cluster_id]] <- cluster_results
  }
  
  return(enrichment_results)
}



# Function to convert named integer vector to list of gene clusters
convert_to_clusters <- function(named_clusters) {
  # Get unique cluster numbers
  cluster_nums <- sort(unique(named_clusters))
  
  # Create list where each element contains genes for one cluster
  clusters_list <- lapply(cluster_nums, function(cluster_num) {
    gsub("\\..*", "", names(named_clusters)[named_clusters == cluster_num])
  })
  
  # Name the list elements as "cluster1", "cluster2", etc.
  names(clusters_list) <- paste0("cluster", cluster_nums)
  
  return(clusters_list)
}


####Atempting to plot enrichment plots next to clusters

library(ggplot2)
library(gridExtra)
library(cowplot)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
clusters_with_enrichment<-c()
#Load the data

#####################

plot_enrichment_results <- function(enrichment_results) {
  # Loop through each cluster
  list_plots<-c()
  for (cluster_id in names(enrichment_results)) {
    cluster_results <- enrichment_results[[cluster_id]]
    
    # Loop through each ontology category
    
    # Get results for this cluster and ontology
    result <- cluster_results[["BP"]]
    ego <- result$ego
    ego_sum <- result$ego_sum
    output_path <- result$output_path
    # Check if there are significant results
    if (!is.null(ego) && nrow(ego@result) > 0 && sum(ego@result$p.adjust < 0.05) > 0) {
      # Generate and save plot
      
      # Custom dotplot
      pdf(file.path(output_path, paste0("GO_BP_dotplot_custom1.pdf")),
          width = 5.5, height = 5)
      
      m<-dotplot_custom(ego_sum, showCategory = 10)+plot_layout(nrow = 1, heights = unit(1, "null"))
      print(m)
      dev.off()
      
      list_plots[[cluster_id]]<-dotplot_custom(ego_sum, showCategory = 10)
      
      
    }
    
  }
  return(list_plots)
}

library(vroom)
prfx
prfx_data = "X/Ole/"
gene_annotation <- vroom(paste0(prfx_data,"tissues/v10/gencode.v39.annotation.bed"))
prfx = "X/"
tissue = "Uterus"
#FOR GTEx v10
counts_tot <- readRDS(paste0(prfx,"Ole/tissues/v10/", tissue, "/counts.rds"))
metadata <- readRDS(paste0(prfx,"Ole/tissues/v10/", tissue, "/metadata.rds"))
prfx2 = "X/Laura/12.Tissues_substructures/DEA/DEA_v10/"
 ## loading the DEGs to subset the clustering to them
dea_res = readRDS(paste0(prfx2, tissue,"/all_cov_no_inter/", tissue, "_all_cov_no_inter_AGE_covariates_and_traits.results.rds"))
dea_res_age = dea_res$Age[dea_res$Age$adj.P.Val<0.05,]

##for loess smoothing
span_t = 0.75

#### analysis 
output_dir = paste0("X/Laura/Figure_plots/", tissue,"_span",span_t, "_manual12")

common_samples = intersect(metadata$Sample, colnames(counts_tot))

##perform clustering
results <- gene_trajectory_analysis(
  counts_tot[intersect(rownames(counts_tot), rownames(dea_res_age)),common_samples], 
  metadata[metadata$Sample %in% common_samples,], 
  span_t = span_t,
  output_dir = output_dir,
  cluster_selection_method = 'decision_done',
  clustering_method = 'complete'
)


plot_gene_clusters(results$loess_curves, results$clusters, results$num_clusters, tissue, single_row = FALSE)
plot_gene_clusters_one_by_one(results$loess_curves, results$clusters, results$num_clusters, tissue, single_row = FALSE)

background_genes <- gsub("\\..*", "",rownames(counts_tot))

gene_clusters = convert_to_clusters(results$clusters)
enrichment_results <- perform_enrichment_analysis(gene_clusters, background_genes, output_dir)
l_plots<-plot_enrichment_results(enrichment_results)
