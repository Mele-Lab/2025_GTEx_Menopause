#Fig 3E-I
library(ComplexHeatmap)
library("MOFA2")
library("MOFAdata")
library(dplyr)
library(tidyr)
library(psych)
library(limma)
library(ggplot2)
library(edgeR)
library(readxl)


input<-"X"
options<- c("all_donors", "only_matching_donors")
option<-options[2]
tissue<-"Ovary" #choose uterus and vagina
n_factors=10
#Read MOFA results
MOFAobject_trained<-readRDS(paste0(input, "/Laura/MOFA/MOFA_", tissue, "_",option,"_","10_ensg_CORRECTED.rds"))

#For ovary, including follicles
MOFAobject_trained<-readRDS(paste0(input, "/Laura/Paper_material/Post_review_analyses/MOFA/MOFA_", tissue, "_",option,"_",n_factors,"_ensg_CORRECTED.rds"))

#Read and filter metadata
metadata_tot <- readRDS(paste0(input, "/Laura/00.Data/v10/", tissue, "/metadata.rds"))
filter<-readRDS(paste0(input, "/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
metadata_tot$Ancestry<- NULL

#Add continuous ancestry
ancestry_file<- read.table(paste0(input,"/Laura/00.Data/GTEx_phenotypes/admixture_inferred_ancestry.txt"))
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
variance<-get_variance_explained(MOFAobject_trained)

#Obtain r values for the correlations with age
covs_to_inspect<-c("Age")

r_df<-as.data.frame(correlate_factors_with_covariates(MOFAobject_trained,
                                                      covariates = c(covs_to_inspect),
                                                      #covariates = colnames(metadata_MOFA_nonas),
                                                      plot="r",
                                                      return_data = TRUE))
r_df <- r_df %>%
  tibble::rownames_to_column(var = "Factor") %>%
  pivot_longer(cols = covs_to_inspect,  # Stack all feature columns
               names_to = "Traits",
               values_to = "r")

log_df<-as.data.frame(correlate_factors_with_covariates(MOFAobject_trained,
                                                        covariates = c(covs_to_inspect),
                                                        #covariates = colnames(metadata_MOFA_nonas),
                                                        plot="log_pval",
                                                        return_data = TRUE))

log_df <- log_df %>%
  tibble::rownames_to_column(var = "Factor")%>%
  pivot_longer(cols = covs_to_inspect,  # Stack all feature columns
               names_to = "Traits",
               values_to = "log_pval")

log_df$pval<-10^(-log_df$log_pval)
log_df$as<-ifelse(log_df$pval<0.05 & log_df$pval>0.005, "*",
                  ifelse(log_df$pval<0.005 & log_df$pval>0.0005, "**",
                         ifelse(log_df$pval<0.0005, "***", "")))

#Join both dataframes with r and log values
long_df<-cbind(r_df, log_df[,c("pval","as")])
p2<-plot_variance_explained(MOFAobject_trained, max_r2=max(variance$r2_per_factor$group1))
p3<-plot_variance_explained(MOFAobject_trained, plot_total = T)[[2]]

# Define your color palette (example: using a custom color palette)
color_palette <-  brewer.pal(9, "BuPu")
df_to_plot<-long_df %>%
  pivot_wider(id_cols = "Factor",
              names_from = "Traits",
              values_from ="r")

df_to_plot<-as.data.frame(df_to_plot)
rownames(df_to_plot)<-df_to_plot$Factor
df_to_plot$Factor<-NULL
df_to_plot<-t(df_to_plot)

# Plot the heatmap

c<-c("#63B8FF", "#8B0000")
# Prepare the asterisk matrix
asterisks_mat <- long_df %>%
  select(Factor, Traits, as) %>%
  pivot_wider(id_cols = "Factor", names_from = "Traits", values_from = "as")

# Set rownames
asterisks_mat <- as.data.frame(asterisks_mat)
rownames(asterisks_mat) <- asterisks_mat$Factor
asterisks_mat$Factor <- NULL
asterisks_mat <- t(asterisks_mat)  # Transpose
ht<-Heatmap(df_to_plot,
            name = "r*log10(adj p-val)",
            # col = colorRampPalette(c("white", "#B22222"))(100),
            col = colorRampPalette(c)(100),
            
            cluster_columns = FALSE,
            cluster_rows = FALSE,
            
            column_names_rot = 25,
            column_names_centered = FALSE,
            show_column_names = FALSE,  # Hides column names
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 14),
            rect_gp = gpar(col = "grey50", lwd = 0.5),
            # column_names_gp = gpar(fontsize = 14, just= "center"),
            #
            height  = unit(length(covs_to_inspect), "cm"),
            
            show_row_names = TRUE,
            heatmap_legend_param = list(
              direction = "horizontal",
              legend_title_position = "topright",
              labels_gp = gpar(fontsize = 12),
              title = "r",
              legend_height = unit(3, "cm"),
              legend_width = unit(3, "cm"),
              grid_width = unit(0.5, "cm")
            ),
            # Add asterisks inside cells
            cell_fun = function(j, i, x, y, width, height, fill) {
              if (asterisks_mat[i, j] != "") {  # If not empty
                grid.text(
                  asterisks_mat[i, j],
                  x, y,
                  gp = gpar(fontsize = 20, col = "black")
                )
              }
            }
)# Set horizontal legend

# draw(ht, heatmap_legend_side = "right")

#################
heatmap_data <- p2$data  # Your original data frame with factor, view, and value

# Reshape the data with 'factor' as rows and 'view' as columns
heatmap_data<-heatmap_data[,1:3]
heatmap_matrix <- reshape(heatmap_data,
                          idvar = "factor",
                          timevar = "view",
                          direction = "wide")
colnames(heatmap_matrix) <- gsub("value.", "", colnames(heatmap_matrix))
rownames(heatmap_matrix) <- heatmap_matrix$factor

# Remove the 'factor' column as it's now the row names
heatmap_matrix <- heatmap_matrix[, -1]
heatmap_matrix<-t(heatmap_matrix)


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)
library(grid)

row_annotation <- rowAnnotation(
  "Var (%)" = anno_barplot(p3$data$R2,
                           gp = gpar(fill = pal, col = pal),  # Fill colors per tissue type
                           border = FALSE,
                           bar_width = 1,
                           axis_param = list(
                             side = "bottom",  # Ensures labels appear at the bottom
                             gp = gpar(fontsize = 10)  # Increase x-axis tick label size
                           )
  ),
  width = unit(1.5, "cm"),  # 👈 Make it wider here
  gap = unit(0.25, "cm"),
  show_legend = TRUE,
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 14)
)
ht2<-Heatmap(heatmap_matrix,
             col = color_palette,
             cluster_columns = FALSE,
             cluster_rows = FALSE,
             
             column_names_rot = 25,
             column_names_centered = FALSE,
             row_names_side = "left",
             row_names_gp = gpar(fontsize = 14),
             rect_gp = gpar(col = "grey50", lwd = 0.5),
             column_names_gp = gpar(fontsize = 14),
             #
             height = unit(length(levels(p2$data$view)), "cm"),
             
             show_row_names = TRUE,
             heatmap_legend_param = list(
               direction = "horizontal",
               legend_title_position = "topright",
               labels_gp = gpar(fontsize = 12),
               title = "Var (%)",
               legend_height = unit(3, "cm"),
               legend_width = unit(3, "cm"),
               grid_width = unit(0.5, "cm")
             ),
             right_annotation = row_annotation
)

# draw(ht2, heatmap_legend_side = "right")

# Draw the heatmaps vertically aligned
ht_list <- ht %v% ht2  # Merge them vertically

# Render the combined heatmap
# draw(ht_list, heatmap_legend_side = "bottom", merge_legend = TRUE)

if(n_factors==10){
  width<- 7
}else{
  width<-8.5
}


pdf(paste0("X/Paper_material/Post_review_analyses/MOFA/",tissue, "_MOFA_",n_factors,".pdf"), width = width, height = unit(length(levels(p2$data$view)), "cm"))  # Adjust width and height as needed
draw(ht_list, heatmap_legend_side = "bottom", merge_legend = TRUE)
dev.off()


########## MOFA enrichments
library(clusterProfiler)
library(org.Hs.eg.db)
library(WebGestaltR)
library(AnnotationDbi)
library(ReactomePA)
library(GSEA)

tissue<- tissue_of_interest # choose uterus factor 4, vagina fator 5 and ovary factor 7

#################for GO database
gene_annotation  <- read.csv(paste0(input,"/Laura/00.Data/gencode.v39.annotation.bed"))

# Retrieve entrez id
entrez.id <- idMapping(organism = "hsapiens",
                       inputGene = sub("\\..*", "", gene_annotation$ensembl.id),
                       sourceIdType = "ensembl_gene_id",
                       targetIdType = "entrezgene")$mapped

# Add gene entrez_id to gene_annotation data
gene_annotation$entrez.id <- sapply(sub("\\..*", "", gene_annotation$ensembl.id), function(gene)
  ifelse(gene %in% entrez.id$userId,
         entrez.id[entrez.id$userId==gene, "entrezgene"],
         NA)
)
# MOFAobject_trained<-readRDS(paste0(input, "/Laura/MOFA/MOFA_", tissue, "_",option,"_","10_ensg_CORRECTED.rds"))

#Obtain the weights for the factor of interst for our genes
p<-plot_top_weights(MOFAobject_trained,
                    view = "gene_expression",
                    factor = 5,
                    nfeatures = 5000,     # Top number of features to highlight
                    scale = T           # Scale weights from -1 to 1
)


degs_scores<-p$data[,c("feature","value", "sign")]
degs_scores$value <- ifelse(degs_scores$sign == "-", -degs_scores$value, degs_scores$value)
degs_scores$sign <- NULL
gene_info<- gene_annotation[sub("\\..*", "", gene_annotation$ensembl.id) %in% degs_scores$feature, c("ensembl.id", "gene.name.x")]
gene_info$ensembl.id<-sub("\\..*", "", gene_info$ensembl.id)
degs_scores<- merge(gene_info, degs_scores, by.x= "ensembl.id", by.y="feature")
degs_scores<- degs_scores[order(degs_scores$value, decreasing= TRUE),]

n<-gene_annotation[gene_annotation$gene.name.x %in% degs_scores$gene.name.x, c("gene.name.x", "entrez.id")]
degs<-merge(degs_scores, n, by= "gene.name.x")
degs<- degs[order(degs$value, decreasing = TRUE),]

gene_ranks <- setNames(degs$value, degs$entrez.id)
gene_ranks <- gene_ranks[!is.na(names(gene_ranks))]
gene_ranks <- gene_ranks[!duplicated(names(gene_ranks))]  # Remove duplicates

gsea_result <- gseGO(geneList = gene_ranks,
                     OrgDb = org.Hs.eg.db,
                     ont = "ALL",          # Can be "BP", "CC", or "MF"
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 0.05,
                     verbose = TRUE)

gsea_simplified <- simplify(gsea_result, 
                            cutoff = 0.7, 
                            measure = "Wang")


# View the simplified result
r<-gsea_simplified@result[gsea_simplified@result$p.adjust<0.05, c("Description","setSize", "enrichmentScore", "NES", "p.adjust")]

#For uterus and vagina, we select the downregulated terms:
r_down<-r[r$NES<0,]
r_down<-r_down[order(r_down$p.adjust),]
r_down<-r_down[r_down$setSize>40,]

r_t<- r_down

#For ovary, we select the upregulated terms
r_up<-r[r$NES>0,]
r_up<-r_up[order(r_up$p.adjust),]
r_up<-r_up[r_up$setSize>40,]
r_up<-r_up[1:12,]
r_t<-r_up

r_t<-rbind(r_down, r_up)
library(ggtext)
library(gridExtra)
library(patchwork)
library(stringr)
library(cowplot)
library(reshape2)
#PLOT
plot_enrich <- ggplot(r_t, aes(x = NES, y = reorder(Description, NES), color = p.adjust)) +
  geom_point(size=4) +
  scale_color_gradient(low = "#8B0000", high = "#63B8FF") +
  labs(x = "NES (Normalized Enrichment Score)", title ="", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_markdown(size = 14, color = "black"), 
        axis.text.x = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 14, color = "black"), 
        plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 14, color = "black"),
        panel.grid = element_blank(),                    # Remove all grid lines
        # axis.line.x = element_line(size = 0.5, color = "grey80"), # Add x-axis line
        # axis.line.y = element_line(size = 0.5, color = "grey80"),
        panel.border = element_rect(color = "grey80", size = 0.5, fill = NA)
        # Add full border
  )+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 35)
  )  


plot_enrich
factor_of_interest <- "5"
fac<-factor_of_interest # choose the factor
pdf(paste0("X/Figures/enrichment",tissue, "factor", fac,".pdf"), width = 5.5, height = 6) 

plot_enrich
dev.off()


