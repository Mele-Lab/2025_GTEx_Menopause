#GSEA for MOFA factors

library(clusterProfiler)
library(org.Hs.eg.db)
library(WebGestaltR)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(ggtext)
library("MOFA2")
library("MOFAdata")



gene_annotation  <- read.csv("~/Desktop/TFM/bsc83671/GTEx_v8//Laura/00.Data/gencode.v39.annotation.bed")

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


tissue<-"Ovary"

MOFAobject_trained<-readRDS(paste0(input, "/Laura/MOFA/MOFA_", tissue, "_only_matching_donors_10_ensg_CORRECTED.rds"))

p<-plot_top_weights(MOFAobject_trained,
                    view = "gene_expression",
                    #view = "medulla",
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
saveRDS(r, paste0("Desktop/TFM/bsc83671/GTEx_v8/Laura/12.Tissues_substructures/Enrichment/Enrichments_GSEA/gsea_", tissue, "_goterms_MOFA.rds"))
tail(r)
tissue
r_down<-r[r$NES<0,]

# Define terms to bold (modify as needed)
bold_terms <- c("extracellular matrix structural constituent", "Z disc", "synapse")  
bold_terms <- c("lipid metabolic process", "epithelium development")  

# Format y-axis labels: Make selected terms bold
r_down$Description <- ifelse(r_down$Description %in% bold_terms, 
                             paste0("<b>", r_down$Description, "</b>"), 
                             r_down$Description)

#
r_down<-r[1:20,]
###Uterus filtering
r_down<-r_down[r_down$setSize>40,]
plot_upregulated <- ggplot(r_down, aes(x = NES, y = reorder(Description, NES), color = p.adjust)) +
  geom_point(size=4) +
  scale_color_gradient(low = "#8B0000", high = "#63B8FF") +
  labs(x = "NES (Normalized Enrichment Score)", title ="", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_markdown(size = 14, color = "black"), 
        axis.text.x = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, color = "black"), 
        plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 14, color = "black"),
        panel.grid = element_blank(),                    # Remove all grid lines
        # axis.line.x = element_line(size = 0.5, color = "grey80"), # Add x-axis line
        # axis.line.y = element_line(size = 0.5, color = "grey80"),
        panel.border = element_rect(color = "grey80", size = 0.5, fill = NA)
        # Add full border
  )

pdf(paste0("~/X/Laura/MOFA/Plots/",tissue, "_enrichment_factor4_down_PAPER.pdf"), width = 8, height = 4)  # Adjust width and height as needed
plot_upregulated
dev.off()
