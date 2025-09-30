library(clusterProfiler)
library(org.Hs.eg.db)
library(WebGestaltR)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(ggplot2)
library(ggtext)
library(writexl)


# === FILE LISTS ===
 dea_files <- c(
  # Age model   #these were ran initially but at the end menopausal model is better
  # "/X/SC_ANALYSES/DEA_results/Endothelium_DEA_results_Age.rds",
  # "/X/SC_ANALYSES/DEA_results/Epithelial_like_cells_DEA_results_Age.rds",
  #  "/X/SC_ANALYSES/DEA_results/Fibroblasts_DEA_results_Age.rds",
  #  "/X/SC_ANALYSES/DEA_results/LEC_DEA_results_Age.rds",
  # "/X/SC_ANALYSES/DEA_results/Lymphoid_cells_DEA_results_Age.rds",
  # "/X/SC_ANALYSES/DEA_results/Mast_cells_DEA_results_Age.rds",
  #  "/X/SC_ANALYSES/DEA_results/Myeloid_cells_DEA_results_Age.rds",
  #  "/X/SC_ANALYSES/DEA_results/PNS_DEA_results_Age.rds",
  #  "/X/SC_ANALYSES/DEA_results/Perivascular_DEA_results_Age.rds",
  # "/X/SC_ANALYSES/DEA_results/SMC_DEA_results_Age.rds",
  
  # Menopause model

  "/X/SC_ANALYSES/DEA_results/Perivascular_DEA_results_Menopausal.status_Post_vs_Peri.rds",
  "/X/SC_ANALYSES/DEA_results/SMC_DEA_results_Menopausal.status_Post_vs_Peri.rds",
  "/X/SC_ANALYSES/DEA_results/Endothelium_DEA_results_Menopausal.status_Post_vs_Peri.rds",
  "/X/SC_ANALYSES/DEA_results/Epithelial_like_cells_DEA_results_Menopausal.status_Post_vs_Peri.rds",
  "/X/SC_ANALYSES/DEA_results/Fibroblasts_DEA_results_Menopausal.status_Post_vs_Peri.rds",
  "/X/SC_ANALYSES/DEA_results/LEC_DEA_results_Menopausal.status_Post_vs_Peri.rds",
  "/X/SC_ANALYSES/DEA_results/Lymphoid_cells_DEA_results_Menopausal.status_Post_vs_Peri.rds",
  "/X/SC_ANALYSES/DEA_results/Mast_cells_DEA_results_Menopausal.status_Post_vs_Peri.rds",
  "/X/SC_ANALYSES/DEA_results/Myeloid_cells_DEA_results_Menopausal.status_Post_vs_Peri.rds",
  "/X/SC_ANALYSES/DEA_results/PNS_DEA_results_Menopausal.status_Post_vs_Peri.rds"
)


# === GENE ANNOTATION ===
gene_annotation <- read.csv("/X/SC_ANALYSES/hpc/projects/bsc83/Projects/GTEx_v8/Laura/00.Data/gencode.v39.annotation.bed")

# Retrieve Entrez IDs
entrez.id <- idMapping(organism = "hsapiens",
                       inputGene = sub("\\..*", "", gene_annotation$ensembl.id),
                       sourceIdType = "ensembl_gene_id",
                       targetIdType = "entrezgene")$mapped

gene_annotation$entrez.id <- sapply(
  sub("\\..*", "", gene_annotation$ensembl.id),
  function(gene) ifelse(
    gene %in% entrez.id$userId,
    entrez.id[entrez.id$userId == gene, "entrezgene"],
    NA
  )
)

# === MAIN LOOP ===
for (file in dea_files) {
  
  fname <- basename(file)
  fname_noext <- sub("\\.rds$", "", fname)
  cell_type <- sub("_DEA_results.*", "", fname_noext)
  model <- sub(".*DEA_results_", "", fname_noext)
  
  message("Processing: ", cell_type, " - ", model)
  
  dea <- readRDS(file)
  dea$gene.name.x <- row.names(dea)
  
  n <- gene_annotation[gene_annotation$gene.name.x %in% dea$gene.name.x,
                       c("gene.name.x", "entrez.id")]
  
  degs <- merge(dea, n, by = "gene.name.x")

  # Signed ranking: sign(logFC) * -log10(P.Value)
  gene_ranks <- with(degs, sign(logFC) * -log10(P.Value))
  names(gene_ranks) <- degs$entrez.id
  gene_ranks <- gene_ranks[!duplicated(names(gene_ranks))]
  gene_ranks <- gene_ranks[!is.na(names(gene_ranks))]
  gene_ranks <- sort(gene_ranks, decreasing = TRUE)
  
  gsea_result <- gseGO(
    geneList = gene_ranks,
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    keyType = "ENTREZID",
    verbose = FALSE
  )
  
  gsea_simplified <- simplify(gsea_result, cutoff = 0.5, measure = "Wang")
  r <- gsea_simplified@result[gsea_simplified@result$p.adjust < 0.05,
                              c("Description", "setSize", "enrichmentScore", "NES", "p.adjust")]
  
  saveRDS(r, paste0("/X/SC_ANALYSES/GSEA_results/GSEA_", 
                    cell_type, "_", model, "_goterms.rds"))
  


