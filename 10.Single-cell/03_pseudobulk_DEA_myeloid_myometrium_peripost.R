#### this is just illustrative of ONE CELL TYPE pseudobulk + DEA example, in this case it takes the age model BUT
# the model that was taken forward for the GSEA is the following, comparing postmenopausal vs premenopausal groups:
## form <- ~ Menopausal.status + ever_given_birth + (1 | Patient)

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(edgeR)
  library(limma)
  library(variancePartition)
  library(BiocParallel)
  library(dplyr)
  library(Seurat)
  library(Matrix)
})

## Demo test for fibroblasts

so_subset <- readRDS("/X/Data/scRNAseq/PunzonGimenez2024_Myometrium/seurat_objects/Myeloid cells_AllZones_SeuratObject.rds")

Group_order.vec <- list(Menopausal.status = c('Peri', 'Post'))

#so_subset@meta.data$Previous.pregnancy <- as.numeric(so_subset@meta.data$Previous.pregnancy)
so_subset@meta.data$Live.births <- as.numeric(so_subset@meta.data$Live.births)
so_subset@meta.data$ever_given_birth <- ifelse(so_subset@meta.data$Live.births > 0, 1, 0)

so_subset@meta.data$Menopausal.status <- factor(so_subset@meta.data$Menopausal.status,
                                                levels = Group_order.vec[['Menopausal.status']])

so_subset@meta.data$Age <- as.numeric(so_subset@meta.data$Age)

# Optional metadata cleaning (Surgery.indication removed from model)
so_subset@meta.data$Surgery.indication <- gsub("Hysterectomy due to prolapse", "Hysterectomy", so_subset@meta.data$Surgery.indication)
so_subset@meta.data$Surgery.indication <- gsub("Uterine removal from organ donor programme", "Uterine_removal", so_subset@meta.data$Surgery.indication)

sce_subset <- as.SingleCellExperiment(so_subset, assay = "RNA")

# Ensure colnames consistent
colnames(sce_subset) <- colnames(so_subset)
stopifnot(identical(colnames(sce_subset), rownames(colData(sce_subset))))

# Create combined sample identifier
sce_subset$Patient_Zone <- paste0(sce_subset$Patient, "_", sce_subset$Zone)

# Aggregate pseudobulk counts matrix
mat <- counts(sce_subset, assay = "counts")
if (!inherits(mat, "dgCMatrix")) {
  mat <- as(mat, "dgCMatrix")
}
grouping <- sce_subset$Patient_Zone
sample_matrix <- sparse.model.matrix(~ 0 + grouping)
pseudobulk_counts <- as.matrix(mat %*% sample_matrix)
colnames(pseudobulk_counts) <- gsub("^grouping", "", colnames(pseudobulk_counts))

# Extract sample metadata matching pseudobulk samples
sample_md <- as.data.frame(colData(sce_subset)) %>%
  dplyr::group_by(Patient_Zone) %>%
  dplyr::summarise(
    #    Menopausal.status = unique(Menopausal.status),
    #    Previous.pregnancy = unique(Previous.pregnancy),
    Age = unique(Age),
    ever_given_birth = unique(ever_given_birth),
    Patient = unique(Patient),
    .groups = "drop"
  ) %>%
  as.data.frame()

rownames(sample_md) <- sample_md$Patient_Zone
sample_md <- sample_md[colnames(pseudobulk_counts), ]

# Filter lowly expressed genes
keep <- rowSums(edgeR::cpm(pseudobulk_counts) > 1) >= 2
pb_filtered <- pseudobulk_counts[keep, ]


# --------- SUBSET 20 WELL-EXPRESSED GENES FOR QUICK TEST --------- NB if any sample has 0 counts for those 20 genes execution will be halted
#this i only did it for the first run with fibroblasts, but whwn i saw it takes approx 20mins i skip it

#set.seed(123)
#expressed_genes <- rownames(pb_filtered)[rowSums(pb_filtered > 0) >= ncol(pb_filtered) / 2]
#test_genes <- sample(expressed_genes, 20)
#pb_test <- pb_filtered[test_genes, ]
#pb_test <- pb_test[, colSums(pb_test) > 0]

# ---------------------------------------------------

# DGEList for test
#dge_test <- DGEList(counts = pb_test)
#dge_test <- calcNormFactors(dge_test)
###################################################################################################
# DEA model
#form <- ~ Menopausal.status + ever_given_birth + (1 | Patient)

form <- ~ Age + ever_given_birth + (1 | Patient)

# Parallel backend
param <- SnowParam(4, "SOCK", progressbar = TRUE)
register(param)

# voom + dream on test
#vobj_test <- voomWithDreamWeights(dge_test, form, sample_md, BPPARAM = param)
#fit_test <- dream(vobj_test, form, sample_md)
#res_test <- topTable(fit_test, coef = "Menopausal.statusPost", number = Inf)
#saveRDS(res_test, file = "/X/results_DEA/SMC_TEST_DEA_results_Menopausal.status_Post_vs_Peri.rds")

# --- Full run on all genes (remove subset) ---
dge <- DGEList(counts = pb_filtered)
dge <- calcNormFactors(dge)
vobj <- voomWithDreamWeights(dge, form, sample_md, BPPARAM = param)
fit <- dream(vobj, form, sample_md)
res <- topTable(fit, coef = "Age", number = Inf)

# Save final result
saveRDS(res, file = "/X/results_DEA/Myeloid_cells_DEA_results_Age.rds")

