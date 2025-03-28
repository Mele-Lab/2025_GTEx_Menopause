#!/usr/bin/env Rscript

#' Differential Expression - Sliding Window ANalysis

Sys.setenv(TZ="Europe/Madrid")
# ---------------------- #
start_time <- Sys.time()
# ---------------------- #
# Libraries ####
suppressMessages(library(edgeR))
suppressMessages(library(limma))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(parallel))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(parallel))
suppressMessages(library(data.table))
suppressMessages(library(gplots))

# Command line arguments ####
args <- commandArgs(trailingOnly=TRUE)

if (getwd()== "/X/Laura/sliding_window_DEA" ){
  input<-"/X/"
  args <- commandArgs(trailingOnly = TRUE)
  tissue<- args[1]
  sex<- args[2]
  outpath<- paste0("/X/Laura/sliding_window_DEA/All_tissues/", tissue, "/")
  if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}
  
}else{
  input<-"/X//"
  tissue<-"Uterus"
  sex<- "female"
  outpath<- paste0("/X//Laura/sliding_window_DEA/All_tissues/", tissue, "/")
  if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}
  
}

# Functions ####
source(paste0(input, "/Laura/01.DEA/DEA_scripts/DEA_and_DSA.R_functions.R"))

#FOR GTEX V10
gene_annotation  <- read.csv(paste0(input, "/Laura/00.Data/gencode.v39.annotation.bed"))
Y_genes <- gene_annotation[gene_annotation$chr=="chrY",]$ensembl.id

# Tissue ----
tissues<-c("Uterus", "Ovary", "Vagina", "BreastMammaryTissue", "CervixEndocervix", "CervixEctocervix", "FallopianTube")
both_tissues<-c("BreastMammaryTissue" ,"ArteryAorta"   ,"ArteryTibial","ArteryCoronary"  )
tissues_with_filter<- c("Uterus", "Ovary", "Vagina", "BreastMammaryTissue", "CervixEndocervix", "CervixEctocervix", "FallopianTube","ArteryAorta"   ,"ArteryTibial","ArteryCoronary"  )
all_tissues <- basename(list.dirs(paste0(input, "Laura/00.Data/v10/"), full.names = TRUE, recursive = FALSE))
tissues_without_filter<- all_tissues[!all_tissues %in% tissues_with_filter]
  
# 1.1 Read in data ----
#FOR GTEx v10
if(tissue == "MuscleSkeletal"){
  metadata_tot <- readRDS(paste0(input, "Laura/00.Data/1.Final_tissues_to_use/", tissue, "/metadata.rds"))
  counts_tot <- readRDS(paste0(input, "Laura/00.Data/1.Final_tissues_to_use/", tissue, "/counts.rds"))
  tpm_tot <- readRDS(paste0(input, "Laura/00.Data/1.Final_tissues_to_use/", tissue, "/tpm.rds"))
}else{
  metadata_tot <- readRDS(paste0(input, "Laura/00.Data/v10/", tissue, "/metadata.rds"))
  counts_tot <- readRDS(paste0(input, "Laura/00.Data/v10/", tissue, "/counts.rds"))
  tpm_tot <- readRDS(paste0(input, "Laura/00.Data/v10/", tissue, "/tpm.rds"))
  
}

metadata_tot$Ancestry<- NULL
ancestry_file<- read.table(paste0(input, "/Laura/00.Data/admixture_inferred_ancestry.txt"))
ances<-ancestry_file[,c(1,3)]
colnames(ances)<- c("Donor", "Ancestry")
metadata_tot<-merge(metadata_tot, ances, by= "Donor")

if (tissue %in% tissues_with_filter){
   if(tissue %in% both_tissues){
      if(sex == "male"){
        filter<-readRDS(paste0(input, "/Laura/03.Image_processing/Second_filtering_images/",tissue, "_male_filtered_images.rds"))
        metadata<- metadata_tot[metadata_tot$Donor %in% filter$Subject.ID,]
        # metadata$Ancestry <- relevel(metadata$Ancestry, ref="EUR")
        
      }else{
        filter<-readRDS(paste0(input, "/Laura/03.Image_processing/Second_filtering_images/",tissue, "_female_filtered_images.rds"))
        metadata<- metadata_tot[metadata_tot$Donor %in% filter$Subject.ID,]
        # metadata$Ancestry <- relevel(metadata$Ancestry, ref="EUR")
      }
  }else if (tissue %in% tissues){
    filter<-readRDS(paste0(input, "/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
    metadata<- metadata_tot[(metadata_tot$Donor %in% filter$Subject.ID),]   
  }
  
}else{
  if (tissue %in% c("Prostate", "Testis")){
    metadata<-metadata_tot
  }else {
    if (sex=="male"){
    metadata<-metadata_tot[metadata_tot$Sex == "1",]
    
  }else
   metadata<-metadata_tot[metadata_tot$Sex == "2",]
  }
}


  
counts <- counts_tot[, metadata$Sample]
tpm<-tpm_tot[, metadata$Sample]


# Gene annotation ----

# 1.2 Genes expressed per tissue ----
# 1.2.1 TPM>=0.1 in at least 20% of the tissue samples
exprs_genes.tpm <- rownames(tpm)[apply(tpm, 1, function(x) sum(x>=1) ) >= 0.2*ncol(tpm)  ]  

# 1.2.2 Count >=6 in at least 20% of the tissue samples
exprs_genes.counts <- rownames(counts)[ apply(counts, 1, function(x) sum(x>=10) ) >= 0.2*ncol(counts)  ]  

# 1.2.3. Intersect gene lists
exprs_genes <- intersect(exprs_genes.tpm,
                         exprs_genes.counts) 

# Exclude chrY genes in female-only tissues
if(sex == "female"){
  exprs_genes <- exprs_genes[!exprs_genes %in% Y_genes]
}

# 2.1 Variables ----

covariates<- colnames(metadata)[colnames(metadata) %in% c("HardyScale","IschemicTime","RIN","ExonicRate","Cohort", "NucIsoBatch", "PEER1", "PEER2")]

#Subset by window
window.center = seq(30,60,1)
buckets.size = 10
interval.size = 0

# window.center = seq(0,1,0.05)
# buckets.size = 0.1
# interval.size = 0
# Count samples in each window
sample_counts <- sapply(window.center, function(center) {
  sum(metadata$Age >= (center - buckets.size - interval.size) & 
        (metadata$Age <= (center + buckets.size + interval.size)) &
        !(metadata$Age > (center - interval.size) & 
            metadata$Age < (center + interval.size)))  
})


# Convert to DataFrame for plotting
df_counts <- data.frame(Window = window.center, Samples = sample_counts)

plot<-ggplot(df_counts, aes(x = factor(Window), y = Samples)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = Samples), vjust = -0.5, color = "black") +
  labs(title = paste0("Samples per Test Window - ", tissue, " - ", sex), x = "Window Center (Age)", y = "Number of Samples") +
  theme_minimal()

ggsave( paste0(outpath, "/sliding_window_DEA_", tissue, "_", sex, "_", buckets.size, "_", interval.size, "_sampling.pdf"), plot = plot, width = 10, height = 6)

#####################################
# 4. Differential expression analysis ####
print("# ---- Running differential expression analysis ---- #")

# 4.1 limma fit : expression ~ covariates + traits ----
my_data <- list()
individual_traits <- c("Age_bin","Ancestry", "BMI")

# Define function for each window
dea_res<-list()
process_window <- function(k) {
  # Subset by window
  qt.tmp  <-  rep(NA, length(metadata$Age))
  # qt.tmp[which(metadata$Age < window.center[k] & metadata$Age >= (window.center[k] - buckets.size))] <- 0
  # qt.tmp[which(metadata$Age > window.center[k] & metadata$Age <= (window.center[k] + buckets.size))] <- 1
  qt.tmp[which(metadata$Age < (window.center[k]-interval.size) & metadata$Age >= (window.center[k] - buckets.size - interval.size))] <- 0
  qt.tmp[which(metadata$Age > (window.center[k]+ interval.size) & metadata$Age <= (window.center[k] + buckets.size + interval.size))] <- 1
  qt.tmp <- factor(qt.tmp)
  metadata$Age_bin<-qt.tmp
  metadata_filtered <- metadata[!is.na(metadata$Age_bin), ]
  metadata_filtered$HardyScale<-as.factor(as.numeric(metadata_filtered$HardyScale))
  
  
  if (length(unique(metadata_filtered$Age_bin)) < 2) {
    message(paste("Skipping window.center:", window.center[k], "- Not enough levels in Age_bin"))
    return(NULL)
  }else if(length(unique(metadata_filtered$HardyScale)) <2){
    covariates<- covariates[!covariates == "HardyScale"]
    
  }
  # 2.2 Create DGEList object ----
  dge <- DGEList(counts[exprs_genes,metadata_filtered$Sample])
  
  # 2.3 Calculate normalization factors (does not do the normalization, only computes the factors) ----
  dge <- calcNormFactors(dge)
  
  # 2.4 Voom ----
  v <- voom(dge, design = NULL, normalize="quantile", save.plot=F, plot = F) # samples are treated as replicates
  
  print(paste0("# ---- Processing window.center: ", window.center[k], " ---- #"))
  
  # Fit model
  fml_args_mod <- paste(c(individual_traits, covariates), collapse = " + ")
  mod <- model.matrix( as.formula(paste(" ~  ", paste(fml_args_mod,collapse = " "))), data =  metadata_filtered)
  model<-as.formula(paste(" ~  ", paste0(fml_args_mod,collapse = " ")))
  print(paste0("# ----Model: ",model, " ---- #" ))
  
  fit <- lmFit(v, mod)
  # Run Limma test
  k_res <- limma_lm(fit, "Age_bin", metadata_filtered, interaction = NULL)
  
  # Compute average TPM, median, variance
  avrg_TPM <- apply(tpm, 1, function(x) mean(log2(x+1)))
  median_TPM <- apply(tpm, 1, function(x) median(log2(x+1)))
  var_TPM <- apply(tpm, 1, function(x) var(x))

  if ((!is.null(k_res) || length(k_res) > 0 ) && !is.null(names(k_res))) {
    k_res$AvgTPM <- avrg_TPM[rownames(k_res)]
    k_res$MedianTPM <- median_TPM[rownames(k_res)]
    k_res$VarTPM <- var_TPM[rownames(k_res)]
    
    gene_names <- sapply(rownames(k_res), function(gene) {
      gene_annotation[gene_annotation$ensembl.id == gene, "gene.name.x"]
    })
    names(gene_names) <- NULL
    k_res[["gene.name.x"]] <- gene_names
    
    return(k_res)  # Return only k_res, not a nested list
    
  } else {
    message(paste("Skipping window.center:", window.center[k], "- No valid DE results"))
    return(NULL)
  }

}
# Run in parallel
ncores <- detectCores() - 1  # Use all available cores except one
results_list <- mclapply(seq_along(window.center), process_window, mc.cores = ncores)
r<-results_list
saveRDS(results_list, paste0(outpath,"/sliding_window_DEA_", tissue, "_", sex,"_", buckets.size,"_", interval.size, "_results.rds"))

tissue<-"Vagina"
outpath<- paste0("/X//Laura/sliding_window_DEA/All_tissues/", tissue, "/")

r<- readRDS(paste0(outpath,"/sliding_window_DEA_", tissue, "_", sex,"_", buckets.size,"_", interval.size, "_results.rds"))

#Now plot the results
pvals<-c("P.Value", "adj.P.Val")

#We choose depending on the results
pval<-pvals[2]

####1. Extract the significant genes in at least one comparison

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
nrow(final_significant_genes)
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
count_significant <- function(p_df, thresholds = c(0.1, 0.05, 0.01, 0.001)) {
  return(sapply(thresholds, function(thresh) colSums(p_df < thresh, na.rm = TRUE)))
}

res.signif <- count_significant(p_value_df)
rownames(res.signif) <- window.center
colnames(res.signif) <- c(0.1, 0.05, 0.01, 0.001)

toPlot <- t(res.signif)
x <- as.numeric(colnames(toPlot))

### 4. Generate Plot
pdf(paste0(input, "/Laura/sliding_window_DEA/All_tissues/", tissue, "/curve_", tissue, "_", sex, "_", pval, "_", buckets.size, "_", interval.size, ".pdf"), width = 7, height = 7)

plot(1, type = "n", xlim = range(x, na.rm = TRUE), ylim = range(toPlot, na.rm = TRUE), 
     ylab = "# significant", xlab = "Age windows")

for (i in seq_len(nrow(toPlot))) {
  lines(x, toPlot[i, ], type = 'l', lwd = i)
}

legend("topleft", legend = paste("p<", rownames(toPlot), sep = ""), lwd = c(1, 2, 3))

dev.off()

### Plot
threshold_interest<- "0.05"
thresholds<- c("0.05")

ovary_counts<-df_counts
uterus_counts<-df_counts
vagina_counts<-df_counts

ovary_counts$Tissue<-"Ovary"
uterus_counts$Tissue<-"Uterus"
vagina_counts$Tissue<-"Vagina"

t_t<-rbind(ovary_counts, uterus_counts, vagina_counts)


p <- ggplot(uterus_counts, aes(x = Window, y = Genes, group = Tissue, color = Tissue)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.5, size = 1) + 
  scale_color_manual(values = pal) +  # Adjust colors as needed
  scale_size(range = c(2, 6)) +
  labs(title = "", 
       x = "Window Center (Age)", 
       y = paste0("#DEGs (adj p-val<", threshold_interest, ")"), 
       color = "Tissue", size = "Number of Samples") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 13, color = "black"), 
        axis.text.x = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13, color = "black"), 
        plot.title = element_text(size = 13, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        panel.grid = element_blank(),                    # Remove all grid lines
        # axis.line.x = element_line(size = 0.5, color = "grey80"), # Add x-axis line
        # axis.line.y = element_line(size = 0.5, color = "grey80"),
        panel.border = element_rect(color = "grey80", size = 0.5, fill = NA) # Add full border
  )
p
pdf(paste0(input, "/Laura/sliding_window_DEA/All_tissues/curves_uterus_ovary_vagina_", pval, "_", threshold_interest, "_", buckets.size,"_", interval.size, "_normalized.pdf"), width = 4.5, height = 4.3)
print(p)
dev.off()

if(!getwd()== "/X/Laura/sliding_window_DEA"){

  # list_plots_adj<-list()
  tissue<- "Uterus"
    if(tissue %in% tissues){
      s<- "female"
    }else if (tissue %in% c("Prostate", "Testis")){
      s<-c("male")
    }else{
      s<-c("female", "male")
    }
  s<-"female"
  
  r<-readRDS(paste0("/X//Laura/sliding_window_DEA/All_tissues/", tissue, "/sliding_window_DEA_", tissue, "_", s, "_results.rds"))
  tissue<- "Uterus"
  
  buckets.size<- 20
  interval.size<-0
  r<-readRDS(paste0("/X//Laura/sliding_window_DEA/All_tissues/", tissue,"/sliding_window_DEA_", tissue, "_", s,"_", buckets.size,"_", interval.size, "_results.rds"))
  
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
  nrow(final_significant_genes)
  # saveRDS(final_significant_genes, paste0("/X//Laura/sliding_window_DEA/All_tissues/", tissue, "/GENES_sliding_window_DEA_", tissue, "_", s,"_", buckets.size,"_", interval.size, "_results.rds" ))

  nrow(readRDS(paste0("/X//Laura/sliding_window_DEA/All_tissues/", tissue, "/GENES_sliding_window_DEA_", tissue, "_", s,"_", buckets.size,"_", interval.size, "_results.rds" )))
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
  if(tissue == "MuscleSkeletal"){
    metadata_tot <- readRDS(paste0(input, "Laura/00.Data/1.Final_tissues_to_use/", tissue, "/metadata.rds"))
    counts_tot <- readRDS(paste0(input, "Laura/00.Data/1.Final_tissues_to_use/", tissue, "/counts.rds"))
    tpm_tot <- readRDS(paste0(input, "Laura/00.Data/1.Final_tissues_to_use/", tissue, "/tpm.rds"))
  }else{
    metadata_tot <- readRDS(paste0(input, "Laura/00.Data/v10/", tissue, "/metadata.rds"))
    counts_tot <- readRDS(paste0(input, "Laura/00.Data/v10/", tissue, "/counts.rds"))
    tpm_tot <- readRDS(paste0(input, "Laura/00.Data/v10/", tissue, "/tpm.rds"))
    
  }
  
  metadata_tot$Ancestry<- NULL
  ancestry_file<- read.table(paste0(input, "/Laura/00.Data/admixture_inferred_ancestry.txt"))
  ances<-ancestry_file[,c(1,3)]
  colnames(ances)<- c("Donor", "Ancestry")
  metadata_tot<-merge(metadata_tot, ances, by= "Donor")

  
  
  if (tissue %in% tissues_with_filter){
    if(tissue %in% both_tissues){
      if(sex == "male"){
        filter<-readRDS(paste0(input, "/Laura/03.Image_processing/Second_filtering_images/",tissue, "_male_filtered_images.rds"))
        metadata<- metadata_tot[metadata_tot$Donor %in% filter$Subject.ID,]
        # metadata$Ancestry <- relevel(metadata$Ancestry, ref="EUR")
        
      }else{
        filter<-readRDS(paste0(input, "/Laura/03.Image_processing/Second_filtering_images/",tissue, "_female_filtered_images.rds"))
        metadata<- metadata_tot[metadata_tot$Donor %in% filter$Subject.ID,]
        # metadata$Ancestry <- relevel(metadata$Ancestry, ref="EUR")
      }
    }else if (tissue %in% tissues){
      filter<-readRDS(paste0(input, "/Laura/03.Image_processing/Second_filtering_images/",tissue, "_final_filtered_images.rds"))
      metadata<- metadata_tot[(metadata_tot$Donor %in% filter$Subject.ID),]   
    }

  }else{
    if (tissue %in% c("Prostate", "Testis")){
      metadata<-metadata_tot
    }else {
      if (sex=="male"){
        metadata<-metadata_tot[metadata_tot$Sex == "1",]
        
      }else
        metadata<-metadata_tot[metadata_tot$Sex == "2",]
    }
  }

    z.score <- function(ft.matrix){
      t(scale(t(ft.matrix), center=TRUE, scale=TRUE))
    }
    
    counts_tot <- readRDS(paste0("~//X//Laura/00.Data/v10/", tissue, "/counts.rds"))
    tpm_tot <- readRDS(paste0("~//X//Laura/00.Data/v10/", tissue, "/tpm.rds"))
    counts <- counts_tot[, metadata$Sample]
    tpm<-tpm_tot[, metadata$Sample]
    
    ######RESIDUALS
    library(dplyr)
    residuals<-readRDS(paste0("~//X//Laura/12.Tissues_substructures/DEA/DEA_v10/",tissue,"/no_subs_no_age/", tissue,"_no_subs_no_age_residuals.results.rds"))
    residuals <- as.data.frame(residuals)
    residuals[1:5, 1:5]
    residuals <- residuals %>%
      mutate(`GTEX-R55G-1626` = rowMeans(dplyr::select(., `GTEX-R55G-1626`, `GTEX-R55G-1626.1`), na.rm = TRUE)) %>%
      select(-`GTEX-R55G-1626.1`)
    
    residuals<- as.data.frame(t(residuals))
    residuals$Age<- metadata$Age
    final_significant_genes<-final_significant_genes[rownames(final_significant_genes) %in% colnames(residuals), , drop = FALSE]
    
    selected_cols <- intersect(c(rownames(final_significant_genes), "Age"), colnames(residuals))
    genes_selected <- residuals[, selected_cols, drop = FALSE]
    colnames(genes_selected)<-c(final_significant_genes$gene.name.x, "Age")
    rownames(genes_selected)<-as.factor(gsub("^(\\w+-\\w+).*", "\\1", rownames(genes_selected)))
    
    
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
    # 
     ncol(genes_selected)
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
    max(z_scores_matrix)
    min(z_scores_matrix)
    
    min<--5
    max<-5

    pairs.breaks <- seq(min, max, by = 0.01)
    
    # Generate the color palette
    library(gplots)
    mycol <- gplots::colorpanel(n = length(pairs.breaks) - 1, low = "cyan", mid = "black", high = "yellow")
    # Plot the heatmap
    pdf(paste0("/X//Laura/sliding_window_DEA/All_tissues/",tissue,"/heatmap_", tissue, "_", sex, "_", buckets.size, "_", interval.size, "_adj_z-score.pdf"), width = 7, height = 6)
    
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
    
    pdf(paste0(input, "/Laura/sliding_window_DEA/All_tissues/", tissue, "/", buckets.size, "_",interval.size, "_adj_zscore_legend_plot.pdf"), width = 3, height = 2)  # Save as PNG (adjust size)
    library(gplots)
    
    mycol <- colorpanel(n = 100, low = "cyan", mid = "black", high = "yellow")
    breaks <- seq(min, max, length.out = 100)
    labels <- c(min, 0, max)  # Only these three labels
    plot(1, type = "n", xlab = "", ylab = "", xlim = c(-6, 6), ylim = c(0, 1), axes = FALSE)
    image(x = breaks, y = 0, z = matrix(breaks, ncol = 1), col = mycol, axes = FALSE, xlab = "", ylab = "")
    box()
    axis(side = 1, at = labels, labels = labels, las = 1)
    title(main = "z-score", line = 1)
    
    dev.off()
    
    
    
    #4.2. With p values
    
    heatmap_data <- list()
    names(r)<- window.center
    
    # Loop through the list of window centers
    for (win_idx in names(r)) {
      window_data <- r[[win_idx]]  # Extract the data for this window
      if (!is.null(window_data) && nrow(window_data) > 0) {
        # Calculate the Z-scores (you can skip if already calculated)
        window_data[[pval]] <- window_data[[pval]]  # Assuming Z-score is the same as 't'
        
        # Multiply the sign of logFC by Z-scores
        window_data$heatmap_value <- sign(window_data$logFC) * -log10(window_data[[pval]])
        
        # Store the data for this window
        heatmap_data[[win_idx]] <- data.frame(
          gene = window_data$gene.name.x,  # Use gene names
          heatmap_value = window_data$heatmap_value
        )
      }
    }
    
    

    # Filter out NULL or invalid entries
    valid_heatmap_data <- heatmap_data[!sapply(heatmap_data, is.null)]
    
    # Convert list to a single data frame with an additional 'window.center' column
    all_windows <- bind_rows(lapply(names(valid_heatmap_data), function(i) {
      cbind(window.center = i, valid_heatmap_data[[i]])
    }))
    all_windows <- all_windows %>%
      group_by(gene, window.center) %>%
      summarise(heatmap_value = mean(heatmap_value, na.rm = TRUE), .groups = "drop")
    
    all_windows<-all_windows[all_windows$gene %in% final_significant_genes$gene.name.x,]
    max(all_windows$heatmap_value)
    min(all_windows$heatmap_value)
    min<-min(all_windows$heatmap_value)
    max<-max(all_windows$heatmap_value)

    # Spread the data into a wide format
    heatmap_matrix <- all_windows %>%
      pivot_wider(names_from = window.center, values_from = heatmap_value) %>%
      column_to_rownames("gene")
    heatmap_matrix<-heatmap_matrix[order(rownames(z_scores_matrix)),]
    pairs.breaks <- seq( min,max, by = 0.01)  # Adjust range for Z-scores
    mycol <- gplots::colorpanel(length(pairs.breaks) - 1, low = "cyan", mid = "black", high = "yellow")
    
    library(gplots)
    
    # Plot the heatmap
    selected_ages <- c(20,30, 40,50, 60, 70)
    
    # Replace all other column names with empty strings
    colnames(heatmap_matrix) <- ifelse(as.numeric(colnames(heatmap_matrix)) %in% selected_ages, 
                                        colnames(heatmap_matrix), "")
    
    
    pdf(paste0(input, "/Laura/sliding_window_DEA/All_tissues/", tissue, "/heatmap_", tissue, "_", sex, "_", buckets.size, "_", interval.size, "_adj_pval.pdf"), width = 30, height = 20)  # Adjust width and height as needed
    par(mar = c(5, 5, 2, 2))  # Adjust margins (bottom, left, top, right)
    heatmap.2(as.matrix(heatmap_matrix),
              cexRow = 0.7, cexCol = 3,  # Adjust font sizes
              trace = "none",  # No trace lines
               dendrogram = "none",  # Cluster rows (genes)
              breaks = pairs.breaks, 
              col = mycol, 
              Rowv = FALSE,  # Cluster rows
              Colv = FALSE, # Do not cluster columns (ages stay in order)
              lhei = c(0.2, 10), 
              lwid = c(0.2, 3),
              key = TRUE,  # Display legend
              labRow = NA,  # Removes row names
              notecex = 0,  # Ensure no annotations appear
              sepcolor = NA,  # Disable border color
              sepwidth = c(0, 0),
              srtCol = 0) 
    dev.off()

    pdf(paste0(input, "/Laura/sliding_window_DEA/All_tissues/", tissue, "/", buckets.size, "_", interval.size, "_adj_pval_legend_plot.pdf"), width = 3, height = 2)  # Save as PNG (adjust size)
    library(gplots)
    
    mycol <- colorpanel(n = 100, low = "cyan", mid = "black", high = "yellow")
    breaks <- seq(min, max, length.out = 100)
    labels <- c(min, 0, max)  # Only these three labels
    plot(1, type = "n", xlab = "", ylab = "", xlim = c(-6, 6), ylim = c(0, 1), axes = FALSE)
    image(x = breaks, y = 0, z = matrix(breaks, ncol = 1), col = mycol, axes = FALSE, xlab = "", ylab = "")
    box()
    axis(side = 1, at = labels, labels = labels, las = 1)
    title(main = "adj p val", line = 1)
    
    dev.off()
    
    
}