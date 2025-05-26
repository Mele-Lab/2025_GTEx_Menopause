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

  tissue<-"Uterus"
  sex<- "female"
  outpath<- paste0("./")
  
# Functions ####
source(paste0("DEA_and_DSA.R_functions.R"))

#FOR GTEX V10
gene_annotation  <- read.csv("../00.Data/gencode.v39.annotation.bed")
Y_genes <- gene_annotation[gene_annotation$chr=="chrY",]$ensembl.id

# Tissue ----
tissues<-c("Uterus", "Ovary", "Vagina", "BreastMammaryTissue", "CervixEndocervix", "CervixEctocervix", "FallopianTube")

# Read in data ----
#FOR GTEx v10
  metadata_tot <- readRDS(paste0("../00.Data/", tissue, "/SIMULATED_metadata.rds"))
  counts_tot <- readRDS(paste0("../00.Data/", tissue, "/SIMULATED_counts.rds"))
  tpm_tot <- readRDS(paste0("../00.Data/", tissue, "/SIMULATED_tpm.rds"))
  

metadata_tot$Ancestry<- NULL
ancestry_file<- read.table(paste0("../00.Data/SIMULATED_admixture_inferred_ancestry.txt"))
metadata_tot<-merge(metadata_tot, ancestry_file, by= "Donor")

filter<-readRDS(paste0("../00.Data/",tissue, "_final_filtered_images.rds"))
metadata<- metadata_tot[(metadata_tot$Donor %in% filter$Subject.ID),]   


counts <- counts_tot[, metadata$Sample]
tpm<-tpm_tot[, metadata$Sample]


# Gene annotation ----

# Genes expressed per tissue ----
#  TPM>=0.1 in at least 20% of the tissue samples
exprs_genes.tpm <- rownames(tpm)[apply(tpm, 1, function(x) sum(x>=1) ) >= 0.2*ncol(tpm)  ]  

# Count >=6 in at least 20% of the tissue samples
exprs_genes.counts <- rownames(counts)[ apply(counts, 1, function(x) sum(x>=10) ) >= 0.2*ncol(counts)  ]  

# Intersect gene lists
exprs_genes <- intersect(exprs_genes.tpm,
                         exprs_genes.counts) 

# Exclude chrY genes in female-only tissues
if(sex == "female"){
  exprs_genes <- exprs_genes[!exprs_genes %in% Y_genes]
}

# Variables ----

covariates<- colnames(metadata)[colnames(metadata) %in% c("HardyScale","IschemicTime","RIN","ExonicRate","Cohort", "NucIsoBatch", "PEER1", "PEER2")]

#Subset by window
window.center = seq(21,70,1)
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



# Differential expression analysis ####
print("# ---- Running differential expression analysis ---- #")

# limma fit : expression ~ covariates + traits ----
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

saveRDS(results_list, paste0(outpath,"/sliding_window_DEA_", tissue, "_", sex,"_", buckets.size,"_", interval.size, "_results.rds"))

