#DEA SUBSTRUCTURES

Sys.setenv(TZ="Europe/Madrid")
# ---------------------- #
start_time <- Sys.time()

# Libraries ####
suppressMessages(library(edgeR))
suppressMessages(library(limma))
suppressMessages(library(parallel))
suppressMessages(library(car))

# Command line arguments ####
args <- commandArgs(trailingOnly=TRUE)

# Functions ####
source("./DEA_and_DSA.R_functions.R")

gene_annotation  <- read.csv("../00.Data/gencode.v39.annotation.bed")
Y_genes <- gene_annotation[gene_annotation$chr=="chrY",]$ensembl.id

# Tissue ----
tissues<-c("Uterus", "Ovary", "Vagina", "BreastMammaryTissue", "CervixEndocervix", "CervixEctocervix", "FallopianTube")
possible_analysis<- c("ALL", "SEPARATING_SEX")

#There are 4 possible analyses:
# 1. Include all the covariates (structures) in the model and test the interaction with age of one of interest
# 2. Include one structure of interest in the model and test its interaction with age
# 3. Include all the structures in the model and do not test the interaction with age
# 4. Include one structure of interest in the model and do not test its interaction with age

tests<-c("all_cov_inter", "one_cov_inter", "all_cov_no_inter", "one_cov_no_inter")
sex<- "female"
analyses<-c("age", "sub", "inter")
both_tissues<-c("BreastMammaryTissue" ,"ArteryAorta"   ,"ArteryTibial","ArteryCoronary"  )
analysis<- possible_analysis[2]
all_tissues_results<-list()
all_metrics_results<-list()
total_results<-list()
total_metrics<-list()
tissue<-"Uterus"
tissues<-c("CervixEndocervix", "CervixEctocervix", "FallopianTube")

for (tissue in tissues){
  tests<-"all_cov_no_inter"
  test<-tests
  test<-"_subs_no_age"
for (test in tests){
  outpath <-paste0("./", tissue, "/", test, "/")
  print(paste0("---------------------Analyzing ", tissue, "---------------------"))
  print(paste0("## ", test, " ##"))
  if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}
  
    # 1.1 Read in data ----
    #FOR GTEx v10
    counts_tot <- readRDS(paste0("../00.Data/", tissue, "/SIMULATED_counts.rds"))
    tpm_tot <- readRDS(paste0("../00.Data/", tissue, "/SIMULATED_tpm.rds"))
    metadata_tot <- readRDS(paste0("../00.Data/", tissue, "/SIMULATED_metadata.rds"))
    metadata_tot$Ancestry<- NULL
    ancestry_file<- read.table("../00.Data/SIMULATED_admixture_inferred_ancestry.txt")
    metadata_tot<-merge(metadata_tot, ancestry_file, by= "Donor")

    sex <- "female"
    filter<-readRDS(paste0("../00.Data/",tissue, "_final_filtered_images.rds"))
    metadata<- metadata_tot[(metadata_tot$Donor %in% filter$Subject.ID),]
  
    #Read substructures' proportions
    #Example data provided only for uterus!!
    if(tissue =="BreastMammaryTissue"){
      cell_prop<-read.csv(paste0("~/X/derived_proportions_Craig/", tissue, "/", tissue, "_pivot.csv"))
      substructures<- c("adipocyte", "lobule", "duct", "stroma", "nerve", "gynecomastoid_hyperplasia")
      
    }else if(tissue == "Ovary"){
      cell_prop<-read.csv(paste0("~/X/12.Tissues_substructures/", tissue, "_no_follicles_pivot.csv"))
      substructures<- c("cortex", "corpora", "medulla", "vessels")
      
    }else if (tissue =="Vagina"){
      cell_prop<-read.csv(paste0("~/X/12.Tissues_substructures/", tissue, "_pivot.csv"))
      substructures<- c("epithelium", "lamina_propria", "vessels", "stroma")
      
    }else if (tissue =="Uterus"){
      cell_prop<-read.csv(paste0("../00.Data/", tissue, "_pivot.csv"))
      substructures<- c("myometrium", "endometrium", "vessels")
      
    }else if (tissue =="CervixEndocervix"){
      cell_prop<-read.csv(paste0("~/X/12.Tissues_substructures/Endocervix_pivot.csv"))
      substructures<- c("vessels", "glandular_epithelium", "stroma")
      
    }else if (tissue =="CervixEctocervix"){
      cell_prop<-read.csv(paste0("~/X/12.Tissues_substructures/Ectocervix_pivot.csv"))
      substructures<- c("epithelium", "stroma", "vessels", "glands")

    }else if (tissue == "FallopianTube"){
      cell_prop<-read.csv(paste0("~/X/12.Tissues_substructures/FallopianTube_pivot.csv"))
      substructures<- c("vessels", "lumen", "smooth_muscle", "epithelium", "stroma") ##always control for adipose, since it is contamination
      
    }
    
    cell_prop$Donor<-  gsub("^(\\w+-\\w+).*", "\\1", cell_prop$slide_id)
    metadata<- merge(metadata, cell_prop, by="Donor")

    counts <- counts_tot[, metadata$Sample]
    tpm<-tpm_tot[, metadata$Sample]

    # 1.2 Genes expressed per tissue ----
    # 1.2.1 TPM>=1 in at least 20% of the tissue samples
    exprs_genes.tpm <- rownames(tpm)[apply(tpm, 1, function(x) sum(x>=1) ) >= 0.2*ncol(tpm)  ]  
    
    # 1.2.2 Count >=10 in at least 20% of the tissue samples
    exprs_genes.counts <- rownames(counts)[ apply(counts, 1, function(x) sum(x>=10) ) >= 0.2*ncol(counts)  ]  
    
    # 1.2.3. Intersect gene lists
    exprs_genes <- intersect(exprs_genes.tpm,
                             exprs_genes.counts) 
    
    # Exclude chrY genes in female-only tissues
    exprs_genes <- exprs_genes[!exprs_genes %in% Y_genes]

    # 2 Variables ----

    if (tissue %in% c("Uterus", "Ovary", "Vagina", "BreastMammaryTissue")){
      peers<- c("PEER1", "PEER2")
    }else{
      peers<-c()
    }
    if (tissue =="FallopianTube"){
      tec<-c(colnames(metadata)[colnames(metadata) %in% c("HardyScale","IschemicTime","RIN","ExonicRate","Cohort", "NucIsoBatch")], "adipose")
    }else{
      tec<- colnames(metadata)[colnames(metadata) %in% c("HardyScale","IschemicTime","RIN","ExonicRate","Cohort", "NucIsoBatch")]
      
    }
    covariates <- c(substructures, peers, tec)
    
    # 2.2 Create DGEList object ----
    dge <- DGEList(counts[exprs_genes,])
    
    # 2.3 Calculate normalization factors (does not do the normalization, only computes the factors) ----
    dge <- calcNormFactors(dge)
    
    # 2.4 Voom ----
    v <- voom(dge, design = NULL, normalize="quantile", save.plot=F, plot = F) # samples are treated as replicates
  
  
    # 4. Differential expression analysis ####
    print("# ---- Running differential expression analysis ---- #")
    
    # 4.1 limma fit : expression ~ covariates + traits ----
    my_data <- list()
    resu<- list()
    summary_results <- list()
    for (sub in substructures){
       individual_traits <- c("Age","Ancestry", "BMI")
        if (test %in% c("all_cov_inter", "one_cov_inter")){
           interaction_terms<-paste(c("Age", sub), collapse= ":")
        }else{
          interaction_terms<-c()
        }
        if(test %in% c("all_cov_inter", "all_cov_no_inter")){
          covariates <- c(substructures, peers, tec)
        }else{
          covariates <- c(sub, peers, tec)
        }
      
      fml_args_mod <- paste(c(covariates, individual_traits, interaction_terms), collapse = " + ")
      mod <- model.matrix( as.formula(paste(" ~  ", paste(fml_args_mod,collapse = " "))), data =  metadata)
      model<-as.formula(paste(" ~  ", paste0(fml_args_mod,collapse = " ")))
      print(paste0("# ----Model: ",model, " ---- #" ))
      
      # 4.2 Limma fit ----
      fit <- lmFit(v, mod)
      #If we want to save residuals:
      
      # exprs_residuals <- resid(fit, v)
      # saveRDS(exprs_residuals, paste0(outpath,tissue,"_", test, "_residuals.results.rds"))


      # Add objects to data
      my_data[["dge"]] <- dge
      my_data[["v"]] <- v
      my_data[["fit"]] <- fit

      # 4.3 Limma test with interaction term ----
      dea_res <- list()
      
      if (test =="all_cov_inter"){

       trait_res <- lapply(c(individual_traits, substructures), function(phenotype) limma_lm(fit, phenotype, metadata))
       names(trait_res) <- c(individual_traits, substructures)
       dea_res <- c(dea_res, trait_res)
       interaction_res <- limma_lm(fit, covariate = NULL, covariate_data = metadata, interaction = paste0(sub, ":Age"))
       name_inter<-paste0("Interaction_Age_", sub)
       dea_res[[name_inter]] <- interaction_res
       new_traits<- c(individual_traits, substructures, name_inter)


      }else if(test=="one_cov_inter"){
        trait_res <- lapply(c(individual_traits, sub), function(phenotype) limma_lm(fit, phenotype, metadata))
        names(trait_res) <- c(individual_traits, sub)
        dea_res <- c(dea_res, trait_res)
        interaction_res <- limma_lm(fit, covariate = NULL, covariate_data = metadata, interaction = paste0(sub, ":Age"))
        name_inter<-paste0("Interaction_Age_", sub)
        dea_res[[name_inter]] <- interaction_res
        new_traits<- c(individual_traits, sub, name_inter)

      }else if(test =="all_cov_no_inter"){
        trait_res <- lapply(c(individual_traits, substructures), function(phenotype) limma_lm(fit, phenotype, metadata))
        names(trait_res) <- c(individual_traits, substructures)
        dea_res <- c(dea_res, trait_res)
        new_traits<- c(individual_traits, substructures)


      }else if(test =="one_cov_no_inter"){
        trait_res <- lapply(c(individual_traits, sub), function(phenotype) limma_lm(fit, phenotype, metadata))
        names(trait_res) <- c(individual_traits, sub)
        dea_res <- c(dea_res, trait_res)
        new_traits<- c(individual_traits, sub)

      }


    # 5. Computing avrg TPM and exprs var ####
      print("# ---- Calculating avrg TPM and var ---- #")

      # 5.1 Compute average TPM expression and variance for each event
      avrg_TPM <- apply(tpm, 1, function(x) mean(log2(x+1)))
      median_TPM <- apply(tpm, 1, function(x) median(log2(x+1)))
      var_TPM <- apply(tpm, 1, function(x) var(x))


      # Add AvgExprs & ExprsVar and order data.frame
      for(trait in new_traits){
        print(trait)
        if(length(dea_res[[trait]])>1){
          # Add average expression TPM and variance
          dea_res[[trait]]$AvgTPM <- sapply(rownames(dea_res[[trait]]), function(gene) avrg_TPM[gene])
          dea_res[[trait]]$MedianTPM <- sapply(rownames(dea_res[[trait]]), function(gene) median_TPM[gene])
          dea_res[[trait]]$VarTPM <- sapply(rownames(dea_res[[trait]]), function(gene) var_TPM[gene])
          gene_names <- sapply(rownames(dea_res[[trait]]), function(gene) gene_annotation[gene_annotation$ensembl.id==gene, "gene.name.x"])
          names(gene_names) <- NULL
          dea_res[[trait]][["gene.name.x"]] <- gene_names
        }
      }

      if (test =="all_cov_no_inter"){
        saveRDS(dea_res,
                paste0(outpath,tissue,"_", test, "_AGE_covariates_and_traits.results.rds"))
        print(paste0("# ---- saved file: ",tissue,"_", test, "_AGE_covariates_and_traits.results.rds", " ---- #") )
      }else{
        saveRDS(dea_res,
          paste0(outpath,tissue,"_", sub, "_", test, "_covariates_and_traits.results.rds"))

        print(paste0("# ---- saved file: ",tissue,"_", sub, "_",test, "_covariates_and_traits.results.rds", " ---- #") )

      }
      }
}
}
   
