suppressMessages({
  library(data.table)
  library(stringr)
  library(VennDiagram)
library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(reshape2)
library(patchwork)
library(tibble)
library(pheatmap)
library(gridExtra)
})


# Define input path
input_path <- "/x/"


##GWASes studies included

##Reproductive aging traits
### 1. ANM_ruth_tables are the supplementary tables of PMID: 34349265

# Read and process the ANM table
temp_gwas_anm <- read_excel(file.path("../00.Data/GWAS_datasets/ANM_ruth_tables.xlsx"),
                            sheet = "ST2", skip = 5) %>%
  separate_rows(Consensus, sep = " / ") %>%
  distinct(Consensus, .keep_all = TRUE) %>%
  mutate(Trait = "ANM", Genes = Consensus)

gwas_anm <- dplyr::select(temp_gwas_anm, Trait, Genes)


##2. AAM GWAS 

AAM_EUR_GWAS <- read_excel("../00.Data/GWAS_datasets/41588_2024_1798_MOESM4_ESM.xlsx", sheet = 19)
gwas_AAM_EUR <- data.frame(                                                                                                              
  Trait = "AAM",
  Genes = unique(AAM_EUR_GWAS$gene),
  stringsAsFactors = FALSE
)
gwas_AAM_EUR <- distinct(gwas_AAM_EUR, Genes, .keep_all = TRUE)

##Clinical traits

##3. Ovarian cysts
fem_repro <- fread(file.path("../00.Data/GWAS_datasets/femrepro_natalia_ST8.tsv"), sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE, skip = 1)

cysts <- fem_repro %>%
  filter(grepl("Ovarian cyst", `Phenotype name`))

cysts_long <- cysts %>%
  pivot_longer(cols = c(`Top 3 Candidate Genes (OpenTargetsGenetics v8)`, V5, V6), 
               names_to = "Source", values_to = "Genes") %>%
  mutate(
    Trait = "Ovarian cysts",
    Genes = trimws(Genes)  # Removes leading/trailing whitespace
  ) %>%  
  select(Trait, Genes) %>%
  filter(!is.na(Genes), Genes != "")

cysts<-cysts_long
head(cysts)

##4. Polyps of female genital tract PMID: 39986329

polyps_FGT <- data.frame(
  Trait = "Female genital tract polyps",
  Genes = c("NFIA","EXO1","EEFSEC","ZBTB38","IGF2BP2","ARHGAP26","TERT-CLPTM1L","TRPS1","POLE3","NOC3L","PLCE1","PSMD13","ODF3","PRIM1","ACTL9","TTC28","CHEK2"),
  stringsAsFactors = FALSE
)


##5. Polycystic ovary syndrome (consortium work, preprint 2024)
preprint_pcos <- read_excel(file.path("../00.Data/GWAS_datasets/supptables_Loes_PCOSpreprint_media-2.xlsx"),
                            sheet = "ST2 Identified sig. in the age ", skip = 3) %>%
  separate_rows(`Consensus Gene (based on biological and literature evidence)`, sep = " / ") %>%
  distinct(`Consensus Gene (based on biological and literature evidence)`, .keep_all = TRUE) %>%
  mutate(Trait = "PCOS (preprint)", Genes = `Consensus Gene (based on biological and literature evidence)`)

preprint_pcos<-dplyr::select(preprint_pcos,Trait, Genes)

length(preprint_pcos$Genes)

##6. Endometriosis PMID: 36914876

endometriosis_paper_genes <- data.frame(
  Genes = c("WNT4", "GREB1", "ETAA1", "KDR", "ID4", "SYNE1", "CDKN2B-AS1", "FSHB", "VEZT", 
            "NGF", "SLC19A2", "DNM3", "BMPR2", "BSN", "PDLIM5", "EBF1", "CD109", "HEY2", 
            "FAM120B", "HOXA10", "KCTD9", "GDAP1", "VPS13B", "ASTN2", "ABO", "MLLT10", 
            "RNLS", "WT1", "PTPRO", "HOXC10", "IGF1", "DLEU1", "RIN3", "SRP14-AS1", 
            "SKAP1", "CEP112", "ACTL9", "TEX11", "FRMD7", "LINC00629"),
  Trait = "endometriosis_paper_genes"
)

##7. Pelvic organ prolapse PMID: 35739095

prolapse_genes_data <- data.frame(
  Trait = "all_genes_prolapse", # Assign "all_genes_prolapse" to all rows
  Genes = c(
    "CDC42", "WNT4","GDF7", "LDAH", "RP11-130L8.1", "EFEMP1", "FAT4", "IMPDH1", "TBX3", "TBX5", "TBX5-AS1", "HNRNPA1P48", "SALL1",
    "HOXD13", "HOXD3", "HOXD8", "SEC61A1", "EPHA5", "ADAMTS16", "CTC-228N24.3", "FBN2", "SLC12A2", "LACTB2", "LACTB2-AS1", "XKR9", "TRAM1",   
    "ADAMTSL1", "KAT6B", "VCL","ADM", "SBF2", "CCDC73", "EIF3M", "WT1", "WT1-AS", "CHRD2", "LIPT2", "POLD3", "DUSP16", "LOH12CR1", "LOH12CR2",
    "FAN1", "KLF13",  "FMN1", "GREM1", "SCG5", "LOXL1", "LOXL1-AS1", "TBC1D21", "CRISPLD2", "ACADVL", "DLG4", "DVL2", "ADAMTS1", "ADAMTS8", 
    "MAFF", "PLA2G6", "TMEM184B"
  )
)


##8-12 Uterine fibroids and Hormonal traits (genes extracted from the GWAS catalog)

#From GWASCatalog we filter out reproductive studies by PUBMEDID

#in hpc: awk -F'\t' 'NR==1 || $2 ~ /^(PUBMED|34211149|36320039|38815977|36914876|31649266|30202859|33239672|28346442|35739095|30566500|32042192|36726022|28836065)$/' gwas_catalog_v1.0-associations_e113_r2025-01-30.tsv > repro_traits_melinda_updated_amh.tsv

reprotraits<-fread(paste0("../00.Data/GWAS_datasets/repro_traits_melinda_updated_amh.tsv"))
dim(reprotraits) 
filtered_df <- reprotraits[reprotraits$PUBMEDID %in% c(38815977, 32042192, 31649266), ]
dim(filtered_df)

subgwas <- filtered_df[, .(INTERGENIC, SNP_GENE_IDS, DOWNSTREAM_GENE_ID, UPSTREAM_GENE_ID, `DISEASE/TRAIT`)]
subgwas[INTERGENIC==0, newid:=strsplit(SNP_GENE_IDS, ", ")]
subgwas[INTERGENIC==1, newid:=strsplit(paste(UPSTREAM_GENE_ID, DOWNSTREAM_GENE_ID, sep=","), ",")]

subgwas[, .(newid, `DISEASE/TRAIT`)]


# Transform `newid` to a list of vectors per unique trait
ensg_list <- subgwas[, .(newid = unlist(newid)), by = "DISEASE/TRAIT"]  # Flatten `newid`
ensg_list <- ensg_list[, .(ENSG_IDs = list(newid)), by = "DISEASE/TRAIT"]  # Aggregate into list


ensg_list[, parsed_ensembls := lapply(ENSG_IDs, function(ids) {
  # Convert to a character string if it's a list
  ids <- as.character(ids)
  
  # Remove c(...) wrappers and any extra quotes
  ids <- gsub("c\\(|\\)|\"", "", ids)
  
  # Split by commas to flatten the list
  unlist(strsplit(ids, ","))
})]

ensg_list[, parsed_ensembls:=gsub(" ", "", parsed_ensembls)]
ensg_list[, parsed_ensembls:=gsub("\n", "", parsed_ensembls)]
new <-ensg_list[, .(gene = unlist(strsplit(parsed_ensembls, ","))), by = `DISEASE/TRAIT`]
new[, geneid:=gsub("c\\(|\\)|\"", "", gene)]
new[, gene:=NULL]
colnames(new) <- c("Trait", "Genes")

#important line to keep with unique gene-trait appearences
new_unique<-unique(new)
new_unique


gene_annotation  <- fread(paste0("../00.Data/gencode.v39.annotation.bed"))
head(new_unique)
table(new_unique$Genes)

# Remove version numbers (anything after the dot) from ensembl.id
gene_annotation <- gene_annotation %>%
  mutate(ensembl.id.clean = str_replace(ensembl.id, "\\.\\d+$", ""))

# Merge with new$Genes on the cleaned ensembl.id
reprotraits_parsed <- new_unique %>%
  left_join(gene_annotation, by = c("Genes" = "ensembl.id.clean")) 

reprotraits_parsed <- dplyr::select(reprotraits_parsed, Trait, gene.name.x)
names(reprotraits_parsed) <- c("Trait", "Genes")

duplicates_check <- reprotraits_parsed %>%
  group_by(Trait, Genes) %>%
  filter(n() > 1)

rev_reprotraits_parsed <- reprotraits_parsed %>%
  filter(Trait %in% c("Anti-Mullerian hormone levels in pre-menopausal women", 
                      "Bioavailable testosterone levels", 
                      "Endometriosis", 
                      "Estradiol levels", 
                      "Uterine leiomyomata", 
                      "Sex hormone-binding globulin levels")) 

all_studies_fem <- rbind(gwas_anm,gwas_AAM_EUR,rev_reprotraits_parsed,preprint_pcos,polyps_FGT,cysts,prolapse_genes_data)

table(all_studies_fem$Trait)

all_studies_fem <- all_studies_fem %>%
  mutate(Trait = case_when(
    Trait == "ANM" ~ "Age at Menopause",
    Trait == "AAM" ~ "Age at Menarche",
    Trait == "PCOS (preprint)" ~ "Polycystic Ovary Syndrome",
    Trait == "Uterine leiomyomata" ~ "Uterine fibroids",
    Trait == "Anti-Mullerian hormone levels in pre-menopausal women" ~ "Anti-Müllerian hormone levels",
    Trait == "all_genes_prolapse" ~ "Pelvic organ prolapse",
    TRUE ~ Trait  # Keep other traits unchanged
  ))


###################################################### reading in the DEGs in uterus, ovary, vagina, breast and myometrium

# Define input path
input_path <- "./"

# Define tissues and file paths
tissues <- c("Uterus", "Ovary", "Vagina", "Breast", "Myometrium")

file_paths <- c(paste0("../07.DEA_and_enrichments/Uterus_all_cov_no_inter_AGE_covariates_and_traits.results.rds"),
                paste0("../07.DEA_and_enrichments/Ovary_all_cov_no_inter_AGE_covariates_and_traits.results.rds"),
                paste0("../07.DEA_and_enrichments/Vagina_all_cov_no_inter_AGE_covariates_and_traits.results.rds"),
                paste0("../07.DEA_and_enrichments/BreastMammaryTissue_all_cov_no_inter_AGE_covariates_and_traits.results.rds"),
                paste0("../07.DEA_and_enrichments/Myometrium_all_cov_no_inter_AGE_covariates_and_traits.results.rds"),
)

# Read RDS files into a named list
data_list <- setNames(lapply(file_paths, readRDS), tissues)

# Extract the 'Age' subset for each dataset
age_data <- lapply(data_list, function(x) x$Age)

# Get unique gene names for each tissue
unique_genes <- lapply(age_data, function(x) unique(x$gene.name.x))

# Assign unique genes dynamically as variables (total_uterus, total_ovary, etc.)
list2env(setNames(unique_genes, paste0("total_", tissues)), envir = .GlobalEnv)

# Compute intersections with table_anm$Consensus
intersect_counts <- sapply(unique_genes, function(genes) {
  length(intersect(unique(gwas_anm$Genes), genes))
})

# Assign intersections dynamically as variables (intersect_uterus, intersect_ovary, etc.)
list2env(setNames(as.list(intersect_counts), paste0("intersect_", tissues)), envir = .GlobalEnv)

# Print the intersection counts
sapply(unique_genes, length)
print(intersect_counts)

# Extract 'Age' subset and filter differentially expressed genes (p-adj < 0.05)
degs_list <- lapply(data_list, function(x) {
  age_subset <- x$Age
  degs <- age_subset[age_subset$adj.P.Val < 0.05, ]
  degs$gene <- gsub("\\.\\d+$", "", rownames(degs))
  degs$value <- degs$gene.name.x
  return(degs$value)
})

# Assign DEGs dynamically as variables
total_degs <- setNames(degs_list, paste0("degs_", tissues))
list2env(total_degs, envir = .GlobalEnv)

# Print DEGs for each tissue
str(total_degs)

degs_counts <- sapply(total_degs, length)
print(degs_counts)

# If you want a more readable output:
cat("\nTotal DEGs count in each tissue:\n")
for (tissue in names(degs_counts)) {
  cat(tissue, ":", degs_counts[tissue], "\n")
}

table(all_studies_fem$Trait)

##################### Preparing data for the Fisher test

# Define datasets and tissues to iterate over
datasets <- list("studies_fem_repro")
tissues <- c("Uterus", "Ovary", "Vagina", "Breast", "Myometrium")

# Get the list of unique Traits in the dataset
unique_traits <- unique(all_studies_fem$Trait)
length(unique_traits)


# Initialize a dataframe to store results
fisher_results_df <- data.frame(Trait = character(),
                                Tissue = character(),
                                Dataset = character(),
                                GWAS_Genes_a = numeric(), 
                                Non_GWAS_b = numeric(), 
                                Non_DEGs_in_GWAS_c = numeric(), 
                                Non_DEGs_and_Non_GWAS_d = numeric(),
                                P_Value = numeric(),
                                Odds_Ratio = numeric(),
                                CI_Lower = numeric(),
                                CI_Upper = numeric(),
                                stringsAsFactors = FALSE)

# Initialize an empty list to store contingency tables
contingency_tables_list <- list()

# Loop through each Trait (Disease/Condition)
for (trait in unique_traits) {
  # Subset the GWAS dataset for this specific trait
  trait_gwas_genes <- all_studies_fem %>%
    filter(Trait == trait) %>%
    pull(Genes) # Extract the genes for this trait
  
  # Loop through each tissue type
  for (tissue in tissues) {
    
    # Dynamically generate total and DEGs variable per tissue
    total_Tissue <- get(paste0("total_", tissue))  # All genes in tissue
    degs_Tissue <- get(paste0("degs_", tissue))    # Differentially expressed genes (DEGs)  
    
    head(total_Uterus)
    head(degs_Uterus)
    
    # Compute overlap of GWAS genes with all genes in the tissue
    overlap_gwas_Tissue <- intersect(total_Tissue, trait_gwas_genes)
    
    # Compute values for contingency table
    a <- length(intersect(degs_Tissue, overlap_gwas_Tissue))  # Overlap between DEGs and GWAS genes
    b <- length(degs_Tissue) - a   # DEGs that are not GWAS genes
    c <- length(overlap_gwas_Tissue) - a    # GWAS genes that are not DEGs
    d <- length(total_Tissue) - (a + b + c)  # Other genes in tissue that are neither DEGs nor GWAS genes
    
    # Create contingency table
    contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                                dimnames = list(c("DEGs", "Non-DEGs"), c("GWAS", "Non-GWAS")))
    
    # Save contingency table
    contingency_tables_list[[paste0(trait, "_", tissue)]] <- contingency_table
    
    # Perform Fisher's Exact Test
    fisher_result <- fisher.test(contingency_table)
    
    # Extract test results
    p_value <- fisher_result$p.value
    odds_ratio <- fisher_result$estimate
    ci_lower <- fisher_result$conf.int[1]
    ci_upper <- fisher_result$conf.int[2]
    
    # Save results to dataframe
    fisher_results_df <- rbind(fisher_results_df, 
                               data.frame(Trait = trait,
                                          Tissue = tissue,
                                          Dataset = "studies_fem",
                                          GWAS_Genes_a = a, 
                                          Non_GWAS_b = b, 
                                          Non_DEGs_in_GWAS_c = c, 
                                          Non_DEGs_and_Non_GWAS_d = d,
                                          P_Value = p_value,
                                          Odds_Ratio = odds_ratio,
                                          CI_Lower = ci_lower,
                                          CI_Upper = ci_upper))
  }
  
}

ppp_sign<-filter(fisher_results_df, P_Value < 0.05)
fisher_results_df_1<-filter(fisher_results_df,GWAS_Genes_a!=0)

# Adjust p-values using Benjamini-Hochberg correction 
fisher_results_df <- fisher_results_df_1 %>%
  group_by(Trait) %>%  # Group by Trait
  mutate(P_Adjusted_FDR = p.adjust(P_Value, method = "fdr")) %>%  # Apply FDR correction per Trait
  ungroup() 

# Filter significant results after FDR correction
p_significant_fdr <- fisher_results_df %>% filter(P_Adjusted_FDR < 0.05)
View(p_significant_fdr)

############################## Plotting

##bubble plot
fisher_results_df <- fisher_results_df %>%
  mutate(
    log2_enrichment = log2(Odds_Ratio),   # Log transform Odds Ratio
    log10_pvalue = -log10(P_Value)        # Log transform p-value
  )

fisher_results_df$Trait[fisher_results_df$Trait == "all_genes_prolapse"] <- "Pelvic organ prolapse"


tissue_order <- c("Uterus", "Myometrium", "Ovary")
trait_order<-c("Pelvic organ prolapse","Sex hormone-binding globulin levels","Age at Menarche")

# Filter and set factor levels for ordered plotting
df_filtered <- fisher_results_df %>%
  filter(Trait %in% trait_order, Tissue %in% tissue_order) %>%
  mutate(
    Trait = factor(Trait, levels = trait_order),
    Tissue = factor(Tissue, levels = tissue_order)
  )


plot1<- ggplot(df_filtered, aes(x = Tissue, y = Trait)) +
  # Main circles
  geom_point(aes(size = log2_enrichment, fill = log10_pvalue), 
             color = "black", shape = 21, stroke = 0.5, alpha = 0.8) +
  scale_size_continuous(range = c(2, 10)) +
  scale_fill_gradient(low = "lightyellow", high = "#8B0000") +
  
  # Add stars for significant associations
  geom_text(data = df_filtered %>% filter(P_Adjusted_FDR < 0.05), 
            aes(label = "*"), size = 6, color = "black", fontface = "bold", hjust = 0.5, vjust = 0.8) +
  
  # Custom publication-ready theme
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
    axis.title = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 14, color = "black"),
    legend.position = "bottom",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey80", size = 0.5, fill = NA)
  ) +
  labs(
    title = "Enrichment of age-DEGs with reproductive traits",
    x = "Tissue",
    y = "Trait",
    fill = "-log10(p-value)",
    size = "log2 Enrichment"
  )


ggsave(
  filename = "/X/Desktop/bubble_plot_GWASenrichment.pdf",
  plot = plot1,
  width = 10, height = 4, dpi = 300
)
# Save the combined plot with the heatmap
plot2<-heatmap_plot2 #Obtained in Fig4I_4SC.R

ggsave(
  filename = "/X/Desktop/Final_combined_plots_heatmap_bubble_plot.pdf",
  plot = arrangeGrob(plot1, plot2$gtable,  ncol = 2),
  width = 10, height = 4, dpi = 300
)