
#Enrichment script

##############################Sections of the code##############################
#1. Select the DEGs from the files
#2. See the common DEGs between tissues
#3. Functional enrichment analysis
################################################################################

library(WebGestaltR)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(DOSE)
library(vroom)
library(tidyr)
library(msigdbr)
library(xCell)
library(GSEABase)
library(gplots)
library(hgu95av2.db)
library(AnnotationDbi)
library(org.Hs.eg.db) 
library(GOstats)
library(reshape2)
library(vroom)
library(ggplot2)
library(rrvgo)

source("Funcional_enrichments.R")

#To run a toy example:
#- Use the example file provided: Uterus_all_cov_no_inter_AGE_covariates_and_traits.results.rds
#- analyses <- c("age")
#- tissues<- "Uterus"

###1. Select the DEGs

#These sets of genes are all the genes, we have to select now the DEGs

sexes<- c("female")
tests<- c("all_cov_no_inter", "all_cov_inter", "one_cov_no_inter", "one_cov_inter")
tissues<-c("Uterus", "Ovary", "Vagina", "BreastMammaryTissue", "CervixEndocervix", "CervixEctocervix", "FallopianTube")
analyses<-c("age","sub", "inter")

gene_annotation  <- read.csv("../00.Data/gencode.v39.annotation.bed")

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

#Establish the function to calculate the log odds ratio
calc_log_odds_ratio <- function(db) {
  if(nrow(db)>0){
    db$GeneRatio_numeric <- sapply(strsplit(as.character(db$GeneRatio), "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
    db$BgRatio_numeric <- sapply(strsplit(as.character(db$BgRatio), "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
    log_odds <- log((db$GeneRatio_numeric / db$BgRatio_numeric) / ((1 - db$GeneRatio_numeric) / (1 - db$BgRatio_numeric)))
    db$log_odds_ratio <- log_odds
  }
  return(db)
}



for (tissue in tissues){
    outpath <-"./"
    if(tissue =="BreastMammaryTissue"){
      substructures<- c("adipocyte", "lobule", "duct", "stroma", "nerve", "gynecomastoid_hyperplasia")
      
    }else if(tissue == "Ovary"){
      substructures<- c("cortex", "corpora", "medulla", "vessels")
      
    }else if (tissue =="Vagina"){
      substructures<- c("epithelium", "lamina_propria", "vessels", "stroma")
      
    }else if (tissue =="Uterus"){
      substructures<- c("myometrium", "endometrium", "vessels")
      
    }else if (tissue =="CervixEndocervix"){
      substructures<- c("vessels", "glandular_epithelium", "stroma")

      
    }else if (tissue =="CervixEctocervix"){
      substructures<- c("epithelium", "stroma", "vessels", "glands")
      
    }else if (tissue == "FallopianTube"){
      substructures<- c("vessels", "lumen", "smooth_muscle", "epithelium", "stroma") 

      }
    
    #For age we only use the model with all the covariates that change with age
    if(analy =="age" & test =="all_cov_no_inter"){
        all<- readRDS( paste0("./", tissue, "_", test, "_AGE_covariates_and_traits.results.rds"))
        all<-all$Age
        print(paste0("#-------------- ", tissue, " - ", test, " - ", analy, " --------------#"))

        if (!tissue %in% c("CervixEndocervix","CervixEctocervix", "FallopianTube")){
          
          degs<- all[all$adj.P.Val<0.05,]
        }else{
          degs<- all[all$P.Value<0.01,]
          
        }
        
        neg_ut<- degs[degs$logFC<0,]
        pos_ut<-degs[degs$logFC>0,]

        neg_value <- ifelse(is.null(neg), 0, neg)
        pos_value <- ifelse(is.null(pos), 0, pos)

        degs<- intersect(rownames(pos_ut), rownames(pos_ov))
          
        if ((neg_value + pos_value) > 3) {
            
        ### 3.Enrichment
        dea_res_all<-all
        dea_res <- rownames(all) ##all the genes that appear in the tissue
            
        # Lists of DEGs ----
        de.genes <- rownames(degs)
        #Gene annotation
        
        # Functional enrichments 1 ----
        dbs <- c("GO:BP", "GO:MF","GO:CC", "KEGG", "ReactomePA","human_phenotype","DisGeNET","OMIM", "gwas")#, "cell_marker", "xCell")#"H", "DO","cell_marker", "xCell"
        
        # input -> variable-DEGs
        # bg -> tissue expressed genes
        
        #Create two sets of both tissue_exprs_genes and variable_DEGs, one selecting the 
        #entrez.id for those databases that do not use the gene.name and another one selecting the
        #gene.name, since with this latter approach we remove less genes.
        
        tissue_exprs_genes_entrez.id <- sapply(dea_res, function(gene) gene_annotation[gene_annotation$ensembl.id==gene, "entrez.id"])
        tissue_exprs_genes_entrez.id <- unlist(tissue_exprs_genes_entrez.id[!is.na(tissue_exprs_genes_entrez.id)])
        
        tissue_exprs_genes_symbols <- sapply(dea_res, function(gene) gene_annotation[gene_annotation$ensembl.id==gene, "gene.name.x"])
        tissue_exprs_genes_symbols <- unlist(tissue_exprs_genes_symbols[!is.na(tissue_exprs_genes_symbols)])
        
        
        variable_DEGs_entrez.id <-  sapply(de.genes, function(gene) gene_annotation[gene_annotation$ensembl.id==gene, "entrez.id"])
        variable_DEGs_entrez.id<- unlist(variable_DEGs_entrez.id[!is.na(variable_DEGs_entrez.id)])
        
        variable_DEGs_symbols <-  sapply(de.genes, function(gene) gene_annotation[gene_annotation$ensembl.id==gene, "gene.name.x"])
        variable_DEGs_symbols<- variable_DEGs_symbols[!is.na(variable_DEGs_symbols)]
        
        
        print("Number of variable-DEGs with entrez.id")
        print(length(variable_DEGs_entrez.id))
        print("Number of variable-DEGs with symbols")
        print(length(variable_DEGs_symbols))
        
        
        #Select the db that uses entrez or symbols
        dbs_entrez<- c("KEGG", "ReactomePA","DO","DisGeNET","OMIM", "human_phenotype","H","cell_marker")
        dbs_symbols<- dbs[!dbs %in% dbs_entrez]
        
        #Perform the enrichments
        
        # Functional enrichments 2 ----
        # input -> variable-DEGs-upregulated(downregulated)
        # bg -> tissue expressed genes
        
        directions<- c("Upregulated", "Downregulated")
        
        if (tissue %in% c("Vagina", "BreastMammaryTissue", "Uterus", "Ovary")){
          get_degs <- function(db,direction){
          if(direction == "Upregulated" & db %in% dbs_symbols){
            sapply(rownames(dea_res_all[dea_res_all$adj.P.Val < 0.05 &
                                          dea_res_all$logFC > 0,]), function(gene)
                                            gene_annotation[gene_annotation$ensembl.id==gene, "gene.name.x"]
            )
          }else if(direction == "Downregulated" & db %in% dbs_symbols){
            sapply(rownames(dea_res_all[dea_res_all$adj.P.Val < 0.05 &
                                          dea_res_all$logFC < 0,]), function(gene)
                                            gene_annotation[gene_annotation$ensembl.id==gene, "gene.name.x"]
            )
          }else if(direction == "Upregulated" & db %in% dbs_entrez){
            sapply(rownames(dea_res_all[dea_res_all$adj.P.Val < 0.05 &
                                          dea_res_all$logFC > 0,]), function(gene)
                                            gene_annotation[gene_annotation$ensembl.id==gene,  "entrez.id"])
          }else if(direction == "Downregulated" & db %in% dbs_entrez){
            sapply(rownames(dea_res_all[dea_res_all$adj.P.Val < 0.05 &
                                          dea_res_all$logFC < 0,]), function(gene)
                                            gene_annotation[gene_annotation$ensembl.id==gene, "entrez.id"])

          }
            }
            }else{
              get_degs <- function(db,direction){
                if(direction == "Upregulated" & db %in% dbs_symbols){
                  sapply(rownames(dea_res_all[dea_res_all$P.Value < 0.01 &
                                                dea_res_all$logFC > 0,]), function(gene)
                                                  gene_annotation[gene_annotation$ensembl.id==gene, "gene.name.x"]
                  )
                }else if(direction == "Downregulated" & db %in% dbs_symbols){
                  sapply(rownames(dea_res_all[dea_res_all$P.Value < 0.01 &
                                                dea_res_all$logFC < 0,]), function(gene)
                                                  gene_annotation[gene_annotation$ensembl.id==gene, "gene.name.x"]
                  )
                }else if(direction == "Upregulated" & db %in% dbs_entrez){
                  sapply(rownames(dea_res_all[dea_res_all$P.Value < 0.01 &
                                                dea_res_all$logFC > 0,]), function(gene)
                                                  gene_annotation[gene_annotation$ensembl.id==gene,  "entrez.id"])
                }else if(direction == "Downregulated" & db %in% dbs_entrez){
                  sapply(rownames(dea_res_all[dea_res_all$P.Value < 0.01 &
                                                dea_res_all$logFC < 0,]), function(gene)
                                                  gene_annotation[gene_annotation$ensembl.id==gene, "entrez.id"])
                }
              }
            }
            
        
            degs_up_down <- lapply(directions, function(dir) {
              sapply(dbs, function(i) get_degs(i, dir))
            })
            
            
            names(degs_up_down)<- directions
            
            degs_up_down$Upregulated<- degs_up_down$Upregulated[rownames(degs_up_down$Upregulated) %in% de.genes,]
            degs_up_down$Downregulated<- degs_up_down$Downregulated[rownames(degs_up_down$Downregulated) %in% de.genes,]
            
            print("Number of upregulated and downregulated variable-DEGs")
            print(sapply(directions, function(dir) nrow(degs_up_down[[dir]])))
            
            degs_up_down <- lapply(degs_up_down, as.data.frame)
            
            
            ora_results_2 <- lapply(directions, function(i) {
              lapply(dbs, function(db) {
                print(paste0(tissue, ": ", i ," - ", db))
                if (db %in% dbs_symbols) {
                  as.data.frame(ora.fun_symbols(
                    gene.list = degs_up_down[[i]][[db]],
                    bg.list = tissue_exprs_genes_symbols,
                    db = db,
                    pvalueCutoff = 0.05
                  ))
                } else if (db %in% dbs_entrez) {
                  as.data.frame(ora.fun_entrez(
                    gene.list = degs_up_down[[i]][[db]],
                    bg.list = tissue_exprs_genes_entrez.id,
                    db = db,
                    pvalueCutoff = 0.05
                  ))
                }
              })
            })
            
            for (i in seq_along(ora_results_2)) {
              # Assign names to the first level
              names(ora_results_2)[i] <- directions[i]
              
              # Iterate through the second level of the nested list
              for (j in seq_along(ora_results_2[[i]])) {
                # Assign names to the sublevels
                names(ora_results_2[[i]])[j] <- dbs[j]
              }
            }
            
            
            ora_results_2<- lapply(ora_results_2, function(inner_list){
              inner_list_log_calculated<- lapply(inner_list, calc_log_odds_ratio)})
            
            if(analy == "age"){
              saveRDS(ora_results_2,
                      paste0(outpath,"/", tissue, "_AGE_DEGs.ORA.up_down_variable_DEGs.variable_DEGs.rds"))
              print(paste0("--Saving: /",tissue, "_AGE_DEGs.ORA.up_down_variable_DEGs.variable_DEGs.rds"))

            }else if(analy== "sub"){
              saveRDS(ora_results_2,
                      paste0(outpath,"/", tissue, "_", sub, "_DEGs.ORA.up_down_variable_DEGs.variable_DEGs.rds"))
              print(paste0("--Saving: /",tissue, "_", sub, "_DEGs.ORA.up_down_variable_DEGs.variable_DEGs.rds"))
              
            }else{
              saveRDS(ora_results_2,
                      paste0(outpath,"/", tissue, "_", sub, "_interaction_DEGs.ORA.up_down_variable_DEGs.variable_DEGs.rds"))
              print(paste0("--Saving: /",tissue, "_", sub, "interaction_DEGs.ORA.up_down_variable_DEGs.variable_DEGs.rds"))
              
                  
                }
              }
        }
    
      for (test in tests){
        for (sub in substructures){
          for (analy in c("sub", "inter")){
              if(analy== "sub" & test =="all_cov_no_inter"){
                all<- readRDS( paste0(inpath, "/", tissue, "_covariates_and_traits.results.rds"))
                all<-all[[sub]]
                print(paste0("#-------------- ", tissue, " - ", test, " - ", analy, " - ", sub, " --------------#"))
                
              }else if(analy== "inter" & test =="one_cov_inter"){
                all<- readRDS( paste0(inpath, test, "/", tissue, "_",sub, "_", test, "_covariates_and_traits.results.rds"))
                all<- all[[paste0("Interaction_Age_", sub)]]
                print(paste0("#-------------- ", tissue, " - ", test, " - ", analy, " - ", sub, " --------------#"))
                
              }else if(analy== "sub" & test =="one_cov_no_inter"){
                all<- readRDS( paste0(inpath, test, "/", tissue, "_",sub, "_", test, "_covariates_and_traits.results.rds"))
                all<- all[[sub]]
                print(paste0("#-------------- ", tissue, " - ", test, " - ", analy, " - ", sub, " --------------#"))
                
              }else{
                 break
               }
          
          if (!tissue %in% c("CervixEndocervix","CervixEctocervix", "FallopianTube")){
            degs<- all[all$adj.P.Val<0.05,]
          }else{
            degs<- all[all$P.Value<0.01,]
          }
            
          neg<-degs[degs$logFC<0,]
          pos<-degs[degs$logFC>0,]
          neg<-nrow(degs[degs$logFC<0,])
          pos<-nrow(degs[degs$logFC>0,])
          neg_value <- ifelse(is.null(neg), 0, neg)
          pos_value <- ifelse(is.null(pos), 0, pos)
            
          if ((neg_value + pos_value) > 3) {
          ### 3.Enrichment
          
          dea_res_all<-all
          dea_res <- rownames(all) ##all the genes that appear in the tissue
          
          
          # Lists of DEGs ----
          de.genes <- rownames(degs)
          #Gene annotation
  
          # Functional enrichments 1 ----
          dbs <- c("GO:BP", "GO:MF","GO:CC", "KEGG", "ReactomePA","human_phenotype","DisGeNET","OMIM", "gwas")
          
          # input -> variable-DEGs
          # bg -> tissue expressed genes
          
          tissue_exprs_genes_entrez.id <- sapply(dea_res, function(gene) gene_annotation[gene_annotation$ensembl.id==gene, "entrez.id"])
          tissue_exprs_genes_entrez.id <- unlist(tissue_exprs_genes_entrez.id[!is.na(tissue_exprs_genes_entrez.id)])
          
          tissue_exprs_genes_symbols <- sapply(dea_res, function(gene) gene_annotation[gene_annotation$ensembl.id==gene, "gene.name.x"])
          tissue_exprs_genes_symbols <- unlist(tissue_exprs_genes_symbols[!is.na(tissue_exprs_genes_symbols)])
          
          
          variable_DEGs_entrez.id <-  sapply(de.genes, function(gene) gene_annotation[gene_annotation$ensembl.id==gene, "entrez.id"])
          variable_DEGs_entrez.id<- unlist(variable_DEGs_entrez.id[!is.na(variable_DEGs_entrez.id)])
          
          variable_DEGs_symbols <-  sapply(de.genes, function(gene) gene_annotation[gene_annotation$ensembl.id==gene, "gene.name.x"])
          variable_DEGs_symbols<- variable_DEGs_symbols[!is.na(variable_DEGs_symbols)]
          
          print("Number of variable-DEGs with entrez.id")
          print(length(variable_DEGs_entrez.id))
          print("Number of variable-DEGs with symbols")
          print(length(variable_DEGs_symbols))
          
          
          #Select the db that uses entrez or symbols
          dbs_entrez<- c("KEGG", "ReactomePA","DO","DisGeNET","OMIM", "human_phenotype","H","cell_marker")
          dbs_symbols<- dbs[!dbs %in% dbs_entrez]
          
          if(test =="all_cov_no_inter" | test =="one_cov_no_inter"){
            
            directions<- c("Upregulated", "Downregulated")
            
            if (tissue %in% c("Vagina", "BreastMammaryTissue", "Uterus", "Ovary")){
              get_degs <- function(db,direction){
                if(direction == "Upregulated" & db %in% dbs_symbols){
                  sapply(rownames(dea_res_all[dea_res_all$adj.P.Val < 0.05 &
                                                dea_res_all$logFC > 0,]), function(gene)
                                                  gene_annotation[gene_annotation$ensembl.id==gene, "gene.name.x"]
                  )
                }else if(direction == "Downregulated" & db %in% dbs_symbols){
                  sapply(rownames(dea_res_all[dea_res_all$adj.P.Val < 0.05 &
                                                dea_res_all$logFC < 0,]), function(gene)
                                                  gene_annotation[gene_annotation$ensembl.id==gene, "gene.name.x"]
                  )
                }else if(direction == "Upregulated" & db %in% dbs_entrez){
                  sapply(rownames(dea_res_all[dea_res_all$adj.P.Val < 0.05 &
                                                dea_res_all$logFC > 0,]), function(gene)
                                                  gene_annotation[gene_annotation$ensembl.id==gene,  "entrez.id"])
                }else if(direction == "Downregulated" & db %in% dbs_entrez){
                  sapply(rownames(dea_res_all[dea_res_all$adj.P.Val < 0.05 &
                                                dea_res_all$logFC < 0,]), function(gene)
                                                  gene_annotation[gene_annotation$ensembl.id==gene, "entrez.id"])
                  
                }
              }
            }else{
              get_degs <- function(db,direction){
                if(direction == "Upregulated" & db %in% dbs_symbols){
                  sapply(rownames(dea_res_all[dea_res_all$P.Value < 0.01 &
                                                dea_res_all$logFC > 0,]), function(gene)
                                                  gene_annotation[gene_annotation$ensembl.id==gene, "gene.name.x"]
                  )
                }else if(direction == "Downregulated" & db %in% dbs_symbols){
                  sapply(rownames(dea_res_all[dea_res_all$P.Value < 0.01 &
                                                dea_res_all$logFC < 0,]), function(gene)
                                                  gene_annotation[gene_annotation$ensembl.id==gene, "gene.name.x"]
                  )
                }else if(direction == "Upregulated" & db %in% dbs_entrez){
                  sapply(rownames(dea_res_all[dea_res_all$P.Value < 0.01 &
                                                dea_res_all$logFC > 0,]), function(gene)
                                                  gene_annotation[gene_annotation$ensembl.id==gene,  "entrez.id"])
                }else if(direction == "Downregulated" & db %in% dbs_entrez){
                  sapply(rownames(dea_res_all[dea_res_all$P.Value < 0.01 &
                                                dea_res_all$logFC < 0,]), function(gene)
                                                  gene_annotation[gene_annotation$ensembl.id==gene, "entrez.id"])
                }
              }
            }
            
            
            degs_up_down <- lapply(directions, function(dir) {
              sapply(dbs, function(i) get_degs(i, dir))
            })
            
    
            names(degs_up_down)<- directions
            
            degs_up_down$Upregulated<- degs_up_down$Upregulated[rownames(degs_up_down$Upregulated) %in% de.genes,]
            degs_up_down$Downregulated<- degs_up_down$Downregulated[rownames(degs_up_down$Downregulated) %in% de.genes,]
            
            print("Number of upregulated and downregulated variable-DEGs")
            print(sapply(directions, function(dir) nrow(degs_up_down[[dir]])))
            
            degs_up_down <- lapply(degs_up_down, as.data.frame)
            
            
            ora_results_2 <- lapply(directions, function(i) {
              lapply(dbs, function(db) {
                print(paste0(tissue, ": ", i ," - ", db))
                if (db %in% dbs_symbols) {
                  as.data.frame(ora.fun_symbols(
                    gene.list = degs_up_down[[i]][[db]],
                    bg.list = tissue_exprs_genes_symbols,
                    db = db,
                    pvalueCutoff = 0.05
                  ))
                } else if (db %in% dbs_entrez) {
                  as.data.frame(ora.fun_entrez(
                    gene.list = degs_up_down[[i]][[db]],
                    bg.list = tissue_exprs_genes_entrez.id,
                    db = db,
                    pvalueCutoff = 0.05
                  ))
                }
              })
            })
            
            for (i in seq_along(ora_results_2)) {
              # Assign names to the first level
              names(ora_results_2)[i] <- directions[i]
              
              # Iterate through the second level of the nested list
              for (j in seq_along(ora_results_2[[i]])) {
                # Assign names to the sublevels
                names(ora_results_2[[i]])[j] <- dbs[j]
              }
            }
            
    
            ora_results_2<- lapply(ora_results_2, function(inner_list){
              inner_list_log_calculated<- lapply(inner_list, calc_log_odds_ratio)})
            
          }else{
            
            if (tissue %in% c("Vagina", "BreastMammaryTissue", "Uterus", "Ovary")){
              
              get_degs <- function(db){
                if(db %in% dbs_symbols){
                  sapply(rownames(dea_res_all[dea_res_all$adj.P.Val < 0.05,]), function(gene)
                    gene_annotation[gene_annotation$ensembl.id==gene, "gene.name.x"]
                  )
                }else if( db %in% dbs_entrez){
                  sapply(rownames(dea_res_all[dea_res_all$adj.P.Val < 0.05,]), function(gene)
                    gene_annotation[gene_annotation$ensembl.id==gene,  "entrez.id"])
                }
              }
            }else{
              get_degs <- function(db){
                if(db %in% dbs_symbols){
                  sapply(rownames(dea_res_all[dea_res_all$P.Value < 0.01,]), function(gene)
                    gene_annotation[gene_annotation$ensembl.id==gene, "gene.name.x"]
                  )
                }else if( db %in% dbs_entrez){
                  sapply(rownames(dea_res_all[dea_res_all$P.Value < 0.01,]), function(gene)
                    gene_annotation[gene_annotation$ensembl.id==gene,  "entrez.id"])
                }
              }
            }
            
            
            degs_all <- sapply(dbs, function(i) get_degs(i))
            degs_all<-as.data.frame(degs_all)
            
            
            ora_results_2 <- lapply(dbs, function(db) {
              print(paste0(tissue, ": ", db))
              if (db %in% dbs_symbols) {
                as.data.frame(ora.fun_symbols(
                  gene.list = degs_all[[db]],
                  bg.list = tissue_exprs_genes_symbols,
                  db = db,
                  pvalueCutoff = 0.05
                ))
              } else if (db %in% dbs_entrez) {
                as.data.frame(ora.fun_entrez(
                  gene.list = degs_all[[db]],
                  bg.list = tissue_exprs_genes_entrez.id,
                  db = db,
                  pvalueCutoff = 0.05
                ))
              }
            })
            
            names(ora_results_2) <- dbs
            
            ora_results_2<- lapply(ora_results_2, calc_log_odds_ratio)
            
          }
          
          if(analy == "age"){
            saveRDS(ora_results_2,
                    paste0(outpath,"/", tissue, "_AGE_DEGs.ORA.up_down_variable_DEGs.variable_DEGs.rds"))
            print(paste0("--Saving: /",tissue, "_AGE_DEGs.ORA.up_down_variable_DEGs.variable_DEGs.rds"))
            
          }else if(analy== "sub"){
            saveRDS(ora_results_2,
                    paste0(outpath,"/", tissue, "_", sub, "_DEGs.ORA.up_down_variable_DEGs.variable_DEGs.rds"))
            print(paste0("--Saving: /",tissue, "_", sub, "_DEGs.ORA.up_down_variable_DEGs.variable_DEGs.rds"))
            
          }else if (analy =="inter"){
            saveRDS(ora_results_2,
                    paste0(outpath,"/", tissue, "_", sub, "_interaction_DEGs.ORA.up_down_variable_DEGs.variable_DEGs.rds"))
            print(paste0("--Saving: /",tissue, "_", sub, "_interaction_DEGs.ORA.up_down_variable_DEGs.variable_DEGs.rds"))
            
            }
            }
          }

      }
      }
    }
 


  