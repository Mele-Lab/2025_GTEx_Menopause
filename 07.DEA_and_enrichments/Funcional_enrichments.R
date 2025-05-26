# Functional enrichment ####

#1. Retrieve databases
#2. Establish functional enrichment functions

library(WebGestaltR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(DOSE)
library(vroom)
library(tidyr)
library(msigdbr)
library(xCell)
library(GSEABase)
library(reshape2)
library(readxl)
library(dplyr)
library(tidyr)


input_path<- "../00.Data/enrichments_dbs_files/"
OMIM.term2gene<- read.csv(paste0(input_path,"/OMIM_term2gene.csv"))[2:3]
human_pheno.term2gene<-read.csv(paste0(input_path,"/human_pheno_term2gene.csv"))[2:3]
cells<-read.csv(paste0(input_path,"/cell_marker_file.csv"))[2:3]
xCell.signatures<-read.csv(paste0(input_path,"/xcell_file.csv"))[2:3]
gwas.term2gene<-read.csv(paste0(input_path,"/gwas_catalog_v1.0-associations_e113_r2025-01-30.tsv"), sep="\t")
gwas.term2gene<-gwas.term2gene[,c("DISEASE.TRAIT", "MAPPED_GENE")]
gwas.term2gene <- gwas.term2gene %>%
  mutate(Genes = gsub(" - |,|;", ",", MAPPED_GENE)) %>%  # unify separators to comma
  separate_rows(Genes, sep = ",\\s*")              # split into separate rows

H<-read.csv(paste0(input_path,"/hallmark_file.csv"))[2:3]



# Funcional enrichments ----

#########################################Functions#############################################
#Functions for dbs that only support entrez.id
# Funcional enrichments ----
ora.fun_entrez <- function(gene.list, 
                           bg.list, 
                           db = "GO:BP", 
                           pvalueCutoff = 0.05, 
                           qvalueCutoff = 0.2){
  if(db=="KEGG"){
    ora.kegg <- enrichKEGG(gene     = unlist(gene.list), # only supports entrez.id
                           universe =    bg.list,
                           keyType = "kegg", 
                           organism = "hsa",
                           minGSSize    = 10,
                           maxGSSize    = 500,
                           pAdjustMethod = "BH",
                           pvalueCutoff = pvalueCutoff,
                           qvalueCutoff = qvalueCutoff,
    )
    return(ora.kegg)
  }else if(db=="ReactomePA"){
    ora.reactome <- enrichPathway(gene     = gene.list, # only supports entrez.id
                                  universe =    bg.list,
                                  organism = "human",
                                  minGSSize    = 10,
                                  maxGSSize    = 500,
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = pvalueCutoff,
                                  qvalueCutoff = qvalueCutoff,
                                  readable = F)
    return(ora.reactome)
  }else if(db=="DO"){
    ora.DO <- enrichDO(gene          = gene.list,  # only supports entrez.id
                       universe      = bg.list,
                       ont           = "DO",
                       pAdjustMethod = "BH",
                       minGSSize     = 10,
                       maxGSSize     = 500,
                       pvalueCutoff = pvalueCutoff,
                       qvalueCutoff = qvalueCutoff,
                       readable      = FALSE)
    return(ora.DO)
  }else if(db=="DisGeNET"){
    ora.DisGeNET <- enrichDGN(gene          = gene.list,  # only supports entrez.id
                              universe      = bg.list,
                              pAdjustMethod = "BH",
                              minGSSize     = 10,
                              maxGSSize     = 500,
                              pvalueCutoff = pvalueCutoff,
                              qvalueCutoff = qvalueCutoff,
                              readable      = F)
    return(ora.DisGeNET)
  }else if(db=="OMIM"){
    ora.OMIM <- enricher(gene          = gene.list, # only supports entrez.id
                         universe      = bg.list,
                         TERM2GENE     = OMIM.term2gene,
                         minGSSize    = 10,
                         maxGSSize    = 500,
                         pAdjustMethod = "BH",
                         pvalueCutoff = pvalueCutoff,
                         qvalueCutoff = qvalueCutoff)
    return(ora.OMIM)
  }else if(db=="human_phenotype"){
    ora.human_pheno <- enricher(gene          = gene.list, # only supports entrez.id
                                universe      = bg.list,
                                TERM2GENE     = human_pheno.term2gene,
                                minGSSize    = 10,
                                maxGSSize    = 500,
                                pAdjustMethod = "BH",
                                pvalueCutoff = pvalueCutoff,
                                qvalueCutoff = qvalueCutoff)
    return(ora.human_pheno)
  }else if(db=="cell_marker"){
    ora.cell_marker <- enricher(gene       = gene.list,  # only supports gene symbol
                                universe      = bg.list,
                                TERM2GENE    = cells,
                                minGSSize    = 5,
                                maxGSSize    = 500,
                                pAdjustMethod = "BH",
                                pvalueCutoff = pvalueCutoff,
                                qvalueCutoff = qvalueCutoff)
    return(ora.cell_marker)
  }else if(db=="H"){
    ora.H <- enricher(gene          = gene.list, # only supports gene symbol
                      universe      = bg.list,
                      TERM2GENE     = H,
                      minGSSize    = 10,
                      maxGSSize    = 500,
                      pAdjustMethod = "BH",
                      pvalueCutoff = pvalueCutoff,
                      qvalueCutoff = qvalueCutoff)
    return(ora.H)
  }
}

#Functions for dbs that support gene symbols
ora.fun_symbols <- function(gene.list, 
                            bg.list, 
                            db = "GO:BP", 
                            pvalueCutoff = 0.05, 
                            qvalueCutoff = 0.2){
  
  if(db=="GO:BP"){
    ora.go <- enrichGO(gene     = gene.list,
                       universe =    bg.list,
                       keyType = "SYMBOL", #"ENSEMBL" "SYMBOL" "ENTREZID"
                       OrgDb        = org.Hs.eg.db,
                       ont          = "BP",
                       minGSSize    = 10,
                       maxGSSize    = 500,
                       pAdjustMethod = "BH",
                       pvalueCutoff = pvalueCutoff,
                       qvalueCutoff = qvalueCutoff,
                       readable = F)
    return(ora.go)
  }else if(db=="GO:MF"){
    ora.go <- enrichGO(gene     = gene.list,
                       universe =    bg.list,
                       keyType = "SYMBOL", #"ENSEMBL" "SYMBOL" "ENTREZID"
                       OrgDb        = org.Hs.eg.db,
                       ont          = "MF",
                       minGSSize    = 10,
                       maxGSSize    = 500,
                       pAdjustMethod = "BH",
                       pvalueCutoff = pvalueCutoff,
                       qvalueCutoff = qvalueCutoff,
                       readable = F)
    return(ora.go)
  }else if(db=="GO:CC"){
    ora.go <- enrichGO(gene     = gene.list,
                       universe =    bg.list,
                       keyType = "SYMBOL", #"ENSEMBL" "SYMBOL" "ENTREZID"
                       OrgDb        = org.Hs.eg.db,
                       ont          = "CC",
                       minGSSize    = 10,
                       maxGSSize    = 500,
                       pAdjustMethod = "BH",
                       pvalueCutoff = pvalueCutoff,
                       qvalueCutoff = qvalueCutoff,
                       readable = F)
    return(ora.go)
  }else if(db=="xCell"){
    ora.xCell <- enricher(gene          = gene.list,
                          universe      = bg.list,
                          TERM2GENE     = xCell.signatures,
                          minGSSize    = 10,
                          maxGSSize    = 500,
                          pAdjustMethod = "BH",
                          pvalueCutoff = pvalueCutoff,
                          qvalueCutoff = qvalueCutoff)
    return(ora.xCell)
  }else if(db=="gwas"){
    ora.gwas <- enricher(gene          = gene.list, # SYMBOL
                         universe      = bg.list, # SYMBOL
                         TERM2GENE     = gwas.term2gene,
                         minGSSize    = 5,
                         maxGSSize    = 500,
                         pAdjustMethod = "BH",
                         pvalueCutoff = pvalueCutoff,
                         qvalueCutoff = qvalueCutoff)
    
    return(ora.gwas)
  }else{
    return(NA)
  }
  
}