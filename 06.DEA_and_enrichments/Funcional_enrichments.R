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

n<-2 #CHANGE TO 0 IF YOU WANT TO REDO THE FILE WITH THE ANNOTATIONS FOR EACH DATABASE

if(n<1){
    # cell marker data ----
  print("---- Preparing cell marker data ----")
  cell_marker_data <- read_excel("~/X/Laura/00.Data/Cell_marker_Human.xlsx")
  
  ## instead of `cell_name`, users can use other features (e.g. `cancerType`)
  cells <- cell_marker_data %>%
    dplyr::select(cell_name, GeneID) %>%
    dplyr::mutate(GeneID = strsplit(GeneID, ', ')) %>%
    tidyr::unnest(cols = c(GeneID))
  #cells <- cells[!cells$cell_name %in% c(unique(grep("Cancer", unique(cells$cell_name), value = T)),
  #                                      unique(grep("cancer", unique(cells$cell_name), value = T))),]
  cells <- cells[cells$cell_name %in% names(table(cells$cell_name))[table(cells$cell_name) >= 5 &
                                                                  table(cells$cell_name) <= 50 ],]
  cells <- as.data.frame(cells)
  cells <- cells[!is.na(cells$GeneID),]
  #cells <- cells[-which(cells$GeneID=="NA"),]
  #cells <- cells[-grep("family", cells$GeneID),]
  cells$GeneID <- gsub("^\\[", "", cells$GeneID, perl = T)
  cells$GeneID <- gsub("\\]$", "", cells$GeneID, perl = T)
  colnames(cells) <- c("term", "gene")
  
  
  write.csv(cells, "~/X/Laura/00.Data/cell_marker_file.csv")
  
  # msigDB ----
  print("---- Preparing hallmark data ----")
  
  #Retrieve human H (hallmark gene sets)
  msigdbr_H <- msigdbr(species = "human", category = "H")
  H <- msigdbr_H %>%
    dplyr::select(gs_name, entrez_gene)
  H$entrez_gene <- as.character(H$entrez_gene)
  colnames(H) <- c("term", "gene")
  H <- as.data.frame(H)
  
  
  write.csv(H, "~/X/Laura/00.Data/hallmark_file.csv")
  
  # xCell cell type signatures ----
  print("---- Preparing xCell data ----")
  cell_signatures <- geneIds(xCell.data$signatures)
  xCell.signatures <- as.data.frame(unlist(cell_signatures))
  genesym<-xCell.signatures$`unlist(cell_signatures)`
  xCell.signatures[,1]<- rownames(xCell.signatures)
  xCell.signatures[,2]<- genesym
  
  xCell.signatures <- xCell.signatures[,c(1,2)]
  rownames(xCell.signatures) <- NULL
  colnames(xCell.signatures) <- c("term", "gene")
  xCell.signatures <- xCell.signatures[!is.na(xCell.signatures$gene),]
  
  write.csv(xCell.signatures, "~/X/Laura/00.Data/xcell_file.csv")

  # GWAS 2019 ----
  
  print("---- Preparing GWAS data ----")
  sample_text <- "~/X/Laura/09.Functional_enrichment/GWAS_Catalog_2019.txt"
  read_the_text <- scan(file = sample_text, # if your data is in a text file, use the file argument
                        what = character(),
                        sep = "\n")
  
  split_each_line_by_tab <- strsplit(x = read_the_text,
                                     split = "\t")
  get_element_names <- lapply(X = split_each_line_by_tab,
                              FUN = `[[`,
                              i = 1)
  get_element_values <- lapply(X = split_each_line_by_tab,
                               FUN = `[`,
                               i = c(-1, -2))
  gwas_db <- setNames(object = get_element_values,
                      nm = get_element_names)
  gwas.term2gene <- do.call(rbind.data.frame,
          lapply(1:length(gwas_db), function(i)
    cbind.data.frame("gwas_trait" = rep(names(gwas_db)[[i]], nrow(melt(gwas_db[[i]])) ),
                    "gene_name" = as.data.frame(melt(gwas_db[[i]]) ))
    )
  )
  
  gwas.term2gene <- gwas.term2gene[,c(1,2)]
  colnames(gwas.term2gene) <- c("term", "gene")
  gwas.term2gene <- gwas.term2gene[!is.na(gwas.term2gene$gene),]
  gwas.term2gene <- gwas.term2gene[gwas.term2gene$term %in% names(table(gwas.term2gene$term)[table(gwas.term2gene$term)>=5 &
                                                                                              table(gwas.term2gene$term)<=50 ]), ]

  write.csv(gwas.term2gene, "~/X/Laura/00.Data/gwas_term2gene_2019Catalog.csv")
  
  
  ##NEW CATALOG
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  m<-read.csv("~/X/Laura/full_gwas_catalog.parsed4_enrichments.geneids.tsv", sep= "\t")
  colnames(m)<-c("gwas_trait","value")
  gene_mapping <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),  # Columns to retrieve
    filters = "ensembl_gene_id",                      # Filter by Ensembl IDs
    values = m$value,                             # Ensembl IDs to map
    mart = ensembl
  )
  colnames(gene_mapping)<-c("ensembl_gene_id", "symbol")
  colnames(m)<-c("gwas_trait", "ensembl_gene_id")
  
  data_with_symbols <- merge(m, gene_mapping, by = "ensembl_gene_id", all.x = FALSE)
  gwas2term<-  data_with_symbols %>% distinct()
  gwas2term<-gwas2term[,2:3]
  colnames(gwas2term)<-c("gwas_trait","value")
  write.csv(gwas2term, "~/X/Laura/00.Data/gwas_term2gene_CatalogPau.csv")
  nrow(gwas2term)
  
  # Human phenotype ----
  print("---- Preparing human phenotype data ----")
  human_pheno_db <- loadGeneSet(
    organism = "hsapiens",
    enrichDatabase =  "phenotype_Human_Phenotype_Ontology",
    enrichDatabaseFile = NULL,
    enrichDatabaseType = NULL,
    enrichDatabaseDescriptionFile = NULL,
    cache = NULL,
    hostName = "http://www.webgestalt.org/"
  )
  human_pheno.term2gene <- human_pheno_db$geneSet
  human_pheno.term2name <- as.data.frame(human_pheno_db$geneSetDes)
  human_pheno.term2gene$term <- sapply(human_pheno.term2gene$geneSet, function(gs) human_pheno.term2name[human_pheno.term2name$geneSet==gs, "description"])
  human_pheno.term2gene <- human_pheno.term2gene[,c(4,3)]
  colnames(human_pheno.term2gene) <- c("term", "gene")
  human_pheno.term2gene <- human_pheno.term2gene[-which(human_pheno.term2gene$term=="All"),]
  
  write.csv(human_pheno.term2gene, "~/X/Laura/00.Data/human_pheno_term2gene.csv")
  
  # OMIM ---
  print("---- Preparing OMIM data ----")
  OMIM_db <- loadGeneSet(
    organism = "hsapiens",
    enrichDatabase =  "disease_OMIM",
    enrichDatabaseFile = NULL,
    enrichDatabaseType = NULL,
    enrichDatabaseDescriptionFile = NULL,
    cache = NULL,
    hostName = "http://www.webgestalt.org/"
  )
  OMIM.term2gene <- OMIM_db$geneSet
  OMIM.term2name <- as.data.frame(OMIM_db$geneSetDes)
  OMIM.term2gene$term <- sapply(OMIM.term2gene$geneSet, function(gs) OMIM.term2name[OMIM.term2name$geneSet==gs, "description"])
  OMIM.term2gene <- OMIM.term2gene[,c(4,3)]
  colnames(OMIM.term2gene) <- c("term", "gene")
  
  write.csv(OMIM.term2gene, "~/X/Laura/00.Data/OMIM_term2gene.csv")
  

}else{
  input_path<- "~/X/Laura/00.Data/enrichments_dbs_files/"
  OMIM.term2gene<- read.csv(paste0(input_path,"/OMIM_term2gene.csv"))[2:3]
  human_pheno.term2gene<-read.csv(paste0(input_path,"/human_pheno_term2gene.csv"))[2:3]
  cells<-read.csv(paste0(input_path,"/cell_marker_file.csv"))[2:3]
  xCell.signatures<-read.csv(paste0(input_path,"/xcell_file.csv"))[2:3]
  gwas.term2gene<-read.csv(paste0(input_path,"/gwas_term2gene_CatalogPau.csv"))[2:3]
  H<-read.csv(paste0(input_path,"/hallmark_file.csv"))[2:3]
}


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