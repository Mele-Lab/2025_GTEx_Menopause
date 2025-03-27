#Enrichment plots
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(grid)
library(rrvgo)

for (sex in sexes){
  for (tissue in tissues){
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
      substructures<- c("vessels", "lumen", "smooth_muscle", "epithelium", "stroma") ##always control for adipose, since it is contamination
    }
    
  
  for (test in tests){
    if (test == "one_cov_inter" | test=="one_cov_no_inter"){
      inpath <-paste0("~/X/Laura/12.Tissues_substructures/Enrichment/Enrichment_substructures_v10//", tissue, "/one_cov/")
      outpath<- paste0("~/X/Laura/12.Tissues_substructures/Figures/Final_v10/Enrichments/one_cov_enrichments/")
    }else{
      inpath <-paste0("~/X/Laura/12.Tissues_substructures/Enrichment/Enrichment_substructures_v10//", tissue, "/")
      outpath<- paste0("~/X/Laura/12.Tissues_substructures/Figures/Final_v10/Enrichments/")
    }
    if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}
    
    for (analy in analyses){
      if(analy == "age" & test =="all_cov_no_inter"){
            if(file.exists(paste0(inpath, tissue, "_AGE_DEGs.ORA.up_down_variable_DEGs.variable_DEGs.rds"))) {
              print(paste0("-----",tissue, "------", analy," ", test, "------"))
              results2_prob<- readRDS(paste0(inpath, tissue, "_AGE_DEGs.ORA.up_down_variable_DEGs.variable_DEGs.rds"))
        }else{
          next
        }
      dbs <- c("GO:BP", "GO:MF","GO:CC", "KEGG", "ReactomePA","human_phenotype","H", "DO","DisGeNET", "gwas")

      #############################################################

      
      ##REDUCE TERMS##
      #Extract the GO terms and p-values from your enrichment results
      go_terms <- results2_prob$Downregulated$`GO:BP`$ID  # GO IDs
      pvalues <- results2_prob$Downregulated$`GO:BP`$p.adjust  # Adjusted p-values

      # # Calculate semantic similarity between GO terms
      simMatrix <- calculateSimMatrix(go_terms,
                                      orgdb = "org.Hs.eg.db",
                                      ont = "BP",
                                      method = "Wang")  # "Rel" is one method; others include "Resnik", "Lin", etc.
      # Reduce the terms based on similarity
      reduced_terms <- reduceSimMatrix(simMatrix,
                                       scores = setNames(-log10(pvalues), go_terms),  # Use -log10 of p-values as the score
                                       threshold = 0.7,  # Similarity threshold: lower = stricter grouping
                                       orgdb = "org.Hs.eg.db")

      # # # View reduced results
        head(reduced_terms)

      ###################No clusters###############################
      #
      ora_results<- lapply(results2_prob, function(sublist) {
        lapply(sublist, function(inner_list) {
          inner_list_df <- as.data.frame(inner_list)
          arrange(inner_list_df, inner_list_df[["p.adjust"]])
        })
      })
      top_30<- lapply(ora_results, function(sublist) {
        lapply(sublist, function(inner_list) {
          inner_list_df <- head(inner_list, 10)
        })
      })

      plot_list_up <- list()
      plot_list_down <- list()

      for (db in dbs) {
        for (element in names(top_30)){
          if ( "log_odds_ratio" %in% names(top_30[[element]][[db]])) {
            if (element == "Upregulated"){

              # Plot for Age_Upregulated
              plot_upregulated <- ggplot(top_30$Upregulated[[db]], aes(x = log_odds_ratio, y = reorder(Description, log_odds_ratio), color = p.adjust, size = Count)) +
                geom_point() +
                scale_color_gradient(low = "#8B0000", high = "#63B8FF") +
                labs(x = "Log odds ratio", title = paste0("Upregulated - ", db), y = NULL) +
                theme_minimal(base_size = 12) +
                theme(axis.text.y = element_text(size = 13, color = "black"), 
                      axis.text.x = element_text(size = 13, color = "black"),
                      axis.title = element_text(size = 13, color = "black"), 
                      plot.title = element_text(size = 13, hjust = 0.5, color = "black"),
                      legend.text = element_text(size = 13, color = "black"),
                      panel.grid = element_blank(),                    # Remove all grid lines
                      axis.line.x = element_line(size = 0.5, color = "grey80"), # Add x-axis line
                      axis.line.y = element_line(size = 0.5, color = "grey80"),
                      panel.border = element_rect(color = "grey80", size = 0.5, fill = NA) # Add full border
                )

              plot_list_up[[element]][[db]] <- plot_upregulated
            }else{

              # Plot for Age_Downregulated
              plot_downregulated <- ggplot(top_30$Downregulated[[db]], aes(x = log_odds_ratio, y = reorder(Description, log_odds_ratio), color = p.adjust, size = Count)) +
                geom_point() +
                scale_color_gradient(low = "#8B0000", high = "#63B8FF") +
                labs(x = "Log odds ratio", title = paste0("Downregulated - ", db), y = NULL) +
                theme_minimal(base_size = 12) +
                theme(axis.text.y = element_text(size = 13, color = "black"), 
                      axis.text.x = element_text(size = 13, color = "black"),
                      axis.title = element_text(size = 13, color = "black"), 
                      plot.title = element_text(size = 13, hjust = 0.5, color = "black"),
                      legend.text = element_text(size = 13, color = "black"),
                      panel.grid = element_blank(),                    # Remove all grid lines
                      axis.line.x = element_line(size = 0.5, color = "grey80"), # Add x-axis line
                      axis.line.y = element_line(size = 0.5, color = "grey80"),
                      panel.border = element_rect(color = "grey80", size = 0.5, fill = NA) # Add full border
                )
              plot_list_down[[element]][[db]] <- plot_downregulated
            }
          }
          }
      #

      }
      
      element_to_plot1 <- names(plot_list_up)  # Use a valid index 'n'
      element_to_plot2 <- names(plot_list_down)  # Use a valid index 'n'
      db_to_plot<-dbs[1]
      # Example to use the first database, adjust as needed
      for (db_to_plot in c(dbs[1], dbs[4], dbs[5], dbs[10])){
        if(db_to_plot == dbs[1]){
          d<- "_gobp.pdf"
        }else if(db_to_plot == dbs[4]){
          
          d<- "_KEGG.pdf"
        }else if(db_to_plot == dbs[5]){
          d<- "_ReactPA.pdf"
        }else{
          d<- "_gwas.pdf"
          
        }
        if (!is.null(element_to_plot1) ){
          if (db_to_plot %in% names(plot_list_up[[element_to_plot1]])){
            
           print(paste0("-----Saving up enrichments for ", analy, "------"))
            pdf(paste0(outpath,tissue,"_",analy, "_", test,"_UP_", d), width = 10, height = 6)  # Adjust width and height as needed
            grid.arrange(plot_list_up[[element_to_plot1]][[db_to_plot]])
            dev.off()
          }
        }

        if( !is.null(element_to_plot2)){
          if(db_to_plot %in% names(plot_list_down[[element_to_plot2]])){

            print(paste0("-----Saving down enrichments for ", analy, "------"))
            
          pdf(paste0(outpath,tissue,"_",analy, "_", test,"_DOWN_", d), width = 7.5, height = 5)  # Adjust width and height as needed
          
          grid.arrange(plot_list_down[[element_to_plot2]][[db_to_plot]])
          dev.off()
        }
            
        }

      }
      
      }else{
          for (sub in substructures){
            if(analy== "sub" & test =="one_cov_no_inter"){
                if(file.exists(paste0(inpath, tissue, "_", sub, "_DEGs.ORA.up_down_variable_DEGs.variable_DEGs.rds"))){
                print(paste0("-----", tissue, "------", analy,"-- ", sub, "-- ", test," ", sub, "------"))
                results2_prob<- readRDS(paste0(inpath, "/", tissue, "_", sub, "_DEGs.ORA.up_down_variable_DEGs.variable_DEGs.rds"))

              }else{
                next
              }
            
            
            }else if(analy== "inter" & test =="all_cov_inter"){
                if(file.exists(paste0(inpath, tissue,"_", sub, "_interaction_DEGs.ORA.up_down_variable_DEGs.variable_DEGs.rds"))){
                  print(paste0("-----",tissue, "------", analy,"-- ", sub, "-- ", test," ", sub, "------"))
                  results2_prob<- readRDS(paste0(inpath, "/", tissue,"_", sub, "_interaction_DEGs.ORA.up_down_variable_DEGs.variable_DEGs.rds"))
                
                }else{
                  next
                }
              
            }else{
              next
            }
            
            dbs <- c("GO:BP", "GO:MF","GO:CC", "KEGG", "ReactomePA","human_phenotype","H", "DO","DisGeNET", "gwas")
            if(test=="one_cov_no_inter"){
              ora_results<- lapply(results2_prob, function(sublist) {
              lapply(sublist, function(inner_list) {
                inner_list_df <- as.data.frame(inner_list)
                arrange(inner_list_df, inner_list_df[["p.adjust"]])
              })
            })
            top_30<- lapply(ora_results, function(sublist) {
              lapply(sublist, function(inner_list) {
                inner_list_df <- head(inner_list, 30)
              })
            })
            
            plot_list_up <- list()
            plot_list_down <- list()
            
            for (db in dbs) {
              for (element in names(top_30)){
                if ( "log_odds_ratio" %in% names(top_30[[element]][[db]])) {
                  if (element == "Upregulated"){
                    
                    # Plot for Age_Upregulated
                    plot_upregulated <- ggplot(top_30$Upregulated[[db]], aes(x = log_odds_ratio, y = reorder(Description, log_odds_ratio), color = p.adjust, size = Count)) +
                      geom_point() +
                      scale_color_gradient(low = "#8B0000", high = "#63B8FF") +
                      labs(x = "Log odds ratio", title = paste0("Upregulated - ", db), y = NULL) +
                      theme_minimal(base_size = 12)+
                      theme(axis.text.y = element_text(size = 13), 
                            axis.text.x = element_text(size = 13),
                            axis.title = element_text(size = 13), 
                            plot.title = element_text(size= 13, hjust=0.5),
                            legend.text = element_text(size=13),
                            panel.grid = element_blank(),                    # Remove all grid lines
                            axis.line.x = element_line(size = 0.5, color = "grey80"), # Add x-axis line
                            axis.line.y = element_line(size = 0.5, color = "grey80"),
                            panel.border = element_rect(color = "grey80", size = 0.5, fill = NA) # Add full border# Add y-axis line
                      )
                    plot_list_up[[element]][[db]] <- plot_upregulated
                  }else{
                    
                    # Plot for Age_Downregulated
                    plot_downregulated <- ggplot(top_30$Downregulated[[db]], aes(x = log_odds_ratio, y = reorder(Description, log_odds_ratio), color = p.adjust, size = Count)) +
                      geom_point() +
                      scale_color_gradient(low = "#8B0000", high = "#63B8FF") +
                      labs(x = "Log odds ratio", title = paste0("Downregulated - ", db), y = NULL) +
                      theme_minimal(base_size = 12)+
                      theme(axis.text.y = element_text(size = 13), 
                            axis.text.x = element_text(size = 13),
                            axis.title = element_text(size = 13), 
                            plot.title = element_text(size= 13, hjust=0.5),
                            legend.text = element_text(size=13),
                            panel.grid = element_blank(),                    # Remove all grid lines
                            axis.line.x = element_line(size = 0.5, color = "grey80"), # Add x-axis line
                            axis.line.y = element_line(size = 0.5, color = "grey80"),
                            panel.border = element_rect(color = "grey80", size = 0.5, fill = NA) # Add full border# Add y-axis line
                      )
                    plot_list_down[[element]][[db]] <- plot_downregulated
                  }
                }
                }
              }
              #
              
              
            
            
             element_to_plot1 <- names(plot_list_up)  # Use a valid index 'n'
            element_to_plot2 <- names(plot_list_down)  # Use a valid index 'n'
            
            # Example to use the first database, adjust as needed
            for (db_to_plot in c(dbs[1], dbs[4], dbs[10])){
              if(db_to_plot == dbs[1]){
                d<- "_gobp.pdf"
              }else if(db_to_plot == dbs[4]){
                d<- "_KEGG.pdf"
                
              }else{
                d<- "_gwas.pdf"
                
              }
              

              if (!is.null(element_to_plot1)){
                if(db_to_plot %in% names(plot_list_up[[element_to_plot1]])){
                                  print(paste0("-----Saving up enrichments for ", analy, " - ", sub))
                pdf(paste0(outpath,tissue,"_",analy, "_", test, "_",sub,"_UP_", d), width = 7, height = 5)  # Adjust width and height as needed
                
                grid.arrange(plot_list_up[[element_to_plot1]][[db_to_plot]])
                dev.off()
                }
                

              
              } 
              if( !is.null(element_to_plot2)){
                if( db_to_plot %in% names(plot_list_down[[element_to_plot2]])){
                                  print(paste0("-----Saving down enrichments for ", analy, " - ", sub))
                pdf(paste0(outpath,tissue,"_",analy, "_", test, "_",sub,"_DOWN_", d), width = 10, height = 7)  # Adjust width and height as needed
                
                grid.arrange(plot_list_down[[element_to_plot2]][[db_to_plot]])
                dev.off()
                }

                }
              }
                
        
              
            }
            
            else{
              ora_results<- lapply(results2_prob, function(inner_list) {
                  inner_list_df <- as.data.frame(inner_list)
                  arrange(inner_list_df, inner_list_df[["p.adjust"]])
              })
              top_30<- lapply(ora_results, function(inner_list) {
                  inner_list_df <- head(inner_list, 30)
              })
              
              plot_list<- list()

              for (db in dbs) {
                  if ( "log_odds_ratio" %in% names(top_30[[db]])) {

                      # Plot for Age_Upregulated
                      plot_upregulated <- ggplot(top_30[[db]], aes(x = log_odds_ratio, y = reorder(Description, log_odds_ratio), color = p.adjust, size = Count)) +
                        geom_point() +
                        scale_color_gradient(low = "#8B0000", high = "#63B8FF") +
                        labs(x = "Log odds ratio", title = paste0(db), y = NULL) +
                        theme_minimal(base_size = 12)+
                        theme(axis.text.y = element_text(size = 14), 
                              axis.text.x = element_text(size = 14),
                              axis.title = element_text(size = 14), 
                              plot.title = element_text(size= 14, hjust=0.5),
                              legend.text = element_text(size=14),
                              panel.grid = element_blank(),                    # Remove all grid lines
                              axis.line.x = element_line(size = 0.5, color = "grey80"), # Add x-axis line
                              axis.line.y = element_line(size = 0.5, color = "grey80"),
                              panel.border = element_rect(color = "grey80", size = 0.5, fill = NA) # Add full border# Add y-axis line
                        )
                      plot_list[[db]] <- plot_upregulated
                               
                    
                  }}
                #
              
              element_to_plot1 <- names(plot_list)  # Use a valid index 'n'

              # Example to use the first database, adjust as needed
              for (element_to_plot1 in c(dbs[1], dbs[4], dbs[10])){
                if(element_to_plot1 == dbs[1]){
                  d<- "_gobp.pdf"
                }else if(db_to_plot == dbs[4]){
                  d<- "_KEGG.pdf"
                  
                }else{
                  d<- "_gwas.pdf"
                  
                }
                
                if (length(plot_list) > 0) {
                  if(element_to_plot1 %in% names(plot_list)){
                    print(paste0("-----Saving enrichments for ", analy, " - ", sub))
                    pdf(paste0(outpath,tissue,"_",analy, "_", test, "_",sub,"_", d), width = 10, height = 7)  # Adjust width and height as needed
                    
                    grid.arrange(plot_list[[element_to_plot1]], ncol=1)
                    dev.off()
                    
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}




genes<- unlist(strsplit(results2_prob$Upregulated$`GO:BP`$geneID[1], "/"))
# Add a title
title <- textGrob("Ovary", gp = gpar(fontsize = 16, fontface = "bold"))
title_with_bg <- gTree(children = gList(rectGrob(gp = gpar(fill = "white", col = NA)), title))

# Combine the title and the plot
final_plot <- arrangeGrob(title_with_bg, combined_plot, heights = c(0.1, 1))
# Combine the title and the plot
plot(final_plot)
