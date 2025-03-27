#Fig 4AB
library(ggplot2)
library(reshape2) 
library(gridExtra)  
library(colorBlindness)
library(tidyr)
library(ComplexHeatmap)
library(cowplot)
library(patchwork)  # For combining the plots
library(RColorBrewer)
library(dplyr)
library(paletteer)
palette <- paletteer_d("colorBlindness::paletteMartin")


#Samples sizes for images and RNA data
pal <- c("Uterus" = "#6DB6FFFF", "Ovary" = "#009292FF", "Vagina" = "#B66DFFFF", "BreastMammaryTissue" = "#490092FF", "CervixEndocervix"= "#FF6DB6FF", "CervixEctocervix"="#FFB6DBFF", "FallopianTube"="#DB6D00FF")
tissues<-c("Uterus", "Ovary", "Vagina", "BreastMammaryTissue", "CervixEndocervix", "CervixEctocervix", "FallopianTube")
samples_list<-list()
probs_list<-list()
age_file<-read.csv("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/00.Data/GTEx_Subject_Phenotypes.GRU.csv", sep="\t", header=TRUE)
all_tissues_results<-readRDS("Desktop/TFM/bsc83671/GTEx_v8/Laura/12.Tissues_substructures/DEA/DEA_v10/Table_all_substructures_no_inter_results.rds")

###only age

df_age<- data.frame(
  Tissue= c("Uterus", "Myometrium", "Ovary", "Vagina", "Breast"),
  `Upregulated DEGs` =c(2313, 589, 968, 29, 8),
  `Downregulated DEGs` =c(1843,579,1448, 18, 2),
  Sample_size= c(130,61, 169, 146, 157)
)
colnames(df_age)<- c("Tissue", "Up-DEGs","Down-DEGs","Sample_size" )
pal <- c("Uterus" = "#6DB6FFFF", "Myometrium"= "#6DB6FFFF", "Ovary" = "#009292FF", "Vagina" = "#B66DFFFF", "Breast" = "#490092FF",  "Endocervix"= "#FF6DB6FF", "Ectocervix"="#FFB6DBFF", "Fallopian Tube"="#DB6D00FF")

row_annotation <- rowAnnotation(
  "Sample size" = anno_barplot(
    sapply(df_age$Tissue, function(tissue) {
      sample_size <- df_age[df_age$Tissue == tissue, "Sample_size"]
      if (length(sample_size) == 0) return(0)
      return(sample_size)
    }),
    gp = gpar(fill = pal, col = pal),  # Fill colors per tissue type
    border = FALSE,
    bar_width = 1,
    axis_param = list(
      side = "bottom",  # Ensures labels appear at the bottom
      gp = gpar(fontsize = 10)  # Increase x-axis tick label size
    )
  ),
  gap = unit(0.25, "cm"),
  show_legend = TRUE, 
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 12)
)
column_annotation <- columnAnnotation(
  Direction =  c("Up-DEGs", "Down-DEGs"),
  col= list(Direction=c("Up-DEGs" = "#63B8FF", 
                        "Down-DEGs" = "#8B0000")),
  # annotation_legend_param = list(title = "Direction", labels=c("Downregulated DEGs", "Upregulated DEGs" ), 
  #                                labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 14) ),
  show_annotation_name = FALSE,
  show_legend=FALSE
)


heatmap_data<-df_age[,2:3]
rownames(heatmap_data)<- df_age$Tissue
colnames(heatmap_data)<-c("Up-DEGs", "Down-DEGs")

color_palette<-c("Uterus"= "#6DB6FFFF", "Ovary" = "#009292FF", "Vagina" = "#B66DFFFF", "BreastMammaryTissue" = "#490092FF", "CervixEndocervix"= "#FF6DB6FF", "CervixEctocervix"="#FFB6DBFF", "FallopianTube"="#DB6D00FF")

heatmap <- Heatmap(
  heatmap_data,
  name = "#DEGs",
  col = brewer.pal(9, "BuPu"),
  bottom_annotation =column_annotation,
  
  na_col = "white",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_rot = 0,
  column_names_centered = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 14),
  column_names_gp = gpar(fontsize = 14, just= "center"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(prettyNum(heatmap_data[i, j], big.mark = ","), x, y, gp = gpar(fontsize = 14))
  },
  height = unit(8, "cm"),
  
  show_row_names = TRUE,
  heatmap_legend_param = list(
    direction = "horizontal",
    legend_title_position = "topright",
    labels_gp = gpar(fontsize = 14),
    title = "#DEGs",
    legend_height = unit(3, "cm"),
    legend_width = unit(4, "cm"),
    grid_width = unit(0.5, "cm")
  ),
  left_annotation = row_annotation
)
draw(heatmap, heatmap_legend_side = "top")

# Draw the heatmap
heatmaptot<- grid.grabExpr(draw(heatmap, heatmap_legend_side = "top"))
pdf("Desktop/Figures/Fig4A.pdf", width = 4.5, height = 5) 
draw(heatmap, heatmap_legend_side = "top")
dev.off()

svg("Desktop/Figures/Fig4A.svg", width = 4.5, height = 5, pointsize = 12)
draw(heatmap, heatmap_legend_side = "top")
dev.off()

pdf("Desktop/TFM/bsc83671/GTEx_v8/Laura/Figure_plots/Fig4A.pdf", width = 4.5, height = 5) 
draw(heatmap, heatmap_legend_side = "top")
dev.off()

svg("Desktop/TFM/bsc83671/GTEx_v8/Laura/Figure_plots/Fig4A.svg", width = 4.5, height = 5, pointsize = 12)
draw(heatmap, heatmap_legend_side = "top")
dev.off()



####Fig 4B --upset plot
genes<- list()
test<- "all_cov_no_inter"
###uterus, vagina, breast, ovary
tissues<-c("Uterus", "Ovary", "Vagina", "BreastMammaryTissue", "CervixEndocervix", "CervixEctocervix", "FallopianTube")
#tissues<-c("CervixEndocervix", "CervixEctocervix", "FallopianTube")
for (tissue in tissues[1:3]){
  inpath<-paste0("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/12.Tissues_substructures/DEA/DEA_v10/", tissue, "/")
  all<- readRDS( paste0(inpath, test, "/", tissue, "_", test, "_AGE_covariates_and_traits.results.rds"))
  all<-all$Age
  print(paste0("#-------------- ", tissue, " - ", test, " - ", analy, " --------------#"))
  degs<- all[all$adj.P.Val<0.05,]
  neg<-(degs[degs$logFC<0,])
  pos<-(degs[degs$logFC>0,])
  
  genes[[tissue]]$Up<-rownames(pos)
  genes[[tissue]]$Down<-rownames(neg)
}

library(UpSetR)
gene_sets <- list(
  UterusUp = genes$Uterus$Up,
  UterusDown = genes$Uterus$Down,
  OvaryUp = genes$Ovary$Up,
  OvaryDown = genes$Ovary$Down,
  VaginaUp = genes$Vagina$Up,
  VaginaDown = genes$Vagina$Down
  # BreastMammaryTissueUp = genes$BreastMammaryTissue$Up,
  # BreastMammaryTissueDown = genes$BreastMammaryTissue$Down
  # # EndocervixUp = genes$CervixEndocervix$Up,
  # # EndocervixDown = genes$CervixEndocervix$Down,
  # # EctocervixUp = genes$CervixEctocervix$Up,
  # # EctocervixDown = genes$CervixEctocervix$Down,
  # # FallopianTubeUp = genes$FallopianTube$Up,
  # # FallopianTubeDown = genes$FallopianTube$Down
)
top<-intersect(intersect(genes$Uterus$Up, genes$Ovary$Up), genes$Vagina$Up)
all_genes <- unique(unlist(gene_sets))

binary_matrix <- sapply(gene_sets, function(set) all_genes %in% set)

# Convert TRUE/FALSE to 1/0
binary_matrix <- as.data.frame(apply(binary_matrix, 2, as.integer))
# Set row names to be the gene names
rownames(binary_matrix) <- all_genes
# Create the upset plot
theme_update(
  axis.text = element_text(color = "black", size = 14),  # Change axis text color
  axis.title = element_text(color = "black", size = 14)  # Change axis title color
)
# Generate UpSet plot
m<-upset(binary_matrix, 
         main.bar.color = "black", 
         sets.bar.color = c("#6DB6FFFF", "#6DB6FFFF", "#009292FF","#009292FF","#B66DFFFF","#B66DFFFF"), 
         order.by = "freq", 
         sets = colnames(binary_matrix),
         keep.order = TRUE,
         text.scale = 1.5)

library(ggplot2)
library(dplyr)

pdf("Desktop/Figures/Fig4B.pdf", width = 5.5, height = 4.5) 
m
dev.off()

svg("Desktop/Figures/Fig4B.svg", width = 5.5, height = 4.5, pointsize = 12)
m
dev.off()

pdf("Desktop/TFM/bsc83671/GTEx_v8/Laura/Figure_plots/Fig4.Bpdf", width = 5.5, height = 4.5) 
m
dev.off()

svg("Desktop/TFM/bsc83671/GTEx_v8/Laura/Figure_plots/Fig4B.svg", width = 5.5, height = 4.5, pointsize = 12)
m
dev.off()


