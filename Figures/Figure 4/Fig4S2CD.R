#Fig 4S CD

#C 
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
df_age_ut<-data.frame(
  Tissue="Uterus",
  `Upregulated DEGs` =c(2313, 16, 589 ),
  `Downregulated DEGs` =c(1843,17, 579)
)
rownames(df_age_ut)<- c("Uterus-Age", "Myometrium","Myometrium-Age*")
colnames(df_age_ut)<- c("Tissue", "Upregulated DEGs","Downregulated DEGs")
df_age_ut$Sample_size<- c(130, 130, 61)


df_age_ov<-data.frame(
  Tissue="Ovary",
  `Upregulated DEGs`  =c(968, 196, 43),
  `Downregulated DEGs` =c(1448, 176, 2)
)

rownames(df_age_ov)<- c("Ovary-Age", "Cortex","Cortex    Age")
colnames(df_age_ov)<- c("Tissue", "Upregulated DEGs","Downregulated DEGs")
df_age_ov$Sample_size<- c(169)

df_age_vg<-data.frame(
  Tissue="Vagina",
  `Upregulated DEGs`  =c(29, 54, 0 ),
  `Downregulated DEGs` = c(18, 50,0)
)

rownames(df_age_vg)<- c("Vagina-Age", "Epithelium","Epithelium    Age")
colnames(df_age_vg)<- c("Tissue", "Upregulated DEGs","Downregulated DEGs")
df_age_vg$Sample_size<- c(146)

df_age_br<-data.frame(
  Tissue="Breast",
  `Upregulated DEGs`  =c(8, 52, 136, 62),
  `Downregulated DEGs` = c(2, 9, 180, 37)
)

rownames(df_age_br)<- c("Breast-Age", "Adipocytes", "Lobules",  "Stroma")
colnames(df_age_br)<- c("Tissue", "Upregulated DEGs","Downregulated DEGs")
df_age_br$Sample_size<- c(157)

total_heatmap<-rbind(df_age_ut, df_age_ov, df_age_vg, df_age_br)
df_age<-total_heatmap
df_age$Color<-ifelse(df_age$Tissue== "Uterus",  "#6DB6FFFF", 
                     ifelse(df_age$Tissue== "Ovary",  "#009292FF", 
                            ifelse(df_age$Tissue== "Vagina",  "#B66DFFFF", "#490092FF")))

# Create a named vector for colors
sample_colors <- setNames(df_age$Color, df_age$Tissue)

row_annotation <- rowAnnotation(
  "Sample size" = anno_barplot(
    df_age$Sample_size,
    gp = gpar(fill = sample_colors, col = sample_colors),  # Dynamic color filling
    border = FALSE,
    bar_width = 1 
  ),
  gap = unit(0.25, "cm"),
  show_legend = TRUE, 
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 10)
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
pdf("Desktop/Figures/fig4S2_C.pdf", width = 4.5, height = 5 )
draw(heatmap, heatmap_legend_side = "top")
dev.off()

pdf("Desktop/TFM/bsc83671/GTEx_v8/Laura/Figure_plots/fig4S2_C.pdf", width = 4.5, height = 5) 
draw(heatmap, heatmap_legend_side = "top")
dev.off()



##D- upset
genes<- list()
test<- "all_cov_no_inter"
###uterus, vagina, breast, ovary
tissues<-c("Uterus", "Ovary", "Vagina", "BreastMammaryTissue")
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
  VaginaDown = genes$Vagina$Down,
  BreastMammaryTissueUp = genes$BreastMammaryTissue$Up,
  BreastMammaryTissueDown = genes$BreastMammaryTissue$Down
  # # EndocervixUp = genes$CervixEndocervix$Up,
  # # EndocervixDown = genes$CervixEndocervix$Down,
  # # EctocervixUp = genes$CervixEctocervix$Up,
  # # EctocervixDown = genes$CervixEctocervix$Down,
  # # FallopianTubeUp = genes$FallopianTube$Up,
  # # FallopianTubeDown = genes$FallopianTube$Down
)
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
         sets.bar.color = c("#6DB6FFFF", "#6DB6FFFF", "#009292FF","#009292FF","#B66DFFFF","#B66DFFFF", "#490092FF", "#490092FF"), 
         order.by = "freq", 
         sets = colnames(binary_matrix),
         keep.order = TRUE,
         text.scale = 1.5)


pdf("Desktop/Figures/fig4S2_D.pdf", width =6.5, height = 4 )
m
dev.off()

pdf("Desktop/TFM/bsc83671/GTEx_v8/Laura/Figure_plots/fig4S2_D.pdf", width = 5, height = 5) 
m
dev.off()

