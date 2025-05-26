.libPaths(c("/x/bscxxxx/Rlibs", .libPaths()))

library(limma)
library(argparse)
library(variancePartition)
library(ggplot2)
library(vroom)
library(parallel)
library(dplyr)
library(stringr)


tissue <- "Uterus"
##in/out directories
tile_features =vroom(paste0("../00.Data/2000_subset_balanced_tile_features_",tissue, "_1000.csv"))

# Subsampling the dataframe
downsampling <- function(tile_features, n_tiles) {
  tile_features %>%
    group_by(slide_name, label) %>%
    slice_sample(n = n_tiles) %>%
    ungroup()
}

#Downsample for uterus 
#tile_features =vroom('/X/Ole/paper/for_plotting/Uterusdownsampled_uterus.csv')


mean_features <- function(tile_features){
  # Compute the mean of the features grouped by slide_name and label
  tile_features %>%
  select(-tile_number) %>% # Remove the tile_number column
  group_by(slide_name, label) %>% # Group by slide_name and label
  summarise(across(starts_with("feature"), mean, na.rm = TRUE)) # Compute mean for feature columns

}

median_features <- function(tile_features){
  # Compute the mean of the features grouped by slide_name and label
  tile_features %>%
  select(-tile_number) %>% # Remove the tile_number column
  group_by(slide_name, label) %>% # Group by slide_name and label
  summarise(across(starts_with("feature"), median, na.rm = TRUE)) # Compute mean for feature columns

}


tile_features <- mean_features(tile_features)
# View the result
dim(tile_features)
endometrium_slides <- tile_features$slide_name[tile_features$label == "endometrium"]

#SUbset the images to the only_myometrium ones in case we need it
only_myometrium_id <- readRDS('/Uterus_filtered_myometrium_images.rds')

only_myometrium_id <- only_myometrium_id$Subject.ID

#Remove observations that have endometrium for donors with only myometrium
#tile_features <- tile_features[!(str_extract(tile_features$slide_name, "GTEX-\\w+") %in% only_myometrium_id & tile_features$label == "endometrium"), ]

#Keep only donors that have both endometrium & myometrium
tile_features <- tile_features[!(str_extract(tile_features$slide_name, "GTEX-\\w+") %in% only_myometrium_id), ]


plot_dir <- paste0("/variance_partition/plots/downsampled_uterus","/", tissue)

# Create the directory if it does not exist
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)  # recursive=TRUE ensures parent directories are also created
}

results_dir <- paste0("/variance_partition/results/downsampled_uterus","/", tissue)

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)  # recursive=TRUE ensures parent directories are also created
}



GTEx_Subject_Phenotypes.GRU <- read.csv('../00.Data/SIMULATED_Subject_Phenotypes.GRU.csv', sep = '\t', header = TRUE)

#metadata files, as well as Subject_Phenotypes.GRU.csv contains Donor's metadata

if (tissue %in% c("Ectocervix", "Endocervix")) {
  metadata <- readRDS(paste0("../00.Data/Cervix", tissue, "/metadata.rds"))
} else {
  metadata <- readRDS(paste0("../00.Data/SIMULATED_", tissue, "_metadata.rds"))
}


#tile_features <- vroom(paste0("/X/Allal/RNAPath/tile_features_balanced_classes_100/balanced_tile_features_Ovary_no_follicles.csv"))


# function for GTEx identifiers
GTEX.namematching <- function(NamesRows){
  snames <- c()
  for (j in NamesRows){
    if (substr(j, start = 10, stop = 10)=="-"){
      snames<- c(snames, substr(j, start = 1, stop = 9))
    }
    else{
      snames<- c(snames, substr(j, start = 1, stop = 10))
    }
  }
  snames
}

#View(tile_features)
tile_features$gtex = GTEX.namematching(tile_features$slide_name)
#tiles_pca = prcomp(tile_features[,4:387])
tile_features$age = GTEx_Subject_Phenotypes.GRU$AGE[match(tile_features$gtex, GTEx_Subject_Phenotypes.GRU$SUBJID)]
#tile_features$HardyScale = GTEx_Subject_Phenotypes.GRU$DTHHRDY[match(tile_features$gtex, GTEx_Subject_Phenotypes.GRU$SUBJID)]

#smaller validation set 
if (tissue %in% c("Uterus", "Vagina", "Ovary")){
  filtered_images = readRDS(paste0("../00.Data/", tissue,"_final_filtered_images.rds"))
  tile_features = tile_features[tile_features$gtex %in% filtered_images$Subject.ID,]
}

atrophicValidationDonors = c("GTEX-1399S", "GTEX-1PIEJ", "GTEX-1R9JW", "GTEX-QDT8", "GTEX-S7SF")
if (tissue=="BreastMammaryTissue"){
  atrophicValidationDonors= c("GTEX-1399S", "GTEX-1MA7W", "GTEX-1PIEJ", "GTEX-1QPFJ", "GTEX-QDT8",  "GTEX-S7SF",  "GTEX-WRHK")
}

tile_features = tile_features[!tile_features$gtex %in% atrophicValidationDonors,]
tile_features$age = GTEx_Subject_Phenotypes.GRU$AGE[match(tile_features$gtex, GTEx_Subject_Phenotypes.GRU$SUBJID)]
tile_features$HardyScale = GTEx_Subject_Phenotypes.GRU$DTHHRDY[match(tile_features$gtex, GTEx_Subject_Phenotypes.GRU$SUBJID)]
tile_features$IschemicTime = GTEx_Subject_Phenotypes.GRU$TRDNISCH[match(tile_features$gtex, GTEx_Subject_Phenotypes.GRU$SUBJID)]
tile_features$bmi = GTEx_Subject_Phenotypes.GRU$BMI[match(tile_features$gtex, GTEx_Subject_Phenotypes.GRU$SUBJID)]
ancestry_file<- read.table(paste0("SIMULATED_admixture_inferred_ancestry.txt"))
tile_features$ancestry = ancestry_file$Ancestry[match(tile_features$gtex,ancestry_file$Donor)]
subtissues = unique(tile_features$label)
traits = c('bmi', 'ancestry', 'age', 'IschemicTime', 'HardyScale', 'gtex')

param = SnowParam(40, "SOCK", progressbar=TRUE)
register(param)

tile_features <- na.omit(tile_features)

num_cores <- parallel::detectCores() - 1  # Use all available cores minus one for system processes

process_subtissue <- function(subtissue) {
  print(subtissue)
  
  data_all <- tile_features[tile_features$label == subtissue, c(3:ncol(tile_features))] 
  data <- data_all[, 1:384]
  df <- data_all[, traits]
  
  print('create model')
  colnames(df) <- c('bmi', 'ancestry', 'age', 'IschemicTime', 'HardyScale', "Donor")
  for_formula <- c('age', 'bmi', 'ancestry', 'HardyScale', 'IschemicTime')
  metadata_2 <- df
  form <- paste0("~", paste0(for_formula, collapse = "+"))
  print(form)
  
  varPart <- fitExtractVarPartModel(t(data), form, metadata_2)
  vp <- sortCols(varPart)
  pl <- plotVarPart(vp)
  
  saveRDS(vp, file = paste0("/variance_partition/results/downsampled_uterus","/", tissue, "/", subtissue, ".rds"))
  
  ggsave(
    paste0("/variance_partition/plots/downsampled_uterus","/", tissue, "/", subtissue, "_plot.svg"),
    plot = pl,
    device = "svg",
    width = 10,
    height = 8,
    dpi = 300
  )
}

mclapply(subtissues, process_subtissue, mc.cores = num_cores)

