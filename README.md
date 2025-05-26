# 2025 GTEx Menopause Project
Published in https://www.biorxiv.org/content/10.1101/2025.05.16.654406v1
## 🔍 Abstract
Female reproductive aging is a complex process with systemic health implications, yet how human aging unfolds across reproductive organs and tissues remains understudied. 
Here, we integrate deep learning–based analysis of 1,112 histological images with RNA-sequencing data from 659 samples across seven reproductive organs in 304 female donors aged 20–70
years. We reveal asynchronous aging dynamics: while ovary and vagina age gradually, the uterus undergoes an abrupt transition around menopause. 
Tissue segmentation highlights the myometrium as the most age-affected structure, marked by extracellular matrix remodeling and immune activation. 
Vaginal epithelial layers also show sharp menopausal shifts, mirroring transcriptional changes in developmental pathways. 
Transcriptomic analysis identifies thousands of age-associated genes, many with nonlinear trajectories enriched in extracellular matrix, angiogenesis, and oogenesis pathways. 
These findings position menopause as a key inflection point in reproductive aging and provide a multimodal, tissue-resolved map with implications for understanding gynecological
aging and menopause-associated disease risk.
![Git_image](https://github.com/user-attachments/assets/1e58ef44-31c2-4f04-b3a5-8d089cf19405)


## 📝 Table of Contents
### - 00.Data
- enrichments_dbs_files: databases' files for the enrichments performed in  0.7.DEA and enrichments.
- GWAS_datasets: reference datasets to assess traits in 08.GWAS.
- 2000_subset_balanced_tile_features_Uterus_1000.csv: subset of 2000 tiles from uterus with their values for the 384 features extracted.
- SIMULATED_Subject_Phenotypes.csv: simulated (not real!!) donors' data, use it to try the code. For confidentiality reasons, the true data cannot be shared.
- SIMULATED_Uterus_metadata.csv, SIMULATED_Uterus_counts.csv, SIMULATED_Uterus_tpms.csv: simulated (not real!!) donors' data for uterus (example tissue), use it to try the code. For confidentiality reasons, the true data cannot be shared.
- Uterus_final_filtered_images.rds, Vagina_final_filtered_images.rds: selected donors for the analyses of uterus and vagina samples.
- Uterus_pivot.csv: tissue proportions in uterus, converted to pivot coordinates.
- Deconvolution_uterus_pivot_v10.csv: cell-type proportions in uterus samples, converted to pivot coordinates.
- gencode.v39.annotation.bed: GENCODE annotation for genes.

### - 01. CNNs
- CNN_training.py, CNN_predict.py: VGG19-based CNN models' code for histological tissue tiles (training and predicting).
- pretrained_vgg19_bn.pt: pretrained model used.
- data_examples: uterus tiles for train, test, and external validation set + middle age.
      
### - 02. LIME (Local Interpretable Model-Agnostic Explanations)
- LIME.py: adaptation of LIME code to our purpose: identify tile areas that contribute to each classification category using our CNN models. LIME article available at: https://arxiv.org/abs/1606.05386.
- hybrid_tile_256_Vagina_trained_model_vgg19_bn_1024RSM.pt: trained model for vagina tiles.
- GTEX-11P81-2125_126.png: example of a vagina tile to run the code and perform the interpretation.
      
### - 03. Tissue segmentation
- tissue_segmentation.m: MATLAB code used for manually segmenting female reproductive organs into their constituent tissues, and then using them for KNN and later segmentation for all the WSIs. MATLAB version R2024b, https://es.mathworks.com/help/install/ug/install-products-with-internet-connection.html.
- GTEX-14PJ6-1625.svs: image example to segment.
- Segmentation: example of output from the segmentation for vagina.
      
### - 04. CellProfiler postprocessing
- newSecondaryFilteredVagina_Image.csv: epithelium measurements obtained from CellProfiler.
- epithelium_measurments_cellprofiler.R: analysis of epithelium measurements obtained with our CellProfiler pipeline.
CellProfiler v4.2.8, https://cellprofiler.org/
https://github.com/CellProfiler/CellProfilerprocessing of the measurements obtained from CellProfiler.
      
### - 05. Variance partition analysis
- variance_partition.R: analysis of the proportion of variation in image features attributable to demographic traits for each organ and tissue structure, while controlling for demographic variables and batch effects.
      
### - 06. MOFA (Multi-Omics Factor Analysis)
- MOFA.R: implementation of MOFA with two types of data modalities: image features and gene expression data.
- GSEA.R: code for Gene Set Enrichment Analysis for the ranked genes according to their contribution to factors.
MOFA original repository: https://github.com/bioFAM/MOFA2 
      
### - 07. DEA and enrichments 
- DEA.R. DEA_sliding_windows.R: Differential Expression Analysis (DEA) code with GTEx gene expression data (functions coded in DEA_and_DSA.R_functions.R).
- enrichments.R: functional enrichments against several databases (read some dbs files through Functional_enrichments.R).
- *all_cov_no_inter_AGE_covariates_and_traits.results.rds: output examples from the DEA for uterus, myometrium, vagina, ovary and breast, so they can be used for the enrichment analysis.
      
### - 08. GWAS
- gwas.R: GWAS overlap of the differentially expressed genes against different female reproductive GWAS.
      
### - 09. CODA (Compositional Data Analysis)
- coda_pivot_coordinates.R: implementation of CODA to test both cell-type and tissue proportion changes with age.

## System requirements
All the originally generated code has been run in R v4.4.1 or Python v3.11.5, and can be run on any operating system (Linux, macOS, Windows). The specific packages needed for each step are specified in the corresponding scripts.
