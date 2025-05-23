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
### - 01. CNNs
      VGG19-based CNN models' code for histological tissue tiles.
### - 02. LIME (Local Interpretable Model-Agnostic Explanations)
      Adaptation of LIME code to our purpose: identify tile areas that contribute to each classification category using our CNN models.
      LIME article available at: https://arxiv.org/abs/1606.05386 
### - 03. Tissue segmentation
      MATLAB code used for manually segmenting female reproductive organs into their constituent tissues, and then use them for KNN and later segmentation for all the WSIs.
      MATLAB version R2024b, https://es.mathworks.com/help/install/ug/install-products-with-internet-connection.html
### - 04. CellProfiler postprocessing
      Analysis of epithelium measurements obtained with our CellProfiler pipeline.
      CellProfiler v4.2.8, https://cellprofiler.org/
      https://github.com/CellProfiler/CellProfiler
### - 05. Variance partition analysis
      Analysis of the proportion of variation in image features attributable to demographic traits for each organ and tissue structure, while controlling for demographic variables and batch effects.
### - 06. MOFA (Multi-Omics Factor Analysis)
      Implementation of MOFA with two types of data modalities: image features and gene expression data. We share also the code for Gene Set Enrichment Analysis for the ranked genes according to their contribution to factors.
      MOFA original repository: https://github.com/bioFAM/MOFA2 
### - 07. DEA and enrichments 
      Differential Expression Analysis (DEA) code with GTEx gene expression data and functional enrichments against several databases.
### - 08. GWAS
      GWAS overlap of the differentially expressed genes against different female reproductive GWAS.
### - 09. CODA (Compositional Data Analysis)
      Implementation of CODA to test both cell-type and tissue proportion changes with age.

## System requirements
All the originally generated code has been run in R v4.4.1 or Python v3.11.5, and can be run on any operating system (Linux, macOS, Windows). The specific packages needed for each step are specified in the corresponding scripts.
