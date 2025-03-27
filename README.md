# 2025 GTEx Menopause Project
Published in (insert bioarxiv link when submitted)
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
      ####VGG19-based CNN models' code for histological tissue tiles.
### - 02. LIME (Local Interpretable Model-Agnostic Explanations)
      Adaptation of LIME code to our purpose: identify tile areas that contribute to each classification category using our CNN models.
### - 03. Tissue segmentation
      Matlab code used for manually segmenting female reproductive organs into their constituent tissues and then fine-tune ViT-S/8 with them.
### - 04. CellProfiler postprocessing
      Analysis of epithelium measurements obtained with our CellProfiler pipeline.
### - 05. MOFA (Multi-Omics Factor Analysis)
      Implementation of MOFA with two types of data modalities: image features and gene expression data. We share also the code for Gene Set Enrichment Analysis for the ranked genes according to their contribution to factors.
### - 06. DEA and enrichments 
      Differential Expression Analysis (DEA) code with GTEx gene expression data and functional enrichments against several databases.
### - 07. GWAS
      GWAS overlap of the differentially expressed genes against different female reproductive GWAS.
### - 08. CODA (Compositional Data Analysis)
      Implementation of CODA to test both cell-type and tissue proportion changes with age.

## 🔍 Abstract
Female reproductive aging is a complex process with systemic health implications, yet how human aging unfolds across reproductive organs and tissues remains understudied. Here, we integrate deep learning–based analysis of 1,112 histological images with RNA-sequencing data from 659 samples across six reproductive organs in 304 female donors aged 20–70 years. We reveal asynchronous aging dynamics: while ovarian and vaginal tissues age gradually, the uterus undergoes an abrupt transition around menopause. Tissue segmentation highlights the myometrium as the most age-affected structure, marked by extracellular matrix remodeling and immune activation. Vaginal epithelial layers also show sharp menopausal shifts, mirroring transcriptional changes in developmental pathways. Transcriptomic analysis identifies thousands of age-associated genes, many with nonlinear trajectories enriched in extracellular matrix, angiogenesis, and oogenesis pathways. These findings position menopause as a key inflection point in reproductive aging and provide a multimodal, tissue-resolved map with implications for understanding gynecological aging and menopause-associated disease risk.

