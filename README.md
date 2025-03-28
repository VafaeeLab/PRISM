# Multi-Omics Prognostic Marker Discovery and Survival Modeling: A Case Study on Pan-Cancer Survival Analyses in Women's Cancers

Survival analysis is essential for predicting patient outcomes and guiding personalized cancer treatments. While multi-omics data offers valuable insights, its high dimensionality complicates analysis and clinical application. Many studies still rely on the traditional Cox proportional hazards model, with limited exploration of alternative survival algorithms or robust feature selection methods. Few frameworks effectively integrate features across multiple omics modalities. To address these issues, we developed PRISM (PRognostic marker Identification and Survival Modelling through Multi-omics Integration), a comprehensive framework designed to improve survival predictions and identify key prognostic markers. PRISM systematically compares various feature selection methods and survival models, and employs a robust pipeline that selects features from single-omics data, integrating them through feature-level fusion and multi-stage refinement. Applied to TCGA data for Breast Invasive Carcinoma (BRCA), Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma (CESC), Ovarian Serous Cystadenocarcinoma (OV) and Uterine Corpus Endometrial Carcinoma (UCEC), PRISM demonstrated that integrating DNA methylation, miRNA, and copy number variation data provided complementary information, significantly outperforming other modality combinations in three cancers (BRCA: C-index 0.77; CESC: 0.80; UCEC: 0.76). Pan-cancer analysis further revealed shared oncogenic pathways and therapeutic targets. PRISM provides a scalable, generalizable solution for multi-omics integration, advancing cancer research and precision medicine. 

![image](https://github.com/user-attachments/assets/67e2bb8e-19ea-4038-9f6f-5084e87272d1)

## **Computational Environment**  
All analyses were conducted on the **University of New South Wales (UNSW) Katana HPC** and **Gadi (NCI)**, utilizing **16 CPUs and 124 GB of memory per job** to process large-scale omics data efficiently.  

- While optimized for HPC, the pipeline can also run on **local machines** by reducing data size or enabling **parallelization**.  
- The current parallelization code assumes **16 cores**.  

## **Pipeline Overview**  
To ensure reproducibility, we recommend running the scripts in the following order:  

### **1. TCGA Data Extraction & Preprocessing**  
📌 *Extract and preprocess omics data (GE, ME, CNV, DM) for each cancer type.*  
- **Scripts:** `TCGA_BRCA.R`, `TCGA_CESC.R`, `TCGA_OV.R`, `TCGA_UCEC.R`  
- **Note:** Due to API changes from the GDC Portal (as of *March 17, 2025*), the `GDCquery_clinic` function in `TCGAbiolinks` may not work. You may need to manually download clinical data from the GDC Portal. For convenience, we have also provided preprocessed data in [`Data`](./Data) folder of this repository.
- **Output:** `BRCA/`,`CESC/`,`OV/`, `UCEC/` etc. These directories for each cancer type will contain the preprocessed omics data as a csv. For example BRCA would have:
  `BRCA_GE_data.csv`, `BRCA_ME_data.csv`, `BRCA_METH_data.csv`, `BRCA_CNV_data.csv`. Where the first three columns (excluding the index) are the case_id, survival information, and the rest are the features. The rows will be the samples.

### **2. Evaluation of Individual Filter Methods**  
📌 *Evaluate individual feature selection methods.*  
- **Script:** `individual_selection.R`
- **Inputs:** This requires you to have `BRCA/`,`CESC/`,`OV/`, `UCEC/` directories with corresponding `*_GE_data.csv`, `*_ME_data.csv`, `*_METH_data.csv`, `*_CNV_data.csv` omics data inside.
- **Outputs:** Results for c-index & features selection for each modality per cancer type. `*_cindex_results.csv`, `*_feature_results.csv`,  `*_cindex_heatmap.pdf`,  `*_feature_heatmap.pdf`.

### **3. Cross-Validation Feature Selection**  
📌 *Run cross-validation-based feature selection to identify key prognostic features.*  
- **Script:** `CV_method.R`
- **Inputs:** This requires you to have `BRCA/`,`CESC/`,`OV/`, `UCEC/` directories with corrosponding `*_GE_data.csv`, `*_ME_data.csv`, `*_METH_data.csv`, `*_CNV_data.csv` omics data inside.
- **Outputs:** Creates `ME/`,`GE/`,`CNV/`, `METH/` omics subdirectories inside of `BRCA/`,`CESC/`,`OV/`, `UCEC/`. Where we have a csv of the features selected by CV `features_cv.csv`, as well as the results of performance against no feature selected `results_with_fs` `results_without_fs`. 

### **4. Multi-Omics Integration - One-Stage vs. Two-Stage Refinement**  
📌 *Compare one-stage and two-stage refinement for multi-omics feature selection.*  
- **Scripts:** `first_stage_refinement.R`, `second_stage_refinement.R`
- **Inputs:** This requires you to have features selected by CV `features_cv.csv` inside of `BRCA/`,`CESC/`,`OV/`, `UCEC/`.
- **Outputs:** Creates a `LF/` directory inside of `ME/`,`GE/`,`CNV/`, `METH/` for Two-stage, will make `EF/` for First_stage, where every modality combintation is stored as a csv.

### **5. Pan-Cancer Downstream Analysis**  
📌 *Perform pan-cancer analyses to identify shared pathways and therapeutic targets.*  
- **Script:** `pan-cancer.R`
- **Inputs:** This requires you to have every modality combintation is stored in `LF/`, in particular `BRCA/LF/BRCA_METH_ME_CNV_data.csv`, `OV/LF/OV_ME_METH_data.csv`, `CESC/LF/CESC_METH_ME_CNV_data.csv`, `UCEC/LF/UCEC_METH_ME_CNV_data.csv`
- **Outputs:**
  1. Pancancer signature overlap - `upset_plot.png`
  2. miRNA disease-assoications across all cancers - `miRNA_Disease_Associations.pdf`
  3. Disease gene enrichment network - `emapplot_output.pdf`
  4. Go Terms & KEGG Pathway Enrichment plots
  
  


