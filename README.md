# Multi-Omics Prognostic Marker Discovery and Survival Modelling: A Case Study on Multi-cancer Survival Analysis of Women’s Specific Tumours

Survival analysis plays a critical role in predicting patient outcomes and guiding personalized cancer therapies. Although multi-omics data provide rich biological insights, their high dimensionality poses significant challenges for robust analysis and clinical implementation. While many studies rely on the traditional Cox proportional hazards model, few have explored alternative survival algorithms combined with rigorous feature selection to identify low-dimensional, clinically feasible prognostic signatures that retain strong predictive power comparable to models using the full feature set. To address these gaps, we developed PRISM (PRognostic marker Identification and Survival Modelling through Multi-omics Integration), a comprehensive framework aimed at improving survival prediction and discovering minimal yet robust biomarker panels across multiple omics modalities. PRISM systematically evaluates various feature selection methods and survival models through a robust pipeline that selects features within single-omics datasets before integrating them via feature-level fusion and multi-stage refinement. Applied to TCGA cohorts of Breast Invasive Carcinoma (BRCA), Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma (CESC), Ovarian Serous Cystadenocarcinoma (OV), and Uterine Corpus Endometrial Carcinoma (UCEC), PRISM revealed that cancer types benefit from unique combinations of omics modalities reflecting their molecular heterogeneity. Notably, miRNA expression consistently provided complementary prognostic information across all cancers, enhancing integrated model performance (C-index: BRCA 0.698, CESC 0.754, UCEC 0.754, OV 0.618). PRISM advances cancer prognosis by delivering scalable, interpretable multi-omics integration and identifying concise biomarker signatures with performance comparable to full-feature models, promoting clinical feasibility and precision oncology.

![image](https://github.com/user-attachments/assets/67e2bb8e-19ea-4038-9f6f-5084e87272d1)

## **Computational Environment**  
All analyses were conducted on the **University of New South Wales (UNSW) Katana HPC** and **Gadi (NCI)**, utilizing **16 CPUs and 124 GB of memory per job** to process large-scale omics data efficiently.  

- While optimized for HPC, the pipeline can also run on **local machines** by reducing data size or enabling **parallelization**.  
- The current parallelization code assumes **16 cores**.  

## **Pipeline Overview**  
To ensure reproducibility, we recommend running the scripts in the following order:  

### **1. TCGA Data Download & Preprocessing**  
📌 *Extract and preprocess omics data (GE, ME, CNV, DM) for each cancer type.*  
- **Scripts:** `TCGA_download.R`, `TCGA_prepreprocessing.R`, 
- **Note:** Data is downloaded from TCGA Hub from UCSC Xena (https://xenabrowser.net/datapages/)
- **Output:** `BRCA/`,`CESC/`,`OV/`, `UCEC/` etc. These directories for each cancer type will contain the preprocessed omics data as a csv. For example BRCA would have:
  `BRCA_GE_clean.csv`, `BRCA_ME_clean.csv`, `BRCA_DM_clean.csv`, `BRCA_CNV_clean.csv`. Where the first three columns are the sample_ID, survival information, and the rest are the features. The rows will be the samples.

### **2. Evaluation of Individual Filter Methods**  
📌 *Evaluate individual feature selection methods.*  
- **Script:** `individual_selection.R`
- **Inputs:** This requires you to have `BRCA/`,`CESC/`,`OV/`, `UCEC/` directories with corresponding `*_GE_clean.csv`, `*_ME_clean.csv`, `*_DM_clean.csv`, `*_CNV_clean.csv` omics data inside.
- **Outputs:** Results for c-index & features selection for each modality per cancer type, for example `CESC_CNV_plots/`, would have files:
   `*_cindex_results.csv`, `*_feature_results.csv`,  `*_cindex_heatmap.pdf`,  `*_feature_heatmap.pdf`.  `*_summary_table.csv`

### **3. Cross-Validation & Bootstrapping Feature Selection**  
📌 *Run both feature selection pipelines to identify key prognostic features.*  
- **Script:** `feature_selection_pipelines.R`
- **Inputs:** This requires you to have `BRCA/`,`CESC/`,`OV/`, `UCEC/` directories with corrosponding `*_GE_clean.csv`, `*_ME_clean.csv`, `*_DM_clean.csv`, `*_CNV_clean.csv`  omics data inside.
- **Outputs:** Creates `ME/`,`GE/`,`CNV/`, `DM/` omics subdirectories inside of `BRCA/`,`CESC/`,`OV/`, `UCEC/`. Where we have a csv of the features selected by CV `features_cv` & BS `features_bs`, as well as the results of performance against no feature selected `results_with_fs` `results_without_fs`. 

### **4. Multi-Omics Integration - One-Stage vs. Two-Stage Refinement**  
📌 *Compare one-stage and two-stage refinement for multi-omics feature selection.*  
- **Scripts:** `first_stage_refinement.R`, `second_stage_refinement.R`
- **Inputs:** This requires you to have features selected by CV `features_cv.csv` inside of `BRCA/`,`CESC/`,`OV/`, `UCEC/`.
- **Outputs:** Creates a `2S/` directory inside of `ME/`,`GE/`,`CNV/`, `DM/` for Two-stage, will make `1S/` for First_stage, where every modality combintation is stored as a csv.

### **5. Comparative-Cancer Downstream Analysis**  
📌 *Perform comparative-cancer analyses to identify shared pathways and therapeutic targets.*  
- **Script:** `pcomparative_cancer_analysis.R`
- **Inputs:** This requires you to have every modality combintation is stored in `1S/` or `2S`.
- **Outputs:**
  1. Gene Target signature overlap - `pan-cancer_gene_target_overlaps.pdf`
  2. miRNA disease-assoications across all cancers - `miRNA_Disease_Associations.pdf`
  3. Disease gene enrichment network - `emapplot_output_miRNA.pdf`
  4. Go Terms & KEGG Pathway Enrichment plots
  
  


