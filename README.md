# Cell-Type-Resolved-Analysis-of-Sex-Dependent-Transcriptional-Alterations-in-Schizophrenia
-------------------------------------------------------------------------------------------
## Master's Thesis in Bioinformatics - Universitat de València
***************************************************************
This repository contains all the code and resources used for the single nuclei RNA sequencing (snRNA-seq) analysis of schizophrenia (SCZ) tissue and healthy control (CTRL) tissue samples, focusing on sex differences and using the workflow of `Seurat`, from R, and `Scanpy`, from Python. 


____________________________________________________________________________________________
### Metadata
The raw datasets are available in the following repositories:
- **[GSE255469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254569)**: From GEO public repository.
    - CTRL_Female: 22
    - CTRL_Male: 25
    - SCZ_Female: 13
    - SCZ_Male: 13
    - **Total nuclei**: 637.154
      
- **[BrainScope](https://brainscope.gersteinlab.org/data/snrna_expr_matrices_zip/SZBDMulti-Seq.zip)**: From BrainScope repository.
    - CTRL_Female: 12
    - CTRL_Male: 12
    - SCZ_Female: 12
    - SCZ_Male: 12
    - **Total nuclei**: 420.903

- **[Zenodo6921620](https://zenodo.org/records/6921620)**: Data related to [_Upper cortical layer-driven network impairment in schizophrenia_](https://www.science.org/doi/10.1126/sciadv.abn8367).
    - CTRL_Female: 5
    - CTRL_Male: 8
    - SCZ_Female: 5
    - SCZ_Male: 6
    - **Total nuclei**: 209.053

*******************************************************************************************
### Packages

#### R (4.4.2)
| Package | Version |
| :--- | :--- |
| `dplyr` | 1.2.1 |
| `ggplot2` | 4.0.3 | 
| `harmony`| 1.3.0 |
| `Matrix`| 1.7.5 |
| `reticulate`| 1.46.6 |
| `Seurat` | 5.5.0 |
| `SeuratIntegrate` | 0.4.1 |
| `SeuratWrappers` | 0.1.0 |
 
#### Python (3.11.4)

| **Package** | **Version** |
| :--- | :--- |
| `anndata` | 0.12.6 |
| `gseapy` | 1.2.1 | 
| `harmonypy`| 2.0.0 |
| `matplotlib-base` | 3.10.8 |
| `networkx` | 3.6.1 |
| `numpy` | 2.2.5 |
| `pandas` | 3.0.1 |
| `scanpy` | 1.11.5 |
| `scib-metrics` | 0.5.9 |
| `scipy` | 1.17.1 |
| `scvi-tools` | 1.4.2 |
| `scikit-learn` | 1.8.0 |
| `seaborn` | 0.13.2 |


