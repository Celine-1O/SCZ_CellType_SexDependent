#' ############################################################################
#' File: 2_Preprocessing.R
#' Author: Celine I. Garcia Rodriguez
#' Date: 27/05/2026
#' 
#' Description:
#'  This script preprocesses Seurat Object which contents the three datasets. 
#'  Quality Control, Normalization, Feature Selection, Scale Data
#' 
#' Packages: Seurat, qs2
#' 
#' ############################################################################

#| 1. Packages
library(Seurat)
library(qs2)

SCZ_Seurat <- qs_read("/home/cgarcia/SCZ_Seurat.qs2")

#| 2. Quality Control
SCZ_Seurat[["percent.mt"]] <- PercentageFeatureSet(SCZ_Seurat, pattern = "^MT-") 
SCZ_Seurat <- subset(SCZ_Seurat, subset = nFeature_RNA > 200 & percent.mt < 5)

#| 3. Normalization
counts_median <- median(SCZ_Seurat$nCount_RNA)
SCZ_Seurat <- NormalizeData(SCZ_Seurat, 
                            normalization.method = "LogNormalize", 
                            scale.factor = counts_median
)


#| 4. HVGs (Highly Variable Genes)
SCZ_Seurat <- FindVariableFeatures(SCZ_Seurat, 
                                    selection.method = "vst", 
                                    nfeatures = 2000
)
                                    
#| 5. Scaling
SCZ_Seurat <- ScaleData(SCZ_Seurat, features = VariableFeatures(SCZ_Seurat))

#| 6. Saving file
qs_save(SCZ_Seurat, "/home/cgarcia/SCZ_Seurat_preprocessed.qs2")

