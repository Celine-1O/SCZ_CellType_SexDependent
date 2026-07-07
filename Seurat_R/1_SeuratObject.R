#' ############################################################################
#' File: 1_SeuratObject.R
#' Author: Celine I. Garcia Rodriguez
#' Date: 27/05/2026
#' 
#' Description:
#'  This script inputs and converts three datasets into a Seurat Object.
#'  Metadata across all datasets is strictly unified into: 
#'                orig.ident, Sex, Disease, Database.
#' 
#' Packages: Seurat, qs2, dplyr, reticulate
#' 
#' Functions: 
#'  - Brainscope_s1: Creates a Seurat Object from a different TXT files, 
#'                   combining matrices with their associated metadata.  
#'  - Zenodo_s1: Gets a large list of RDS matrices and generates a Seurat Object
#'              Metadata comes from an Excel file. 
#' ############################################################################


#| 1. Packages

library(Seurat)
library(qs2)
library(dplyr)

#| 2. SeuratObjects

#############################################
#                 BRAINSCOPE                #
#############################################

#' There is one TXT file per person. Previously, they have been divided into 
#' four different folders based on the Sex and Disease being:
#' - CON_F: Control females
#' - CON_M: Control males
#' - SCZ_F: Schizophrenia females
#' - SCZ_M: Schizophrenia males


BrainScope_s1 <- function(dir_path, sex, disease) {
  files_BS <- list.files(
    path = dir_path,
    pattern = "\\.txt$",
    full.names = TRUE
  )
  
  object_list <- lapply(files_BS, function(f) {
    m <- read.table(f, header = TRUE, sep = "\t", row.names = 1)
       
    # ID sample
    sample_id <- tools::file_path_sans_ext(basename(f))
    sample_id <- sub("-annotated_matrix$", "", sample_id)

    seu <- CreateSeuratObject(counts = m, assay = "RNA", min.cells = 3, project = sample_id)
    
    # Metadata
    seu$Sex <- as.factor(sex)
    seu$Disease <- as.factor(disease)
    seu$Database <- "BrainScope"
    
    return(seu)
  })
  
  return(object_list)
}

scz_f_list <- BrainScope_s1("/home/cgarcia/files/SCZ_descompr/SCZ_F", "Female", "Schizophrenia")
scz_m_list <- BrainScope_s1("/home/cgarcia/files/SCZ_descompr/SCZ_M", "Male", "Schizophrenia")
con_f_list <- BrainScope_s1("/home/cgarcia/files/SCZ_descompr/CON_F", "Female", "Control")
con_m_list <- BrainScope_s1("/home/cgarcia/files/SCZ_descompr/CON_M", "Male", "Control")

bs_objects <- c(scz_f_list, scz_m_list, con_f_list, con_m_list)
common_genes <- Reduce(intersect, lapply(bs_objects, function(x) rownames(x)))

# Filtering
bs_objects_filtered <- lapply(bs_objects, function(s) s[common_genes, ])

ids <- sapply(bs_objects_filtered, function(x) unique(x$orig.ident))

# Fusion
SCZ_Bs <- merge(
  x = bs_objects_filtered[[1]],
  y = bs_objects_filtered[-1],
  add.cell.ids = ids,
  project = "SCZ_BrainScope"
)
SCZ_Bs <- JoinLayers(SCZ_Bs, merge.assays = TRUE)

qs_save(SCZ_Bs,"SCZ_Bs.qs2")

#############################################
#                   ZENODO                  #
#############################################

# Matrix in RDS file and metadata in XLSX

Zenodo_s1 <- function(matrix_list, metadata) {
  
  sample_names <- names(matrix_list)
  sample_list <- lapply(sample_names, function(id) {
    
    matriz <- matrix_list[[id]] # Matrix extraction
    seu <- CreateSeuratObject(counts = matriz, min.cells = 3, project = id )
    
    # Metadata
    info_MB <- metadata[[id]]

    seu$Sex <- as.factor(info_MB[1])
    seu$Disease <- as.factor(info_MB[2])
    seu$Database <- "Zenodo"
    
    return(seu)
  })
  
  names(sample_list) <- sample_names
  return(sample_list)
}

dict_zenodo <- list(
  "MB6"    = c("Male", "Schizophrenia"),
  "MB7"    = c("Female", "Control"),
  "MB8"    = c("Male", "Schizophrenia"), 
  "MB8-2"  = c("Male", "Schizophrenia"),
  "MB9"    = c("Female", "Control"),
  "MB10"   = c("Female", "Schizophrenia"),
  "MB11"   = c("Male", "Control"),
  "MB12"   = c("Female", "Schizophrenia"),
  "MB13"   = c("Male", "Control"),
  "MB14"   = c("Female", "Schizophrenia"),
  "MB15"   = c("Female", "Control"),
  "MB16"   = c("Male", "Control"),
  "MB17"   = c("Male", "Control"),
  "MB18-2" = c("Male", "Control"),
  "MB19"   = c("Male", "Control"),
  "MB21"   = c("Female", "Control"),
  "MB22"   = c("Female", "Schizophrenia"),
  "MB23"   = c("Male", "Schizophrenia"),
  "MB50"   = c("Female", "Schizophrenia"), 
  "MB51"   = c("Male", "Control"),
  "MB52"   = c("Female", "Schizophrenia"), 
  "MB53"   = c("Female", "Control"),
  "MB54"   = c("Female", "Schizophrenia"),
  "MB55"   = c("Male", "Control"),
  "MB56"   = c("Male", "Schizophrenia"),
  "MB57"   = c("Female", "Control")
)

sRNA <- readRDS("snRNA-seq_raw_countmatrices.RDS")

zenodo_list <- Zenodo_s1(sRNA, dict_zenodo)

# Fusion
SCZ_Zen <- merge(
  x = zenodo_list[[1]],
  y = zenodo_list[-1],
  project = "SCZ_Zenodo"
)
SCZ_Zen <- JoinLayers(SCZ_Zen)

qs_save(SCZ_Zen, "SCZ_Zen_Obj.qs2")

#############################################
#                     GEO                   #
#############################################

#' GEO file is an Anndata Object (H5AD). Seurat v.5 is not natively programmed 
#' to transform this object into a Seurat format, nor are other packages such 
#' as ZellKonventer or AnndataR working due to environment issues. This is the 
#' alternative the author follows to resolve the problem related to Seurat 
#' and the Python version of the cluster.


# FROM LINUX CONSOLE

#| - module purge
#| - module load Python/3.9.6-GCCcore-11.2.0
#| - module load h5py/3.2.1-foss-2021b

#| module load R --> R:
#| - library(reticulate)
#| - library(Matrix)
#| - library(Seurat)
#| - library(qs2)
#| - use_python(Sys.which("python"), required = TRUE) 
#|        --> Indicate the prefered Python version. 
#| - py_config()
#| - py_module_available("h5py") 
#|         --> TRUE
#| - py_module_available("anndata") 
#|         --> TRUE
#| - anndata <- import("anndata")
#| - ad <- anndata$read_h5ad("/home/cgarcia/GEO/GSE254569_adata_RNA.h5ad") 
#|         --> class()= "anndata._core.anndata.AnnData" "python.builtin.object" 
#| - counts <- t(py_to_r(ad$X))
#| - counts <- as(counts, "CsparseMatrix")
#| - counts <- as(counts, "dgCMatrix")
#| - meta <- py_to_r(ad$obs)
#| 
#| - cell_names <- py_to_r(ad$obs_names$to_list())
#| - gene_names <- py_to_r(ad$var_names$to_list())
#| 
#| - rownames(counts) <- gene_names
#| - colnames(counts) <- cell_names
#| - rownames(meta) <- cell_names
#| 
#| - seurat_obj <- CreateSeuratObject(counts = counts,meta.data = meta,
#|                                    min.cells = 3, project = "SCZ_GEO")
#| - qs_save(seurat_obj, "/home/cgarcia/GEO/SCZ_GEO.qs2")

#| In Rstudio
SCZ_GEO <- qs_read("/home/cgarcia/GEO/SCZ_GEO.qs2")

SCZ_GEO <- subset(SCZ_GEO, Classification %in% c("Schizophrenia", "Control" )) 

SCZ_GEO$Disease <- as.factor(SCZ_GEO$Disease_Status)
SCZ_GEO$orig.ident <- SCZ_GEO$Donor
SCZ_GEO$Database <- "GEO"
metadata_GEO <- c("orig.ident", "nCount_RNA", "nFeature_RNA","Sex", "Disease", "Database")
SCZ_GEO@meta.data <- SCZ_GEO@meta.data[, metadata_GEO]

qs_save(SCZ_GEO, "/home/cgarcia/GEO/SCZ_GEO.qs2")


#| 3. Final Seurat Object

SCZ_Seurat <- merge(x = SCZ_Bs, y = c(SCZ_GEO, SCZ_Zen), project = "Schizophrenia_3")
#' This Seurat Object has three layers, one per each dataset.


#| 4. Saving Seurat Object
qs_save(SCZ_Seurat, "/home/cgarcia/SCZ_Seurat.qs2")




